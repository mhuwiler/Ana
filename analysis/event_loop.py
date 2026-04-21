"""Unified event loop: histogram booking, running, and materializing."""

import sys
import time
from dataclasses import dataclass

# Acceptance baseline cut applied to both data and MC before trigger selection.
# Exported so scripts can compute data_acc = data_df.Filter(BASE_CUT).
BASE_CUT = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

from analysis.cutflow import book_cutflow_actions, extract_and_print_cutflow
from analysis.histograms import (
    build_groups, book_mc_denom, book_mc_per_trigger,
    book_2d_histograms, book_gen_histograms, book_exclusive_histograms,
    materialize_hists, materialize_gen_hists,
)
from utils.data import book_lumi_actions, extract_lumi


@dataclass
class EventLoopResult:
    """All histogram results produced by book_and_run."""
    mc_hists_by_trig: dict
    mc_denom_hists: dict
    mc_items_by_trig: dict
    lumi: float
    gen_hists_by_trig: dict
    mc_hists_excl_by_trig: dict
    cutflow_tables: dict
    trig_selections: dict   # kept for Phase 4/5 data histogram booking
    sig_indices: list
    bkg_indices: list
    h2d_book: dict          # lazy 2D histogram pointers (for plotting)
    excl_sig_mHH_ptrs: dict # {excl_trig_name: [TH1 ptrs]} for exclusive trigger efficiency


def book_and_run(
    mc, data_df, args,
    ana_cfg,
    trig_list, excl_trig_list,
    plot_vars, gen_plot_vars, plot_vars_2d,
    brilcalc_default, lumi,
    sig_samples, dy_samples, tt_samples, qcd_samples,
    plot_dir,
    cuts_by_trig=None,
):
    """Book all analysis actions, run one unified event loop, and return results.

    Parameters
    ----------
    mc : dict[str, RDataFrame]
        MC RDataFrames with kinematics columns and "w" weight defined.
    data_df : RDataFrame or None
        Data RDataFrame with kinematics columns, or None for MC-only.
    args : argparse.Namespace
        Runtime flags (max_mc_files, no_data, etc.).
    ana_cfg : AnalysisConfig
        Decay modes, bkg_modes, sig_modes.
    trig_list : list of (name, expr)
        Inclusive trigger definitions.
    excl_trig_list : list of (name, expr)
        Exclusive trigger definitions.
    plot_vars : list
        1D plot variable definitions (PlotVar or tuple).
    gen_plot_vars : list
        Gen-level plot variable definitions.
    plot_vars_2d : list
        2D histogram definitions.
    brilcalc_default : dict
        Brilcalc data for lumi extraction.
    lumi : float or None
        Pre-computed lumi (from --no-data path); None = extract from data.
    sig_samples, dy_samples, tt_samples, qcd_samples : list[str]
        Sample name lists.
    plot_dir : str
        Output plot directory (used by cutflow extraction).

    Returns
    -------
    EventLoopResult
    """
    import ROOT
    import gc

    decay_modes = ana_cfg.decay_modes
    bkg_modes = ana_cfg.bkg_modes
    sig_modes = ana_cfg.sig_modes
    if cuts_by_trig is None:
        cuts_by_trig = {}
    active_modes = bkg_modes + sig_modes

    # mc here is PRE-CUT — cuts are applied per-trigger below
    mc_base = {name: df.Filter(BASE_CUT) for name, df in mc.items()}
    data_base = data_df.Filter(BASE_CUT) if data_df is not None else None

    proc_to_samples = {
        "DY": dy_samples, "TT": tt_samples,
        "HHbbtt": sig_samples, "QCD": qcd_samples,
    }

    # ── Book all actions ──
    unified_ptrs = []
    trig_selections = {}

    if data_df is not None:
        lumi_run_take, lumi_ls_take = book_lumi_actions(data_df)
        unified_ptrs.extend([lumi_run_take, lumi_ls_take])
    else:
        lumi_run_take = lumi_ls_take = None

    for trig_name, trig in trig_list:
        print(f"Booking actions for {trig_name} ...")

        # Apply trigger filter
        
        if trig is not None:
            data_sel = data_base.Filter(trig) if data_base is not None else None
            mc_sel = {name: df.Filter(trig) for name, df in mc_base.items()}
        else:
            data_sel = data_base
            mc_sel = dict(mc_base)

        for name, df in mc_sel.items(): 
            mc_sel[name] = mc_sel[name].Define("prescaleWeight", "prescaleWeight(L1_HTT200er, L1_HTT255er, L1_HTT280er, L1_HTT320er, L1_HTT360er, L1_HTT400er, L1_HTT450er, L1_ETT2000, L1_SingleJet180, L1_SingleJet200, L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)").Define("finalWeight", "w*prescaleWeight")

        # Apply per-trigger cuts from cuts_by_trig
        trig_cuts = cuts_by_trig.get(trig_name, [])
        for _, cut_expr in trig_cuts:
            if data_sel is not None:
                data_sel = data_sel.Filter(cut_expr)
            mc_sel = {name: df.Filter(cut_expr) for name, df in mc_sel.items()}

        sig_mHH_ptrs = [
            mc_sel[s].Histo1D(
                (f"h_mHH_{s}_{trig_name}",
                 "m_{HH} gen-level;m_{HH} [GeV];Events", 24, 0, 1200),
                "gen_mHH", "finalWeight")
            for s in sig_samples
        ]
        unified_ptrs.extend(sig_mHH_ptrs)

        trig_selections[trig_name] = {
            "trig": trig, "data_sel": data_sel, "mc_sel": mc_sel,
            "sig_mHH_ptrs": sig_mHH_ptrs,
        }

    # Book cutflow actions (per-trigger, using pre-cut mc with progressive cuts)
    phase1_data = None
    has_cutflow = any(len(cuts_by_trig.get(tn, [])) > 0 for tn, _ in trig_list)
    if has_cutflow:
        phase1_data = book_cutflow_actions(
            trig_list, cuts_by_trig, mc, data_df, mc_base, data_base,
            decay_modes, active_modes, sig_samples, dy_samples, tt_samples,
            qcd_samples, unified_ptrs)

    mc_denom_groups, _ = build_groups(active_modes, decay_modes, proc_to_samples, mc_base)
    denom_book = book_mc_denom(mc_base, mc_denom_groups, plot_vars, unified_ptrs)

    mc_items_by_trig, histo_books, sig_indices, bkg_indices, _ = \
        book_mc_per_trigger(trig_selections, plot_vars, decay_modes, active_modes,
                            proc_to_samples, sig_samples, sig_modes, bkg_modes,
                            unified_ptrs)

    h2d_book = book_2d_histograms(trig_selections, plot_vars_2d, sig_samples, unified_ptrs)
    gen_histo_book = book_gen_histograms(trig_selections, gen_plot_vars,
                                          sig_samples, sig_modes, unified_ptrs)
    excl_histo_books = book_exclusive_histograms(
        excl_trig_list, mc_base, plot_vars, active_modes, decay_modes,
        proc_to_samples, unified_ptrs)

    # Book exclusive gen_mHH histograms for signal (trigger efficiency studies)
    excl_sig_mHH_ptrs = {}
    for trig_name, trig in excl_trig_list:
        if trig is None:
            continue
        ptrs = [
            mc_base[s].Filter(trig).Define("prescaleWeight", "prescaleWeight(L1_HTT200er, L1_HTT255er, L1_HTT280er, L1_HTT320er, L1_HTT360er, L1_HTT400er, L1_HTT450er, L1_ETT2000, L1_SingleJet180, L1_SingleJet200, L1_DoubleJet30er2p5_Mass_Min250_dEta_Max1p5, L1_DoubleJet30er2p5_Mass_Min300_dEta_Max1p5, L1_DoubleJet30er2p5_Mass_Min330_dEta_Max1p5)").Define("finalWeight", "w*prescaleWeight").Histo1D(
                (f"h_mHH_{s}_{trig_name}",
                 "m_{HH} gen-level;m_{HH} [GeV];Events", 24, 0, 1200),
                "gen_mHH", "finalWeight")
            for s in sig_samples
        ]
        unified_ptrs.extend(ptrs)
        excl_sig_mHH_ptrs[trig_name] = ptrs

    # ── Run unified event loop ──
    print(f"\nBooked {len(unified_ptrs)} total actions")
    print("Running unified event loop ...", flush=True)
    sys.stdout.flush()
    t0 = time.time()
    ROOT.RDF.RunGraphs(unified_ptrs)
    print(f"Unified event loop done in {time.time() - t0:.1f}s")

    if lumi is None and lumi_run_take is not None:
        print("  Extracting lumi from Take vectors ...", flush=True)
        t0_lumi = time.time()
        lumi, data_runs, n_missing = extract_lumi(
            lumi_run_take, lumi_ls_take, brilcalc_default)
        print(f"  Lumi: {len(data_runs)} runs, {n_missing} not in brilcalc, "
              f"extracted in {time.time() - t0_lumi:.1f}s")
        print(f"  Data lumi (LS-matched) = {lumi:.4f} fb^-1")

    if phase1_data is not None:
        cutflow_tables = extract_and_print_cutflow(
            trig_list, phase1_data, decay_modes, active_modes,
            dy_samples, tt_samples, qcd_samples, sig_samples,
            lumi, plot_dir)
    else:
        print("[--skip-cutflow] Skipping cutflow extraction")
        cutflow_tables = {}

    # Materialize histograms
    mc_denom_hists = materialize_hists(denom_book, len(mc_denom_groups), lumi)
    mc_hists_by_trig = {}
    for trig_name in trig_selections:
        mc_hists_by_trig[trig_name] = materialize_hists(
            histo_books[trig_name], len(mc_items_by_trig[trig_name]), lumi)
    gen_hists_by_trig = materialize_gen_hists(gen_histo_book, sig_modes, lumi)
    mc_hists_excl_by_trig = {}
    for trig_name, histo_book in excl_histo_books.items():
        if not histo_book:
            continue
        n_groups = len(next(iter(histo_book.values())))
        mc_hists_excl_by_trig[trig_name] = materialize_hists(histo_book, n_groups, lumi)

    # Materialize 2D histograms (unscaled; lumi scaling happens at plot time)
    h2d_hists = {
        pname: [(tn, s, ptr.GetValue().Clone()) for tn, s, ptr in entries]
        for pname, entries in h2d_book.items()
    }

    del unified_ptrs, histo_books, excl_histo_books, denom_book, gen_histo_book, phase1_data
    gc.collect()

    return EventLoopResult(
        mc_hists_by_trig=mc_hists_by_trig,
        mc_denom_hists=mc_denom_hists,
        mc_items_by_trig=mc_items_by_trig,
        lumi=lumi,
        gen_hists_by_trig=gen_hists_by_trig,
        mc_hists_excl_by_trig=mc_hists_excl_by_trig,
        cutflow_tables=cutflow_tables,
        trig_selections=trig_selections,
        sig_indices=sig_indices,
        bkg_indices=bkg_indices,
        h2d_book=h2d_hists,
        excl_sig_mHH_ptrs=excl_sig_mHH_ptrs,
    )
