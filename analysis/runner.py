"""Framework core: setup() and load_and_run() for ultra-thin user scripts.

Usage in a user script::

    from analysis.runner import setup, load_and_run
    from analysis.variables import PlotVar

    SAMPLES  = ["HHbbtt", "DY", "TT", "QCD"]
    DATA     = True
    NTHREADS = 4

    PLOT_VARS = [PlotVar("ak4_pt0", r"Leading jet $p_T$ [GeV]", 50, 0, 250), ...]

    ctx = setup(__file__, samples=SAMPLES, data=DATA, nthreads=NTHREADS)
    result = load_and_run(ctx, PLOT_VARS)
    result.plot.stacked(ratio="significance")
"""

import os
from dataclasses import dataclass, field
from types import SimpleNamespace

from analysis.setup import init_root, load_macros
from utils.logging import init_logging


# ═══════════════════════════════════════════════════════════════════════════════
#  AnalysisContext — returned by setup()
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class AnalysisContext:
    """All configuration and state from setup()."""
    ana_cfg: object               # AnalysisConfig
    trig_list: list               # [(name, expr)] filtered to user's TRIGGERS
    excl_trig_list: list          # [(name, expr)] exclusive triggers
    brilcalc_default: dict        # parsed brilcalc data for default trigger
    plot_dir: str                 # output directory for plots
    samples: list                 # user's SAMPLES list (group names)
    data: bool
    nthreads: int
    nplot_workers: int
    theme: str
    overwrite: bool
    max_mc_files: int             # 0 = use MAX_EVENTS limit
    max_data_files: int           # 0 = all data files
    max_events: object            # None = use default from samples.yaml
    cuts: object                  # None = use cuts.yaml, [] = no cuts, [...] = custom
    skim: bool                    # True = auto-skim to slim files, False = always use raw EOS
    ana_dir: str                  # repo root directory


# ═══════════════════════════════════════════════════════════════════════════════
#  setup()
# ═══════════════════════════════════════════════════════════════════════════════

def setup(script_file, samples, triggers=None, data=True,
          nthreads=4, nplot_workers=8, theme="light", overwrite=False,
          max_mc_files=0, max_data_files=0, max_events=None, cuts=None,
          skim=True):
    """One-line framework init. Returns AnalysisContext.

    Parameters
    ----------
    script_file : str
        Pass ``__file__`` from the user script.
    samples : list[str]
        Group names to load (e.g. ["HHbbtt", "DY", "TT", "QCD"]).
    triggers : list[str] or None
        Trigger names to use (e.g. ["DST_JetHT"]). None = all triggers.
    data : bool
        Whether to load scouting data.
    nthreads : int
        Number of ImplicitMT threads.
    theme : str
        "light" or "dark" plot theme.
    overwrite : bool
        Overwrite existing plot files.
    max_mc_files : int
        Max files per MC sample (0 = use MAX_EVENTS limit from config).
    max_events : int or None
        Max events cap (None = use default from samples.yaml). When both
        max_mc_files and max_events are set, whichever gives fewer files wins.
    cuts : list[str] or None
        C++ filter expressions applied to all samples. None = use cuts.yaml,
        empty list [] = no cuts.
    skim : bool
        If True (default), auto-skim to slim ROOT files for faster re-runs.
        If False, always load from raw EOS and run Define chains (no disk writes).
    """
    init_logging()
    init_root(nthreads=nthreads)
    load_macros()

    from analysis.config import load_analysis_config, load_triggers
    from utils.data import parse_brilcalc
    from utils.plotting import setup_style

    ana_cfg = load_analysis_config()
    trig_list_all, excl_trig_list_all, brilcalc_files, brilcalc_fallback = load_triggers()

    # Extract triggers — from CUTS tuple format or explicit TRIGGERS list
    if triggers is not None:
        _trigger_set = set(triggers)
    elif isinstance(cuts, list) and cuts and isinstance(cuts[0], tuple):
        # Extract unique trigger names from CUTS tuples
        _trigger_set = set(trig for trig, _ in cuts)
    else:
        _trigger_set = None

    if _trigger_set is not None:
        trig_list = [(n, e) for n, e in trig_list_all if n in _trigger_set]
        excl_trig_list = [(n, e) for n, e in excl_trig_list_all
                          if n.replace("_excl", "") in _trigger_set]
    else:
        trig_list = trig_list_all
        excl_trig_list = excl_trig_list_all

    # Load brilcalc
    brilcalc = {}
    for trig_key, path in brilcalc_files.items():
        if os.path.exists(path):
            brilcalc[trig_key] = parse_brilcalc(path)
            total = sum(v["lumi"] for v in brilcalc[trig_key].values())
            print(f"Loaded brilcalc [{trig_key}]: {len(brilcalc[trig_key])} runs, "
                  f"total = {total:.3f} fb^-1")
        else:
            print(f"WARNING: brilcalc file not found for {trig_key}: {path}")
    for fb_key, fb_src in brilcalc_fallback.items():
        if fb_key not in brilcalc and fb_src in brilcalc:
            brilcalc[fb_key] = brilcalc[fb_src]
            print(f"  (using {fb_src} brilcalc as fallback for {fb_key})")

    brilcalc_default = brilcalc.get("DST_JetHT", {})

    # Plot directory
    ana_dir = os.path.dirname(os.path.abspath(script_file))
    script_stem = os.path.splitext(os.path.basename(script_file))[0]
    plot_dir = os.path.join(ana_dir, "output", "plots", script_stem)
    os.makedirs(plot_dir, exist_ok=True)
    print(f"Plots will be saved to: {plot_dir}")

    setup_style(dark=(theme == "dark"))

    return AnalysisContext(
        ana_cfg=ana_cfg,
        trig_list=trig_list,
        excl_trig_list=excl_trig_list,
        brilcalc_default=brilcalc_default,
        plot_dir=plot_dir,
        samples=samples,
        data=data,
        nthreads=nthreads,
        nplot_workers=nplot_workers,
        theme=theme,
        overwrite=overwrite,
        max_mc_files=max_mc_files,
        max_data_files=max_data_files,
        max_events=max_events,
        cuts=cuts,
        skim=skim,
        ana_dir=ana_dir,
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  load_and_run()
# ═══════════════════════════════════════════════════════════════════════════════

def load_and_run(ctx, plot_vars, vars_2d=None, cm_2d=None):
    """Load samples, book histograms, run event loop, return PlotResult.

    Parameters
    ----------
    ctx : AnalysisContext
        From setup().
    plot_vars : list[PlotVar]
        1D plot variable definitions.
    vars_2d : list[tuple] or None
        2D histogram definitions. Each entry is a tuple of two tuples:
        ``((x_name, x_label, nx, x0, x1), (y_name, y_label, ny, y0, y1))``.
    cm_2d : list[tuple] or None
        Confusion matrix definitions. Same format as vars_2d — each entry is
        a tuple of two tuples. The two variables are booked as 1D histograms
        internally, and the confusion matrix is built from their bin counts.

    Returns
    -------
    PlotResult
    """
    from analysis.variables import filter_vars, PlotVar2D
    from analysis.data_loading import load_mc_samples, load_data
    from analysis.event_loop import book_and_run, BASE_CUT
    from analysis.plots import book_data_histograms
    from utils.data import (
        ls_nanoaod_files_groups, BASE, GROUPS, XSEC, MAX_EVENTS,
        DATA_YEARS, DATA_RUNS,
    )
    from utils.skim import rdf_exprs

    # Filter GROUPS to user's requested samples
    filtered_groups = {k: v for k, v in GROUPS.items() if k in ctx.samples}
    if not filtered_groups:
        print(f"WARNING: No matching sample groups for {ctx.samples}")
        print(f"  Available groups: {list(GROUPS.keys())}")

    # Build an args-like namespace for load_mc_samples / load_data compatibility
    args = SimpleNamespace(
        max_mc_files=ctx.max_mc_files,
        all_events=False,
        no_data=not ctx.data,
        max_data_files=ctx.max_data_files,
    )

    # Load MC
    _, group_files_by_sample = ls_nanoaod_files_groups(
        base_dir=BASE, groups=filtered_groups, verbose=0)
    max_events = ctx.max_events if ctx.max_events is not None else MAX_EVENTS
    print(f"\nPreparing MC weights  MAX_EVENTS = {max_events}")
    mc, mc_files_map, dy_samples, tt_samples, sig_samples, qcd_samples = load_mc_samples(
        group_files_by_sample, XSEC, max_events, args, rdf_exprs=rdf_exprs,
        skim=ctx.skim)

    # Load data
    data_df, lumi_precomputed = load_data(
        args, brilcalc_default=ctx.brilcalc_default,
        data_years=DATA_YEARS, data_runs=DATA_RUNS)

    # Convert vars_2d tuples to PlotVar2D for backwards compat with book_2d_histograms
    plot_vars_2d = []
    if vars_2d:
        for x_def, y_def in vars_2d:
            x_name, x_label, nx, x0, x1 = x_def
            y_name, y_label, ny, y0, y1 = y_def
            plot_vars_2d.append(PlotVar2D(
                name=f"{x_name}_vs_{y_name}",
                xvar=x_name, yvar=y_name,
                xlabel=x_label, ylabel=y_label,
                nx=nx, x0=x0, x1=x1,
                ny=ny, y0=y0, y1=y1,
            ))

    # Add confusion matrix variables to plot_vars for 1D booking
    from analysis.variables import PlotVar
    cm_var_names = set()
    if cm_2d:
        for entry in cm_2d:
            _cm_name, x_def, y_def = entry
            for var_def in (x_def, y_def):
                vname = var_def[0]
                if vname not in cm_var_names and vname not in {v.name for v in plot_vars}:
                    cm_var_names.add(vname)
                    plot_vars = list(plot_vars) + [PlotVar(vname, var_def[1], var_def[2], var_def[3], var_def[4])]

    # Build per-trigger cut map from CUTS
    from analysis.data_loading import resolve_cuts

    cuts_by_trig = {}  # {trig_name: [(cut_name, expr), ...]}
    if isinstance(ctx.cuts, list) and ctx.cuts and isinstance(ctx.cuts[0], tuple):
        # New tuple format: [("NoTrigger", ["common", "tauhtauh"]), ...]
        for trig_name, cut_cards in ctx.cuts:
            named_cuts, _src = resolve_cuts(cut_cards)
            cuts_by_trig[trig_name] = named_cuts
            print(f"\n[cuts] {trig_name}: {len(named_cuts)} cuts from {_src}")
            for cut_name, cut_expr in named_cuts:
                tag = f"{cut_name}: " if cut_name else ""
                print(f"  \u2192 {tag}{cut_expr}")
    else:
        # Old flat format — same cuts for all triggers
        named_cuts, _src = resolve_cuts(ctx.cuts)
        if named_cuts:
            print(f"\n[cuts] Applying {len(named_cuts)} cuts from {_src}:")
            for cut_name, cut_expr in named_cuts:
                tag = f"{cut_name}: " if cut_name else ""
                print(f"  \u2192 {tag}{cut_expr}")
        for trig_name, _ in ctx.trig_list:
            cuts_by_trig[trig_name] = named_cuts

    # Run event loop
    result = book_and_run(
        mc, data_df, args,
        ctx.ana_cfg,
        ctx.trig_list, ctx.excl_trig_list,
        plot_vars, [], plot_vars_2d,
        ctx.brilcalc_default, lumi_precomputed,
        sig_samples, dy_samples, tt_samples, qcd_samples,
        ctx.plot_dir,
        cuts_by_trig=cuts_by_trig,
    )

    # Book data histograms if data is loaded
    data_hists_by_trig = None
    data_denom_hists = None
    if ctx.data and data_df is not None:
        data_acc = data_df.Filter(BASE_CUT)
        data_hists_by_trig, data_denom_hists = book_data_histograms(
            data_acc, result.trig_selections, plot_vars)

    return PlotResult(
        event_result=result,
        ctx=ctx,
        plot_vars=plot_vars,
        vars_2d=vars_2d or [],
        plot_vars_2d=plot_vars_2d,
        cm_2d=cm_2d or [],
        data_hists_by_trig=data_hists_by_trig,
        data_denom_hists=data_denom_hists,
        sig_samples=sig_samples,
        dy_samples=dy_samples,
        tt_samples=tt_samples,
        qcd_samples=qcd_samples,
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  PlotResult + Plotter
# ═══════════════════════════════════════════════════════════════════════════════

class Plotter:
    """Plot accessor attached to PlotResult as ``.plot``."""

    def __init__(self, result):
        self._r = result

    def stacked(self, ratio=None, vars=None, triggers=None):
        """Stacked MC histogram for each variable x trigger.

        Parameters
        ----------
        ratio : str or None
            "efficiency", "significance", or "cuml_significance" — adds a ratio panel below.
        vars : list[str] or None
            Variable names to plot. None = all.
        triggers : list[str] or None
            Trigger names to plot. None = all.
        """
        from analysis.plots import plot_mc_distributions
        r = self._r
        plot_mc_distributions(
            r.trig_selections, r.mc_hists_by_trig, r.mc_denom_hists,
            r.mc_items_by_trig, r.gen_hists_by_trig,
            r.plot_vars, [],
            r.sig_indices, r.bkg_indices,
            r.ctx.ana_cfg.sig_modes, r.ctx.ana_cfg.decay_modes,
            r.lumi, r.ctx.plot_dir, r.ctx.overwrite,
            do_stacked=True, do_shape=False,
            do_efficiency=(ratio == "efficiency"),
            do_significance=(ratio == "significance"),
            do_cuml_significance=(ratio == "cuml_significance"),
            n_workers=r.ctx.nplot_workers)

    def shapes(self, vars=None, triggers=None):
        """Shape overlay (normalized to unity) for each variable x trigger."""
        from analysis.plots import plot_mc_distributions
        r = self._r
        plot_mc_distributions(
            r.trig_selections, r.mc_hists_by_trig, r.mc_denom_hists,
            r.mc_items_by_trig, r.gen_hists_by_trig,
            r.plot_vars, [],
            r.sig_indices, r.bkg_indices,
            r.ctx.ana_cfg.sig_modes, r.ctx.ana_cfg.decay_modes,
            r.lumi, r.ctx.plot_dir, r.ctx.overwrite,
            do_stacked=False, do_shape=True,
            n_workers=r.ctx.nplot_workers)

    def data(self, ratio=None, vars=None, triggers=None):
        """Stacked MC + data overlay for each variable x trigger.

        Parameters
        ----------
        ratio : str or None
            "data_mc" — adds data/MC ratio panel below.
        """
        from analysis.plots import plot_data_mc
        r = self._r
        if r.data_hists_by_trig is None:
            print("[plot.data] No data loaded — skipping")
            return
        plot_data_mc(
            r.trig_selections, r.mc_hists_by_trig, r.mc_denom_hists,
            r.data_hists_by_trig, r.data_denom_hists,
            r.mc_items_by_trig, r.plot_vars,
            r.sig_indices, r.bkg_indices,
            r.lumi, r.ctx.plot_dir, r.ctx.overwrite,
            do_stacked=True, do_shape=True,
            do_efficiency=(ratio == "efficiency"),
            do_significance=(ratio == "significance" or ratio == "data_mc"),
            do_cuml_significance=(ratio == "cuml_significance"))

    def trigger_overlays(self, vars=None):
        """Same variable across triggers on one plot."""
        from analysis.plots import plot_trigger_overlays
        r = self._r
        plot_trigger_overlays(
            r.trig_selections, r.mc_hists_by_trig, r.mc_hists_excl_by_trig,
            r.plot_vars, r.sig_indices,
            r.lumi, r.ctx.plot_dir, r.ctx.overwrite, do_shape=True)

    def hist2d(self):
        """2D histograms from PLOT_VARS_2D."""
        from analysis.plots import plot_gen_overlays
        r = self._r
        if not r.plot_vars_2d:
            return
        plot_gen_overlays(
            r.trig_selections, r.gen_hists_by_trig,
            [], r.ctx.ana_cfg.sig_modes, r.ctx.ana_cfg.decay_modes,
            r.h2d_book, r.plot_vars_2d,
            r.lumi, r.ctx.plot_dir, r.ctx.overwrite)

    def confusion_matrix(self):
        """Confusion matrices from PLOT_CM_2D. Each entry produces a 2x2 seaborn heatmap."""
        import numpy as np
        import matplotlib.pyplot as plt
        from utils.plotting import plot_confusion_matrix, th1_to_np
        r = self._r

        if not r.cm_2d:
            return

        first_trig = next(iter(r.mc_hists_by_trig))
        h_by_var = r.mc_hists_by_trig[first_trig]

        def _get_correct_frac(var_name):
            """Get fraction of correct predictions (bin[1]) from 1D histogram."""
            if var_name not in h_by_var:
                print(f"[plot.confusion_matrix] {var_name} not booked — skipping")
                return None
            h_list = h_by_var[var_name]
            h_sum = h_list[0].Clone("_tmp_cm")
            for h in h_list[1:]:
                h_sum.Add(h)
            _, vals, _ = th1_to_np(h_sum)
            total = vals.sum()
            if total <= 0:
                return 0.5
            correct = vals[1] if len(vals) >= 2 else 0
            return correct / total

        for entry in r.cm_2d:
            cm_name, x_def, y_def = entry
            x_name, x_label = x_def[0], x_def[1]
            y_name, y_label = y_def[0], y_def[1]

            x_correct = _get_correct_frac(x_name)
            y_correct = _get_correct_frac(y_name)
            if x_correct is None or y_correct is None:
                continue

            matrix = np.array([
                [x_correct, 1 - x_correct],
                [1 - y_correct, y_correct],
            ])

            row_labels = [f"Gen {x_label}", f"Gen {y_label}"]
            col_labels = [f"Tagger {x_label}", f"Tagger {y_label}"]

            cm_dir = os.path.join(r.ctx.plot_dir, "mc", "cm")
            os.makedirs(cm_dir, exist_ok=True)
            outpath = os.path.join(cm_dir, f"{cm_name}.png")
            if not os.path.exists(outpath) or r.ctx.overwrite:
                fig = plot_confusion_matrix(matrix, row_labels, col_labels, lumi=r.lumi)
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

    def all(self, ratio_stacked="significance", ratio_data="data_mc"):
        """Convenience: produce all standard plot types."""
        self.stacked(ratio=ratio_stacked)
        self.shapes()
        self.trigger_overlays()
        self.hist2d()
        self.confusion_matrix()
        if self._r.data_hists_by_trig is not None:
            self.data(ratio=ratio_data)


class PlotResult:
    """Returned by load_and_run(). Holds histograms + ``.plot`` accessor."""

    def __init__(self, event_result, ctx, plot_vars, vars_2d, plot_vars_2d,
                 cm_2d=None, data_hists_by_trig=None, data_denom_hists=None,
                 sig_samples=None, dy_samples=None, tt_samples=None, qcd_samples=None):
        # From EventLoopResult
        self.mc_hists_by_trig = event_result.mc_hists_by_trig
        self.mc_denom_hists = event_result.mc_denom_hists
        self.mc_items_by_trig = event_result.mc_items_by_trig
        self.lumi = event_result.lumi
        self.gen_hists_by_trig = event_result.gen_hists_by_trig
        self.mc_hists_excl_by_trig = event_result.mc_hists_excl_by_trig
        self.cutflow_tables = event_result.cutflow_tables
        self.trig_selections = event_result.trig_selections
        self.sig_indices = event_result.sig_indices
        self.bkg_indices = event_result.bkg_indices
        self.h2d_book = event_result.h2d_book
        self.excl_sig_mHH_ptrs = event_result.excl_sig_mHH_ptrs

        # Data histograms (None if --no-data)
        self.data_hists_by_trig = data_hists_by_trig
        self.data_denom_hists = data_denom_hists

        # Context and config
        self.ctx = ctx
        self.plot_vars = plot_vars
        self.vars_2d = vars_2d
        self.plot_vars_2d = plot_vars_2d
        self.cm_2d = cm_2d or []

        # Sample lists
        self.sig_samples = sig_samples or []
        self.dy_samples = dy_samples or []
        self.tt_samples = tt_samples or []
        self.qcd_samples = qcd_samples or []

        # Plot accessor
        self.plot = Plotter(self)
