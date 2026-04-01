"""Plot orchestration: MC distributions, trigger overlays, MC+data plots."""

import os
import time

import numpy as np
import matplotlib.pyplot as plt

from utils.plotting import (
    cms_label,
    plot_stacked,
    plot_shape,
    plot_stacked_with_efficiency,
    plot_stacked_with_significance,
    plot_stacked_with_cumulative_significance,
    plot_trigger_shape_overlay,
    plot_2d_hist,
    th1_to_np,
)
from analysis.variables import var_attrs, var2d_attrs, is_sentinel


# ═══════════════════════════════════════════════════════════════════════════════
#  Internal helpers
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_path(plot_dir, category, plot_type, trig_name=None, var_name=""):
    """Build plot output path and ensure directory exists."""
    if trig_name:
        d = os.path.join(plot_dir, category, plot_type, trig_name)
    else:
        d = os.path.join(plot_dir, category, plot_type)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{var_name}.png")


def _should_save(outpath, overwrite):
    """Return True if the plot should be (re)generated."""
    return not os.path.exists(outpath) or overwrite


def _save_fig(fig, xlabel_ax, xlabel, lumi, outpath, overwrite, cms_ax=None):
    """Apply CMS style, save figure, and close. Returns True if saved.

    Parameters
    ----------
    cms_ax : Axes or None
        Axes to draw the CMS label on. Defaults to xlabel_ax.
        For two-panel plots, pass the top panel so the label isn't on the ratio.
    """
    if not _should_save(outpath, overwrite):
        plt.close(fig)
        return False
    xlabel_ax.set_xlabel(xlabel)
    cms_label(cms_ax or xlabel_ax, lumi=lumi)
    fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    print(f"Saved: {outpath}")
    plt.close(fig)
    return True


def _combine_hists(h_list):
    """Clone first histogram and add the rest."""
    h = h_list[0].Clone()
    for hi in h_list[1:]:
        h.Add(hi)
    return h


# ═══════════════════════════════════════════════════════════════════════════════
#  Per-variable plot generation (shared between MC-only and MC+Data)
# ═══════════════════════════════════════════════════════════════════════════════

def _plot_var(var_name, xlabel, h_mc_list, mc_items, mc_denom_hists,
              sig_indices, bkg_indices, trig_has_expr,
              lumi, plot_dir, category, trig_name, overwrite,
              do_stacked, do_shape, do_efficiency=False, do_significance=False,
              do_cum_significance=False,
              h_data=None, data_denom_hists=None):
    """Generate plot types for one variable x one trigger."""

    if do_stacked:
        fig, ax = plot_stacked(h_mc_list, mc_items, h_data=h_data, logy=True)
        _save_fig(fig, ax, xlabel, lumi,
                  _plot_path(plot_dir, category, "stacked", trig_name, var_name), overwrite)

    if do_shape:
        fig, ax = plot_shape(h_mc_list, mc_items, h_data=h_data)
        _save_fig(fig, ax, xlabel, lumi,
                  _plot_path(plot_dir, category, "shape", trig_name, var_name), overwrite)

    if do_efficiency and trig_has_expr:
        h_denom_list = mc_denom_hists.get(var_name)
        h_d_denom = data_denom_hists.get(var_name) if data_denom_hists else None
        fig, ax, rax = plot_stacked_with_efficiency(
            h_mc_list, mc_items, h_data=h_data,
            h_mc_denom_list=h_denom_list, h_data_denom=h_d_denom, logy=True)
        _save_fig(fig, rax, xlabel, lumi,
                  _plot_path(plot_dir, category, "eff", trig_name, var_name), overwrite,
                  cms_ax=ax)

    if do_significance:
        fig, ax, rax = plot_stacked_with_significance(
            h_mc_list, mc_items, h_data=h_data,
            sig_indices=sig_indices, bkg_indices=bkg_indices, logy=True)
        _save_fig(fig, rax, xlabel, lumi,
                  _plot_path(plot_dir, category, "sig", trig_name, var_name), overwrite,
                  cms_ax=ax)

    if do_cum_significance:
        fig, ax, rax = plot_stacked_with_cumulative_significance(
            h_mc_list, mc_items, h_data=h_data,
            sig_indices=sig_indices, bkg_indices=bkg_indices, logy=True)
        _save_fig(fig, rax, xlabel, lumi,
                  _plot_path(plot_dir, category, "cum_sig", trig_name, var_name), overwrite,
                  cms_ax=ax)


# ═══════════════════════════════════════════════════════════════════════════════
#  Trigger overlay helper
# ═══════════════════════════════════════════════════════════════════════════════

def _save_trigger_overlay(h_by_trig, xlabel, lumi, outpath, overwrite, title=None):
    """Plot trigger shape overlay and save."""
    if not h_by_trig:
        return
    fig, ax = plot_trigger_shape_overlay(h_by_trig, title=title)
    _save_fig(fig, ax, xlabel, lumi, outpath, overwrite)


# ═══════════════════════════════════════════════════════════════════════════════
#  Public API
# ═══════════════════════════════════════════════════════════════════════════════

def plot_mHH_comparison(sel, h_by_var, sig_indices, outpath, lumi, overwrite,
                        h_data_m4j=None):
    """Plot gen mHH vs reco m4j (and optionally data m4j)."""
    if not _should_save(outpath, overwrite):
        return
    if not sel["sig_mHH_ptrs"] or "m4j" not in h_by_var or not sig_indices:
        return

    fig, ax = plt.subplots(figsize=(8, 6))

    def _norm_step(h, color, label, **kw):
        edges, vals, _ = th1_to_np(h)
        area = float(np.sum(vals * np.diff(edges)))
        if area > 0:
            vals = vals / area
        ax.step(edges, np.r_[vals, vals[-1]], where="post", linewidth=2,
                color=color, label=label, **kw)

    h_gen = _combine_hists([p.GetPtr() for p in sel["sig_mHH_ptrs"]])
    _norm_step(h_gen, "tab:blue", r"$m_{HH}$ gen-level (signal)")

    h_sig_m4j = _combine_hists([h_by_var["m4j"][si] for si in sig_indices])
    _norm_step(h_sig_m4j, "tab:red", r"$m_{4j}$ reco (signal)")

    if h_data_m4j is not None:
        _norm_step(h_data_m4j, "black", r"$m_{4j}$ reco (data)", linestyle="--")

    ax.set_ylabel("Normalised to unity")
    ax.legend()
    ax.grid(True, alpha=0.3)
    _save_fig(fig, ax, r"Mass [GeV]", lumi, outpath, overwrite)


def plot_mc_distributions(trig_selections, mc_hists_by_trig, mc_denom_hists,
                          mc_items_by_trig, gen_hists_by_trig,
                          plot_vars, gen_plot_vars, sig_indices, bkg_indices,
                          sig_modes, decay_modes,
                          lumi, plot_dir, overwrite, do_stacked, do_shape,
                          do_efficiency=False, do_significance=False,
                          do_cum_significance=False,
                          n_workers=1):
    """MC-only plots: stacked, shape, and optionally efficiency/significance.

    When n_workers > 1, renders plots in parallel using ProcessPoolExecutor.
    """
    from utils.plotting import prepare_mc_task, render_task

    # Phase 1: collect tasks (sequential, fast — converts ROOT → numpy)
    tasks = []
    for trig_name in trig_selections:
        sel = trig_selections[trig_name]
        h_by_var = mc_hists_by_trig[trig_name]
        mc_items = mc_items_by_trig[trig_name]

        # mHH comparison (not parallelized — uses custom logic)
        plot_mHH_comparison(sel, h_by_var, sig_indices,
                            _plot_path(plot_dir, "mc", "shape", trig_name, "mHH_m4j"),
                            lumi, overwrite)

        for iv, v in enumerate(plot_vars):
            var_name, xlabel, *_ = var_attrs(v)
            h_mc_list = h_by_var[var_name]
            edges, mc_components, _ = prepare_mc_task(h_mc_list, mc_items,
                                                      warn=(iv == 0))

            base = {"edges": edges, "mc_components": mc_components,
                    "xlabel": xlabel, "lumi": lumi, "overwrite": overwrite,
                    "sig_indices": sig_indices, "bkg_indices": bkg_indices}

            if do_stacked:
                tasks.append({**base, "type": "stacked",
                              "outpath": _plot_path(plot_dir, "mc", "stacked", trig_name, var_name)})
            if do_shape:
                tasks.append({**base, "type": "shape",
                              "outpath": _plot_path(plot_dir, "mc", "shape", trig_name, var_name)})
            if do_significance:
                tasks.append({**base, "type": "sig",
                              "outpath": _plot_path(plot_dir, "mc", "sig", trig_name, var_name)})
            if do_cum_significance:
                tasks.append({**base, "type": "cum_sig",
                              "outpath": _plot_path(plot_dir, "mc", "cum_sig", trig_name, var_name)})
            if do_efficiency and sel["trig"] is not None:
                # Efficiency still uses sequential path (needs ROOT TEfficiency)
                _plot_var(var_name, xlabel, h_mc_list, mc_items, mc_denom_hists,
                          sig_indices, bkg_indices, True,
                          lumi, plot_dir, "mc", trig_name, overwrite,
                          do_stacked=False, do_shape=False, do_efficiency=True)

    # Phase 2: render in parallel
    if tasks:
        if n_workers > 1:
            import multiprocessing
            from concurrent.futures import ProcessPoolExecutor
            ctx = multiprocessing.get_context("spawn")
            print(f"[plots] Rendering {len(tasks)} plots with {n_workers} workers...")
            with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as pool:
                results = list(pool.map(render_task, tasks))
            print(f"[plots] Done — {len(results)} plots saved")
        else:
            for task in tasks:
                render_task(task)
                print(f"Saved: {task['outpath']}")

        # Gen-level shape plots per decay channel
        if trig_name in gen_hists_by_trig:
            for v in gen_plot_vars:
                var_name, xlabel, nbins, vmin, vmax = var_attrs(v)
                gen_dir = os.path.join(plot_dir, "mc", "shape", trig_name, "gen", "HBBTT_SM")
                os.makedirs(gen_dir, exist_ok=True)
                outpath = os.path.join(gen_dir, f"{var_name}.png")
                if not _should_save(outpath, overwrite):
                    continue
                h_by_mode = gen_hists_by_trig[trig_name][var_name]
                fig, ax = plt.subplots(figsize=(8, 6))
                has_entries = False
                for mode in sig_modes:
                    h = h_by_mode.get(mode)
                    if h is None or h.GetEntries() == 0:
                        continue
                    has_entries = True
                    edges, vals, errs = th1_to_np(h)
                    area = float(np.sum(vals * np.diff(edges)))
                    if area > 0:
                        vals, errs = vals / area, errs / area
                    color = decay_modes[mode]["color"]
                    label = decay_modes[mode]["label"]
                    ax.step(edges, np.r_[vals, vals[-1]], where="post",
                            linewidth=2, color=color, label=label)
                    ax.fill_between(edges,
                                    np.r_[vals - errs, (vals - errs)[-1]],
                                    np.r_[vals + errs, (vals + errs)[-1]],
                                    step="post", alpha=0.2, color=color)
                if not has_entries:
                    plt.close(fig)
                    continue
                ax.set_ylabel("Normalised to unity")
                ax.set_xlim(vmin, vmax)
                ax.legend()
                ax.grid(True, alpha=0.3)
                _save_fig(fig, ax, xlabel, lumi, outpath, overwrite=True)


def plot_trigger_overlays(trig_selections, mc_hists_by_trig, mc_hists_excl_by_trig,
                          plot_vars, sig_indices,
                          lumi, plot_dir, overwrite, do_shape):
    """Trigger shape overlays: inclusive and exclusive, per variable and per signal channel."""
    if not do_shape:
        return

    sig_channels = []
    if len(sig_indices) >= 3:
        sig_channels = [
            (sig_indices[0], "hh", r"$HH\to bb\tau_h\tau_h$"),
            (sig_indices[1], "hm", r"$HH\to bb\tau_\mu\tau_h$"),
            (sig_indices[2], "he", r"$HH\to bb\tau_e\tau_h$"),
        ]

    for v in plot_vars:
        var_name, xlabel, *_ = var_attrs(v)

        # All-MC trigger overlay
        h_total_by_trig = {tn: _combine_hists(mc_hists_by_trig[tn][var_name])
                           for tn in trig_selections}
        _save_trigger_overlay(h_total_by_trig, xlabel, lumi,
                              _plot_path(plot_dir, "mc", "overlay", var_name=var_name), overwrite)

        # Per-signal-channel trigger overlay
        for gi, ch_key, ch_label in sig_channels:
            h_by_trig = {tn: mc_hists_by_trig[tn][var_name][gi] for tn in trig_selections}
            _save_trigger_overlay(h_by_trig, xlabel, lumi,
                                  _plot_path(plot_dir, "mc", "overlay", var_name=f"{ch_key}_{var_name}"),
                                  overwrite, title=ch_label)

    print("Trigger overlay plots done")

    # Exclusive-trigger overlays
    if mc_hists_excl_by_trig:
        for v in plot_vars:
            var_name, xlabel, *_ = var_attrs(v)

            h_by_trig = {}
            for tn, h_by_var in mc_hists_excl_by_trig.items():
                h_list = h_by_var.get(var_name, [])
                if h_list:
                    h_by_trig[tn] = _combine_hists(h_list)
            _save_trigger_overlay(h_by_trig, xlabel, lumi,
                                  _plot_path(plot_dir, "mc", "overlay_exclusive", var_name=var_name),
                                  overwrite)

            for gi, ch_key, ch_label in sig_channels:
                h_by_trig = {}
                for tn, h_by_var in mc_hists_excl_by_trig.items():
                    h_list = h_by_var.get(var_name, [])
                    if gi < len(h_list):
                        h_by_trig[tn] = h_list[gi]
                _save_trigger_overlay(h_by_trig, xlabel, lumi,
                                      _plot_path(plot_dir, "mc", "overlay_exclusive",
                                                 var_name=f"{ch_key}_{var_name}"),
                                      overwrite, title=ch_label)

        print("Exclusive-trigger overlay plots done")


def plot_gen_overlays(trig_selections, gen_hists_by_trig,
                      gen_plot_vars, sig_modes, decay_modes,
                      h2d_book, plot_vars_2d,
                      lumi, plot_dir, overwrite):
    """Gen-level trigger overlays (combined + per-channel) and 2D gen-matched plots."""
    for v in gen_plot_vars:
        var_name, xlabel, *_ = var_attrs(v)
        if not gen_hists_by_trig or var_name not in next(iter(gen_hists_by_trig.values()), {}):
            continue

        # Combined overlay
        h_total_by_trig = {}
        for tn in trig_selections:
            vals = list(gen_hists_by_trig.get(tn, {}).get(var_name, {}).values())
            if vals:
                h_total_by_trig[tn] = _combine_hists(vals)
        _save_trigger_overlay(h_total_by_trig, xlabel, lumi,
                              _plot_path(plot_dir, "mc", "overlay", trig_name="gen", var_name=var_name),
                              overwrite)

        # Per-channel overlay
        for mode in sig_modes:
            ch_key = {20: "hh", 21: "hm", 22: "he"}.get(mode, str(mode))
            h_by_trig = {}
            for tn in trig_selections:
                h = gen_hists_by_trig.get(tn, {}).get(var_name, {}).get(mode)
                if h and h.GetEntries() > 0:
                    h_by_trig[tn] = h
            _save_trigger_overlay(h_by_trig, xlabel, lumi,
                                  _plot_path(plot_dir, "mc", "overlay", trig_name="gen",
                                             var_name=f"{ch_key}_{var_name}"),
                                  overwrite, title=decay_modes[mode]["label"])

    print("Gen-level overlay plots done")

    # 2D gen-matched plots
    if h2d_book and plot_vars_2d:
        for v2d in plot_vars_2d:
            pname, _xvar, _yvar, xlabel, ylabel, *_ = var2d_attrs(v2d)
            trigs_seen = {}
            for trig_name, s, ptr in h2d_book.get(pname, []):
                h2 = ptr.Clone()
                h2.Scale(lumi)
                if trig_name not in trigs_seen:
                    trigs_seen[trig_name] = h2
                else:
                    trigs_seen[trig_name].Add(h2)
            for trig_name, h2_combined in trigs_seen.items():
                outpath = _plot_path(plot_dir, "mc", "gen2d", trig_name, var_name=pname)
                if _should_save(outpath, overwrite):
                    fig = plot_2d_hist(h2_combined, xlabel=xlabel, ylabel=ylabel, lumi=lumi)
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)
        print("2D gen-matched plots done")


def plot_data_mc(trig_selections, mc_hists_by_trig, mc_denom_hists,
                 data_hists_by_trig, data_denom_hists,
                 mc_items_by_trig, plot_vars, sig_indices, bkg_indices,
                 lumi, plot_dir, overwrite, do_stacked, do_shape,
                 do_efficiency=False, do_significance=False,
                 do_cum_significance=False):
    """MC+Data overlay plots."""
    for trig_name in trig_selections:
        h_by_var = mc_hists_by_trig[trig_name]
        h_data_var = data_hists_by_trig[trig_name]
        mc_items = mc_items_by_trig[trig_name]

        for v in plot_vars:
            var_name, xlabel, *_ = var_attrs(v)
            h_data = h_data_var.get(var_name)
            if h_data is None:
                continue
            _plot_var(var_name, xlabel, h_by_var[var_name], mc_items, mc_denom_hists,
                      sig_indices, bkg_indices, trig_selections[trig_name]["trig"] is not None,
                      lumi, plot_dir, "data", trig_name, overwrite, do_stacked, do_shape,
                      do_efficiency=do_efficiency, do_significance=do_significance,
                      do_cum_significance=do_cum_significance,
                      h_data=h_data, data_denom_hists=data_denom_hists)

        # mHH comparison with data
        sel = trig_selections[trig_name]
        plot_mHH_comparison(sel, h_by_var, sig_indices,
                            _plot_path(plot_dir, "data", "shape", trig_name, "mHH_m4j"),
                            lumi, overwrite, h_data_m4j=h_data_var.get("m4j"))

        print(f"[{trig_name}] MC+Data plots done")


def book_data_histograms(data_acc, trig_selections, plot_vars):
    """Book and run data histograms.

    Returns
    -------
    data_hists_by_trig : dict {trig_name: {var_name: TH1}}
    data_denom_hists : dict {var_name: TH1}
    """
    import ROOT

    data_hists_by_trig = {}
    data_denom_hists = {}

    data_denom_ptrs = []
    data_denom_book = {}
    for v in plot_vars:
        var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
        if var_name == "gen_mHH":
            continue
        uid_d = f"data_denom_{var_name}"
        df_v = data_acc.Filter(f"{var_name} > -0.5f") if is_sentinel(v) else data_acc
        d_ptr = df_v.Histo1D(
            (f"h_{uid_d}", f";{var_name};Events", nbins, vmin, vmax), var_name)
        data_denom_book[var_name] = d_ptr
        data_denom_ptrs.append(d_ptr)

    for trig_name in trig_selections:
        data_sel = trig_selections[trig_name]["data_sel"]
        data_histo_ptrs = []
        data_histo_book = {}
        for v in plot_vars:
            var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
            if var_name == "gen_mHH":
                continue
            uid_d = f"{trig_name}_{var_name}_data"
            df_v = data_sel.Filter(f"{var_name} > -0.5f") if is_sentinel(v) else data_sel
            d_ptr = df_v.Histo1D(
                (f"h_{uid_d}", f";{var_name};Events", nbins, vmin, vmax), var_name)
            data_histo_book[var_name] = d_ptr
            data_histo_ptrs.append(d_ptr)

        if not data_denom_hists:
            data_histo_ptrs.extend(data_denom_ptrs)

        print(f"\n[{trig_name}] Booked {len(data_histo_ptrs)} data histograms")
        print(f"[{trig_name}] Running data event loop (XCache) ...", flush=True)
        t0 = time.time()
        ROOT.RDF.RunGraphs(data_histo_ptrs)
        print(f"[{trig_name}] Data event loop done in {time.time() - t0:.1f}s")

        if not data_denom_hists:
            for var_name, d_ptr in data_denom_book.items():
                data_denom_hists[var_name] = d_ptr.GetValue()

        data_hists_by_trig[trig_name] = {
            var_name: d_ptr.GetValue() for var_name, d_ptr in data_histo_book.items()
        }

    return data_hists_by_trig, data_denom_hists
