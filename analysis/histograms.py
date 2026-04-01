"""Histogram booking for MC and data."""

import ROOT

from analysis.variables import var_attrs, is_sentinel, var2d_attrs


def _sentinel_vars(plot_vars):
    """Return set of variable names that use -1 sentinel (AK8/gen-matched)."""
    return {var_attrs(v)[0] for v in plot_vars if is_sentinel(v)}


def build_groups(active_modes, decay_modes, proc_to_samples, sample_dfs):
    """Build mc_groups list from active decay modes using decayMode filter.

    Returns:
        groups: list of (dfs, label, color) — one entry per active mode that
                has at least one sample.
        built_modes: list of mode ints that were actually built.
    """
    groups = []
    built_modes = []
    for mode in active_modes:
        info = decay_modes[mode]
        proc = info["process"]
        slist = proc_to_samples.get(proc, [])
        if not slist:
            continue
        built_modes.append(mode)
        if proc == "QCD":
            dfs = [sample_dfs[s] for s in slist]
        else:
            dfs = [sample_dfs[s].Filter(f"decayMode=={mode}") for s in slist]
        groups.append((dfs, info["label"], info["color"]))
    return groups, built_modes


def book_mc_denom(mc_acc, mc_denom_groups, plot_vars, unified_ptrs):
    """Book denominator histograms (no trigger, acceptance only).

    Returns:
        denom_book: dict {var_name: [[ptrs per sample] per group]}
    """
    ak8_vars = _sentinel_vars(plot_vars)
    denom_book = {}
    for v in plot_vars:
        var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
        denom_book[var_name] = []
        for gi, (dfs, _lbl, _col) in enumerate(mc_denom_groups):
            group_ptrs = []
            for si, df in enumerate(dfs):
                uid = f"denom_{var_name}_{gi}_{si}"
                df_v = df.Filter(f"{var_name} > -0.5f") if var_name in ak8_vars else df
                ptr = df_v.Histo1D(
                    (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                    var_name, "w")
                group_ptrs.append(ptr)
                unified_ptrs.append(ptr)
            denom_book[var_name].append(group_ptrs)
    return denom_book


def book_mc_per_trigger(trig_selections, plot_vars, decay_modes, active_modes,
                        proc_to_samples, sig_samples, sig_modes, bkg_modes,
                        unified_ptrs):
    """Book per-trigger MC histograms.

    Returns:
        mc_items_by_trig: dict {trig: [{"label":, "color":}]}
        histo_books: dict {trig: {var: [[ptrs per sample] per group]}}
        sig_indices: list of int (indices into groups for signal)
        bkg_indices: list of int (indices into groups for background)
        built_modes: list of mode ints that were actually built
    """
    ak8_vars = _sentinel_vars(plot_vars)
    mc_items_by_trig = {}
    histo_books = {}

    # Determine which modes have samples
    built_modes = [m for m in active_modes
                   if proc_to_samples.get(decay_modes[m]["process"], [])]

    sig_indices = [i for i, m in enumerate(built_modes) if m in sig_modes]
    bkg_indices = [i for i, m in enumerate(built_modes) if m in bkg_modes]

    for trig_name in trig_selections:
        mc_sel = trig_selections[trig_name]["mc_sel"]
        mc_groups, _ = build_groups(active_modes, decay_modes, proc_to_samples, mc_sel)
        mc_items = [{"label": lbl, "color": col} for _, lbl, col in mc_groups]
        mc_items_by_trig[trig_name] = mc_items

        histo_book = {}
        for v in plot_vars:
            var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
            histo_book[var_name] = []
            for gi, (dfs, _lbl, _col) in enumerate(mc_groups):
                group_ptrs = []
                for si, df in enumerate(dfs):
                    uid = f"{trig_name}_{var_name}_{gi}_{si}"
                    df_v = df.Filter(f"{var_name} > -0.5f") if var_name in ak8_vars else df
                    ptr = df_v.Histo1D(
                        (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                        var_name, "w")
                    group_ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                histo_book[var_name].append(group_ptrs)
        histo_books[trig_name] = histo_book

    return mc_items_by_trig, histo_books, sig_indices, bkg_indices, built_modes


def book_2d_histograms(trig_selections, plot_vars_2d, sig_samples, unified_ptrs):
    """Book 2D histograms (signal only, gen-matched).

    Returns:
        h2d_book: dict {pname: [(trig, sample, ptr)]}
    """
    h2d_book = {}
    for v2d in plot_vars_2d:
        pname, xvar, yvar, _xl, _yl, nx, x0, x1, ny, y0, y1 = var2d_attrs(v2d)
        h2d_book[pname] = []
        for trig_name in trig_selections:
            mc_sel_t = trig_selections[trig_name]["mc_sel"]
            for s in sig_samples:
                df_v = mc_sel_t[s].Filter(f"{xvar} > -0.5f && {yvar} > -0.5f")
                uid = f"h2d_{pname}_{s}_{trig_name}"
                ptr = df_v.Histo2D(
                    ROOT.RDF.TH2DModel(uid, f";{xvar};{yvar}",
                                       nx, x0, x1, ny, y0, y1),
                    xvar, yvar, "w")
                h2d_book[pname].append((trig_name, s, ptr))
                unified_ptrs.append(ptr)
    return h2d_book


def book_gen_histograms(trig_selections, gen_plot_vars, sig_samples, sig_modes,
                        unified_ptrs):
    """Book gen-level histograms (signal only, per decay channel).

    Returns:
        gen_histo_book: dict {trig: {var: {mode: [ptrs per sig sample]}}}
    """
    gen_histo_book = {}
    for trig_name in trig_selections:
        mc_sel_t = trig_selections[trig_name]["mc_sel"]
        gen_histo_book[trig_name] = {}
        for v in gen_plot_vars:
            var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
            gen_histo_book[trig_name][var_name] = {}
            for mode in sig_modes:
                ptrs = []
                for s in sig_samples:
                    uid = f"gen_{trig_name}_{var_name}_{s}_m{mode}"
                    ptr = (mc_sel_t[s]
                           .Filter(f"decayMode=={mode}")
                           .Histo1D((uid, f";{var_name};Events", nbins, vmin, vmax),
                                    var_name, "w"))
                    ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                gen_histo_book[trig_name][var_name][mode] = ptrs
    return gen_histo_book


def book_exclusive_histograms(excl_trig_list, mc_acc, plot_vars,
                              active_modes, decay_modes, proc_to_samples,
                              unified_ptrs):
    """Book exclusive-trigger MC histograms.

    Returns:
        excl_histo_books: dict {trig: {var: [[ptrs per sample] per group]}}
    """
    ak8_vars = _sentinel_vars(plot_vars)
    excl_histo_books = {}
    for trig_name, trig in excl_trig_list:
        if trig is not None:
            mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}
        else:
            mc_sel = dict(mc_acc)
        mc_groups, _ = build_groups(active_modes, decay_modes, proc_to_samples, mc_sel)
        histo_book = {}
        for v in plot_vars:
            var_name, _xlabel, nbins, vmin, vmax = var_attrs(v)
            histo_book[var_name] = []
            for gi, (dfs, _lbl, _col) in enumerate(mc_groups):
                group_ptrs = []
                for si, df in enumerate(dfs):
                    uid = f"excl_{trig_name}_{var_name}_{gi}_{si}"
                    df_v = df.Filter(f"{var_name} > -0.5f") if var_name in ak8_vars else df
                    ptr = df_v.Histo1D(
                        (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                        var_name, "w")
                    group_ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                histo_book[var_name].append(group_ptrs)
        excl_histo_books[trig_name] = histo_book
    return excl_histo_books


def materialize_hists(book, n_groups, lumi):
    """Materialize booked histograms: combine samples per group, scale by lumi.

    Args:
        book: dict {var_name: [[ptrs per sample] per group]}
        n_groups: number of groups
        lumi: luminosity in fb^-1

    Returns:
        hists: dict {var_name: [TH1 per group]}
    """
    hists = {}
    for var_name in book:
        h_list = []
        for gi in range(n_groups):
            ptrs = book[var_name][gi]
            combined = ptrs[0].GetValue().Clone()
            for ptr in ptrs[1:]:
                combined.Add(ptr.GetValue())
            combined.Scale(lumi)
            h_list.append(combined)
        hists[var_name] = h_list
    return hists


def materialize_gen_hists(gen_histo_book, sig_modes, lumi):
    """Materialize gen-level histograms.

    Returns:
        gen_hists_by_trig: dict {trig: {var: {mode: TH1}}}
    """
    gen_hists_by_trig = {}
    for trig_name in gen_histo_book:
        h_by_var = {}
        for var_name in gen_histo_book[trig_name]:
            h_by_mode = {}
            for mode in sig_modes:
                ptrs = gen_histo_book[trig_name][var_name][mode]
                combined = ptrs[0].GetValue().Clone(f"gen_{trig_name}_{var_name}_m{mode}")
                for ptr in ptrs[1:]:
                    combined.Add(ptr.GetValue())
                combined.Scale(lumi)
                h_by_mode[mode] = combined
            h_by_var[var_name] = h_by_mode
        gen_hists_by_trig[trig_name] = h_by_var
    return gen_hists_by_trig
