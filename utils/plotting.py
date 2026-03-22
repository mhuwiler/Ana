"""
Plotting utilities for scouting analysis.

Contains:
  - ROOT TH1 / TEfficiency -> numpy converters
  - Trigger efficiency overlay plot
  - Yield overlay plot
  - Stacked MC + data plot
  - CMS style setup via mplhep
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D

import ROOT


# ---------------------------------------------------------------------------
# Style setup
# ---------------------------------------------------------------------------

_DARK_MODE = False


def setup_style(dark=False):
    """Apply CMS style via mplhep and set standard rcParams.

    Parameters:
        dark : bool – if True, use a dark background theme.
    """
    global _DARK_MODE
    _DARK_MODE = dark

    plt.style.use(hep.style.CMS)
    common = {
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.labelsize": 13,
        "axes.titlesize": 13,
        "xtick.labelsize": 13,
        "ytick.labelsize": 13,
        "legend.fontsize": 13,
        "legend.title_fontsize": 13,
        "axes.linewidth": 1.2,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 6,
        "ytick.major.size": 6,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
    }
    mpl.rcParams.update(common)

    if dark:
        mpl.rcParams.update({
            "figure.facecolor": "#1e1e1e",
            "axes.facecolor": "#2b2b2b",
            "savefig.facecolor": "#1e1e1e",
            "axes.edgecolor": "#cccccc",
            "axes.labelcolor": "#cccccc",
            "xtick.color": "#cccccc",
            "ytick.color": "#cccccc",
            "text.color": "#cccccc",
            "legend.facecolor": "#2b2b2b",
            "legend.edgecolor": "#555555",
            "grid.color": "#444444",
        })


def cms_label(ax, lumi=None, label="Work in Progress"):
    """Add standard CMS label to an axes.

    Draws 'CMS' (bold) + label (italic) top-left, and lumi + energy top-right.

    Parameters:
        ax   : matplotlib Axes
        lumi : float or None – integrated luminosity in fb^-1 (shown top-right)
        label: str – text after 'CMS' (e.g. 'Work in Progress', 'Preliminary')
    """
    hep.cms.label(label, ax=ax, data=True,
                  lumi=f"{lumi:.4g}" if lumi is not None else None,
                  year="2024", com=13.6,
                  fontsize=16)


# ---------------------------------------------------------------------------
# ROOT -> numpy converters
# ---------------------------------------------------------------------------

def th1_to_np(h):
    """Convert a ROOT TH1 to (edges, values, errors) numpy arrays.

    Underflow is folded into the first visible bin, overflow into the last.
    Errors are added in quadrature.
    """
    nb = h.GetNbinsX()
    edges = np.array(
        [h.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)]
        + [h.GetXaxis().GetBinUpEdge(nb)],
        dtype=float,
    )
    vals = np.array([h.GetBinContent(1 + i) for i in range(nb)], dtype=float)
    errs = np.array([h.GetBinError(1 + i) for i in range(nb)], dtype=float)

    # Fold underflow (bin 0) into first visible bin
    vals[0]  += h.GetBinContent(0)
    errs[0]   = np.sqrt(errs[0]**2 + h.GetBinError(0)**2)

    # Fold overflow (bin nb+1) into last visible bin
    vals[-1] += h.GetBinContent(nb + 1)
    errs[-1]  = np.sqrt(errs[-1]**2 + h.GetBinError(nb + 1)**2)

    return edges, vals, errs


def teff_to_np(teff):
    """Convert a ROOT TEfficiency to (x, y, xerr, yerr_lo, yerr_hi) numpy arrays."""
    htot = teff.GetTotalHistogram()
    nb = htot.GetNbinsX()
    x, y, yerr_lo, yerr_hi, xerr = [], [], [], [], []
    for i in range(1, nb + 1):
        tot = htot.GetBinContent(i)
        if tot <= 0:
            continue
        xc = htot.GetXaxis().GetBinCenter(i)
        hw = 0.5 * htot.GetXaxis().GetBinWidth(i)
        eff = float(np.clip(teff.GetEfficiency(i), 0.0, 1.0))
        elo = float(np.clip(teff.GetEfficiencyErrorLow(i), 0.0, 1.0))
        ehi = float(np.clip(teff.GetEfficiencyErrorUp(i), 0.0, 1.0))

        x.append(xc)
        xerr.append(hw)
        y.append(eff)
        yerr_lo.append(elo)
        yerr_hi.append(ehi)

    return (np.array(x), np.array(y), np.array(xerr),
            np.array(yerr_lo), np.array(yerr_hi))


# ---------------------------------------------------------------------------
# Trigger efficiency overlay (multi-channel, multi-trigger)
# ---------------------------------------------------------------------------

def plot_trigger_eff_overlay_channels(
    channels, triggers, *,
    var="gen_mHH", weight="w", pre_expr=None,
    nbins=20, xmin=200.0, xmax=1800.0,
    cms_year="2024", cms_com="13.6", cms_label_text="",
    figsize=(10, 10),
    pre_label="Preselection", trigger_colors=None,
    pre_key="pre",
    pre_color=None,
    xlabel=None,
):
    """
    Overlay trigger efficiency for multiple decay channels and triggers.

    channels: list of dicts like
      [
        {"key":"hh", "label":"HH->bbthth",  "df":sig_df_hh, "ls":"-",  "marker":"o"},
        {"key":"hm", "label":"HH->bbteth",  "df":sig_df_hm, "ls":"--", "marker":"s"},
        {"key":"he", "label":"HH->bbtmth",  "df":sig_df_he, "ls":":",  "marker":"^"},
      ]

    triggers: list of dicts like
      [
        {"key":"all", "expr": DST_ALL_expr, "label":"Data Scouting", "lw":2.5},
        {"key":"nom", "expr": NOM_expr,     "label":"Nominal",       "lw":2.5},
        {"key":"park","expr": PARK_expr,    "label":"Data Parking",  "lw":2.5},
      ]
    """
    if trigger_colors is None:
        trigger_colors = {}

    color_map = {
        pre_key: trigger_colors.get(pre_key, "tab:blue"),
    }
    for t in triggers:
        k = t["key"]
        color_map[k] = trigger_colors.get(k, {
            "nom": "tab:orange",
            "park": "tab:green",
            "all": "tab:red",
        }.get(k, "black"))

    # Book all histogram actions
    h_pre_ptrs = []
    h_num_ptrs = {t["key"]: [] for t in triggers}

    for i, ch in enumerate(channels):
        df = ch["df"]
        df_pre = df if (pre_expr is None or str(pre_expr).strip() in ["", "1"]) else df.Filter(pre_expr)

        h_pre_ptrs.append(
            df_pre.Histo1D(
                (f"h_pre_{i}", f";{var} [GeV];Events", nbins, xmin, xmax),
                var, weight,
            )
        )
        for t in triggers:
            k = t["key"]
            h_num_ptrs[k].append(
                df_pre.Filter(t["expr"]).Histo1D(
                    (f"h_{k}_{i}", f";{var} [GeV];Events", nbins, xmin, xmax),
                    var, weight,
                )
            )

    # Run once
    all_ptrs = list(h_pre_ptrs)
    for k in h_num_ptrs:
        all_ptrs.extend(h_num_ptrs[k])
    try:
        ROOT.RDF.RunGraphs(all_ptrs)
    except Exception:
        pass

    # Materialize + TEfficiency + numpy per channel
    per_ch = []
    for i, ch in enumerate(channels):
        Hpre = h_pre_ptrs[i].GetValue()
        Hnums = {t["key"]: h_num_ptrs[t["key"]][i].GetValue() for t in triggers}

        effs = {}
        for t in triggers:
            k = t["key"]
            te = ROOT.TEfficiency(Hnums[k], Hpre)
            te.SetUseWeightedEvents(True)
            effs[k] = te

        edges, pre_vals, pre_errs = th1_to_np(Hpre)
        num_np = {k: th1_to_np(Hnums[k])[1:] for k in Hnums}
        eff_np = {k: teff_to_np(effs[k]) for k in effs}

        per_ch.append({
            "key": ch.get("key", f"ch{i}"),
            "label": ch.get("label", f"ch{i}"),
            "ls": ch.get("ls", "-"),
            "marker": ch.get("marker", "o"),
            "edges": edges,
            "pre": (pre_vals, pre_errs),
            "num": num_np,
            "eff_np": eff_np,
            "Hpre": Hpre,
            "Hnums": Hnums,
            "effs": effs,
        })

    # Plot: one figure, overlay everything
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    # TOP: preselection + trigger shapes
    for ch in per_ch:
        edges = ch["edges"]
        pre_vals, _ = ch["pre"]
        ax.step(
            edges, np.r_[pre_vals, pre_vals[-1]],
            where="post", linewidth=2.5,
            linestyle=ch["ls"], color=color_map[pre_key],
        )
        for t in triggers:
            k = t["key"]
            vals, _errs = ch["num"][k]
            ax.step(
                edges, np.r_[vals, vals[-1]],
                where="post", linewidth=float(t.get("lw", 2.5)),
                linestyle=ch["ls"], color=color_map[k],
            )

    ax.set_ylabel("Yields")
    ax.tick_params(labelbottom=False)
    hep.cms.label(cms_label_text, data=False, ax=ax, year=str(cms_year), com=str(cms_com))

    # BOTTOM: efficiencies
    for ch in per_ch:
        for t in triggers:
            k = t["key"]
            x, y, xerr_arr, ylo, yhi = ch["eff_np"][k]
            rax.errorbar(
                x, y,
                xerr=xerr_arr,
                yerr=np.vstack([ylo, yhi]),
                fmt=ch["marker"],
                markersize=float(t.get("markersize", 5.0)),
                capsize=2, linewidth=1.2, linestyle="none",
                color=color_map[k],
            )

    rax.set_xlabel(xlabel if xlabel is not None else rf"{var} [GeV]")
    rax.set_ylabel("Trigger Efficiency")
    rax.set_ylim(0.0, 1.1)
    rax.grid(True, axis="y", alpha=0.3)

    # Legends
    channel_handles = [
        Line2D([0], [0], color="black", linestyle=ch.get("ls", "-"),
               marker=ch.get("marker", "o"), markersize=6, linewidth=2,
               label=ch.get("label", ch.get("key", "ch")))
        for ch in channels
    ]
    trigger_handles = [
        Line2D([0], [0], color=color_map[pre_key], linestyle="-",
               linewidth=2.5, label=pre_label)
    ]
    trigger_handles += [
        Line2D([0], [0], color=color_map[t["key"]], linestyle="-",
               linewidth=float(t.get("lw", 2.5)), label=t.get("label", t["key"]))
        for t in triggers
    ]

    leg1 = ax.legend(handles=channel_handles, title="Decay channels",
                     loc="upper left", frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=trigger_handles, loc="upper right", frameon=True)

    out = {"per_channel": per_ch, "color_map": color_map}
    return fig, ax, rax, out


# ---------------------------------------------------------------------------
# Yield overlay (multiple samples, single plot)
# ---------------------------------------------------------------------------

def plot_yield_overlay(
    dfs, labels, *,
    var=None,
    vars=None,
    weight="w",
    pre_expr=None,
    nbins=20, xmin=200.0, xmax=1800.0,
    cms_year="2024", cms_com="13.6", cms_label_text="Simulation",
    figsize=(10, 7),
    xlabel=None,
    ylabel="Weighted yields",
    density=False,
    logy=False,
):
    if len(dfs) != len(labels):
        raise ValueError("dfs and labels must have the same length")
    if (var is None) == (vars is None):
        raise ValueError("Provide exactly one of: var='...' (common) OR vars=[...] (per-DF).")
    if vars is not None and len(vars) != len(dfs):
        raise ValueError("vars must have the same length as dfs")

    vlist = [var] * len(dfs) if vars is None else list(vars)

    h_ptrs = []
    for i, df in enumerate(dfs):
        df_sel = df if (pre_expr is None or str(pre_expr).strip() in ["", "1"]) else df.Filter(pre_expr)
        vname = vlist[i]
        cols = set(str(c) for c in df_sel.GetColumnNames())
        if vname not in cols:
            sample = sorted(list(cols))[:50]
            raise RuntimeError(
                f'Unknown column "{vname}" for dataset "{labels[i]}". '
                f"First ~50 columns: {sample}"
            )
        h_ptrs.append(
            df_sel.Histo1D(
                (f"h_{i}", f";{vname};Yields", int(nbins), float(xmin), float(xmax)),
                vname, weight,
            )
        )

    try:
        ROOT.RDF.RunGraphs(h_ptrs)
    except Exception:
        pass

    fig, ax = plt.subplots(figsize=figsize)
    line_objs = []
    for i, hptr in enumerate(h_ptrs):
        H = hptr.GetValue()
        edges, vals, errs = th1_to_np(H)
        if density:
            bw = np.diff(edges)
            area = float(np.sum(vals * bw))
            if area > 0:
                vals = vals / area
                errs = errs / area
        (ln,) = ax.step(edges, np.r_[vals, vals[-1]], where="post", linewidth=2.5)
        line_objs.append(ln)

    ax.set_xlabel(xlabel if xlabel is not None else rf"{(var if var is not None else 'x')} [arb]")
    ax.set_ylabel("Density" if density else ylabel)
    if logy:
        ax.set_yscale("log")

    hep.cms.label(cms_label_text, data=False, ax=ax, year=str(cms_year), com=str(cms_com))
    ax.legend(line_objs, labels, loc="best", frameon=True)
    return fig, ax


# ---------------------------------------------------------------------------
# Stacked MC + data overlay  (helpers + main function)
# ---------------------------------------------------------------------------

def _get_dfs(item):
    """Get list of DataFrames from an mc_item (supports 'df' or 'dfs')."""
    if "dfs" in item:
        return list(item["dfs"])
    return [item["df"]]


def _book_histograms(data_df, mc_items, var, weight, nbins, xmin, xmax):
    """Book RDataFrame histogram actions for data and MC, run the event loop once.

    Each mc_item may contain either a single "df" or a list "dfs".
    When "dfs" is used, per-sub-df histograms are booked and then summed
    so that the caller still gets one TH1 per mc_item.
    data_df may be None (MC-only mode).
    """
    uid = ROOT.TUUID().AsString().replace("-", "_")

    # Data (optional)
    h_data_ptr = None
    if data_df is not None:
        h_data_ptr = data_df.Histo1D(
            (f"h_data_{uid}", f";{var};Events", nbins, xmin, xmax), var,
        )

    # MC — book one histogram per sub-df
    all_mc_ptrs = []       # flat list of RResultPtrs
    mc_item_slices = []    # (start, count) per mc_item
    idx = 0
    for i, item in enumerate(mc_items):
        dfs = _get_dfs(item)
        mc_item_slices.append((idx, len(dfs)))
        for j, df in enumerate(dfs):
            ptr = df.Histo1D(
                (f"h_mc_{i}_{j}_{uid}", f";{var};Events", nbins, xmin, xmax),
                var, weight,
            )
            all_mc_ptrs.append(ptr)
            idx += 1

    # Run all graphs in one go
    graphs = list(all_mc_ptrs)
    if h_data_ptr is not None:
        graphs.append(h_data_ptr)
    ROOT.RDF.RunGraphs(graphs)

    # Combine sub-df histograms within each mc_item
    h_mc_combined = []
    for start, count in mc_item_slices:
        combined = all_mc_ptrs[start].GetValue().Clone()
        for k in range(1, count):
            combined.Add(all_mc_ptrs[start + k].GetValue())
        h_mc_combined.append(combined)

    return h_data_ptr, h_mc_combined


def _build_mc_components(h_mc_list, mc_items, sort_by_yield=True):
    """Convert MC histograms (TH1 objects) to numpy arrays (vals, errs), optionally sorted."""
    components = []
    for h, item in zip(h_mc_list, mc_items):
        _, vals_raw, errs_raw = th1_to_np(h)
        if np.any(vals_raw < 0):
            neg_count = np.sum(vals_raw < 0)
            print(
                f"  -> WARNING: Plotting '{item['label']}' - found {neg_count} "
                f"bin(s) with negative yields. Clipping to 0 for the stack."
            )
        vals_clipped = np.maximum(vals_raw, 0.0)
        components.append({
            "item": item,
            "vals": vals_clipped,
            "errs": errs_raw,
            "yield": float(np.sum(vals_clipped)),
        })

    if sort_by_yield:
        components.sort(key=lambda x: x["yield"])

    return components


def _draw_mc_stack(ax, centers, widths, mc_components):
    """Draw the stacked MC bars and return the total stack height per bin."""
    bottom = np.zeros_like(centers, dtype=float)
    edge_col = "#cccccc" if _DARK_MODE else "black"
    for comp in mc_components:
        item = comp["item"]
        total_yield = comp["yield"]
        label = f"{item['label']} ({total_yield:.2g})"
        ax.bar(
            centers, comp["vals"], width=widths, bottom=bottom, align="center",
            label=label, color=item["color"], alpha=1.0,
            edgecolor=edge_col, linewidth=0.2,
        )
        bottom += comp["vals"]
    return bottom


def _draw_mc_stat_unc(ax, edges, stack_total, mc_components):
    """Draw hatched MC statistical uncertainty band: sigma = sqrt(sum w^2)."""
    mc_err2 = np.zeros(len(mc_components[0]["errs"]), dtype=float)
    for comp in mc_components:
        mc_err2 += comp["errs"] ** 2
    mc_err_total = np.sqrt(mc_err2)

    band_lo = np.r_[stack_total - mc_err_total, (stack_total - mc_err_total)[-1]]
    band_hi = np.r_[stack_total + mc_err_total, (stack_total + mc_err_total)[-1]]
    hatch_col = "#aaaaaa" if _DARK_MODE else "gray"
    ax.fill_between(
        edges, band_lo, band_hi,
        step="post", facecolor="none", edgecolor=hatch_col,
        hatch="////", linewidth=0, label="MC stat. unc.", zorder=5,
    )
    return mc_err_total


def _draw_data(ax, centers, data_vals, data_errs, label="Data"):
    """Draw data points with Poisson error bars."""
    total_yield = float(np.sum(data_vals))
    label_with_yield = f"{label} ({total_yield:.2g})"
    data_col = "white" if _DARK_MODE else "black"
    ax.errorbar(
        centers, data_vals, yerr=data_errs, fmt="o", color=data_col,
        label=label_with_yield, ms=4, capsize=2, linewidth=1, zorder=10,
    )


def plot_stacked_all_mc(
    data_df=None, mc_items=None, *,
    var="ak4_pt0", weight="w",
    nbins=20, xmin=0, xmax=1800,
    logy=True,
    sort_mc_by_yield=True,
    title=None,
):
    """
    Stacked MC histogram with optional data overlay and MC stat. uncertainty band.

    Parameters:
        data_df  : RDataFrame for data (unweighted), or None for MC-only
        mc_items : list of dicts with keys ("df" | "dfs"), "label", "color"
                   Use "df" for a single RDataFrame, or "dfs" for a list of
                   RDataFrames whose histograms are summed (per-sample normalisation).
        var      : branch name to plot
        weight   : weight column name (MC only)
        nbins, xmin, xmax : histogram binning
        logy     : logarithmic y-axis
        sort_mc_by_yield : stack rare processes at the bottom
    """
    # 1. Book & run histograms
    h_data_ptr, h_mc_list = _book_histograms(
        data_df, mc_items, var, weight, nbins, xmin, xmax,
    )

    # 2. Extract numpy arrays
    # Use MC histogram to get bin edges (works even without data)
    edges_ref = h_mc_list[0] if h_mc_list else None
    if h_data_ptr is not None:
        edges, data_vals, data_errs = th1_to_np(h_data_ptr.GetValue())
    else:
        edges, _, _ = th1_to_np(edges_ref)
        data_vals = None
        data_errs = None

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    mc_components = _build_mc_components(h_mc_list, mc_items, sort_by_yield=sort_mc_by_yield)

    # 3. Draw
    fig, ax = plt.subplots(figsize=(8, 6))

    stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
    _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
    if data_vals is not None:
        _draw_data(ax, centers, data_vals, data_errs)

    # 4. Cosmetics
    ax.set_xlabel(var)
    ax.set_ylabel("Events")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2)

    if logy:
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        if data_vals is not None:
            ymax = max(ymax, float(np.max(data_vals)))
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))

    return fig, ax


# ---------------------------------------------------------------------------
# Shape overlay (normalised to unity)
# ---------------------------------------------------------------------------

def plot_shape_overlay(
    data_df=None, mc_items=None, *,
    var="ak4_pt0", weight="w",
    nbins=20, xmin=0, xmax=1800,
    logy=False,
    title=None,
):
    """
    Overlay normalised-to-unity shapes for each MC group (and optionally data).

    Same interface as plot_stacked_all_mc but draws step histograms
    instead of a stack, each normalised so that its integral equals 1.
    """
    h_data_ptr, h_mc_list = _book_histograms(
        data_df, mc_items, var, weight, nbins, xmin, xmax,
    )

    edges = th1_to_np(h_mc_list[0])[0]

    fig, ax = plt.subplots(figsize=(8, 6))

    for h, item in zip(h_mc_list, mc_items):
        _, vals, _ = th1_to_np(h)
        area = float(np.sum(vals * np.diff(edges)))
        if area > 0:
            vals = vals / area
        ax.step(edges, np.r_[vals, vals[-1]], where="post",
                linewidth=2, color=item["color"], label=item["label"])

    if h_data_ptr is not None:
        _, data_vals, _ = th1_to_np(h_data_ptr.GetValue())
        area = float(np.sum(data_vals * np.diff(edges)))
        if area > 0:
            data_vals = data_vals / area
        data_col = "white" if _DARK_MODE else "black"
        ax.step(edges, np.r_[data_vals, data_vals[-1]], where="post",
                linewidth=2, color=data_col, linestyle="--", label="Data")

    ax.set_xlabel(var)
    ax.set_ylabel("Normalised to unity")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2)

    if logy:
        ax.set_yscale("log")

    return fig, ax


# ---------------------------------------------------------------------------
# Stacked MC + data with trigger efficiency ratio panel
# ---------------------------------------------------------------------------

def plot_stacked_with_efficiency(
    h_mc_list, mc_items, *,
    h_data=None,
    h_mc_denom_list=None,
    h_data_denom=None,
    logy=True,
    sort_mc_by_yield=True,
    title=None,
):
    """Stacked MC plot with a trigger efficiency ratio panel below.

    Top panel: stacked MC + optional data overlay (same as plot_stacked_from_hists).
    Bottom panel: trigger efficiency = numerator / denominator for combined MC and data.

    Parameters
    ----------
    h_mc_list : list of TH1
        MC histograms WITH trigger applied (numerator), one per group.
    mc_items : list of dict
        Label/color info per MC group.
    h_data : TH1 or None
        Data histogram WITH trigger applied (numerator).
    h_mc_denom_list : list of TH1 or None
        MC histograms WITHOUT trigger (denominator), same grouping as h_mc_list.
    h_data_denom : TH1 or None
        Data histogram WITHOUT trigger (denominator).
    """
    import ROOT as _ROOT

    if h_data is not None:
        edges, data_vals, data_errs = th1_to_np(h_data)
    else:
        edges, _, _ = th1_to_np(h_mc_list[0])
        data_vals = None
        data_errs = None

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    mc_components = _build_mc_components(h_mc_list, mc_items, sort_by_yield=sort_mc_by_yield)

    # ── Figure with two panels ──
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    # ── Top panel: stacked MC + data ──
    stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
    _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
    if data_vals is not None:
        _draw_data(ax, centers, data_vals, data_errs)

    ax.set_ylabel("Events")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=10)
    ax.tick_params(labelbottom=False)

    if logy:
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        if data_vals is not None:
            ymax = max(ymax, float(np.max(data_vals)))
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))

    # ── Bottom panel: trigger efficiency ──
    eff_col_mc = "tab:blue"
    eff_col_data = "white" if _DARK_MODE else "black"

    # MC efficiency: sum all groups for numerator and denominator
    if h_mc_denom_list is not None and len(h_mc_denom_list) > 0:
        h_mc_num_combined = h_mc_list[0].Clone("_eff_mc_num")
        for h in h_mc_list[1:]:
            h_mc_num_combined.Add(h)

        h_mc_den_combined = h_mc_denom_list[0].Clone("_eff_mc_den")
        for h in h_mc_denom_list[1:]:
            h_mc_den_combined.Add(h)

        teff_mc = _ROOT.TEfficiency(h_mc_num_combined, h_mc_den_combined)
        teff_mc.SetUseWeightedEvents(True)
        x_mc, y_mc, xerr_mc, ylo_mc, yhi_mc = teff_to_np(teff_mc)

        rax.errorbar(x_mc, y_mc, xerr=xerr_mc,
                     yerr=np.vstack([ylo_mc, yhi_mc]),
                     fmt="s", color=eff_col_mc, ms=4, capsize=2,
                     linewidth=1, label="MC")

    # Data efficiency
    if h_data is not None and h_data_denom is not None:
        teff_data = _ROOT.TEfficiency(h_data, h_data_denom)
        x_d, y_d, xerr_d, ylo_d, yhi_d = teff_to_np(teff_data)

        rax.errorbar(x_d, y_d, xerr=xerr_d,
                     yerr=np.vstack([ylo_d, yhi_d]),
                     fmt="o", color=eff_col_data, ms=4, capsize=2,
                     linewidth=1, label="Data")

    rax.set_ylabel("Trigger Eff.")
    rax.set_ylim(0.0, 1.15)
    rax.grid(True, axis="y", alpha=0.3)
    rax.legend(loc="lower right", fontsize=10)

    return fig, ax, rax


# ---------------------------------------------------------------------------
# Stacked MC + data with S/sqrt(B) ratio panel
# ---------------------------------------------------------------------------

def plot_stacked_with_significance(
    h_mc_list, mc_items, *,
    h_data=None,
    sig_indices=None,
    bkg_indices=None,
    logy=True,
    sort_mc_by_yield=True,
    title=None,
):
    """Stacked MC plot with a per-bin S/sqrt(B) panel below.

    Top panel: stacked MC + optional data overlay.
    Bottom panel: S/sqrt(B) per bin for each signal channel individually.

    Parameters
    ----------
    h_mc_list : list of TH1
        MC histograms (one per group), already scaled to luminosity.
    mc_items : list of dict
        Label/color info per MC group.
    h_data : TH1 or None
        Data histogram.
    sig_indices : list of int
        Indices into h_mc_list that are signal channels (plotted individually).
    bkg_indices : list of int
        Indices into h_mc_list that are background (summed for B).
    """
    if sig_indices is None or bkg_indices is None:
        raise ValueError("sig_indices and bkg_indices must be provided")

    if h_data is not None:
        edges, data_vals, data_errs = th1_to_np(h_data)
    else:
        edges, _, _ = th1_to_np(h_mc_list[0])
        data_vals = None
        data_errs = None

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    mc_components = _build_mc_components(h_mc_list, mc_items, sort_by_yield=sort_mc_by_yield)

    # ── Figure with two panels ──
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    # ── Top panel: stacked MC + data ──
    stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
    _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
    if data_vals is not None:
        _draw_data(ax, centers, data_vals, data_errs)

    ax.set_ylabel("Events")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=10)
    ax.tick_params(labelbottom=False)

    if logy:
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        if data_vals is not None:
            ymax = max(ymax, float(np.max(data_vals)))
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))

    # ── Bottom panel: S / sqrt(B) per bin ──
    # Sum background bins
    bkg_vals = np.zeros(len(centers), dtype=float)
    for bi in bkg_indices:
        _, bv, _ = th1_to_np(h_mc_list[bi])
        bkg_vals += np.maximum(bv, 0.0)

    # Plot each signal channel
    for si in sig_indices:
        _, sv, _ = th1_to_np(h_mc_list[si])
        sv = np.maximum(sv, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            s_over_sqrtb = np.where(bkg_vals > 0, sv / np.sqrt(bkg_vals), 0.0)
        item = mc_items[si]
        rax.step(centers, s_over_sqrtb, where="mid",
                 linewidth=1.5, color=item["color"], label=item["label"])

    rax.set_ylabel(r"$S\,/\,\sqrt{B}$")
    rax.set_ylim(bottom=0)
    rax.grid(True, axis="y", alpha=0.3)
    rax.legend(loc="upper right", fontsize=9)

    return fig, ax, rax


# ---------------------------------------------------------------------------
# Pre-materialized histogram versions  (no RunGraphs — caller provides TH1s)
# ---------------------------------------------------------------------------

def plot_stacked_from_hists(
    h_mc_list, mc_items, *,
    h_data=None,
    logy=True,
    sort_mc_by_yield=True,
    title=None,
):
    """Stacked MC plot from pre-materialized TH1 objects (one per mc_item group).

    Same visual output as plot_stacked_all_mc but skips histogram booking
    and RunGraphs — the caller is responsible for providing ready TH1s.
    """
    if h_data is not None:
        edges, data_vals, data_errs = th1_to_np(h_data)
    else:
        edges, _, _ = th1_to_np(h_mc_list[0])
        data_vals = None
        data_errs = None

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    mc_components = _build_mc_components(h_mc_list, mc_items, sort_by_yield=sort_mc_by_yield)

    fig, ax = plt.subplots(figsize=(8, 6))
    stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
    _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
    if data_vals is not None:
        _draw_data(ax, centers, data_vals, data_errs)

    ax.set_ylabel("Events")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2)

    if logy:
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        if data_vals is not None:
            ymax = max(ymax, float(np.max(data_vals)))
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))

    return fig, ax


def plot_shape_from_hists(
    h_mc_list, mc_items, *,
    h_data=None,
    logy=False,
    title=None,
):
    """Shape overlay from pre-materialized TH1 objects (one per mc_item group).

    Same visual output as plot_shape_overlay but skips histogram booking
    and RunGraphs — the caller is responsible for providing ready TH1s.
    """
    edges = th1_to_np(h_mc_list[0])[0]

    fig, ax = plt.subplots(figsize=(8, 6))

    for h, item in zip(h_mc_list, mc_items):
        _, vals, _ = th1_to_np(h)
        area = float(np.sum(vals * np.diff(edges)))
        if area > 0:
            vals = vals / area
        ax.step(edges, np.r_[vals, vals[-1]], where="post",
                linewidth=2, color=item["color"], label=item["label"])

    if h_data is not None:
        _, data_vals, _ = th1_to_np(h_data)
        area = float(np.sum(data_vals * np.diff(edges)))
        if area > 0:
            data_vals = data_vals / area
        data_col = "white" if _DARK_MODE else "black"
        ax.step(edges, np.r_[data_vals, data_vals[-1]], where="post",
                linewidth=2, color=data_col, linestyle="--", label="Data")

    ax.set_ylabel("Normalised to unity")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend(ncol=2)

    if logy:
        ax.set_yscale("log")

    return fig, ax


# ---------------------------------------------------------------------------
# Trigger shape overlay  (all triggers on one axes, normalised to unity)
# ---------------------------------------------------------------------------

TRIGGER_STYLES = {
    "NoTrigger":  {"color": "tab:gray",   "linestyle": "-",  "linewidth": 2.5},
    "DST_JetHT":  {"color": "tab:blue",   "linestyle": "--", "linewidth": 2.5},
    "PARKING_HH": {"color": "tab:orange",  "linestyle": ":",  "linewidth": 3.0},
}


def plot_trigger_shape_overlay(h_total_by_trig, *, logy=False, title=None):
    """Shape overlay comparing the same variable across triggers.

    Parameters
    ----------
    h_total_by_trig : dict  {trig_name: TH1}
        One combined (all-MC-groups-summed) TH1 per trigger, already lumi-scaled.
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    for trig_name, h in h_total_by_trig.items():
        style = TRIGGER_STYLES.get(trig_name, {"color": "black", "linestyle": "-", "linewidth": 2})
        edges, vals, _ = th1_to_np(h)
        area = float(np.sum(vals * np.diff(edges)))
        if area > 0:
            vals = vals / area
        ax.step(edges, np.r_[vals, vals[-1]], where="post",
                label=trig_name, **style)

    ax.set_ylabel("Normalised to unity")
    if title is not None:
        ax.set_title(title)
    ax.grid(True, axis="both", alpha=0.25)
    ax.legend()
    if logy:
        ax.set_yscale("log")

    return fig, ax


# ── 2D histogram (heatmap) ──────────────────────────────────────────────────

def th2_to_np(h2):
    """Convert ROOT TH2 to numpy arrays: (xedges, yedges, vals)."""
    nx = h2.GetNbinsX()
    ny = h2.GetNbinsY()
    xedges = np.array([h2.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)])
    yedges = np.array([h2.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)])
    vals = np.zeros((ny, nx))
    for iy in range(1, ny + 1):
        for ix in range(1, nx + 1):
            vals[iy - 1, ix - 1] = h2.GetBinContent(ix, iy)
    return xedges, yedges, vals


def plot_2d_hist(h2, *, xlabel=None, ylabel=None, title=None,
                 log_z=True, cmap="viridis"):
    """Plot a ROOT TH2 as a matplotlib pcolormesh heatmap."""
    from matplotlib.colors import LogNorm

    xedges, yedges, vals = th2_to_np(h2)

    fig, ax = plt.subplots(figsize=(8, 6))
    if log_z and vals.max() > 0 and np.any(vals > 0):
        vmin = max(vals[vals > 0].min(), 1e-1)
        norm = LogNorm(vmin=vmin, vmax=vals.max())
    else:
        norm = None
    mesh = ax.pcolormesh(xedges, yedges, vals, cmap=cmap, norm=norm)
    fig.colorbar(mesh, ax=ax, label="Events")
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    hep.cms.label("Work in Progress", data=False, ax=ax)
    return fig
