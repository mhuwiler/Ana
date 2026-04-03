"""
Plotting utilities for scouting analysis.

Contains:
  - ROOT TH1 / TH2 -> numpy converters
  - Stacked MC + data plot (with optional ratio panels)
  - Shape overlay (normalised to unity)
  - Trigger shape overlay
  - 2D histogram heatmap
  - CMS style setup via mplhep
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep


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
# MC component helpers
# ---------------------------------------------------------------------------

_neg_yield_warned = set()


def _build_mc_components(h_mc_list, mc_items, sort_by_yield=True, warn=True):
    """Convert MC histograms (TH1 objects) to numpy arrays (vals, errs), optionally sorted."""
    components = []
    neg_groups = []
    for idx, (h, item) in enumerate(zip(h_mc_list, mc_items)):
        _, vals_raw, errs_raw = th1_to_np(h)
        if np.any(vals_raw < 0):
            neg_groups.append(item['label'])
        vals_clipped = np.maximum(vals_raw, 0.0)
        components.append({
            "item": item,
            "vals": vals_clipped,
            "errs": errs_raw,
            "yield": float(np.sum(vals_clipped)),
            "orig_idx": idx,
        })

    if warn and neg_groups:
        key = frozenset(neg_groups)
        if key not in _neg_yield_warned:
            _neg_yield_warned.add(key)
            names = ", ".join(neg_groups)
            print(f"  [plots] {len(neg_groups)} group(s) clipped negative bins to 0 (NLO): {names}")

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

    Top panel: stacked MC + optional data overlay (same as plot_stacked).
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
# Two-panel stacked MC + ratio helpers
# ---------------------------------------------------------------------------

def _setup_significance_figure(h_mc_list, mc_items, *, h_data=None,
                               sig_indices, bkg_indices, logy=True,
                               sort_mc_by_yield=True, title=None):
    """Shared setup for significance plots: two-panel figure with stacked top.

    Returns (fig, ax, rax, centers, bkg_vals) — caller draws the ratio panel.
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

    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

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

    # Sum background
    bkg_vals = np.zeros(len(centers), dtype=float)
    for bi in bkg_indices:
        _, bv, _ = th1_to_np(h_mc_list[bi])
        bkg_vals += np.maximum(bv, 0.0)

    return fig, ax, rax, centers, bkg_vals


def plot_stacked_with_significance(
    h_mc_list, mc_items, *,
    h_data=None, sig_indices=None, bkg_indices=None,
    logy=True, sort_mc_by_yield=True, title=None,
):
    """Stacked MC + per-bin S/sqrt(B) ratio panel."""
    fig, ax, rax, centers, bkg_vals = _setup_significance_figure(
        h_mc_list, mc_items, h_data=h_data,
        sig_indices=sig_indices, bkg_indices=bkg_indices,
        logy=logy, sort_mc_by_yield=sort_mc_by_yield, title=title)

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


def plot_stacked_with_cumulative_significance(
    h_mc_list, mc_items, *,
    h_data=None, sig_indices=None, bkg_indices=None,
    logy=True, sort_mc_by_yield=True, title=None,
):
    """Stacked MC + cumulative S/sqrt(B) ratio panel (both cut directions)."""
    fig, ax, rax, centers, bkg_vals = _setup_significance_figure(
        h_mc_list, mc_items, h_data=h_data,
        sig_indices=sig_indices, bkg_indices=bkg_indices,
        logy=logy, sort_mc_by_yield=sort_mc_by_yield, title=title)

    # Per-channel, both directions
    bkg_r = np.cumsum(bkg_vals[::-1])[::-1]
    bkg_l = np.cumsum(bkg_vals)
    for si in sig_indices:
        _, sv, _ = th1_to_np(h_mc_list[si])
        sv = np.maximum(sv, 0.0)
        item = mc_items[si]

        sig_r = np.cumsum(sv[::-1])[::-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            cum_r = np.where(bkg_r > 0, sig_r / np.sqrt(bkg_r), 0.0)
        sig_l = np.cumsum(sv)
        with np.errstate(divide="ignore", invalid="ignore"):
            cum_l = np.where(bkg_l > 0, sig_l / np.sqrt(bkg_l), 0.0)

        rax.step(centers, cum_r, where="mid", linewidth=1.5,
                 color=item["color"], linestyle="-",
                 label=f'{item["label"]} $>x$')
        rax.step(centers, cum_l, where="mid", linewidth=1.5,
                 color=item["color"], linestyle="--",
                 label=f'{item["label"]} $<x$')

    rax.set_ylabel(r"Cumulative $S\,/\,\sqrt{B}$")
    rax.set_ylim(bottom=0)
    rax.grid(True, axis="y", alpha=0.3)
    rax.legend(loc="upper right", fontsize=7, ncol=2)
    return fig, ax, rax


# ---------------------------------------------------------------------------
# Stacked MC + data plot
# ---------------------------------------------------------------------------

def plot_stacked(
    h_mc_list, mc_items, *,
    h_data=None,
    logy=True,
    sort_mc_by_yield=True,
    title=None,
):
    """Stacked MC plot from pre-materialized TH1 objects (one per mc_item group).

    The caller is responsible for providing ready TH1s.
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


def plot_shape(
    h_mc_list, mc_items, *,
    h_data=None,
    logy=False,
    title=None,
):
    """Shape overlay from pre-materialized TH1 objects (one per mc_item group).

    Each MC group is normalised to unity. The caller provides ready TH1s.
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
    "NoTrigger":         {"color": "tab:gray",   "linestyle": "-",  "linewidth": 2.5},
    "DST_JetHT":         {"color": "tab:blue",   "linestyle": "--", "linewidth": 2.5},
    "DST_Muon":          {"color": "tab:green",  "linestyle": "-.", "linewidth": 2.5},
    "DST_Electron":      {"color": "tab:red",    "linestyle": "--", "linewidth": 2.5},
    "PARKING_HH":        {"color": "tab:orange", "linestyle": ":",  "linewidth": 3.0},
    # Exclusive variants (same hue, thinner line)
    "DST_JetHT_excl":    {"color": "tab:blue",   "linestyle": "--", "linewidth": 1.5},
    "DST_Muon_excl":     {"color": "tab:green",  "linestyle": "-.", "linewidth": 1.5},
    "DST_Electron_excl": {"color": "tab:red",    "linestyle": "--", "linewidth": 1.5},
    "PARKING_HH_excl":   {"color": "tab:orange", "linestyle": ":",  "linewidth": 2.0},
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
                 log_z=True, cmap="viridis", lumi=None):
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
    cms_label(ax, lumi=lumi)
    fig.tight_layout()
    return fig


# ── Confusion matrix ────────────────────────────────────────────────────────

def plot_confusion_matrix(matrix, row_labels, col_labels, *,
                          title=None, lumi=None):
    """Plot a confusion matrix as an annotated seaborn heatmap.

    Parameters
    ----------
    matrix : np.ndarray
        2D array of fractions (rows=truth, cols=predicted). Each row sums to 1.
    row_labels, col_labels : list[str]
        Labels for rows (y-axis, truth) and columns (x-axis, predicted).
    """
    import seaborn as sns

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(matrix, annot=True, fmt=".1%", cmap="Blues",
                xticklabels=col_labels, yticklabels=row_labels,
                vmin=0, vmax=1, linewidths=1, linecolor="white",
                cbar_kws={"label": "Fraction"}, ax=ax)
    ax.set_xlabel("Tagger prediction")
    ax.set_ylabel("Gen truth")
    if title:
        ax.set_title(title)
    cms_label(ax, lumi=lumi)
    fig.tight_layout()
    return fig


# ── Parallel rendering ──────────────────────────────────────────────────────

def prepare_mc_task(h_mc_list, mc_items, h_data=None, warn=True):
    """Convert ROOT histograms to picklable numpy data for parallel rendering."""
    edges = th1_to_np(h_mc_list[0])[0]
    mc_components = _build_mc_components(h_mc_list, mc_items, warn=warn)
    data = th1_to_np(h_data) if h_data is not None else None
    return edges, mc_components, data


def render_task(task):
    """Render one plot from pre-converted numpy data. Runs in subprocess."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mplhep as hep

    ptype = task["type"]
    outpath = task["outpath"]
    xlabel = task["xlabel"]
    lumi = task["lumi"]
    edges = task["edges"]
    mc_components = task["mc_components"]
    data_tuple = task.get("data")
    overwrite = task.get("overwrite", True)

    import os
    if not overwrite and os.path.exists(outpath):
        return outpath
    os.makedirs(os.path.dirname(outpath), exist_ok=True)

    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)
    data_vals = data_tuple[1] if data_tuple else None
    data_errs = data_tuple[2] if data_tuple else None

    if ptype == "stacked":
        fig, ax = plt.subplots(figsize=(8, 6))
        stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
        _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
        if data_vals is not None:
            _draw_data(ax, centers, data_vals, data_errs)
        ax.set_ylabel("Events")
        ax.grid(True, axis="both", alpha=0.25)
        ax.legend(ncol=2)
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        if data_vals is not None:
            ymax = max(ymax, float(np.max(data_vals)))
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))
        xlabel_ax = ax

    elif ptype == "shape":
        fig, ax = plt.subplots(figsize=(8, 6))
        for comp in mc_components:
            item = comp["item"]
            vals = comp["vals"]
            area = float(np.sum(vals * widths))
            norm_vals = vals / area if area > 0 else vals
            ax.step(edges, np.r_[norm_vals, norm_vals[-1]], where="post",
                    linewidth=1.5, color=item["color"], label=item["label"])
        if data_vals is not None:
            area_d = float(np.sum(data_vals * widths))
            if area_d > 0:
                ax.errorbar(centers, data_vals / area_d, yerr=data_errs / area_d,
                            fmt="o", color="black", label="Data", ms=4, capsize=2, zorder=10)
        ax.set_ylabel("Normalised to unity")
        ax.grid(True, axis="both", alpha=0.25)
        ax.legend()
        xlabel_ax = ax

    elif ptype in ("sig", "cuml_sig"):
        # ── Shared two-panel top: stacked MC ──
        fig = plt.figure(figsize=(8, 8))
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax = fig.add_subplot(gs[0])
        rax = fig.add_subplot(gs[1], sharex=ax)
        stack_total = _draw_mc_stack(ax, centers, widths, mc_components)
        _draw_mc_stat_unc(ax, edges, stack_total, mc_components)
        if data_vals is not None:
            _draw_data(ax, centers, data_vals, data_errs)
        ax.set_ylabel("Events")
        ax.grid(True, axis="both", alpha=0.25)
        ax.legend(ncol=2, fontsize=10)
        ax.tick_params(labelbottom=False)
        ax.set_yscale("log")
        ymax = float(np.max(stack_total)) if len(stack_total) else 1.0
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))

        # ── Separate S and B from mc_components ──
        # Use orig_idx to match against pre-sort sig/bkg indices
        sig_indices = set(task.get("sig_indices", []))
        bkg_indices = set(task.get("bkg_indices", []))
        sig_comps = []  # [(vals, item), ...] per signal channel
        bkg_total = np.zeros_like(centers)
        for comp in mc_components:
            oi = comp.get("orig_idx", -1)
            if oi in sig_indices:
                sig_comps.append((comp["vals"], comp["item"]))
            elif oi in bkg_indices:
                bkg_total += comp["vals"]

        # ── Bottom panel: bin-by-bin or cumulative S/√B ──
        all_ratio_vals = []
        if ptype == "sig":
            for sv, item in sig_comps:
                with np.errstate(divide="ignore", invalid="ignore"):
                    rv = np.where(bkg_total > 0, sv / np.sqrt(bkg_total), 0.0)
                rax.step(edges, np.r_[rv, rv[-1]], where="post",
                         linewidth=1.5, color=item["color"], label=item["label"])
                all_ratio_vals.append(rv)
            rax.set_ylabel(r"$S/\sqrt{B}$")
            rax.legend(loc="upper right", fontsize=9)
        else:  # cuml_sig — per channel, both directions
            bkg_r = np.cumsum(bkg_total[::-1])[::-1]
            bkg_l = np.cumsum(bkg_total)
            for sv, item in sig_comps:
                sig_r = np.cumsum(sv[::-1])[::-1]
                with np.errstate(divide="ignore", invalid="ignore"):
                    cum_r = np.where(bkg_r > 0, sig_r / np.sqrt(bkg_r), 0.0)
                sig_l = np.cumsum(sv)
                with np.errstate(divide="ignore", invalid="ignore"):
                    cum_l = np.where(bkg_l > 0, sig_l / np.sqrt(bkg_l), 0.0)
                rax.step(edges, np.r_[cum_r, cum_r[-1]], where="post",
                         linewidth=1.5, color=item["color"], linestyle="-",
                         label=f'{item["label"]} $>x$')
                rax.step(edges, np.r_[cum_l, cum_l[-1]], where="post",
                         linewidth=1.5, color=item["color"], linestyle="--",
                         label=f'{item["label"]} $<x$')
                all_ratio_vals.extend([cum_r, cum_l])
            rax.set_ylabel(r"Cumulative $S/\sqrt{B}$")
            rax.legend(loc="upper right", fontsize=7, ncol=2)

        rax.grid(True, axis="both", alpha=0.25)
        peak = max((float(np.max(rv)) for rv in all_ratio_vals), default=0.0)
        rax.set_ylim(0, max(peak * 1.3, 1e-6))
        xlabel_ax = rax
    else:
        return outpath

    xlabel_ax.set_xlabel(xlabel)
    hep.cms.label("Work in Progress", ax=ax, data=True,
                  lumi=f"{lumi:.4g}" if lumi is not None else None,
                  year="2024", com=13.6)
    if ptype in ("sig", "cuml_sig"):
        fig.subplots_adjust(hspace=0.05)
    else:
        fig.tight_layout()
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    return outpath
