#!/usr/bin/env python3
"""Collect analysis plots into a PDF.

Usage:
    python make_pdf.py                                          # all vars, dark theme
    python make_pdf.py --theme light                           # light theme
    python make_pdf.py --theme light --compare --plot-type shape stacked --vars "ak4_pt*"
    python make_pdf.py --theme light --compare --vars HT MHT b0_score
    python make_pdf.py --theme light --compare --triggers DST_JetHT PARKING_HH
    python make_pdf.py --theme dark -o my_plots.pdf

--vars accepts fnmatch patterns: ak4_pt*, b*_score, HT, etc.
--triggers selects which columns appear (default: NoTrigger DST_JetHT PARKING_HH)

Plot directory structure:
    plots/{theme}/
    ├── mc/stacked/{trig}/{var}.png
    ├── mc/shape/{trig}/{var}.png
    ├── mc/eff/{trig}/{var}.png
    ├── mc/sig/{trig}/{var}.png
    ├── mc/overlay/{var}.png, mc/overlay/{ch}_{var}.png
    ├── data/stacked/{trig}/{var}.png
    ├── data/shape/{trig}/{var}.png
    ├── data/eff/{trig}/{var}.png
    └── data/sig/{trig}/{var}.png
"""

import argparse
import fnmatch
import os
import glob

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.image import imread

parser = argparse.ArgumentParser()
parser.add_argument("--theme", choices=["light", "dark"], default="dark")
parser.add_argument("-o", "--output", default=None,
                    help="Output PDF path (default: plots/<theme>/all_plots.pdf or comparison_<type>.pdf)")
parser.add_argument("--compare", action="store_true",
                    help="Generate a comparison PDF with all triggers side by side per variable")
parser.add_argument("--plot-type", nargs="+",
                    choices=["stacked", "shape", "sig", "eff",
                             "data_stacked", "data_shape", "data_eff", "data_sig",
                             "trigger_overlay"],
                    default=["shape"],
                    help="Plot type(s) for comparison PDF — each becomes a row (default: shape)")
parser.add_argument("--vars", nargs="+", default=None, metavar="PATTERN",
                    help="Variable name patterns to include (fnmatch, e.g. 'ak4_pt*' 'HT'). Default: all.")
parser.add_argument("--triggers", nargs="+",
                    default=["NoTrigger", "DST_JetHT", "PARKING_HH"],
                    choices=["NoTrigger", "DST_JetHT", "PARKING_HH"],
                    help="Triggers to show as columns (default: all three).")
args = parser.parse_args()

PLOT_DIR = os.path.join("plots", args.theme)
out_path = args.output or os.path.join(PLOT_DIR, "all_plots.pdf")

# ── Discover all PNGs in subdirectories ──
all_pngs = sorted(glob.glob(os.path.join(PLOT_DIR, "**", "*.png"), recursive=True))
if not all_pngs:
    print(f"No PNGs found in {PLOT_DIR}")
    raise SystemExit(1)

# ── Identify triggers and variables from folder structure ──
TRIGGERS = ["NoTrigger", "DST_JetHT", "PARKING_HH"]

# New directory structure: mc/{plot_type}/{trig}/{var}.png, data/{plot_type}/{trig}/{var}.png
# Overlay: mc/overlay/{var}.png, mc/overlay/{ch}_{var}.png

# Build a lookup: (trigger, var_name) -> {plot_type: filepath}
plot_map = {}

for png in all_pngs:
    rel = os.path.relpath(png, PLOT_DIR)
    parts = rel.split(os.sep)

    # mc/overlay/{var}.png or mc/overlay/{ch}_{var}.png (3 parts)
    if len(parts) == 3 and parts[0] == "mc" and parts[1] == "overlay":
        name = parts[2].replace(".png", "")
        overlay_key = "_overlay"
        var = name
        for ch in ("hh_", "hm_", "he_"):
            if name.startswith(ch):
                overlay_key = f"_overlay_{ch[:-1]}"
                var = name[len(ch):]
                break
        key = (overlay_key, var)
        if key not in plot_map:
            plot_map[key] = {}
        plot_map[key]["trigger_overlay"] = png
        continue

    # mc/{plot_type}/{trig}/{var}.png or data/{plot_type}/{trig}/{var}.png (4 parts)
    if len(parts) != 4:
        continue
    category, ptype, trig, fname = parts
    if category not in ("mc", "data"):
        continue
    if trig not in TRIGGERS:
        continue
    var = fname.replace(".png", "")

    # Map to unified plot_type key (prepend "data_" for data category)
    ptype_key = f"data_{ptype}" if category == "data" else ptype

    key = (trig, var)
    if key not in plot_map:
        plot_map[key] = {}
    plot_map[key][ptype_key] = png

# ── Order variables logically ──
VAR_ORDER = [
    # Jet kinematics
    "ak4_pt0", "ak4_pt1", "ak4_pt2", "ak4_pt3",
    "ak4_eta0", "ak4_eta1", "ak4_eta2", "ak4_eta3",
    "ak4_mass0", "ak4_mass1",
    # Global
    "nJets", "nElectrons", "nMuons", "nLeptons",
    "HT", "MHT", "centrality",
    # Dijet / 4-jet
    "mjj_01", "dR_01", "dEta_01", "m4j",
    # b-tagged jets
    "b0_pt", "b1_pt",
    "b0_score", "b1_score",
    "b0_raw", "b1_raw",
    "mbb", "dR_bb", "ptbb",
    # Higgs candidates
    "mbb_cand", "mtautau_cand",
    # Special
    "gen_mHH", "mHH_m4j",
]

def var_sort_key(var):
    try:
        return VAR_ORDER.index(var)
    except ValueError:
        return len(VAR_ORDER)


# ══════════════════════════════════════════════════════════════════════════════
# Mode 1 (default): all_plots.pdf — per trigger, per variable
# ══════════════════════════════════════════════════════════════════════════════

if not args.compare:
    print(f"Generating PDF from {len(all_pngs)} images in {PLOT_DIR} ...")

    with PdfPages(out_path) as pdf:
        for trig in args.triggers:
            vars_for_trig = sorted(
                [var for (t, var) in plot_map if t == trig
                 and (not args.vars or any(fnmatch.fnmatch(var, p) for p in args.vars))],
                key=var_sort_key,
            )

            for var in vars_for_trig:
                entry = plot_map[(trig, var)]
                has_stacked = "stacked" in entry
                has_shape = "shape" in entry
                ncols = int(has_stacked) + int(has_shape)
                if ncols == 0:
                    continue

                if ncols == 2:
                    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
                    for ax, ptype in zip(axes, ["stacked", "shape"]):
                        img = imread(entry[ptype])
                        ax.imshow(img)
                        ax.set_axis_off()
                    fig.suptitle(f"{trig}  —  {var}", fontsize=14, y=0.98)
                else:
                    fig, ax = plt.subplots(1, 1, figsize=(10, 7))
                    ptype = "stacked" if has_stacked else "shape"
                    img = imread(entry[ptype])
                    ax.imshow(img)
                    ax.set_axis_off()
                    fig.suptitle(f"{trig}  —  {var}  ({ptype})", fontsize=14, y=0.98)

                fig.tight_layout(rect=[0, 0, 1, 0.96])
                pdf.savefig(fig)
                plt.close(fig)

    print(f"Saved: {out_path}  ({os.path.getsize(out_path) / 1e6:.1f} MB)")


# ══════════════════════════════════════════════════════════════════════════════
# Mode 2 (--compare): comparison_{plot_type}.pdf — all triggers side by side
# ══════════════════════════════════════════════════════════════════════════════

else:
    pts = args.plot_type  # list, e.g. ["shape", "stacked"]
    suffix = "_".join(pts)
    var_tag = "_".join(args.vars).replace("*", "X") if args.vars else "all"
    cmp_path = args.output or os.path.join(PLOT_DIR, f"comparison_{suffix}_{var_tag}.pdf")
    print(f"Generating comparison PDF ({', '.join(pts)}) from {PLOT_DIR} ...")

    all_vars = sorted(
        {var for (trig, var) in plot_map},
        key=var_sort_key,
    )

    # ── Filter variables by pattern ──
    if args.vars:
        all_vars = [v for v in all_vars if any(fnmatch.fnmatch(v, pat) for pat in args.vars)]
        if not all_vars:
            print(f"No variables matched patterns: {args.vars}")
            raise SystemExit(1)
        print(f"Variables selected ({len(all_vars)}): {', '.join(all_vars)}")

    with PdfPages(cmp_path) as pdf:
        for var in all_vars:
            # Build grid: rows = plot types, cols = triggers (only those that have the image)
            rows = []
            for pt in pts:
                if pt == "trigger_overlay":
                    # Overlay: total + per-channel (hh, hm, he) as columns
                    ov_imgs, ov_labels = [], []
                    for ov_key, ov_lbl in [("_overlay", "All MC"),
                                           ("_overlay_hh", r"$bb\tau_h\tau_h$"),
                                           ("_overlay_hm", r"$bb\tau_\mu\tau_h$"),
                                           ("_overlay_he", r"$bb\tau_e\tau_h$")]:
                        path = plot_map.get((ov_key, var), {}).get("trigger_overlay")
                        if path:
                            ov_imgs.append(imread(path))
                            ov_labels.append(ov_lbl)
                    if ov_imgs:
                        rows.append((pt, ov_imgs, ov_labels))
                    continue

                row_imgs, row_labels = [], []
                for trig in args.triggers:   # respects --triggers selection
                    path = plot_map.get((trig, var), {}).get(pt)
                    if path:
                        row_imgs.append(imread(path))
                        row_labels.append(trig)
                if row_imgs:
                    rows.append((pt, row_imgs, row_labels))

            if not rows:
                continue

            nrows = len(rows)
            ncols = max(len(r[1]) for r in rows)
            fig, axes = plt.subplots(nrows, ncols, figsize=(7 * ncols, 6 * nrows),
                                     squeeze=False)

            for ri, (pt, row_imgs, row_labels) in enumerate(rows):
                for ci in range(ncols):
                    ax = axes[ri][ci]
                    if ci < len(row_imgs):
                        ax.imshow(row_imgs[ci])
                        if ri == 0:
                            ax.set_title(row_labels[ci], fontsize=13)
                        ax.set_ylabel(pt, fontsize=11) if ci == 0 else None
                    ax.set_axis_off()

            fig.suptitle(var, fontsize=15, y=1.01)
            fig.tight_layout()
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)

    print(f"Saved: {cmp_path}  ({os.path.getsize(cmp_path) / 1e6:.1f} MB)")
