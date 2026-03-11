#!/usr/bin/env python3
"""Collect all analysis plots into a single PDF.

Usage:
    python make_pdf.py                    # default: dark theme
    python make_pdf.py --theme light
    python make_pdf.py --theme dark -o my_plots.pdf
"""

import argparse
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
                    help="Output PDF path (default: plots/<theme>/all_plots.pdf)")
args = parser.parse_args()

PLOT_DIR = os.path.join("plots", args.theme)
out_path = args.output or os.path.join(PLOT_DIR, "all_plots.pdf")

# ── Discover all PNGs ──
all_pngs = sorted(glob.glob(os.path.join(PLOT_DIR, "*.png")))
if not all_pngs:
    print(f"No PNGs found in {PLOT_DIR}")
    raise SystemExit(1)

# ── Identify triggers and variables from new naming: {TRIG}_{type}_{var}.png ──
TRIGGERS = ["DST_JetHT", "PARKING_HH"]

# Build a lookup: (trigger, plot_type, var_name) -> filepath
# New naming: DST_JetHT_stacked_ak4_pt0.png  -> trig=DST_JetHT, type=stacked, var=ak4_pt0
# Old naming: stacked_ak4_pt0_DST_JetHT.png  -> same but different pattern
# Also: DST_JetHT_shape_gen_mHH.png (special)

plot_map = {}  # (trigger, var_name) -> {"stacked": path, "shape": path}

for png in all_pngs:
    fname = os.path.basename(png)
    name = fname.replace(".png", "")

    trig = None
    ptype = None
    var = None

    # Try new naming: {TRIG}_{type}_{var}
    for t in TRIGGERS:
        if name.startswith(t + "_"):
            rest = name[len(t) + 1:]
            for pt in ("stacked", "shape"):
                if rest.startswith(pt + "_"):
                    trig = t
                    ptype = pt
                    var = rest[len(pt) + 1:]
                    break
            break

    # Try old naming: {type}_{var}_{TRIG}
    if trig is None:
        for t in TRIGGERS:
            if name.endswith("_" + t):
                rest = name[: -(len(t) + 1)]
                for pt in ("stacked", "shape"):
                    if rest.startswith(pt + "_"):
                        trig = t
                        ptype = pt
                        var = rest[len(pt) + 1:]
                        break
                break

    if trig is None or ptype is None or var is None:
        continue

    key = (trig, var)
    if key not in plot_map:
        plot_map[key] = {}
    # Prefer new naming convention over old (overwrite if both exist)
    if ptype not in plot_map[key] or name.startswith(trig):
        plot_map[key][ptype] = png

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
    # Special
    "gen_mHH",
]

def var_sort_key(var):
    try:
        return VAR_ORDER.index(var)
    except ValueError:
        return len(VAR_ORDER)


# ── Generate PDF ──
print(f"Generating PDF from {len(all_pngs)} images in {PLOT_DIR} ...")

with PdfPages(out_path) as pdf:
    for trig in TRIGGERS:
        # Get all variables for this trigger, sorted
        vars_for_trig = sorted(
            [var for (t, var) in plot_map if t == trig],
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
