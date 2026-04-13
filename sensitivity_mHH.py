#!/usr/bin/env python3
"""Per-bin Asimov significance (Z_A) vs reconstructed m4j.

Shows which mass regions contribute most to the analysis sensitivity.
Uses full tauhtauh cuts from cuts.yaml. Signal and background both
binned by reconstructed m4j (4-jet invariant mass).
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES       = ["HHbbtt", "DY", "TT"]
DATA          = False
NTHREADS      = 32
MAX_MC_FILES  = 22

# Full tauhtauh cuts
CUTS = [
    ("NoTrigger", ["common", "tauhtauh"]),
]

# m4j for per-bin Z_A
PLOT_VARS = [
    PlotVar("m4j", r"Reconstructed $m_{4j}$ [GeV]", 19, 250, 1200),
]

LUMI = 103.965

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import os
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from analysis.runner import setup, load_and_run
    from utils.plotting import plot_per_bin_significance

    ctx = setup(__file__, samples=SAMPLES, cuts=CUTS, data=DATA,
                nthreads=NTHREADS, max_mc_files=MAX_MC_FILES)
    result = load_and_run(ctx, PLOT_VARS)

    # Get m4j histograms — signal vs background
    trig_name = "NoTrigger"
    h_by_var = result.mc_hists_by_trig[trig_name]
    h_list = h_by_var["m4j"]  # list of TH1, one per MC group

    sig_indices = result.sig_indices
    bkg_indices = result.bkg_indices

    # Sum signal groups
    h_sig = h_list[sig_indices[0]].Clone("h_m4j_sig")
    for gi in sig_indices[1:]:
        h_sig.Add(h_list[gi])

    # Sum background groups
    h_bkg = h_list[bkg_indices[0]].Clone("h_m4j_bkg")
    for gi in bkg_indices[1:]:
        h_bkg.Add(h_list[gi])

    print(f"\nSignal integral:     {h_sig.Integral():.2f}")
    print(f"Background integral: {h_bkg.Integral():.2f}")

    fig, ax = plot_per_bin_significance(
        h_sig, h_bkg,
        xlabel=r"Reconstructed $m_{4j}$ [GeV]",
        lumi=LUMI,
        title=r"$HH \to bb\tau\tau$ (full $\tau_h\tau_h$ cuts)")

    outpath = os.path.join(ctx.plot_dir, "sensitivity_m4j.png")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {outpath}")
