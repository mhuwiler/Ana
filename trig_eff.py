#!/usr/bin/env python3
"""Trigger efficiency vs gen_mHH — CMS-style overlay plot.

No analysis cuts applied — pure trigger efficiency on signal MC.
Produces trigger_efficiency_mHH.png with filled signal shape + efficiency curves.
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES       = ["HHbbtt"]
DATA          = False
NTHREADS      = 32
MAX_MC_FILES  = 22

# All triggers to measure efficiency for
TRIGGERS = ["NoTrigger", "DST_JetHT", "PARKING_HH"]

# No cuts — pure trigger efficiency
CUTS = []

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    import os
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from analysis.runner import setup, load_and_run
    from utils.plotting import plot_trigger_efficiency_overlay

    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, cuts=CUTS,
                data=DATA, nthreads=NTHREADS, max_mc_files=MAX_MC_FILES)
    result = load_and_run(ctx, [])

    trig_sel = result.trig_selections

    # Denominator: NoTrigger signal mHH (no trigger requirement)
    h_denom_ptrs = trig_sel["NoTrigger"]["sig_mHH_ptrs"]
    if not h_denom_ptrs:
        print("ERROR: No signal mHH histograms found for NoTrigger")
        exit(1)

    h_denom = h_denom_ptrs[0].GetPtr().Clone("h_mHH_denom")
    for p in h_denom_ptrs[1:]:
        h_denom.Add(p.GetPtr())

    # Numerator: per trigger (skip NoTrigger)
    h_num_by_trig = {}
    for trig_name in TRIGGERS:
        if trig_name == "NoTrigger":
            continue
        if trig_name not in trig_sel:
            continue
        ptrs = trig_sel[trig_name]["sig_mHH_ptrs"]
        if not ptrs:
            continue
        h_num = ptrs[0].GetPtr().Clone(f"h_mHH_{trig_name}")
        for p in ptrs[1:]:
            h_num.Add(p.GetPtr())
        h_num_by_trig[trig_name] = h_num

    # Exclusive triggers
    for trig_name, ptrs in result.excl_sig_mHH_ptrs.items():
        if not ptrs:
            continue
        h_num = ptrs[0].GetPtr().Clone(f"h_mHH_{trig_name}")
        for p in ptrs[1:]:
            h_num.Add(p.GetPtr())
        h_num_by_trig[trig_name] = h_num

    if not h_num_by_trig:
        print("ERROR: No trigger numerator histograms found")
        exit(1)

    fig, ax = plot_trigger_efficiency_overlay(
        h_denom, h_num_by_trig,
        xlabel=r"Generator-level $m_{HH}$ [GeV]",
        lumi=result.lumi,
        title=r"$HH \to bb\tau\tau$")

    outpath = os.path.join(ctx.plot_dir, "trigger_efficiency_mHH.png")
    os.makedirs(os.path.dirname(outpath), exist_ok=True)
    fig.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Saved: {outpath}")
