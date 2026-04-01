#!/usr/bin/env python3
"""COI tagging study — b and tau candidate scores across all MC samples."""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
TRIGGERS     = ["NoTrigger", "DST_JetHT", "DST_Muon", "DST_Electron"]
DATA         = False
NTHREADS     = 32
NPLOT_WORKERS = 32
THEME        = "light"
OVERWRITE    = True
MAX_MC_FILES = 1
MAX_EVENTS   = None
CUTS         = []

# ── Variables to plot ──────────────────────────────────────────────────────────
PLOT_VARS = [
    # B COI scores
    PlotVar("b_coi0_score",  r"$b_{COI0}$ BvsAll",             50, 0, 1),
    PlotVar("b_coi1_score",  r"$b_{COI1}$ BvsAll",             50, 0, 1),
    PlotVar("b_coi0_pt",     r"$b_{COI0}$ $p_T$ [GeV]",        50, 0, 500),
    PlotVar("b_coi1_pt",     r"$b_{COI1}$ $p_T$ [GeV]",        50, 0, 500),

    # Tau COI scores
    PlotVar("tau_coi0_score", r"$\tau_{COI0}$ TauVsAll",       50, 0, 1),
    PlotVar("tau_coi1_score", r"$\tau_{COI1}$ TauVsAll",       50, 0, 1),
    PlotVar("tau_coi0_taup",  r"$\tau_{COI0}$ $\tau_h^+$ score", 50, 0, 1),
    PlotVar("tau_coi0_taum",  r"$\tau_{COI0}$ $\tau_h^-$ score", 50, 0, 1),
    PlotVar("tau_coi0_pt",    r"$\tau_{COI0}$ $p_T$ [GeV]",    50, 0, 500),
    PlotVar("tau_coi1_pt",    r"$\tau_{COI1}$ $p_T$ [GeV]",    50, 0, 500),

    # Cross-scores (does best b-jet look like tau? vice versa?)
    PlotVar("b_coi0_TauVsAll",  r"Best b-jet: TauVsAll",      50, 0, 1),
    PlotVar("tau_coi0_BvsAll",  r"Best $\tau$-jet: BvsAll",    50, 0, 1),

    # Overlap (same jet picked as both b and tau COI)
    PlotVar("b_tau_overlap_0",   r"$b_{COI0} = \tau_{COI0}$ (same jet)", 2, -0.5, 1.5),
    PlotVar("any_b_tau_overlap", r"Any b-$\tau$ COI overlap",            2, -0.5, 1.5),

    # Di-candidate masses
    PlotVar("mbb_coi",       r"$m_{bb}$ COI [GeV]",           30, 0, 300),
    PlotVar("mtautau_coi",   r"$m_{\tau\tau}$ COI [GeV]",     30, 0, 300),
    PlotVar("dR_bb_coi",     r"$\Delta R(b_{COI0}, b_{COI1})$",       30, 0, 6),
    PlotVar("dR_tautau_coi", r"$\Delta R(\tau_{COI0}, \tau_{COI1})$", 30, 0, 6),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, data=DATA,
                nthreads=NTHREADS, theme=THEME, overwrite=OVERWRITE,
                max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS, cuts=CUTS)
    result = load_and_run(ctx, PLOT_VARS)
    result.plot.stacked()
    result.plot.stacked(ratio="significance")
    result.plot.stacked(ratio="cum_significance")
    # result.plot.stacked(ratio="efficiency")
    result.plot.shapes()
