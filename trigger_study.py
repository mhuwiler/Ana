#!/usr/bin/env python3
"""Trigger efficiency and exclusive trigger studies.

Produces per-trigger cutflow tables, stacked MC plots, and trigger overlay
plots (inclusive + exclusive) for all scouting triggers.
"""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES       = ["HHbbtt", "DY", "TT"]
DATA          = False
NTHREADS      = 32
NPLOT_WORKERS = 32
OVERWRITE     = True
MAX_MC_FILES  = 22

# Per-trigger cuts from cuts.yaml
CUTS = [
    ("NoTrigger",     ["common", "tauhtauh"]),
    ("DST_JetHT",     ["common", "tauhtauh"]),
    ("PARKING_HH",    ["common", "tauhtauh"]),

    # ("DST_Muon",      ["common", "taumutauh"]),
    # ("DST_Electron",  ["common", "tauetauh"]),
]

# ── Variables to plot ──────────────────────────────────────────────────────────
PLOT_VARS = [
    # Tagger scores
    # PlotVar("b_coi0_score",    r"$b_{COI0}$ BvsAll",              50, 0, 1),
    # PlotVar("b_coi1_score",    r"$b_{COI1}$ BvsAll",              50, 0, 1),
    # PlotVar("tau_coi0_score",  r"$\tau_{COI0}$ TauVsAll",         50, 0, 1),
    # PlotVar("tau_coi1_score",  r"$\tau_{COI1}$ TauVsAll",         50, 0, 1),

    # # Di-candidate masses
    # PlotVar("mbb_coi",        r"$m_{bb}$ COI [GeV]",              30, 0, 300),
    # PlotVar("mtautau_coi",    r"$m_{\tau\tau}$ COI [GeV]",        30, 0, 300),

    # # Event-level kinematics
    # PlotVar("HT",             r"$H_T$ [GeV]",                     50, 0, 1000),
    # PlotVar("MET",            r"Scouting MET [GeV]",               50, 0, 200),
    # PlotVar("nJets",          r"Number of AK4 jets",               16, 0, 16),

    # # Angular / QCD rejection
    # PlotVar("dphi_bb_tautau", r"$|\Delta\phi(bb, \tau\tau)|$",     30, 0, 3.15),
    # PlotVar("dR_bb_tautau",   r"$\Delta R(bb, \tau\tau)$",        30, 0, 6),
    # PlotVar("D_zeta",         r"$D_\zeta$ [GeV]",                  50, -200, 100),
    # PlotVar("MT_tau0_MET",    r"$M_T(\tau_0, \mathrm{MET})$ [GeV]", 40, 0, 200),
    # PlotVar("dphi_MET_tau0",  r"$|\Delta\phi(\mathrm{MET}, \tau_0)|$", 30, 0, 3.15),

    # Object pT
    PlotVar("b_coi0_pt",      r"$b_{COI0}$ $p_T$ [GeV]",          50, 0, 500),
    PlotVar("b_coi1_pt",      r"$b_{COI1}$ $p_T$ [GeV]",          50, 0, 500),
    PlotVar("tau_coi0_pt",    r"$\tau_{COI0}$ $p_T$ [GeV]",       50, 0, 500),
    PlotVar("tau_coi1_pt",    r"$\tau_{COI1}$ $p_T$ [GeV]",       50, 0, 500),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, cuts=CUTS, data=DATA,
                nthreads=NTHREADS, nplot_workers=NPLOT_WORKERS, overwrite=OVERWRITE,
                max_mc_files=MAX_MC_FILES)
    result = load_and_run(ctx, PLOT_VARS)
    result.plot.stacked(ratio="significance")
    result.plot.trigger_overlays()
