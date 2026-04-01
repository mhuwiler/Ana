#!/usr/bin/env python3
"""AK8 fat jet tagger scores, candidate kinematics, and gen-matched AK8 studies."""

from analysis.variables import PlotVar, expand_ak8_vars, expand_ak8_candidate_vars

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
DATA         = False
NTHREADS     = 4
NPLOT_WORKERS = 8
THEME        = "light"
OVERWRITE    = False
MAX_MC_FILES = 1
MAX_EVENTS   = None

# ── Variables to plot ──────────────────────────────────────────────────────────

# Per-jet AK8 variables (leading / sub-leading / third, pT > 150 GeV)
PLOT_VARS = expand_ak8_vars(
    count=3,
    labels=["Leading", "Sub-leading", "Third"],
    pt_maxes=[1000, 800, 600],
    vars_per_jet=[
        ("pt",       r"$p_T$ [GeV]",                         50, 150, "{pt_max}"),
        ("eta",      r"$\eta$",                               30,  -5,         5),
        ("mass",     r"mass [GeV]",                           40,   0,       400),
        ("msd",      r"$m_{SD}$ [GeV]",                      40,   0,       300),
        ("Xbb",      r"$X_{bb}$",                             25,   0,         1),
        ("Xtt",      r"$X_{\tau_h\tau_h}$",                   25,   0,         1),
        ("Xtm",      r"$X_{\tau_\mu\tau_h}$",                 25,   0,         1),
        ("Xte",      r"$X_{\tau_e\tau_h}$",                   25,   0,         1),
        ("QCD",      r"QCD",                                  25,   0,         1),
        ("XbbVsQCD", r"$X_{bb}$ vs QCD",                     25,   0,         1),
        ("XttVsQCD", r"$X_{\tau_h\tau_h}$ vs QCD",           25,   0,         1),
        ("XtmVsQCD", r"$X_{\tau_\mu\tau_h}$ vs QCD",         25,   0,         1),
        ("XteVsQCD", r"$X_{\tau_e\tau_h}$ vs QCD",           25,   0,         1),
        ("massCorr", r"regressed mass [GeV]",                 40,   0,       400),
        ("massRes",  r"resonance mass [GeV]",                 40,   0,       400),
    ],
)

# AK8 Higgs candidate variables
PLOT_VARS += expand_ak8_candidate_vars(
    candidates=[
        ("Hbb", r"$H \to bb$"),
        ("Htt", r"$H \to \tau_h\tau_h$"),
        ("Htm", r"$H \to \tau_\mu\tau_h$"),
        ("Hte", r"$H \to \tau_e\tau_h$"),
    ],
    candidate_vars=[
        ("pt",    r"$p_T$ [GeV]",     50, 150, 1000),
        ("eta",   r"$\eta$",           30,  -5,    5),
        ("mass",  r"mass [GeV]",       40,   0,  400),
        ("msd",   r"$m_{SD}$ [GeV]",  40,   0,  300),
        ("score", r"score",            25,   0,    1),
    ],
)

# Standalone AK8 + gen-matched AK8 tagger scores
PLOT_VARS += [
    PlotVar("ak8_dR_Hbb_Htt",           r"$\Delta R(H_{bb}, H_{\tau\tau})$ cand AK8",                30, 0, 6,   sentinel=True),
    PlotVar("ak8_genHbb_match_dR",       r"$\Delta R$(gen $H_{bb}$, matched AK8)",                    25, 0, 1,   sentinel=True),
    PlotVar("ak8_genHbb_match_Xbb",      r"Gen-matched $H_{bb}$ AK8 raw $X_{bb}$",                   25, 0, 1,   sentinel=True),
    PlotVar("ak8_genHbb_match_XbbVsAll", r"Gen-matched $H_{bb}$ AK8 $X_{bb}$ vs All",                25, 0, 1,   sentinel=True),
    PlotVar("ak8_genHtt_match_dR",       r"$\Delta R$(gen $H_{\tau\tau}$, matched AK8)",              25, 0, 1,   sentinel=True),
    PlotVar("ak8_genHtt_match_Xtt",      r"Gen-matched $H_{\tau\tau}$ AK8 raw $X_{\tau_h\tau_h}$",   25, 0, 1,   sentinel=True),
    PlotVar("ak8_genHtt_match_XttVsAll", r"Gen-matched $H_{\tau\tau}$ AK8 $X_{\tau_h\tau_h}$ vs All", 25, 0, 1,  sentinel=True),
]

# Gen-level Higgs kinematics (signal only — backgrounds produce empty histograms)
PLOT_VARS += [
    PlotVar("gen_pt_Hbb",       r"Gen $H \to bb$ $p_T$ [GeV]",          50,    0,   500),
    PlotVar("gen_pt_Htautau",   r"Gen $H \to \tau\tau$ $p_T$ [GeV]",    50,    0,   500),
    PlotVar("gen_eta_Hbb",      r"Gen $H \to bb$ $\eta$",               30,   -5,     5),
    PlotVar("gen_eta_Htautau",  r"Gen $H \to \tau\tau$ $\eta$",         30,   -5,     5),
    PlotVar("gen_dR_H1H2",      r"Gen $\Delta R(H_{bb}, H_{\tau\tau})$", 30,   0,     6),
    PlotVar("gen_mHH",          r"Gen $m_{HH}$ [GeV]",                  50,  200,   800),
]

# 2D histograms: tuple of two (name, label, nbins, min, max) tuples
PLOT_VARS_2D = [
    (("ak8_genHbb_match_dR", r"$\Delta R$(gen $H_{bb}$, AK8)", 25, 0, 1.0),
     ("ak8_genHbb_match_XbbVsAll", r"$X_{bb}$ vs All", 25, 0, 1.0)),
    (("ak8_genHtt_match_dR", r"$\Delta R$(gen $H_{\tau\tau}$, AK8)", 25, 0, 1.0),
     ("ak8_genHtt_match_XttVsAll", r"$X_{\tau\tau}$ vs All", 25, 0, 1.0)),
    (("ak8_genHbb_match_XbbVsAll", r"$X_{bb}$ vs All ($H\to bb$ cand)", 25, 0, 1.0),
     ("ak8_genHtt_match_XttVsAll", r"$X_{\tau\tau}$ vs All ($H\to\tau\tau$ cand)", 25, 0, 1.0)),
    (("gen_pt_Hbb", r"Gen $H \to bb$ $p_T$ [GeV]", 25, 150, 600),
     ("ak8_genHbb_match_XbbVsAll", r"$X_{bb}$ vs All", 25, 0, 1.0)),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, data=DATA, nthreads=NTHREADS, nplot_workers=NPLOT_WORKERS,
                theme=THEME, overwrite=OVERWRITE, max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS)
    result = load_and_run(ctx, PLOT_VARS, vars_2d=PLOT_VARS_2D)
    result.plot.stacked(ratio="significance")
    result.plot.shapes()
    result.plot.trigger_overlays()
    result.plot.hist2d()
    if DATA:
        result.plot.data(ratio="data_mc")
