#!/usr/bin/env python3
"""All MC samples — no cuts, no trigger. Baseline distributions."""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
TRIGGERS     = ["NoTrigger", "DST_JetHT", "DST_Muon", "DST_Electron"]
DATA         = False
NTHREADS     = 16
NPLOT_WORKERS = 8
THEME        = "light"
OVERWRITE    = True
MAX_MC_FILES = 1
MAX_EVENTS   = None
CUTS         = []

# ── Variables to plot ──────────────────────────────────────────────────────────
PLOT_VARS = [
    # AK4 jet pT
    PlotVar("ak4_pt0",   r"Leading jet $p_T$ [GeV]",        50,    0,  250),
    PlotVar("ak4_pt1",   r"Sub-leading jet $p_T$ [GeV]",    50,    0,  250),
    PlotVar("ak4_pt2",   r"3rd jet $p_T$ [GeV]",            50,    0,  250),
    PlotVar("ak4_pt3",   r"4th jet $p_T$ [GeV]",            50,    0,  250),
    # AK4 jet eta/mass
    PlotVar("ak4_eta0",  r"Leading jet $\eta$",              30,   -5,     5),
    PlotVar("ak4_eta1",  r"Sub-leading jet $\eta$",          30,   -5,     5),
    PlotVar("ak4_mass0", r"Leading jet mass [GeV]",          30,    0,   100),
    # Multiplicity
    PlotVar("nJets",      r"Number of AK4 jets",             16,    0,    16),
    PlotVar("nLeptons",   r"Number of leptons ($\mu + e$)",  10,    0,    10),
    PlotVar("nFatJets",   r"Number of AK8 fat jets",          8,    0,     8),
    # Global event
    PlotVar("HT",         r"$H_T$ [GeV]",                  200,    0,  2000),
    PlotVar("MHT",        r"$\slash{H}_T$ [GeV]",           30,    0,  1500),
    PlotVar("centrality", r"Centrality ($H_T / \sum E$)",   25,    0,     1),
    # Dijet
    PlotVar("mjj_01",    r"$m_{jj}$ (leading dijet) [GeV]", 40,    0,  2000),
    PlotVar("dR_01",     r"$\Delta R(j_0, j_1)$",           30,    0,     6),
    PlotVar("m4j",       r"$m_{4j}$ (leading 4 jets) [GeV]", 32,   0,   800),
    # B COI (jets with highest BvsAll score)
    PlotVar("b_coi0_score", r"Highest BvsAll jet score",      25,    0,     1),
    PlotVar("b_coi1_score", r"2nd highest BvsAll jet score",  25,    0,     1),
    PlotVar("b_coi0_pt",    r"Highest BvsAll jet $p_T$ [GeV]", 36,  0,  1000),
    PlotVar("mbb_coi",      r"$m_{bb}$ COI [GeV]",           30,    0,   300),
    PlotVar("dR_bb_coi",    r"$\Delta R(b_{COI0}, b_{COI1})$", 30,  0,     6),
    # Tau COI (jets with highest TauVsAll score)
    PlotVar("tau_coi0_score", r"Highest TauVsAll jet score",     50, 0, 1),
    PlotVar("tau_coi1_score", r"2nd highest TauVsAll jet score", 50, 0, 1),
    PlotVar("tau_coi0_pt",    r"Highest TauVsAll jet $p_T$ [GeV]", 50, 0, 500),
    PlotVar("tau_coi0_taup",  r"Highest TauVsAll jet $\tau_h^+$ score", 50, 0, 1),
    PlotVar("tau_coi0_taum",  r"Highest TauVsAll jet $\tau_h^-$ score", 50, 0, 1),
    PlotVar("mtautau_coi",    r"$m_{\tau\tau}$ COI [GeV]",   30,    0,   300),
    PlotVar("dR_tautau_coi",  r"$\Delta R(\tau_{COI0}, \tau_{COI1})$", 30, 0, 6),
    # AK8 leading jet
    # PlotVar("ak8_pt0",       r"Leading AK8 $p_T$ [GeV]",       50, 150, 1000),
    # PlotVar("ak8_mass0",     r"Leading AK8 mass [GeV]",        40,   0,  400),
    # PlotVar("ak8_msd0",      r"Leading AK8 $m_{SD}$ [GeV]",    40,   0,  300),
    # PlotVar("ak8_XbbVsQCD0", r"Leading AK8 $X_{bb}$ vs QCD",   25,   0,    1),
    # PlotVar("ak8_XttVsQCD0", r"Leading AK8 $X_{\tau\tau}$ vs QCD", 25, 0,  1),
    # # Higgs candidates
    # PlotVar("mbb_cand",     r"$m_{jj}^{H \to bb}$ candidate [GeV]",        25,  50,  200),
    # PlotVar("mtautau_cand", r"$m_{jj}^{H \to \tau\tau}$ candidate [GeV]",  25,   0,  200),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, data=DATA,
                nthreads=NTHREADS, nplot_workers=NPLOT_WORKERS, theme=THEME, overwrite=OVERWRITE,
                max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS, cuts=CUTS)
    result = load_and_run(ctx, PLOT_VARS)
    result.plot.stacked()
    result.plot.shapes()
