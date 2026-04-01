#!/usr/bin/env python3
"""AK4 jet kinematics, event-level variables, and b-tagging plots."""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
DATA         = False
NTHREADS     = 4
NPLOT_WORKERS = 8
THEME        = "light"
OVERWRITE    = False
MAX_MC_FILES = 1    # 0 = all files
MAX_EVENTS   = None # None = default (100000 from samples.yaml)

# ── Variables to plot ──────────────────────────────────────────────────────────
PLOT_VARS = [
    # Jet pT
    PlotVar("ak4_pt0",   r"Leading jet $p_T$ [GeV]",        50,    0,  250),
    PlotVar("ak4_pt1",   r"Sub-leading jet $p_T$ [GeV]",    50,    0,  250),
    PlotVar("ak4_pt2",   r"3rd jet $p_T$ [GeV]",            50,    0,  250),
    PlotVar("ak4_pt3",   r"4th jet $p_T$ [GeV]",            50,    0,  250),
    # Jet eta
    PlotVar("ak4_eta0",  r"Leading jet $\eta$",              30,   -5,     5),
    PlotVar("ak4_eta1",  r"Sub-leading jet $\eta$",          30,   -5,     5),
    PlotVar("ak4_eta2",  r"3rd jet $\eta$",                  30,   -5,     5),
    PlotVar("ak4_eta3",  r"4th jet $\eta$",                  30,   -5,     5),
    # Jet mass
    PlotVar("ak4_mass0", r"Leading jet mass [GeV]",          30,    0,   100),
    PlotVar("ak4_mass1", r"Sub-leading jet mass [GeV]",      30,    0,   100),
    # Multiplicity
    PlotVar("nJets",      r"Number of AK4 jets",             16,    0,    16),
    PlotVar("nLeptons",   r"Number of leptons ($\mu + e$)",  10,    0,    10),
    PlotVar("nMuons",     r"Number of muons",                 6,    0,     6),
    PlotVar("nElectrons", r"Number of electrons",             6,    0,     6),
    PlotVar("nFatJets",   r"Number of AK8 fat jets",          8,    0,     8),
    # Global event
    PlotVar("HT",         r"$H_T$ [GeV]",                  200,    0,  2000),
    PlotVar("MHT",        r"$\slash{H}_T$ [GeV]",           30,    0,  1500),
    PlotVar("centrality", r"Centrality ($H_T / \sum E$)",   25,    0,     1),
    # Dijet (leading pair)
    PlotVar("mjj_01",    r"$m_{jj}$ (leading dijet) [GeV]", 40,    0,  2000),
    PlotVar("dR_01",     r"$\Delta R(j_0, j_1)$",           30,    0,     6),
    PlotVar("dEta_01",   r"$|\Delta\eta(j_0, j_1)|$",       25,    0,     5),
    # 4-jet
    PlotVar("m4j",       r"$m_{4j}$ (leading 4 jets) [GeV]", 32,   0,   800),
    # b-tagged jets
    PlotVar("b0_pt",     r"$b_0$ jet $p_T$ [GeV]",          36,    0,  1000),
    PlotVar("b1_pt",     r"$b_1$ jet $p_T$ [GeV]",          36,    0,   600),
    PlotVar("b0_score",  r"$b_0$ UParT BvsAll",             25,    0,     1),
    PlotVar("b1_score",  r"$b_1$ UParT BvsAll",             25,    0,     1),
    PlotVar("b0_raw",    r"$b_0$ UParT raw prob\_b",        25,    0,     1),
    PlotVar("b1_raw",    r"$b_1$ UParT raw prob\_b",        25,    0,     1),
    PlotVar("mbb",       r"$m_{bb}$ [GeV]",                 30,    0,   300),
    PlotVar("dR_bb",     r"$\Delta R(b_0, b_1)$",           30,    0,     6),
    PlotVar("ptbb",      r"$p_T^{bb}$ [GeV]",               30,    0,  1000),
    # Higgs candidate dijet masses
    PlotVar("mbb_cand",     r"$m_{jj}^{H \to bb}$ candidate [GeV]",        25,  50,  200),
    PlotVar("mtautau_cand", r"$m_{jj}^{H \to \tau\tau}$ candidate [GeV]",  25,   0,  200),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, data=DATA, nthreads=NTHREADS, nplot_workers=NPLOT_WORKERS,
                theme=THEME, overwrite=OVERWRITE, max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS)
    result = load_and_run(ctx, PLOT_VARS)
    result.plot.stacked(ratio="significance")
    result.plot.shapes()
    result.plot.trigger_overlays()
    if DATA:
        result.plot.data(ratio="data_mc")
