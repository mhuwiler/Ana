#!/usr/bin/env python3
"""Gen-level b-jet matching and Higgs kinematics — with cuts from cuts.yaml."""

from analysis.variables import PlotVar

# ── User config ────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt"]
TRIGGERS     = ["NoTrigger"]
DATA         = False
NTHREADS     = 4
NPLOT_WORKERS = 8
THEME        = "light"
OVERWRITE    = False
MAX_MC_FILES = 1
MAX_EVENTS   = None
# CUTS not set — uses cuts.yaml defaults

# ── Variables to plot ──────────────────────────────────────────────────────────
PLOT_VARS = [
    # Gen-matched b-jet scores
    PlotVar("ak4_genb1_BvsAll",  r"Gen-matched $b_1$ UParT BvsAll",  50, 0, 1),
    PlotVar("ak4_genb2_BvsAll",  r"Gen-matched $b_2$ UParT BvsAll",  50, 0, 1),
    PlotVar("ak4_genb1_dR",      r"$\Delta R$(gen $b_1$, AK4)",      50, 0, 0.5),
    PlotVar("ak4_genb2_dR",      r"$\Delta R$(gen $b_2$, AK4)",      50, 0, 0.5),
    # Gen-matched tau scores (charge-identified via PDG ID)
    PlotVar("gen_taum_TaumVsAll", r"Gen $\tau^-$ AK4 $\tau_h^-$ vs All (correct)",  50, 0, 1),
    PlotVar("gen_taum_TaupVsAll", r"Gen $\tau^-$ AK4 $\tau_h^+$ vs All (wrong)",    50, 0, 1),
    PlotVar("gen_taum_TauVsAll",  r"Gen $\tau^-$ AK4 $\tau_h$ vs All (combined)",   50, 0, 1),
    PlotVar("gen_taup_TaupVsAll", r"Gen $\tau^+$ AK4 $\tau_h^+$ vs All (correct)",  50, 0, 1),
    PlotVar("gen_taup_TaumVsAll", r"Gen $\tau^+$ AK4 $\tau_h^-$ vs All (wrong)",    50, 0, 1),
    PlotVar("gen_taup_TauVsAll",  r"Gen $\tau^+$ AK4 $\tau_h$ vs All (combined)",   50, 0, 1),
    # Reco b-tagged jet scores (sorted by BvsAll, not pT)
    PlotVar("b0_score",  r"Highest UParT BvsAll jet",      50, 0, 1),
    PlotVar("b1_score",  r"2nd highest UParT BvsAll jet",   50, 0, 1),
    # Gen-level Higgs kinematics
    PlotVar("gen_pt_Hbb",       r"Gen $H \to bb$ $p_T$ [GeV]",          50, 0, 500),
    PlotVar("gen_pt_Htautau",   r"Gen $H \to \tau\tau$ $p_T$ [GeV]",    50, 0, 500),
    PlotVar("gen_eta_Hbb",      r"Gen $H \to bb$ $\eta$",               30, -5, 5),
    PlotVar("gen_eta_Htautau",  r"Gen $H \to \tau\tau$ $\eta$",         30, -5, 5),
    PlotVar("gen_mass_Hbb",     r"Gen $H \to bb$ mass [GeV]",           50, 100, 150),
    PlotVar("gen_mass_Htautau", r"Gen $H \to \tau\tau$ mass [GeV]",     50, 0, 150),
    PlotVar("gen_dR_H1H2",     r"Gen $\Delta R(H_{bb}, H_{\tau\tau})$", 30, 0, 6),
    PlotVar("gen_pt_HH",       r"Gen $p_T^{HH}$ [GeV]",                50, 0, 500),
    PlotVar("gen_mHH",         r"Gen $m_{HH}$ [GeV]",                   50, 200, 800),
]

# Tau charge-ID confusion matrix (2D heatmaps)
PLOT_VARS_2D = [
    # Gen τ⁻: all pairwise
    (("gen_taum_TaumVsAll", r"Gen $\tau^-$: $\tau_h^-$ vs All", 25, 0, 1),
     ("gen_taum_TaupVsAll", r"Gen $\tau^-$: $\tau_h^+$ vs All", 25, 0, 1)),
    (("gen_taum_TaumVsAll", r"Gen $\tau^-$: $\tau_h^-$ vs All", 25, 0, 1),
     ("gen_taum_TauVsAll",  r"Gen $\tau^-$: $\tau_h$ vs All",   25, 0, 1)),
    (("gen_taum_TaupVsAll", r"Gen $\tau^-$: $\tau_h^+$ vs All", 25, 0, 1),
     ("gen_taum_TauVsAll",  r"Gen $\tau^-$: $\tau_h$ vs All",   25, 0, 1)),
    # Gen τ⁺: all pairwise
    (("gen_taup_TaupVsAll", r"Gen $\tau^+$: $\tau_h^+$ vs All", 25, 0, 1),
     ("gen_taup_TaumVsAll", r"Gen $\tau^+$: $\tau_h^-$ vs All", 25, 0, 1)),
    (("gen_taup_TaupVsAll", r"Gen $\tau^+$: $\tau_h^+$ vs All", 25, 0, 1),
     ("gen_taup_TauVsAll",  r"Gen $\tau^+$: $\tau_h$ vs All",   25, 0, 1)),
    (("gen_taup_TaumVsAll", r"Gen $\tau^+$: $\tau_h^-$ vs All", 25, 0, 1),
     ("gen_taup_TauVsAll",  r"Gen $\tau^+$: $\tau_h$ vs All",   25, 0, 1)),
]

# Confusion matrix: tau charge-ID
PLOT_CM_2D = [
    ("tau_charge",
     ("gen_taum_pred_correct", r"$\tau^-$", 2, -0.5, 1.5),
     ("gen_taup_pred_correct", r"$\tau^+$", 2, -0.5, 1.5)),
]

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup, load_and_run
    ctx = setup(__file__, samples=SAMPLES, triggers=TRIGGERS, data=DATA,
                nthreads=NTHREADS, nplot_workers=NPLOT_WORKERS, theme=THEME, overwrite=OVERWRITE,
                max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS)
    result = load_and_run(ctx, PLOT_VARS, vars_2d=PLOT_VARS_2D, cm_2d=PLOT_CM_2D)
    result.plot.stacked()
    result.plot.shapes()
    result.plot.hist2d()
    result.plot.confusion_matrix()
