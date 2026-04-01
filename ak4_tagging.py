#!/usr/bin/env python3
"""Gen-matched AK4 UParT tagger score studies."""

from analysis.variables import PlotVar

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

def _gen_particle_vars(tag, label):
    """PlotVars for one gen-matched AK4 particle."""
    return [
        PlotVar(f"ak4_gen{tag}_dR",          rf"$\Delta R$(gen ${label}$, AK4)",                       25, 0, 0.5),
        PlotVar(f"ak4_gen{tag}_BvsAll",       rf"Gen-matched ${label}$ AK4 UParT BvsAll",               25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_b",        rf"Gen-matched ${label}$ AK4 raw prob\_b",                25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_c",        rf"Gen-matched ${label}$ AK4 raw prob\_c",                25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_g",        rf"Gen-matched ${label}$ AK4 raw prob\_g",                25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_uds",      rf"Gen-matched ${label}$ AK4 raw prob\_uds",              25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_lepb",     rf"Gen-matched ${label}$ AK4 raw prob\_lepb",             25, 0, 1,   sentinel=True),
    ]

def _gen_tau_vars(tag, label):
    """PlotVars for one gen-matched tau AK4 (adds tau score columns)."""
    return _gen_particle_vars(tag, label) + [
        PlotVar(f"ak4_gen{tag}_TauVsAll",     rf"Gen-matched ${label}$ AK4 UParT TauVsAll",            25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_TaumVsAll",    rf"Gen-matched ${label}$ AK4 UParT $\tau_\mu$ VsAll",    25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_TaupVsAll",    rf"Gen-matched ${label}$ AK4 UParT $\tau_h$ VsAll",      25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_taup",     rf"Gen-matched ${label}$ AK4 raw prob $\tau_h^+$",       25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_taum",     rf"Gen-matched ${label}$ AK4 raw prob $\tau_h^-$",       25, 0, 1,   sentinel=True),
        PlotVar(f"ak4_gen{tag}_raw_tau",      rf"Gen-matched ${label}$ AK4 raw prob $\tau$ (p+m)",     25, 0, 1,   sentinel=True),
    ]

PLOT_VARS = (
    _gen_particle_vars("b1", r"b_1")
    + _gen_particle_vars("b2", r"b_2")
    + _gen_tau_vars("tau1", r"\tau_1")
    + _gen_tau_vars("tau2", r"\tau_2")
)

# Tau charge misidentification
PLOT_VARS += [
    PlotVar("gen_taum_TaumVsAll", r"Gen $\tau^-$ AK4 $\tau_h^-$ vs All (correct sign)",  25, 0, 1),
    PlotVar("gen_taum_TaupVsAll", r"Gen $\tau^-$ AK4 $\tau_h^+$ vs All (wrong sign)",    25, 0, 1),
    PlotVar("gen_taum_TauVsAll",  r"Gen $\tau^-$ AK4 $\tau_h$ vs All (combined)",         25, 0, 1),
    PlotVar("gen_taup_TaupVsAll", r"Gen $\tau^+$ AK4 $\tau_h^+$ vs All (correct sign)",  25, 0, 1),
    PlotVar("gen_taup_TaumVsAll", r"Gen $\tau^+$ AK4 $\tau_h^-$ vs All (wrong sign)",    25, 0, 1),
    PlotVar("gen_taup_TauVsAll",  r"Gen $\tau^+$ AK4 $\tau_h$ vs All (combined)",         25, 0, 1),
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
