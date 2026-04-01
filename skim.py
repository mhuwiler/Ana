#!/usr/bin/env python3
"""Write slim ROOT files to /depot/ — keeps all events, drops unused branches.

This is a data-reduction step, separate from analysis. Run once after
changing the column definitions, then analysis scripts read from slim files.
"""

# ── Config ─────────────────────────────────────────────────────────────────────
SAMPLES      = ["HHbbtt", "DY", "TT", "QCD"]
NTHREADS     = 4
MAX_MC_FILES = 1    # 0 = all files
MAX_EVENTS   = None

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    from analysis.runner import setup
    from utils.skim import rdf_exprs, collect_all_expressions, run_skim_pipeline
    from utils.data import ls_nanoaod_files_groups, BASE, GROUPS, XSEC, MAX_EVENTS
    from analysis.data_loading import load_mc_samples
    from types import SimpleNamespace

    ctx = setup(__file__, samples=SAMPLES, data=False, nthreads=NTHREADS,
                max_mc_files=MAX_MC_FILES, max_events=MAX_EVENTS)

    # Load MC with all columns defined
    filtered_groups = {k: v for k, v in GROUPS.items() if k in SAMPLES}
    _, group_files_by_sample = ls_nanoaod_files_groups(
        base_dir=BASE, groups=filtered_groups, verbose=0)
    args = SimpleNamespace(
        max_mc_files=MAX_MC_FILES, all_events=False,
        no_data=True, max_data_files=0, skip_cutflow=True)
    mc, mc_files_map, *_ = load_mc_samples(
        group_files_by_sample, XSEC, MAX_EVENTS, args, rdf_exprs=rdf_exprs)

    # Collect all expressions from library code (triggers, AK8, GenPart, etc.)
    exprs = collect_all_expressions()
    run_skim_pipeline(mc, mc_files_map, exprs, args)
