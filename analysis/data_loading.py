"""MC and data loading into RDataFrames with weights and column definitions."""

import os
import yaml

from utils.data import limit_files, load_scouting_data
from utils.skim import load_slim_or_eos
from analysis.definitions import define_gen_columns, define_kinematics, define_gen_matched_ak4, define_coi_matching


def _extract_named_cuts(section):
    """Extract cuts from a YAML section (dict or list).

    Named format (dict): {name: expr, ...} → [(name, expr), ...]
    Old format (list): [expr, ...] → [(None, expr), ...]
    None/empty → []
    """
    if section is None:
        return []
    if isinstance(section, dict):
        return [(name, expr) for name, expr in section.items()]
    if isinstance(section, list):
        return [(None, expr) for expr in section]
    return []


def resolve_cuts(cuts):
    """Resolve CUTS config to a list of named cut tuples.

    Parameters
    ----------
    cuts : list, str, or None
        [] → no cuts. ["expr", ...] → explicit. "tauhtauh" → card from cuts.yaml.
        None → backwards compat.

    Returns
    -------
    (list[tuple], str) — ([(name_or_None, expr), ...], source description)
    """
    _ana_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    _cuts_cfg = os.path.join(_ana_dir, "config", "cuts.yaml")

    if isinstance(cuts, list):
        return [(None, c) for c in cuts], "script"

    if isinstance(cuts, str):
        if not os.path.exists(_cuts_cfg):
            print(f"WARNING: cuts.yaml not found, no cuts applied")
            return [], f"cuts.yaml [{cuts}] (missing)"
        with open(_cuts_cfg) as f:
            data = yaml.safe_load(f) or {}
        common = _extract_named_cuts(data.get("common"))
        if cuts == "common":
            return common, "cuts.yaml [common]"
        channel = _extract_named_cuts(data.get(cuts))
        return common + channel, f"cuts.yaml [{cuts}]"

    # None — backwards compat
    if os.path.exists(_cuts_cfg):
        with open(_cuts_cfg) as f:
            data = yaml.safe_load(f) or {}
        return _extract_named_cuts(data.get("cuts")), "config/cuts.yaml"
    return [], "config/cuts.yaml"


def load_mc_samples(group_files_by_sample, xsec, max_events, args, rdf_exprs=None, cuts=None):
    """Build weighted MC RDataFrames and define analysis columns.

    Parameters
    ----------
    group_files_by_sample : dict[str, dict[str, list[str]]]
        From ls_nanoaod_files_groups: {group: {sample: [files]}}.
    xsec : dict[str, float]
        Cross-sections in pb per sample name.
    max_events : int
        MAX_EVENTS cap used to derive file count (0 = no cap).
    args : argparse.Namespace
        Must have: max_mc_files, all_events.
    rdf_exprs : list or None
        If given, MC weight expressions are appended (for skim branch detection).

    Returns
    -------
    mc : dict[str, RDataFrame]
        Per-sample RDataFrames with "w" weight column, kinematics, and gen columns.
    mc_files_map : dict[str, list[str]]
        Actual files used per sample.
    dy_samples, tt_samples, sig_samples, qcd_samples : list[str]
        Sample name lists by physics group.
    """
    import ROOT

    mc_dfs = {}
    sum_genw_ptrs = {}
    mc_files_map = {}

    for group_name, samples_dict in group_files_by_sample.items():
        print(f"\n[{group_name}]")
        for sample_name, files in samples_dict.items():
            good_files = limit_files(
                files, sample_name,
                max_files=args.max_mc_files,
                all_events=args.all_events,
                max_events=max_events,
            )
            mc_files_map[sample_name] = good_files
            df = load_slim_or_eos(sample_name, good_files)
            sum_genw_ptrs[sample_name] = df.Sum("genWeight")
            mc_dfs[sample_name] = df

    print("\nRunning sum(genWeight) for all samples in parallel...")
    ROOT.RDF.RunGraphs(list(sum_genw_ptrs.values()))

    mc = {}
    w_exprs = {}
    for sample_name, df in mc_dfs.items():
        sum_genw = sum_genw_ptrs[sample_name].GetValue()
        xsec_pb = xsec[sample_name]
        scale_nolumi = xsec_pb * 1000.0 / sum_genw
        print(f"  {sample_name[:60]:60s}  xsec={xsec_pb:.4g} pb  "
              f"sum_genw={sum_genw:.4g}  scale(no L)={scale_nolumi:.4e}")
        w_expr = f"genWeight * {scale_nolumi:.10e}"
        w_exprs[sample_name] = w_expr
        if rdf_exprs is not None:
            rdf_exprs.append(w_expr)
        # Define weight column (use Redefine if already exists in slim)
        existing_cols = set(str(c) for c in df.GetColumnNames())
        if "w" in existing_cols:
            mc[sample_name] = df.Redefine("w", w_expr)
        else:
            mc[sample_name] = df.Define("w", w_expr)

    dy_samples  = list(group_files_by_sample.get("DY",     {}).keys())
    tt_samples  = list(group_files_by_sample.get("TT",     {}).keys())
    sig_samples = list(group_files_by_sample.get("HHbbtt", {}).keys())
    qcd_samples = list(group_files_by_sample.get("QCD",    {}).keys())

    # Determine which samples loaded from slim (have derived columns)
    from utils.skim import has_complete_slim_sample, ensure_slim, _slim_path
    slim_samples = set()
    needs_define = set()
    for name in mc:
        n_files = len(mc_files_map.get(name, []))
        slim = _slim_path(name, n_files)
        if os.path.exists(slim) and has_complete_slim_sample(slim):
            slim_samples.add(name)
        else:
            needs_define.add(name)

    if slim_samples:
        print(f"\n[slim] {len(slim_samples)} sample(s) loaded with pre-computed columns")
    if needs_define:
        print(f"[slim] {len(needs_define)} sample(s) need Define chains")

    # Run Define chains only on samples that need them
    if needs_define:
        mc_need = {n: mc[n] for n in needs_define}
        sig_need = [s for s in sig_samples if s in needs_define]
        define_gen_columns(mc_need, sig_need)
        for name in mc_need:
            mc_need[name] = define_kinematics(mc_need[name])
        define_gen_matched_ak4(mc_need, sig_need)
        define_coi_matching(mc_need, sig_need)
        mc.update(mc_need)

        # Snapshot everything to slim
        mc_to_skim = {n: mc[n] for n in needs_define}
        ensure_slim(mc_to_skim, {n: mc_files_map[n] for n in needs_define})

        # Re-define "w" on reloaded slim RDataFrames
        for name in needs_define:
            cols = set(str(c) for c in mc[name].GetColumnNames())
            if "w" not in cols:
                mc[name] = mc[name].Define("w", w_exprs[name])

    # Apply cuts
    named_cuts, _src = resolve_cuts(cuts)
    if named_cuts:
        print(f"\n[cuts] Applying {len(named_cuts)} cuts from {_src}:")
        for cut_name, cut_expr in named_cuts:
            if cut_name:
                print(f"  \u2192 {cut_name}: {cut_expr}")
            else:
                print(f"  \u2192 {cut_expr}")
            for name in mc:
                mc[name] = mc[name].Filter(cut_expr)
    elif isinstance(cuts, (list, str)):
        print(f"\n[cuts] No cuts applied ({_src})")

    return mc, mc_files_map, dy_samples, tt_samples, sig_samples, qcd_samples


def load_data(args, brilcalc_default=None, data_years=None, data_runs=None):
    """Load scouting data into an RDataFrame with kinematics columns defined.

    Parameters
    ----------
    args : argparse.Namespace
        Must have: no_data, all_events, max_data_files.
    brilcalc_default : dict or None
        Parsed brilcalc data (used to estimate lumi when --no-data).
    data_years, data_runs : list[str] or None
        Forwarded to load_scouting_data.

    Returns
    -------
    data_df : RDataFrame or None
        None if args.no_data. Has kinematics columns defined.
    lumi : float or None
        Estimated lumi from brilcalc if --no-data; None otherwise (computed later).
    """
    import ROOT

    if args.no_data:
        if brilcalc_default:
            lumi = sum(v["lumi"] for v in brilcalc_default.values())
        else:
            lumi = 103.965  # fallback fb^-1
        print(f"\n[--no-data] Skipping data; brilcalc lumi = {lumi:.3f} fb^-1")
        return None, lumi

    data_files_all, _ = load_scouting_data(years=data_years, runs=data_runs)
    if not args.all_events and args.max_data_files > 0:
        data_files = data_files_all[:args.max_data_files]
        print(f"  Limited data to {len(data_files)}/{len(data_files_all)} files")
    else:
        data_files = data_files_all

    _prev_err = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    data_df = ROOT.RDataFrame("Events", data_files)
    ROOT.gErrorIgnoreLevel = _prev_err
    print(f"  Loaded {len(data_files)} data files into RDataFrame")

    data_df = define_kinematics(data_df)
    return data_df, None
