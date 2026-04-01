"""
Data loading utilities for scouting NanoAOD analysis.

Handles:
  - Local filesystem NanoAOD file discovery (ls_nanoaod_files_*)
  - XCache / DAS-based scouting data loading
  - Sample definitions and constants
"""

import os
import glob
import json
import subprocess

import yaml
import ROOT

# ---------------------------------------------------------------------------
# YAML configuration loader
# ---------------------------------------------------------------------------

_ANA_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_DEFAULT_CONFIG = os.path.join(_ANA_DIR, "config", "samples.yaml")


def load_config(path=None):
    """Load sample configuration from YAML.

    Returns the parsed dict, or None if the file does not exist.
    """
    p = path or _DEFAULT_CONFIG
    if os.path.isfile(p):
        with open(p) as f:
            cfg = yaml.safe_load(f)
        print(f"[config] Loaded sample config from {p}")
        return cfg
    return None


# ---------------------------------------------------------------------------
# Constants  (populated from YAML if available, else hardcoded fallback)
# ---------------------------------------------------------------------------

_cfg = load_config()

if _cfg is not None:
    # ── From YAML ──────────────────────────────────────────────────────────
    _settings = _cfg.get("settings", {})
    BASE = _settings.get("base_dir", "")
    MAX_EVENTS = _settings.get("max_events", 100_000)
    XCACHE_PREFIX = _settings.get("xcache_prefix", "root://xcache.cms.rcac.purdue.edu/")
    DASGOCLIENT = _settings.get("dasgoclient", "/cvmfs/cms.cern.ch/common/dasgoclient")
    SCOUTING_DATA_JSON = _settings.get("scouting_data_json", "")

    _yaml_groups = _cfg.get("groups", {})
    GROUPS = {}
    XSEC = {}
    DY = []
    TT = []
    SM_SIG = []

    for _gname, _samples in _yaml_groups.items():
        _group_dict = {}
        for _sname, _opts in _samples.items():
            # Pass through all keys except xsec_pb as per-sample overrides
            # (dataset_tag, production_id, etc.)
            _group_dict[_sname] = {
                k: v for k, v in _opts.items() if k != "xsec_pb"
            }
            if "xsec_pb" in _opts:
                XSEC[_sname] = _opts["xsec_pb"]
        GROUPS[_gname] = _group_dict

        # Populate convenience lists
        if _gname == "DY":
            DY = list(_samples.keys())
        elif _gname == "TT":
            TT = list(_samples.keys())
        elif _gname == "HHbbtt":
            SM_SIG = list(_samples.keys())

    _data_cfg = _cfg.get("data", {})
    DATA_YEARS = _data_cfg.get("years", ["2024"])
    DATA_RUNS = _data_cfg.get("runs", None)

else:
    # ── Hardcoded fallback (no config/samples.yaml found) ──────────────────
    MAX_EVENTS = 100_000

    BASE = os.path.join(
        "/eos/purdue/store/user/arghyara",
        "production",
        "Scouting",
        "NanoAODv15Scouting24",
        "mcscouting_2024",
    )

    DY = [
        "DYto2L-2Jets_Bin-2J-MLL-50-PTLL-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "DYto2L-2Jets_Bin-MLL-50-PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
        "DYto2L-2Jets_Bin-MLL-50-PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8",
    ]

    TT = [
        "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8",
        "TTto4Q_TuneCP5_13p6TeV_powheg-pythia8",
        "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8",
    ]

    SM_SIG = [
        "GluGluHHto2B2Tau_Par-c2-0p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",
    ]

    GROUPS = {"DY": DY, "TT": TT, "HHbbtt": SM_SIG}

    _TT_INCLUSIVE_PB = 923.6
    _BR_W_LNU = 0.3258
    _BR_W_QQ = 0.6741
    _BR_HH_BBTAUTAU = 2 * 0.5824 * 0.06272
    _SIGMA_GGHH_FB = {"kl-1p00": 34.30}

    XSEC = {
        "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8":  _TT_INCLUSIVE_PB * _BR_W_LNU**2,
        "TTto4Q_TuneCP5_13p6TeV_powheg-pythia8":     _TT_INCLUSIVE_PB * _BR_W_QQ**2,
        "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8":  _TT_INCLUSIVE_PB * 2*_BR_W_LNU*_BR_W_QQ,
        "DYto2L-2Jets_Bin-2J-MLL-50-PTLL-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 170.5,
        "DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":        107.9,
        "DYto2L-2Jets_Bin-MLL-50-PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":         11.13,
        "DYto2L-2Jets_Bin-MLL-50-PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":          0.5914,
        "DYto2L-2Jets_Bin-MLL-50-PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":          0.08008,
        "GluGluHHto2B2Tau_Par-c2-0p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8":
            _SIGMA_GGHH_FB["kl-1p00"] / 1000.0 * _BR_HH_BBTAUTAU,
    }

    XCACHE_PREFIX = "root://xcache.cms.rcac.purdue.edu/"
    DASGOCLIENT = "/cvmfs/cms.cern.ch/common/dasgoclient"
    SCOUTING_DATA_JSON = "/home/das214/HHtobbtautau/Run3_nano_submission/datasets/Scouting_DATA.json"

    DATA_YEARS = ["2024"]
    DATA_RUNS = None


# ---------------------------------------------------------------------------
# Local filesystem NanoAOD file discovery
# ---------------------------------------------------------------------------

def ls_nanoaod_files_from_fs(
    base_dir,
    sample,
    dataset_tag=None,
    production_id=None,
    pattern="*.root",
    recursive=True,
    verbose=0,
    xrootd_prefix=None,
):
    base_dir = os.path.abspath(base_dir)
    sample_dir = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_dir):
        avail = sorted(
            d for d in os.listdir(base_dir)
            if os.path.isdir(os.path.join(base_dir, d))
        )
        raise FileNotFoundError(
            f"Sample dir not found:\n  {sample_dir}\n"
            f"Available samples under {base_dir}:\n  - " + "\n  - ".join(avail)
        )

    # dataset_tag selection
    dataset_tags = sorted(
        d for d in os.listdir(sample_dir)
        if os.path.isdir(os.path.join(sample_dir, d))
    )
    if not dataset_tags:
        raise FileNotFoundError(f"No dataset_tag dirs under:\n  {sample_dir}")

    if dataset_tag is None:
        if len(dataset_tags) == 1:
            dataset_tag = dataset_tags[0]
        else:
            # Pick the most recently modified dataset_tag directory
            dataset_tag = max(
                dataset_tags,
                key=lambda d: os.path.getmtime(os.path.join(sample_dir, d)),
            )
            import warnings
            warnings.warn(
                f"Multiple dataset_tag dirs under {sample_dir}, "
                f"auto-selected newest: {dataset_tag}",
                stacklevel=2,
            )
    else:
        if dataset_tag not in dataset_tags:
            raise ValueError(
                f"dataset_tag '{dataset_tag}' not found under:\n  {sample_dir}\n"
                f"Available:\n  - " + "\n  - ".join(dataset_tags)
            )

    dataset_dir = os.path.join(sample_dir, dataset_tag)

    # production_id selection
    prod_ids = sorted(
        d for d in os.listdir(dataset_dir)
        if os.path.isdir(os.path.join(dataset_dir, d))
    )
    if not prod_ids:
        raise FileNotFoundError(f"No production_id dirs under:\n  {dataset_dir}")

    if production_id is None:
        production_id = prod_ids[-1]
    else:
        if production_id not in prod_ids:
            raise ValueError(
                f"production_id '{production_id}' not found under:\n  {dataset_dir}\n"
                f"Available:\n  - " + "\n  - ".join(prod_ids)
            )

    prod_dir = os.path.join(dataset_dir, production_id)

    # collect files
    if recursive:
        file_glob = os.path.join(prod_dir, "**", pattern)
        files = sorted(glob.glob(file_glob, recursive=True))
    else:
        file_glob = os.path.join(prod_dir, pattern)
        files = sorted(glob.glob(file_glob))

    files = [f for f in files if os.path.isfile(f)]

    if verbose:
        print(f"[ls_nanoaod_files_from_fs]")
        print(f"  base_dir      = {base_dir}")
        print(f"  sample        = {sample}")
        print(f"  dataset_tag   = {dataset_tag}")
        print(f"  production_id = {production_id}")
        print(f"  glob          = {file_glob}")
        print(f"  nfiles        = {len(files)}")
        if verbose > 1:
            for f in files[:10]:
                print("   ", f)

    if not files:
        raise FileNotFoundError(
            f"No ROOT files matched.\n"
            f"  prod_dir   = {prod_dir}\n"
            f"  pattern    = {pattern}\n"
            f"  recursive  = {recursive}\n"
            f"Try pattern='*.root' (already default) or pattern='nano*.root' "
            f"and keep recursive=True."
        )

    if xrootd_prefix is not None:
        files = [xrootd_prefix + f for f in files]

    return files


def ls_nanoaod_files_many(
    base_dir,
    samples,
    dataset_tag=None,
    production_id=None,
    pattern="*.root",
    recursive=True,
    verbose=0,
    allow_missing=False,
    xrootd_prefix=None,
):
    """
    Bulk collect NanoAOD ROOT files for multiple samples.

    samples can be:
      - list of sample names
      - dict: { sample_name: { 'dataset_tag':..., 'production_id':..., 'pattern':... } }
        (per-sample overrides)
    """
    if isinstance(samples, dict):
        sample_names = list(samples.keys())
        per_sample_opts = samples
    else:
        sample_names = list(samples)
        per_sample_opts = {}

    files_by_sample = {}
    all_files = []

    for s in sample_names:
        opts = per_sample_opts.get(s, {})
        try:
            files = ls_nanoaod_files_from_fs(
                base_dir=opts.get("base_dir", base_dir),
                sample=s,
                dataset_tag=opts.get("dataset_tag", dataset_tag),
                production_id=opts.get("production_id", production_id),
                pattern=opts.get("pattern", pattern),
                recursive=opts.get("recursive", recursive),
                verbose=verbose,
                xrootd_prefix=xrootd_prefix,
            )
        except Exception as e:
            if allow_missing:
                if verbose:
                    print(f"[WARN] skipping sample '{s}': {e}")
                continue
            raise

        files_by_sample[s] = files
        all_files.extend(files)

    return files_by_sample, all_files


def ls_nanoaod_files_groups(
    base_dir,
    groups,
    dataset_tag=None,
    production_id=None,
    pattern="*.root",
    recursive=True,
    verbose=0,
    allow_missing=False,
    xrootd_prefix=None,
):
    """
    Collect per-group file lists.
    Returns:
      group_files: dict[group] -> list[str]
      group_files_by_sample: dict[group] -> dict[sample]->list[str]
    """
    group_files = {}
    group_files_by_sample = {}

    for gname, gspec in groups.items():
        by_sample, all_files = ls_nanoaod_files_many(
            base_dir=base_dir,
            samples=gspec,
            dataset_tag=dataset_tag,
            production_id=production_id,
            pattern=pattern,
            recursive=recursive,
            verbose=verbose,
            allow_missing=allow_missing,
            xrootd_prefix=xrootd_prefix,
        )
        group_files[gname] = all_files
        group_files_by_sample[gname] = by_sample

    return group_files, group_files_by_sample


# ---------------------------------------------------------------------------
# XCache / DAS data loading
# ---------------------------------------------------------------------------

def das_files(dataset, prefix=XCACHE_PREFIX):
    """Query DAS for file list and prepend XCache prefix for remote access."""
    result = subprocess.run(
        [DASGOCLIENT, "-query", f"file dataset={dataset}"],
        capture_output=True, text=True, check=True,
    )
    files = [prefix + f.strip() for f in result.stdout.strip().split("\n") if f.strip()]
    return sorted(files)


def load_scouting_data(years=None, runs=None):
    """
    Load scouting data files via XCache.

    Parameters:
        years : list of str, e.g. ["2024"] or ["2022", "2023", "2024"]
                If None, loads all years.
        runs  : list of str, e.g. ["Run2024C", "Run2024D"]
                If None, loads all runs within the selected years.

    Returns:
        all_files : list of XCache file paths
        summary   : dict of {run_name: n_files}
    """
    with open(SCOUTING_DATA_JSON) as f:
        scouting_datasets = json.load(f)

    if years is None:
        years = list(scouting_datasets.keys())

    all_files = []
    summary = {}
    for yr in years:
        if yr not in scouting_datasets:
            print(f"WARNING: year '{yr}' not in JSON, skipping")
            continue
        for run_name, das_path in scouting_datasets[yr].items():
            if runs is not None and run_name not in runs:
                continue
            print(f"  Querying DAS: {run_name} -> {das_path} ...")
            files = das_files(das_path)
            summary[run_name] = len(files)
            all_files.extend(files)
            print(f"    -> {len(files)} files")

    print(f"\nTotal files: {len(all_files)}")
    return all_files, summary


# ---------------------------------------------------------------------------
# Luminosity helpers
# ---------------------------------------------------------------------------

def parse_brilcalc(filepath):
    """Parse brilcalc lumi output into {run: {"lumi": fb, "ncms": int}} dict.

    Handles the pipe-delimited table format from ``brilcalc lumi``.
    Columns: run:fill | time | ncms | hltpath | delivered | recorded
    If a run appears in multiple rows (different HLT versions), lumi and
    ncms are summed (different versions cover different lumisections).
    """
    result = {}
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line.startswith("|"):
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 8:
                continue
            try:
                run = int(parts[1].split(":")[0])
                ncms = int(parts[3])
                rec = float(parts[6])
                if run not in result:
                    result[run] = {"lumi": 0.0, "ncms": 0}
                result[run]["lumi"] += rec
                result[run]["ncms"] += ncms
            except (ValueError, IndexError):
                continue
    return result


def book_lumi_actions(data_df):
    """Book lazy Take actions for run and luminosityBlock columns.

    Returns two RResultPtrs (run_take, ls_take) that should be included
    in a RunGraphs call.  After RunGraphs completes, pass both to
    ``extract_lumi()``.

    Memory: ~62 MB for 7.8M events (vs 7.7 GB for the old 2D histogram
    with ImplicitMT(32)).
    """
    import ROOT as _ROOT
    run_take = data_df.Take["unsigned int"]("run")
    ls_take  = data_df.Take["unsigned int"]("luminosityBlock")
    return run_take, ls_take


def extract_lumi(run_take, ls_take, brilcalc_data):
    """Extract luminosity from filled Take vectors.

    Parameters
    ----------
    run_take, ls_take : RResultPtr<vector<unsigned int>>
        Filled Take results from ``book_lumi_actions()``.
    brilcalc_data : dict
        ``{run: {"lumi": float, "ncms": int}}`` from ``parse_brilcalc()``.

    Returns (lumi_fb, data_runs, n_missing).
    """
    import numpy as np
    import time as _time
    import sys as _sys

    t0 = _time.time()
    runs = np.asarray(run_take.GetValue(), dtype=np.int64)
    ls   = np.asarray(ls_take.GetValue(),  dtype=np.int64)
    t1 = _time.time()
    print(f"    GetValue + asarray: {len(runs):,} events in {t1-t0:.1f}s",
          flush=True)

    # Encode (run, LS) into single int64 key, find unique pairs
    keys = runs * 100000 + ls
    unique_keys = np.unique(keys)
    u_runs = (unique_keys // 100000).astype(np.int32)

    # Count unique LS per run
    unique_run_vals, counts = np.unique(u_runs, return_counts=True)
    run_ls_count = dict(zip(unique_run_vals.tolist(), counts.tolist()))
    t2 = _time.time()
    print(f"    numpy unique: {len(unique_keys):,} unique (run,LS) pairs "
          f"across {len(run_ls_count)} runs in {t2-t1:.1f}s", flush=True)

    # Scale lumi per run by LS fraction
    lumi = 0.0
    n_missing = 0
    for run, n_loaded_ls in run_ls_count.items():
        if run not in brilcalc_data:
            n_missing += 1
            continue
        total_ls = brilcalc_data[run]["ncms"]
        run_lumi = brilcalc_data[run]["lumi"]
        frac = min(n_loaded_ls / total_ls, 1.0) if total_ls > 0 else 0.0
        lumi += run_lumi * frac

    data_runs = sorted(run_ls_count.keys())
    return lumi, data_runs, n_missing


# ---------------------------------------------------------------------------
# File list helpers
# ---------------------------------------------------------------------------

def limit_files(files, sample_name, max_files=0, all_events=False, max_events=0):
    """Filter missing files and apply a file-count cap.

    Parameters
    ----------
    files : list[str]
        Candidate file paths.
    sample_name : str
        Used in warning/error messages.
    max_files : int
        Hard cap on number of files (0 = no hard cap).
    all_events : bool
        If True, skip the max_events-derived cap.
    max_events : int
        Used to derive a file count when max_files == 0 and not all_events.

    Returns
    -------
    list[str]
        Existing files, capped as requested.
    """
    good_files = [f for f in files if os.path.isfile(f)]
    if len(good_files) < len(files):
        print(f"  WARNING: {len(files) - len(good_files)} missing file(s) "
              f"in {sample_name}, using {len(good_files)}/{len(files)}")
    if not good_files:
        raise FileNotFoundError(f"No valid files for {sample_name}")
    caps = []
    if max_files > 0:
        caps.append(max_files)
    if not all_events and max_events > 0:
        caps.append(max(1, max_events // 10_000))
    if caps:
        good_files = good_files[:min(caps)]
    return good_files
