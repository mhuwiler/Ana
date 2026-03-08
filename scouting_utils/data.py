"""
Data loading utilities for scouting NanoAOD analysis.

Handles:
  - Local filesystem NanoAOD file discovery (ls_nanoaod_files_*)
  - XCache / DAS-based scouting data loading
  - RDataFrame helper functions (dropBranchNames, WriteFile, generalise)
  - Sample definitions and constants
"""

import os
import glob
import json
import subprocess

import ROOT

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FILE_SIZE = 100_000

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

# Some samples have ext1 dataset extensions that have been submitted and can be
# loaded similarly (same physics, additional statistics). On disk both the base
# tag (v2) and ext1-v2 directories exist for: kl-0p00, kl-1p00, kl-2p45,
# kl-5p00. The file discovery picks the most recent tag automatically.

SM_SIG = [
    # SM coupling point: c2=0, kl=1, kt=1
    "GluGluHHto2B2Tau_Par-c2-0p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",
]

# BSM_SIG = [
#     # Varied kl (c2=0, kt=1)
#     "GluGluHHto2B2Tau_Par-c2-0p00-kl-0p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # kl=0
#     "GluGluHHto2B2Tau_Par-c2-0p00-kl-2p45-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # kl=2.45
#     "GluGluHHto2B2Tau_Par-c2-0p00-kl-5p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # kl=5
#     # Varied c2 (kt=1)
#     "GluGluHHto2B2Tau_Par-c2-0p10-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # c2=0.1
#     "GluGluHHto2B2Tau_Par-c2-0p35-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # c2=0.35
#     "GluGluHHto2B2Tau_Par-c2-1p00-kl-0p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # c2=1, kl=0
#     "GluGluHHto2B2Tau_Par-c2-3p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",   # c2=3
#     "GluGluHHto2B2Tau_Par-c2-m2p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8",  # c2=-2
#     # Extreme BSM
#     "GluGluHHto2B2Tau_Par-c2-2p24-kl-m20p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8", # c2=2.24, kl=-20
# ] 

GROUPS = {"DY": DY, "TT": TT, "HHbbtt": SM_SIG}

_TT_INCLUSIVE_PB = 923.6 # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO      
_BR_W_LNU  = 0.3258              
_BR_W_QQ   = 0.6741

_BR_HH_BBTAUTAU    = 2 * 0.5824 * 0.06272

# gg->HH cross sections at 13.6 TeV in fb (NNLO FTapprox, LHCHWG)
_SIGMA_GGHH_FB = {
    # "kl-0p00":  75.76,   # kl=0
    "kl-1p00":  34.30,   # kl=1 (SM)
}

XSEC = {
    # -- TT  (inclusive × BR) ------------------------------------------------
    "TTto2L2Nu_TuneCP5_13p6TeV_powheg-pythia8":  _TT_INCLUSIVE_PB * _BR_W_LNU**2,
    "TTto4Q_TuneCP5_13p6TeV_powheg-pythia8":     _TT_INCLUSIVE_PB * _BR_W_QQ**2,
    "TTtoLNu2Q_TuneCP5_13p6TeV_powheg-pythia8":  _TT_INCLUSIVE_PB * 2*_BR_W_LNU*_BR_W_QQ,
    # -- DY NLO  (amcatnloFXFX, binned in pT_ll) ----------------------------
    "DYto2L-2Jets_Bin-2J-MLL-50-PTLL-40to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8": 170.5,
    "DYto2L-2Jets_Bin-MLL-50-PTLL-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":        107.9,
    "DYto2L-2Jets_Bin-MLL-50-PTLL-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":         11.13,
    "DYto2L-2Jets_Bin-MLL-50-PTLL-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":          0.5914,
    "DYto2L-2Jets_Bin-MLL-50-PTLL-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8":          0.08008,
    # -- Signal  (LHCHWG NNLO FTapprox × BR(HH->bbtautau)) ------------------
    # "GluGluHHto2B2Tau_Par-c2-0p00-kl-0p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8":
    #     _SIGMA_GGHH_FB["kl-0p00"] / 1000.0 * _BR_HH_BBTAUTAU,
    "GluGluHHto2B2Tau_Par-c2-0p00-kl-1p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8":
        _SIGMA_GGHH_FB["kl-1p00"] / 1000.0 * _BR_HH_BBTAUTAU,
}

# XCache / DAS settings
XCACHE_PREFIX = "root://xcache.cms.rcac.purdue.edu/"
DASGOCLIENT = "/cvmfs/cms.cern.ch/common/dasgoclient"
SCOUTING_DATA_JSON = "/home/das214/HHtobbtautau/Run3_nano_submission/datasets/Scouting_DATA.json"


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
                base_dir=base_dir,
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
# RDataFrame helpers
# ---------------------------------------------------------------------------

def dropBranchNames(frame, filename, exclusionlist=None):
    """
    Write column (branch) names from an RDataFrame to a text file,
    excluding branches containing substrings in exclusionlist.
    """
    if exclusionlist is None:
        exclusionlist = []
    from ROOT import Ana
    with open(filename, "w") as file:
        for name in frame.GetColumnNames():
            name = str(name)
            if not any(excluded in name for excluded in exclusionlist):
                file.write(f"{name}\n")


def WriteFile(sample, filename, blacklist, treename="Events"):
    from ROOT import Ana
    sample.Snapshot(treename, filename, Ana.purgeColumns(sample.GetColumnNames(), blacklist))


def generalise(df):
    return ROOT.ROOT.RDF.AsRNode(df)
