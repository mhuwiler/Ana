#!/usr/bin/env python3
"""
Cutflow & Trigger Efficiency analysis script.

Output:
    plots/light/  or  plots/dark/   (depending on --theme)

Usage:
    python cutflow_TrigEff.py                              # MC-only plots, light theme
    python cutflow_TrigEff.py --overwrite --theme dark     # regenerate all, dark theme
    python cutflow_TrigEff.py --max-data-files 50          # load 50 data files (default: 10)
    python cutflow_TrigEff.py --no-data                    # skip data entirely (Phase 4+5)
    python cutflow_TrigEff.py --all-events                 # all MC + all data files
    python cutflow_TrigEff.py --max-mc-files 2                # quick test with 2 MC files per sample
    python cutflow_TrigEff.py --plot-type stacked          # only stacked plots
    python cutflow_TrigEff.py --plot-type shape            # only shape overlay plots
    python cutflow_TrigEff.py --skip-cutflow               # skip cutflow, only plots
    python cutflow_TrigEff.py --recache                    # force re-run event loops
    python cutflow_TrigEff.py --plot-vars "ak4_pt*" HT     # only book/plot selected vars

Phases:
    0    Data loading + exact lumi from brilcalc (XCache, quick)
    1+2  Cutflow + MC histograms (single event loop) [--skip-cutflow skips cutflow]
         Cached to .hist_cache.root — subsequent runs skip to Phase 3.
         Use --recache to force re-running.
    3    Draw MC-only plots (no I/O)
    3.5  Trigger shape overlay plots (no I/O, cross-trigger comparison)
    4    Book & run data histograms (XCache, slow) [skip with --no-data]
    5    Draw MC+Data plots (no I/O, reuses Phase 2+4) [skip with --no-data]

Output files (under plots/{theme}/):
    stacked/{trig}_{var}.png          shape/{trig}_{var}.png
    eff_stacked/{trig}_{var}.png      sig_stacked/{trig}_{var}.png
    overlay/{var}.png                 overlay/{ch}_{var}.png
    data_stacked/{trig}_{var}.png     data_shape/{trig}_{var}.png
"""

import os
import sys
import argparse
import fnmatch
import hashlib
import yaml
import math
import time
import resource

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing plots instead of skipping them")
parser.add_argument("--theme", choices=["light", "dark"], default="light",
                    help="Plot colour theme (default: light)")
parser.add_argument("--max-mc-files", type=int, default=0,
                    help="Max files per MC sample (0 = use MAX_EVENTS limit)")
parser.add_argument("--max-data-files", type=int, default=10,
                    help="Max data files per run (0 = all files, default: 10)")
parser.add_argument("--all-events", action="store_true",
                    help="Process all events (ignore MAX_EVENTS and max-data-files limits)")
parser.add_argument("--plot-type", choices=["stacked", "shape", "both"],
                    default="both",
                    help="Which plot types to produce (default: both)")
parser.add_argument("--no-data", action="store_true",
                    help="Skip data overlay on plots (MC only)")
parser.add_argument("--nthreads", type=int, default=4,
                    help="Number of ImplicitMT threads (default: 4)")
parser.add_argument("--skip-cutflow", action="store_true",
                    help="Skip cutflow computation (Phase 1), only produce plots")
parser.add_argument("--recache", action="store_true",
                    help="Force re-running event loops even if histogram cache exists")
parser.add_argument("--plot-vars", nargs="+", default=None, metavar="PATTERN",
                    help="Only book/plot these variables (fnmatch patterns, e.g. 'ak4_pt*' HT)")
parser.add_argument("--skim", action="store_true",
                    help="Write slim ROOT files to /depot/ (keeps all events, drops unused branches)")
parser.add_argument("--reslim", action="store_true",
                    help="Incrementally add new branches to slim files (or full rebuild with --force)")
parser.add_argument("--force", action="store_true",
                    help="With --reslim: full rebuild instead of incremental friend-tree update")
ARGS = parser.parse_args()

# ── Tee stdout/stderr to a timestamped log file ──
LOG_DIR = os.path.join(
    os.path.expanduser("/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana"),
    "logs")
os.makedirs(LOG_DIR, exist_ok=True)
_log_ts = time.strftime("%Y%m%d_%H%M%S")
LOG_PATH = os.path.join(LOG_DIR, f"run_{_log_ts}.log")

class _Tee:
    """Write to both a file and the original stream."""
    def __init__(self, stream, logfile):
        self._stream = stream
        self._log = logfile
    def write(self, msg):
        self._stream.write(msg)
        self._log.write(msg)
        self._log.flush()
    def flush(self):
        self._stream.flush()
        self._log.flush()

_log_fh = open(LOG_PATH, "w")
_log_fh.write(f"# {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
_log_fh.write(f"# Args: {sys.argv}\n\n")
sys.stdout = _Tee(sys.__stdout__, _log_fh)
sys.stderr = _Tee(sys.__stderr__, _log_fh)
print(f"Logging to: {LOG_PATH}")

# Increase stack size to 64 MB to prevent stack overflow from deep
# RDataFrame Define/Filter chains with large TChains.
STACK_SIZE = 64 * 1024 * 1024  # 64 MB
resource.setrlimit(resource.RLIMIT_STACK, (STACK_SIZE, resource.RLIM_INFINITY))

import ROOT
# anaConfig removed — FileFlow.h (old framework) no longer needed

from utils.data import (
    ls_nanoaod_files_groups,
    load_scouting_data,
    parse_brilcalc,
    book_lumi_actions,
    extract_lumi,
    BASE, GROUPS, XSEC, MAX_EVENTS,
    DATA_YEARS, DATA_RUNS,
)
from utils.triggers import (
    DST_JetHT_expr,
    DST_MU_expr,
    DST_EL_expr,
    PARKING_HH_expr,
)
from utils.skim import r_define, r_filter, rdf_exprs, detect_used_branches, run_skim, load_slim_or_eos
from utils.plotting import (
    setup_style,
    cms_label,
    plot_stacked_from_hists,
    plot_shape_from_hists,
    plot_stacked_with_efficiency,
    plot_stacked_with_significance,
    plot_trigger_shape_overlay,
    plot_2d_hist,
    th1_to_np,
)

import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# Histogram cache helpers
# ──────────────────────────────────────────────────────────────────────────────

_CACHE_VERSION = 2   # bump to invalidate all caches

_INVALIDATION_FILES = [
    "cutflow_TrigEff.py",
    "config/samples.yaml",
    "elements/GenMatching.C",
    "elements/RecoObjects.C",
    "elements/common.h",
    "config/objects.yaml",
    "config/acceptance.yaml",
    "config/regions.yaml",
    "config/cuts.yaml",
]


def _cache_meta(max_files, plot_vars_names):
    """Build a dict of metadata used to validate cache freshness."""
    mtimes = {}
    for f in _INVALIDATION_FILES:
        try:
            mtimes[f] = os.path.getmtime(f)
        except OSError:
            mtimes[f] = 0.0
    vars_hash = hashlib.md5(",".join(sorted(plot_vars_names)).encode()).hexdigest()
    return {
        "version": _CACHE_VERSION,
        "max_files": max_files,
        "vars_hash": vars_hash,
        **{f"mtime_{f.replace('/', '_')}": mtimes[f] for f in _INVALIDATION_FILES},
    }


def save_hist_cache(cache_path, mc_hists_by_trig, mc_denom_hists,
                    mc_items_by_trig, lumi, meta, gen_hists=None, excl_hists=None):
    """Save all materialized TH1 objects to a ROOT file."""
    import json
    tf = ROOT.TFile.Open(cache_path, "RECREATE")

    # Save metadata as a TNamed
    tf.mkdir("meta")
    tf.cd("meta")
    meta["lumi"] = lumi
    # mc_items_by_trig is needed to reconstruct labels/colors
    meta["mc_items"] = {tn: items for tn, items in mc_items_by_trig.items()}
    ROOT.TNamed("cache_meta", json.dumps(meta)).Write()

    # Save denominator histograms
    for var_name, h_list in mc_denom_hists.items():
        dname = f"mc_denom/{var_name}"
        tf.mkdir(dname)
        tf.cd(dname)
        for gi, h in enumerate(h_list):
            h.Write(f"g{gi}")

    # Save per-trigger histograms
    for trig_name, h_by_var in mc_hists_by_trig.items():
        for var_name, h_list in h_by_var.items():
            dname = f"mc/{trig_name}/{var_name}"
            tf.mkdir(dname)
            tf.cd(dname)
            for gi, h in enumerate(h_list):
                h.Write(f"g{gi}")

    # Save gen-level histograms (signal only, per decay channel)
    if gen_hists:
        for trig_name, h_by_var in gen_hists.items():
            for var_name, h_by_mode in h_by_var.items():
                dname = f"gen/{trig_name}/{var_name}"
                tf.mkdir(dname)
                tf.cd(dname)
                for mode, h in h_by_mode.items():
                    h.Write(f"m{mode}")

    # Save exclusive-trigger histograms
    if excl_hists:
        for trig_name, h_by_var in excl_hists.items():
            for var_name, h_list in h_by_var.items():
                dname = f"mc_excl/{trig_name}/{var_name}"
                tf.mkdir(dname)
                tf.cd(dname)
                for gi, h in enumerate(h_list):
                    h.Write(f"g{gi}")

    tf.Close()
    sz = os.path.getsize(cache_path) / 1e6
    print(f"Histogram cache saved: {cache_path} ({sz:.1f} MB)")


def load_hist_cache(cache_path, expected_meta, trig_names):
    """Load cached TH1 objects. Returns (mc_hists_by_trig, mc_denom_hists,
    mc_items_by_trig, lumi) or None if cache is invalid."""
    import json

    if not os.path.exists(cache_path):
        return None

    tf = ROOT.TFile.Open(cache_path, "READ")
    if not tf or tf.IsZombie():
        return None

    # Check metadata
    tf.cd("meta")
    named = tf.Get("meta/cache_meta")
    if not named:
        tf.Close()
        return None

    try:
        stored_meta = json.loads(named.GetTitle())
    except (json.JSONDecodeError, AttributeError):
        tf.Close()
        return None

    # Validate cache freshness
    for key in expected_meta:
        if stored_meta.get(key) != expected_meta[key]:
            print(f"Cache invalidated: {key} changed "
                  f"({stored_meta.get(key)} vs {expected_meta[key]})")
            tf.Close()
            return None

    lumi = stored_meta["lumi"]
    mc_items_by_trig = stored_meta.get("mc_items", {})
    # Convert mc_items keys back to proper format
    mc_items_by_trig = {tn: items for tn, items in mc_items_by_trig.items()}

    # Load denominator histograms
    mc_denom_hists = {}
    denom_dir = tf.Get("mc_denom")
    if denom_dir:
        for var_key in denom_dir.GetListOfKeys():
            var_name = var_key.GetName()
            var_dir = tf.Get(f"mc_denom/{var_name}")
            h_list = []
            for gkey in sorted(var_dir.GetListOfKeys(), key=lambda k: k.GetName()):
                h = gkey.ReadObj().Clone()
                h.SetDirectory(0)
                h_list.append(h)
            mc_denom_hists[var_name] = h_list

    # Load per-trigger histograms
    mc_hists_by_trig = {}
    for trig_name in trig_names:
        trig_dir = tf.Get(f"mc/{trig_name}")
        if not trig_dir:
            tf.Close()
            return None
        h_by_var = {}
        for var_key in trig_dir.GetListOfKeys():
            var_name = var_key.GetName()
            var_dir = tf.Get(f"mc/{trig_name}/{var_name}")
            h_list = []
            for gkey in sorted(var_dir.GetListOfKeys(), key=lambda k: k.GetName()):
                h = gkey.ReadObj().Clone()
                h.SetDirectory(0)
                h_list.append(h)
            h_by_var[var_name] = h_list
        mc_hists_by_trig[trig_name] = h_by_var

    # Load gen-level histograms (signal only, per decay channel)
    gen_hists_by_trig = {}
    gen_top = tf.Get("gen")
    if gen_top:
        for trig_key in gen_top.GetListOfKeys():
            trig_name = trig_key.GetName()
            trig_dir  = tf.Get(f"gen/{trig_name}")
            h_by_var  = {}
            for var_key in trig_dir.GetListOfKeys():
                var_name = var_key.GetName()
                var_dir  = tf.Get(f"gen/{trig_name}/{var_name}")
                h_by_mode = {}
                for mkey in var_dir.GetListOfKeys():
                    mode = int(mkey.GetName()[1:])  # "m20" -> 20
                    h = var_dir.Get(mkey.GetName()).Clone()
                    h.SetDirectory(0)
                    h_by_mode[mode] = h
                h_by_var[var_name] = h_by_mode
            gen_hists_by_trig[trig_name] = h_by_var

    # Load exclusive-trigger histograms (optional — missing in old caches → empty dict)
    mc_hists_excl_by_trig = {}
    excl_top = tf.Get("mc_excl")
    if excl_top:
        for trig_key in excl_top.GetListOfKeys():
            trig_name = trig_key.GetName()
            trig_dir  = tf.Get(f"mc_excl/{trig_name}")
            h_by_var  = {}
            for var_key in trig_dir.GetListOfKeys():
                var_name = var_key.GetName()
                var_dir  = tf.Get(f"mc_excl/{trig_name}/{var_name}")
                h_list   = []
                for gkey in sorted(var_dir.GetListOfKeys(), key=lambda k: k.GetName()):
                    h = gkey.ReadObj().Clone()
                    h.SetDirectory(0)
                    h_list.append(h)
                h_by_var[var_name] = h_list
            mc_hists_excl_by_trig[trig_name] = h_by_var

    tf.Close()
    return (mc_hists_by_trig, mc_denom_hists, mc_items_by_trig, lumi,
            gen_hists_by_trig, mc_hists_excl_by_trig)


# ──────────────────────────────────────────────────────────────────────────────
# ROOT setup
# ──────────────────────────────────────────────────────────────────────────────

ANA_DIR = os.path.expanduser(
    "/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana"
)
os.chdir(ANA_DIR)
sys.path.insert(0, ANA_DIR)

ROOT.gInterpreter.AddIncludePath(ANA_DIR)
# ACLiC puts .so files next to source; .gitignore handles them
ROOT.gROOT.LoadMacro("elements/GenMatching.C+")
ROOT.gROOT.LoadMacro("elements/RecoObjects.C+")

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.ROOT.EnableImplicitMT(ARGS.nthreads)

print("CWD =", os.getcwd())
print("sys.path[0] =", sys.path[0])


BTAG_WP = 0.1  # loose b-tag working point for cutflow

# ── Global decay mode LUT ──
# Maps GlobalDecayMode() return value → process, label, color.
# Gap-based numbering: DY=1-5, TT=10-12, HH=20-22, room for QCD=30+, W=40+.
DECAY_MODES = {
    0:  {"process": "unknown",  "label": "Unknown",                          "color": "gray"},
    1:  {"process": "DY",       "label": r"DY (Z$\to$ee)",                   "color": "gold"},
    2:  {"process": "DY",       "label": r"DY (Z$\to\mu\mu$)",              "color": "orange"},
    3:  {"process": "DY",       "label": r"DY (Z$\to\tau_h\tau_h$)",        "color": "darkorange"},
    4:  {"process": "DY",       "label": r"DY (Z$\to\tau_\mu\tau_h$)",      "color": "coral"},
    5:  {"process": "DY",       "label": r"DY (Z$\to\tau_e\tau_h$)",        "color": "tomato"},
    10: {"process": "TT",       "label": r"TT (had)",                        "color": "limegreen"},
    11: {"process": "TT",       "label": r"TT (semi)",                       "color": "forestgreen"},
    12: {"process": "TT",       "label": r"TT (dilep)",                      "color": "teal"},
    20: {"process": "HHbbtt",   "label": r"$HH\to bb\tau_h\tau_h$",         "color": "tab:red"},
    21: {"process": "HHbbtt",   "label": r"$HH\to bb\tau_\mu\tau_h$",       "color": "tab:pink"},
    22: {"process": "HHbbtt",   "label": r"$HH\to bb\tau_e\tau_h$",         "color": "tab:purple"},
    30: {"process": "QCD",      "label": "QCD multijet",                     "color": "tab:cyan"},
}
BKG_MODES = [1, 2, 3, 4, 5, 10, 11, 12, 30]
SIG_MODES = [20, 21, 22]
ACTIVE_MODES = BKG_MODES + SIG_MODES

# Per-trigger brilcalc data: {run: {"lumi": fb, "ncms": int}}
BRILCALC_FILES = {
    "DST_JetHT": os.path.join(ANA_DIR, "output", "brilcalc", "brilcalc_DST_PFScouting_JetHT.txt"),
    "PARKING_HH": os.path.join(ANA_DIR, "output", "brilcalc", "brilcalc_PARKING_HH.txt"),
}
BRILCALC = {}
for _trig_key, _path in BRILCALC_FILES.items():
    if os.path.exists(_path):
        BRILCALC[_trig_key] = parse_brilcalc(_path)
        _total = sum(v["lumi"] for v in BRILCALC[_trig_key].values())
        print(f"Loaded brilcalc [{_trig_key}]: {len(BRILCALC[_trig_key])} runs, "
              f"total = {_total:.3f} fb^-1")
    else:
        print(f"WARNING: brilcalc file not found for {_trig_key}: {_path}")
# Fallback: if PARKING_HH missing, use DST_JetHT
if "PARKING_HH" not in BRILCALC and "DST_JetHT" in BRILCALC:
    BRILCALC["PARKING_HH"] = BRILCALC["DST_JetHT"]
    print("  (using DST_JetHT brilcalc as fallback for PARKING_HH)")


_, group_files_by_sample = ls_nanoaod_files_groups(
    base_dir=BASE,
    groups=GROUPS,
    verbose=0,
)


def _limit_files(files, sample_name):
    """Filter missing files and apply --max-mc-files / MAX_EVENTS file cap."""
    good_files = [f for f in files if os.path.isfile(f)]
    if len(good_files) < len(files):
        print(f"  WARNING: {len(files) - len(good_files)} missing file(s) "
              f"in {sample_name}, using {len(good_files)}/{len(files)}")
    if not good_files:
        raise FileNotFoundError(f"No valid files for {sample_name}")
    if ARGS.max_mc_files > 0:
        good_files = good_files[:ARGS.max_mc_files]
    elif not ARGS.all_events and MAX_EVENTS > 0:
        n_files = max(1, MAX_EVENTS // 10_000)
        good_files = good_files[:n_files]
    return good_files


# ──────────────────────────────────────────────────────────────────────────────
# Data loading  (before MC so we can compute effective lumi for normalisation)
# ──────────────────────────────────────────────────────────────────────────────

YEARS = DATA_YEARS
RUNS = DATA_RUNS  # Configured in config/samples.yaml; null = all runs

if not ARGS.no_data:
    data_files_all, data_summary = load_scouting_data(years=YEARS, runs=RUNS)
    n_total_data = len(data_files_all)

    # Apply data file limit unless --all-events
    if not ARGS.all_events and ARGS.max_data_files > 0:
        data_files = data_files_all[:ARGS.max_data_files]
        print(f"  Limited data to {len(data_files)}/{n_total_data} files "
              f"(--max-data-files {ARGS.max_data_files}, use --all-events for full dataset)")
    else:
        data_files = data_files_all

    # Suppress XRootD error messages from unavailable files (non-fatal)
    _prev_err_level = ROOT.gErrorIgnoreLevel
    ROOT.gErrorIgnoreLevel = ROOT.kFatal
    data_df = ROOT.RDataFrame("Events", data_files)
    ROOT.gErrorIgnoreLevel = _prev_err_level
    print(f"  Loaded {len(data_files)} data files into RDataFrame")

    # Lumi computation is deferred to after Phase 1 (avoids double XCache read).
    # MC weights are defined WITHOUT LUMI; LUMI is applied post-hoc via Scale().
    _brilcalc_default = BRILCALC.get("DST_JetHT", {})
    LUMI = None  # computed after Phase 1
else:
    data_files = []
    data_df = None
    # Use total brilcalc luminosity for MC normalisation
    _brilcalc_default = BRILCALC.get("DST_JetHT", {})
    LUMI = sum(v["lumi"] for v in _brilcalc_default.values()) if _brilcalc_default else 103.965
    print(f"\n[--no-data] Skipping data loading; using brilcalc lumi = {LUMI:.3f} fb^-1")


# ──────────────────────────────────────────────────────────────────────────────
# MC normalisation  (LUMI-free weights; LUMI applied post-hoc)
# ──────────────────────────────────────────────────────────────────────────────

print(f"\nPreparing MC weights (LUMI deferred)  MAX_EVENTS = {MAX_EVENTS}")
mc = {}           # sample_name -> RDataFrame (with column "w")
mc_dfs = {}       # raw DataFrames before weight definition
sum_genw_ptrs = {}

mc_files_map = {}  # sample_name -> list of EOS file paths (for skim)
for group_name, samples_dict in group_files_by_sample.items():
    print(f"\n[{group_name}]")
    for sample_name, files in samples_dict.items():
        good_files = _limit_files(files, sample_name)
        mc_files_map[sample_name] = good_files
        df = load_slim_or_eos(sample_name, good_files)
        sum_genw_ptrs[sample_name] = df.Sum("genWeight")
        mc_dfs[sample_name] = df

print("\nRunning sum(genWeight) for all samples in parallel...")
ROOT.RDF.RunGraphs(list(sum_genw_ptrs.values()))

for sample_name, df in mc_dfs.items():
    sum_genw = sum_genw_ptrs[sample_name].GetValue()
    xsec_pb = XSEC[sample_name]
    scale_nolumi = xsec_pb * 1000.0 / sum_genw  # LUMI applied post-hoc
    print(f"  {sample_name[:60]:60s}  xsec={xsec_pb:.4g} pb  "
          f"sum_genw={sum_genw:.4g}  scale(no L)={scale_nolumi:.4e}")
    w_expr = f"genWeight * {scale_nolumi:.10e}"
    rdf_exprs.append(w_expr)
    mc[sample_name] = df.Define("w", w_expr)

# Convenience: lists of sample names per group (empty if group not in config)
dy_samples  = list(group_files_by_sample.get("DY", {}).keys())
tt_samples  = list(group_files_by_sample.get("TT", {}).keys())
sig_samples = list(group_files_by_sample.get("HHbbtt", {}).keys())
qcd_samples = list(group_files_by_sample.get("QCD", {}).keys())




# ──────────────────────────────────────────────────────────────────────────────
# Gen-level columns
# ──────────────────────────────────────────────────────────────────────────────

_GP = "GenPart"
decay_mode_expr = f"Ana::GlobalDecayMode({_GP}_pdgId, {_GP}_genPartIdxMother, {_GP}_statusFlags)"
gen_hh_expr     = f"Ana::HHGenMatching({_GP}_pdgId, {_GP}_genPartIdxMother, {_GP}_statusFlags)"

for name in mc:
    mc[name] = mc[name].Define("decayMode", decay_mode_expr)

_AK8_GM = "ScoutingFatPFJetRecluster"
_SGP_GM = f"{_AK8_GM}_scoutGlobalParT"
_GM_PROB_NAMES = ["Xbb", "Xbc", "Xbs", "Xcc", "Xcs", "Xss", "Xud", "Xgg", "Xqq",
                  "Xtauhtauh", "Xtauhtaum", "Xtauhtaue", "QCD"]

def _vs_all_at(idx_var, prob_name):
    """Build XvsAll expression at a specific jet index."""
    num = f"{_SGP_GM}_prob_{prob_name}[{idx_var}]"
    denom = " + ".join(f"{_SGP_GM}_prob_{p}[{idx_var}]" for p in _GM_PROB_NAMES)
    return f"(float)({num} / ({denom}))"

for name in sig_samples:
    mc[name] = (mc[name]
        .Define("gen_HH",           gen_hh_expr)
        .Define("genHbb_idx",       "(int)gen_HH.Htob")
        .Define("genHtautau_idx",   "(int)gen_HH.Htotau")
        .Define("genHbb_p4",        "Ana::getP4(genHbb_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("genHtautau_p4",    "Ana::getP4(genHtautau_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("gen_mHH",          "(genHbb_p4 + genHtautau_p4).M()")
        .Define("gen_pt_Hbb",       "genHbb_p4.Pt()")
        .Define("gen_pt_Htautau",   "genHtautau_p4.Pt()")
        .Define("gen_eta_Hbb",      "genHbb_p4.Eta()")
        .Define("gen_eta_Htautau",  "genHtautau_p4.Eta()")
        .Define("gen_phi_Hbb",      "genHbb_p4.Phi()")
        .Define("gen_phi_Htautau",  "genHtautau_p4.Phi()")
        .Define("gen_E_Hbb",        "genHbb_p4.E()")
        .Define("gen_E_Htautau",    "genHtautau_p4.E()")
        .Define("gen_mass_Hbb",     "genHbb_p4.M()")
        .Define("gen_mass_Htautau", "genHtautau_p4.M()")
        .Define("gen_dR_H1H2",     "genHbb_p4.DeltaR(genHtautau_p4)")
        .Define("gen_HH_p4",       "genHbb_p4 + genHtautau_p4")
        .Define("gen_pt_HH",       "gen_HH_p4.Pt()")
        .Define("gen_eta_HH",      "gen_HH_p4.Eta()")
        .Define("gen_phi_HH",      "gen_HH_p4.Phi()")
        .Define("gen_E_HH",        "gen_HH_p4.E()")
        .Define("gen_tau1_charge",
                "gen_HH.tau1 >= 0 ? (float)(GenPart_pdgId[gen_HH.tau1] > 0 ? -1 : 1) : -99.f")
        .Define("gen_tau2_charge",
                "gen_HH.tau2 >= 0 ? (float)(GenPart_pdgId[gen_HH.tau2] > 0 ? -1 : 1) : -99.f")
    )
    # -- Gen-match AK8 jets to gen H→bb and H→ττ --
    mc[name] = (mc[name]
        # ΔR of each AK8 jet to gen H→bb
        .Define("ak8_dR_to_genHbb",
                f"ROOT::VecOps::DeltaR("
                f"{_AK8_GM}_eta, ROOT::RVecF(n{_AK8_GM}, (float)genHbb_p4.Eta()), "
                f"{_AK8_GM}_phi, ROOT::RVecF(n{_AK8_GM}, (float)genHbb_p4.Phi()))")
        .Define("ak8_genHbb_match_idx",
                f"n{_AK8_GM} > 0 && genHbb_p4.Pt() > 0 ? "
                f"(Min(ak8_dR_to_genHbb) < 0.8f ? (int)ArgMin(ak8_dR_to_genHbb) : -1) : -1")
        .Define("ak8_genHbb_match_dR",
                "ak8_genHbb_match_idx >= 0 ? (float)ak8_dR_to_genHbb[ak8_genHbb_match_idx] : -1.f")
        .Define("ak8_genHbb_match_Xbb",
                f"ak8_genHbb_match_idx >= 0 ? {_SGP_GM}_prob_Xbb[ak8_genHbb_match_idx] : -1.f")
        .Define("ak8_genHbb_match_XbbVsAll",
                f"ak8_genHbb_match_idx >= 0 ? {_vs_all_at('ak8_genHbb_match_idx', 'Xbb')} : -1.f")
        # ΔR of each AK8 jet to gen H→ττ
        .Define("ak8_dR_to_genHtt",
                f"ROOT::VecOps::DeltaR("
                f"{_AK8_GM}_eta, ROOT::RVecF(n{_AK8_GM}, (float)genHtautau_p4.Eta()), "
                f"{_AK8_GM}_phi, ROOT::RVecF(n{_AK8_GM}, (float)genHtautau_p4.Phi()))")
        .Define("ak8_genHtt_match_idx",
                f"n{_AK8_GM} > 0 && genHtautau_p4.Pt() > 0 ? "
                f"(Min(ak8_dR_to_genHtt) < 0.8f ? (int)ArgMin(ak8_dR_to_genHtt) : -1) : -1")
        .Define("ak8_genHtt_match_dR",
                "ak8_genHtt_match_idx >= 0 ? (float)ak8_dR_to_genHtt[ak8_genHtt_match_idx] : -1.f")
        .Define("ak8_genHtt_match_Xtt",
                f"ak8_genHtt_match_idx >= 0 ? {_SGP_GM}_prob_Xtauhtauh[ak8_genHtt_match_idx] : -1.f")
        .Define("ak8_genHtt_match_XttVsAll",
                f"ak8_genHtt_match_idx >= 0 ? {_vs_all_at('ak8_genHtt_match_idx', 'Xtauhtauh')} : -1.f")
    )

# Dummy gen-match columns for non-signal samples (AK8 only here; AK4 added after define_kinematics)
_gen_match_cols_ak8 = ["ak8_genHbb_match_dR", "ak8_genHbb_match_Xbb",
                       "ak8_genHbb_match_XbbVsAll", "ak8_genHtt_match_dR",
                       "ak8_genHtt_match_Xtt", "ak8_genHtt_match_XttVsAll"]
for name in mc:
    if name not in sig_samples:
        for col in _gen_match_cols_ak8:
            mc[name] = mc[name].Define(col, "-1.f")


# ──────────────────────────────────────────────────────────────────────────────
# Kinematic definitions
# ──────────────────────────────────────────────────────────────────────────────

_AK4 = "ScoutingPFJetRecluster2"

def define_kinematics(df):
    """Define derived kinematic columns for jets.

    Intermediate columns are avoided to keep the branch-proxy chain shallow
    and prevent stack overflows with large TChains.
    """
    df = (df
        .Define("ak4_pt0", "ScoutingPFJetRecluster2_pt[0]") # Check this (Could be a bug in the processing)
        .Define("ak4_pt1", "ScoutingPFJetRecluster2_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster2_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster2_pt[3]")
        .Define("HT",  "Sum(ScoutingPFJetRecluster2_pt)")
        .Define("nJets", "nScoutingPFJetRecluster2")
        .Define("nLeptons", "nScoutingMuonVtx + nScoutingElectron")
    )
    # -- AK8 fat jet kinematics + ScoutGlobalParT tagger + XvsQCD (pT > 150 GeV) --
    _AK8 = "ScoutingFatPFJetRecluster"
    _SGP = f"{_AK8}_scoutGlobalParT"
    _AK8_PT_MIN = 150.0
    df = df.Define("nFatJets", f"(int)Sum({_AK8}_pt > {_AK8_PT_MIN}f)")
    for _i in range(3):
        _has = f"(n{_AK8}>={_i+1} && {_AK8}_pt[{_i}]>{_AK8_PT_MIN}f)"
        df = (df
            .Define(f"ak8_pt{_i}",   f"{_has} ? {_AK8}_pt[{_i}] : -1.f")
            .Define(f"ak8_eta{_i}",  f"{_has} ? {_AK8}_eta[{_i}] : -99.f")
            .Define(f"ak8_mass{_i}", f"{_has} ? {_AK8}_mass[{_i}] : -1.f")
            .Define(f"ak8_msd{_i}",  f"{_has} ? {_AK8}_msoftdrop[{_i}] : -1.f")
            .Define(f"ak8_Xbb{_i}",  f"{_has} ? {_SGP}_prob_Xbb[{_i}] : -1.f")
            .Define(f"ak8_Xtt{_i}",  f"{_has} ? {_SGP}_prob_Xtauhtauh[{_i}] : -1.f")
            .Define(f"ak8_Xtm{_i}",  f"{_has} ? {_SGP}_prob_Xtauhtaum[{_i}] : -1.f")
            .Define(f"ak8_Xte{_i}",  f"{_has} ? {_SGP}_prob_Xtauhtaue[{_i}] : -1.f")
            .Define(f"ak8_QCD{_i}",  f"{_has} ? {_SGP}_prob_QCD[{_i}] : -1.f")
            .Define(f"ak8_XbbVsQCD{_i}", f"{_has} ? (float)({_SGP}_prob_Xbb[{_i}] / ({_SGP}_prob_Xbb[{_i}] + {_SGP}_prob_QCD[{_i}])) : -1.f")
            .Define(f"ak8_XttVsQCD{_i}", f"{_has} ? (float)({_SGP}_prob_Xtauhtauh[{_i}] / ({_SGP}_prob_Xtauhtauh[{_i}] + {_SGP}_prob_QCD[{_i}])) : -1.f")
            .Define(f"ak8_XtmVsQCD{_i}", f"{_has} ? (float)({_SGP}_prob_Xtauhtaum[{_i}] / ({_SGP}_prob_Xtauhtaum[{_i}] + {_SGP}_prob_QCD[{_i}])) : -1.f")
            .Define(f"ak8_XteVsQCD{_i}", f"{_has} ? (float)({_SGP}_prob_Xtauhtaue[{_i}] / ({_SGP}_prob_Xtauhtaue[{_i}] + {_SGP}_prob_QCD[{_i}])) : -1.f")
            .Define(f"ak8_massCorr{_i}", f"{_has} ? (float)({_AK8}_mass[{_i}] * {_SGP}_massCorrGeneric[{_i}]) : -1.f")
            .Define(f"ak8_massRes{_i}",  f"{_has} ? (float)({_AK8}_mass[{_i}] * {_SGP}_massCorrResonance[{_i}]) : -1.f")
        )
    # -- AK8 tagger-based candidates (best XvsAll score per hypothesis) --
    _pt_cut = f"{_AK8}_pt > {_AK8_PT_MIN}f"
    # Denominator: sum of ALL ScoutGlobalParT probabilities
    _all_probs = " + ".join(f"{_SGP}_prob_{p}" for p in [
        "Xbb", "Xbc", "Xbs", "Xcc", "Xcs", "Xss", "Xud", "Xgg", "Xqq",
        "Xtauhtauh", "Xtauhtaum", "Xtauhtaue", "QCD",
    ])
    _cand_scores = {
        "Hbb": f"{_SGP}_prob_Xbb / ({_all_probs})",
        "Htt": f"{_SGP}_prob_Xtauhtauh / ({_all_probs})",
        "Htm": f"{_SGP}_prob_Xtauhtaum / ({_all_probs})",
        "Hte": f"{_SGP}_prob_Xtauhtaue / ({_all_probs})",
    }
    for _cand, _score_expr in _cand_scores.items():
        _masked = (f"ROOT::VecOps::Where({_pt_cut}, "
                   f"(ROOT::RVecF)({_score_expr}), "
                   f"ROOT::RVecF(n{_AK8}, -1.f))")
        df = (df
            .Define(f"ak8_{_cand}_idx",   f"nFatJets > 0 ? (int)ArgMax({_masked}) : -1")
            .Define(f"ak8_{_cand}_score", f"ak8_{_cand}_idx >= 0 ? ({_masked})[ak8_{_cand}_idx] : -1.f")
            .Define(f"ak8_{_cand}_pt",    f"ak8_{_cand}_idx >= 0 ? {_AK8}_pt[ak8_{_cand}_idx] : -1.f")
            .Define(f"ak8_{_cand}_eta",   f"ak8_{_cand}_idx >= 0 ? {_AK8}_eta[ak8_{_cand}_idx] : -99.f")
            .Define(f"ak8_{_cand}_mass",  f"ak8_{_cand}_idx >= 0 ? {_AK8}_mass[ak8_{_cand}_idx] : -1.f")
            .Define(f"ak8_{_cand}_msd",   f"ak8_{_cand}_idx >= 0 ? {_AK8}_msoftdrop[ak8_{_cand}_idx] : -1.f")
            .Define(f"ak8_{_cand}_phi",   f"ak8_{_cand}_idx >= 0 ? {_AK8}_phi[ak8_{_cand}_idx] : -99.f")
        )
    # -- ΔR between Hbb and Htt candidates --
    df = df.Define("ak8_dR_Hbb_Htt",
                   "(ak8_Hbb_idx >= 0 && ak8_Htt_idx >= 0) ? "
                   "ROOT::VecOps::DeltaR(ak8_Hbb_eta, ak8_Htt_eta, ak8_Hbb_phi, ak8_Htt_phi) : -1.f")
    df = (df.Define("mjj_01",
                "(float)(ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[0],ScoutingPFJetRecluster2_eta[0],"
                "ScoutingPFJetRecluster2_phi[0],ScoutingPFJetRecluster2_mass[0])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[1],ScoutingPFJetRecluster2_eta[1],"
                "ScoutingPFJetRecluster2_phi[1],ScoutingPFJetRecluster2_mass[1])).M()")
        .Define("dR_01",
                "ROOT::VecOps::DeltaR("
                "ScoutingPFJetRecluster2_eta[0],ScoutingPFJetRecluster2_eta[1],"
                "ScoutingPFJetRecluster2_phi[0],ScoutingPFJetRecluster2_phi[1])")
        .Define("MHT",
                "(float)sqrt("
                "pow(Sum(ScoutingPFJetRecluster2_pt*cos(ScoutingPFJetRecluster2_phi)),2) + "
                "pow(Sum(ScoutingPFJetRecluster2_pt*sin(ScoutingPFJetRecluster2_phi)),2))")
        .Define("m4j",
                "(float)(ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[0],ScoutingPFJetRecluster2_eta[0],"
                "ScoutingPFJetRecluster2_phi[0],ScoutingPFJetRecluster2_mass[0])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[1],ScoutingPFJetRecluster2_eta[1],"
                "ScoutingPFJetRecluster2_phi[1],ScoutingPFJetRecluster2_mass[1])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[2],ScoutingPFJetRecluster2_eta[2],"
                "ScoutingPFJetRecluster2_phi[2],ScoutingPFJetRecluster2_mass[2])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster2_pt[3],ScoutingPFJetRecluster2_eta[3],"
                "ScoutingPFJetRecluster2_phi[3],ScoutingPFJetRecluster2_mass[3])).M()")
    )
    # b-jet selection: ParticleNet discriminators (prob_bb is for AK8 fat jets, not AK4)
    # upart_b_raw  = raw UParT probb score per jet
    # ak4_BvsAll   = probb / (probb + probc + probg + probuds + problepb + probtaum + probtaup)
    # Sort jets by BvsAll in descending order; [0] = highest, [1] = second-highest
    df = (df
        .Define("upart_b_raw",    f"{_AK4}_scoutUParT_probb")
        .Define("upart_taup_raw", f"{_AK4}_scoutUParT_probtaup")
        .Define("upart_taum_raw", f"{_AK4}_scoutUParT_probtaum")
        .Define("upart_tau_raw",  f"{_AK4}_scoutUParT_probtaup + {_AK4}_scoutUParT_probtaum")
        .Define("upart_c_raw",    f"{_AK4}_scoutUParT_probc")
        .Define("upart_g_raw",    f"{_AK4}_scoutUParT_probg")
        .Define("upart_uds_raw",  f"{_AK4}_scoutUParT_probuds")
        .Define("upart_lepb_raw", f"{_AK4}_scoutUParT_problepb")
        .Define("ak4_BvsAll",
                f"{_AK4}_scoutUParT_probb"
                f" / ({_AK4}_scoutUParT_probb"
                f" + {_AK4}_scoutUParT_probc + {_AK4}_scoutUParT_probg"
                f" + {_AK4}_scoutUParT_probuds + {_AK4}_scoutUParT_problepb"
                f" + {_AK4}_scoutUParT_probtaum + {_AK4}_scoutUParT_probtaup)")
        .Define("ak4_TaupVsAll",
                f"{_AK4}_scoutUParT_probtaup"
                f" / ({_AK4}_scoutUParT_probb"
                f" + {_AK4}_scoutUParT_probc + {_AK4}_scoutUParT_probg"
                f" + {_AK4}_scoutUParT_probuds + {_AK4}_scoutUParT_problepb"
                f" + {_AK4}_scoutUParT_probtaum + {_AK4}_scoutUParT_probtaup)")
        .Define("ak4_TaumVsAll",
                f"{_AK4}_scoutUParT_probtaum"
                f" / ({_AK4}_scoutUParT_probb"
                f" + {_AK4}_scoutUParT_probc + {_AK4}_scoutUParT_probg"
                f" + {_AK4}_scoutUParT_probuds + {_AK4}_scoutUParT_problepb"
                f" + {_AK4}_scoutUParT_probtaum + {_AK4}_scoutUParT_probtaup)")
        .Define("ak4_TauVsAll",
                f"({_AK4}_scoutUParT_probtaup + {_AK4}_scoutUParT_probtaum)"
                f" / ({_AK4}_scoutUParT_probb"
                f" + {_AK4}_scoutUParT_probc + {_AK4}_scoutUParT_probg"
                f" + {_AK4}_scoutUParT_probuds + {_AK4}_scoutUParT_problepb"
                f" + {_AK4}_scoutUParT_probtaum + {_AK4}_scoutUParT_probtaup)")
        .Define("bsort_idx",
                "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(ak4_BvsAll))")
        .Define("b0_idx", "(int)bsort_idx[0]")
        .Define("b1_idx", "(int)bsort_idx[1]")
        .Define("b0_pt",    "ScoutingPFJetRecluster2_pt[b0_idx]")
        .Define("b0_eta",   "ScoutingPFJetRecluster2_eta[b0_idx]")
        .Define("b0_phi",   "ScoutingPFJetRecluster2_phi[b0_idx]")
        .Define("b0_mass",  "ScoutingPFJetRecluster2_mass[b0_idx]")
        .Define("b0_score", "ak4_BvsAll[b0_idx]")
        .Define("b1_pt",    "ScoutingPFJetRecluster2_pt[b1_idx]")
        .Define("b1_eta",   "ScoutingPFJetRecluster2_eta[b1_idx]")
        .Define("b1_phi",   "ScoutingPFJetRecluster2_phi[b1_idx]")
        .Define("b1_mass",  "ScoutingPFJetRecluster2_mass[b1_idx]")
        .Define("b1_score", "ak4_BvsAll[b1_idx]")
        .Define("b0_raw",   "upart_b_raw[b0_idx]")
        .Define("b1_raw",   "upart_b_raw[b1_idx]")
        .Define("mbb",
                "(float)(ROOT::Math::PtEtaPhiMVector(b0_pt,b0_eta,b0_phi,b0_mass)"
                " + ROOT::Math::PtEtaPhiMVector(b1_pt,b1_eta,b1_phi,b1_mass)).M()")
        .Define("dR_bb",
                "ROOT::VecOps::DeltaR(b0_eta,b1_eta,b0_phi,b1_phi)")
        .Define("ptbb",
                "(float)(ROOT::Math::PtEtaPhiMVector(b0_pt,b0_eta,b0_phi,b0_mass)"
                " + ROOT::Math::PtEtaPhiMVector(b1_pt,b1_eta,b1_phi,b1_mass)).Pt()")
    )
    df = (df
        .Define("ak4_eta0", "ScoutingPFJetRecluster2_eta[0]")
        .Define("ak4_eta1", "ScoutingPFJetRecluster2_eta[1]")
        .Define("ak4_eta2", "ScoutingPFJetRecluster2_eta[2]")
        .Define("ak4_eta3", "ScoutingPFJetRecluster2_eta[3]")
        .Define("ak4_mass0", "ScoutingPFJetRecluster2_mass[0]")
        .Define("ak4_mass1", "ScoutingPFJetRecluster2_mass[1]")
        .Define("nMuons", "nScoutingMuonVtx")
        .Define("nElectrons", "nScoutingElectron")
        .Define("centrality",
                "(float)(Sum(ScoutingPFJetRecluster2_pt) / "
                "Sum(ScoutingPFJetRecluster2_pt * cosh(ScoutingPFJetRecluster2_eta)))")
        .Define("dEta_01",
                "(float)abs(ScoutingPFJetRecluster2_eta[0] - ScoutingPFJetRecluster2_eta[1])")
    )
    # H→bb candidate: best dijet pair with m_jj ∈ [100, 150] GeV (closest to 125)
    df = (df
        .Define("hbb_pair",
                "Ana::findDijetInWindow("
                "ScoutingPFJetRecluster2_pt, ScoutingPFJetRecluster2_eta,"
                "ScoutingPFJetRecluster2_phi, ScoutingPFJetRecluster2_mass,"
                "100.f, 150.f, 125.f)")
        .Define("has_hbb",      "hbb_pair.i1 >= 0")
        .Define("mbb_cand",     "hbb_pair.mass")
        # H→ττ candidate: best dijet pair from remaining jets, m_jj ∈ [40, 150] GeV
        .Define("htautau_pair",
                "Ana::findDijetInWindow("
                "ScoutingPFJetRecluster2_pt, ScoutingPFJetRecluster2_eta,"
                "ScoutingPFJetRecluster2_phi, ScoutingPFJetRecluster2_mass,"
                "40.f, 150.f, 80.f, {hbb_pair.i1, hbb_pair.i2})")
        .Define("has_htautau",  "htautau_pair.i1 >= 0")
        .Define("mtautau_cand", "htautau_pair.mass")
    )
    return df

if data_df is not None:
    data_df = define_kinematics(data_df)
for name in mc:
    mc[name] = define_kinematics(mc[name])

# ── Gen-matched AK4 jets (signal only, after define_kinematics for ak4_BvsAll/TaupVsAll) ──
_AK4_GM = "ScoutingPFJetRecluster2"
_ak4_gen_particles = [
    ("b1",   "gen_HH.b1",   "BvsAll"),
    ("b2",   "gen_HH.b2",   "BvsAll"),
    ("tau1", "gen_HH.tau1", "TauVsAll"),
    ("tau2", "gen_HH.tau2", "TauVsAll"),
]
_ak4_gen_cols = []
for _gp_label, _gp_expr, _score_name in _ak4_gen_particles:
    _ak4_gen_cols.extend([
        f"ak4_gen{_gp_label}_dR",
        f"ak4_gen{_gp_label}_{_score_name}",
    ])

for name in sig_samples:
    for _gp_label, _gp_expr, _score_name in _ak4_gen_particles:
        _p4 = f"gen_{_gp_label}_p4"
        _dR_vec = f"ak4_dR_to_gen{_gp_label}"
        _idx = f"ak4_gen{_gp_label}_idx"
        mc[name] = (mc[name]
            .Define(_p4, f"Ana::getP4((int){_gp_expr}, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
            .Define(_dR_vec,
                    f"ROOT::VecOps::DeltaR("
                    f"{_AK4_GM}_eta, ROOT::RVecF(n{_AK4_GM}, (float){_p4}.Eta()), "
                    f"{_AK4_GM}_phi, ROOT::RVecF(n{_AK4_GM}, (float){_p4}.Phi()))")
            .Define(_idx,
                    f"n{_AK4_GM} > 0 && {_p4}.Pt() > 0 ? "
                    f"(int)(Min({_dR_vec}) < 0.4f ? ArgMin({_dR_vec}) : -1) : -1")
            .Define(f"ak4_gen{_gp_label}_dR",
                    f"{_idx} >= 0 ? {_dR_vec}[{_idx}] : -1.f")
            .Define(f"ak4_gen{_gp_label}_{_score_name}",
                    f"{_idx} >= 0 ? ak4_{_score_name}[{_idx}] : -1.f")
        )

# Additional tau scores (TaumVsAll, TaupVsAll) using existing gen-matched indices
_ak4_gen_tau_extra = []
for _tl in ["tau1", "tau2"]:
    for _sc in ["TaumVsAll", "TaupVsAll"]:
        _ak4_gen_tau_extra.append(f"ak4_gen{_tl}_{_sc}")
for name in sig_samples:
    for _tl in ["tau1", "tau2"]:
        _idx = f"ak4_gen{_tl}_idx"
        mc[name] = (mc[name]
            .Define(f"ak4_gen{_tl}_TaumVsAll", f"{_idx} >= 0 ? ak4_TaumVsAll[{_idx}] : -1.f")
            .Define(f"ak4_gen{_tl}_TaupVsAll", f"{_idx} >= 0 ? ak4_TaupVsAll[{_idx}] : -1.f")
        )

# Charge misidentification columns: combined across tau1/tau2, keyed by true gen charge.
# τ⁻ (pdgId=15): charge=-1 → charge < 0 && > -90 (excludes sentinel -99)
# τ⁺ (pdgId=-15): charge=+1 → charge > 0 (sentinel -99 is negative, so naturally excluded)
for name in sig_samples:
    def _neg(tl): return (f"gen_{tl}_charge < 0.f && gen_{tl}_charge > -90.f"
                          f" && ak4_gen{tl}_dR > -0.5f")
    def _pos(tl): return f"gen_{tl}_charge > 0.f && ak4_gen{tl}_dR > -0.5f"
    mc[name] = (mc[name]
        # τ⁻ jet: correct-sign score (should peak at 1)
        .Define("gen_taum_TaumVsAll",
                f"({_neg('tau1')}) ? ak4_gentau1_TaumVsAll :"
                f" ({_neg('tau2')}) ? ak4_gentau2_TaumVsAll : -1.f")
        # τ⁻ jet: wrong-sign score (should peak at 0 if charge ID works)
        .Define("gen_taum_TaupVsAll",
                f"({_neg('tau1')}) ? ak4_gentau1_TaupVsAll :"
                f" ({_neg('tau2')}) ? ak4_gentau2_TaupVsAll : -1.f")
        # τ⁻ jet: combined charge-agnostic score (uses ak4_TauVsAll directly)
        .Define("gen_taum_TauVsAll",
                f"({_neg('tau1')}) ? ak4_gentau1_TauVsAll :"
                f" ({_neg('tau2')}) ? ak4_gentau2_TauVsAll : -1.f")
        # τ⁺ jet: correct-sign score (should peak at 1)
        .Define("gen_taup_TaupVsAll",
                f"({_pos('tau1')}) ? ak4_gentau1_TaupVsAll :"
                f" ({_pos('tau2')}) ? ak4_gentau2_TaupVsAll : -1.f")
        # τ⁺ jet: wrong-sign score (should peak at 0 if charge ID works)
        .Define("gen_taup_TaumVsAll",
                f"({_pos('tau1')}) ? ak4_gentau1_TaumVsAll :"
                f" ({_pos('tau2')}) ? ak4_gentau2_TaumVsAll : -1.f")
        # τ⁺ jet: combined charge-agnostic score (uses ak4_TauVsAll directly)
        .Define("gen_taup_TauVsAll",
                f"({_pos('tau1')}) ? ak4_gentau1_TauVsAll :"
                f" ({_pos('tau2')}) ? ak4_gentau2_TauVsAll : -1.f")
    )

# Raw UParT prob scores for gen-matched jets
_raw_scores = [
    ("raw_b",    "upart_b_raw"),
    ("raw_taup", "upart_taup_raw"),
    ("raw_taum", "upart_taum_raw"),
    ("raw_tau",  "upart_tau_raw"),
    ("raw_c",    "upart_c_raw"),
    ("raw_g",    "upart_g_raw"),
    ("raw_uds",  "upart_uds_raw"),
    ("raw_lepb", "upart_lepb_raw"),
]
_ak4_gen_raw_cols = []
for name in sig_samples:
    for _gp_label in ["b1", "b2", "tau1", "tau2"]:
        _idx = f"ak4_gen{_gp_label}_idx"
        for _score_suffix, _vec_name in _raw_scores:
            col = f"ak4_gen{_gp_label}_{_score_suffix}"
            mc[name] = mc[name].Define(col, f"{_idx} >= 0 ? {_vec_name}[{_idx}] : -1.f")
            if col not in _ak4_gen_raw_cols:
                _ak4_gen_raw_cols.append(col)

# Dummy AK4 gen-match columns for non-signal samples
_all_ak4_gen_cols = _ak4_gen_cols + _ak4_gen_tau_extra + _ak4_gen_raw_cols
for name in mc:
    if name not in sig_samples:
        for col in _all_ak4_gen_cols:
            mc[name] = mc[name].Define(col, "-1.f")

# ── Load and apply user cuts from config/cuts.yaml ──
_cuts_cfg = os.path.join(ANA_DIR, "config", "cuts.yaml")
USER_CUTS = []
if os.path.exists(_cuts_cfg):
    with open(_cuts_cfg) as _f:
        _cuts_data = yaml.safe_load(_f) or {}
    USER_CUTS = _cuts_data.get("cuts", []) or []

if USER_CUTS:
    print(f"\n[cuts] Applying {len(USER_CUTS)} user cuts from config/cuts.yaml:")
    for _cut in USER_CUTS:
        print(f"  → {_cut}")
        for name in mc:
            mc[name] = mc[name].Filter(_cut)
        if data_df is not None:
            data_df = data_df.Filter(_cut)


# ──────────────────────────────────────────────────────────────────────────────
# Selection & style
# ──────────────────────────────────────────────────────────────────────────────

base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_acc = data_df.Filter(base_cut) if data_df is not None else None
mc_acc = {name: df.Filter(base_cut) for name, df in mc.items()}

setup_style(dark=(ARGS.theme == "dark"))

PLOT_DIR = os.path.join(ANA_DIR, "plots", ARGS.theme)
os.makedirs(PLOT_DIR, exist_ok=True)


def _plot_path(category, plot_type, trig_name=None, var_name=""):
    """Build plot output path and ensure directory exists.

    category:  "mc" or "data"
    plot_type: "stacked", "shape", "eff", "sig", "overlay", "gen"
    trig_name: trigger name (creates subdirectory), None for overlay/gen
    var_name:  variable name for the filename (without .png)
    """
    if trig_name:
        d = os.path.join(PLOT_DIR, category, plot_type, trig_name)
    else:
        d = os.path.join(PLOT_DIR, category, plot_type)
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, f"{var_name}.png")

TRIG_LIST = [
    ("NoTrigger",    None),
    ("DST_JetHT",    DST_JetHT_expr),
    ("DST_Muon",     DST_MU_expr),
    ("DST_Electron", DST_EL_expr),
    ("PARKING_HH",   PARKING_HH_expr),
]

# Exclusive triggers: events that fired ONLY one trigger (not any other)
_trig_exprs = [e for _, e in TRIG_LIST if e is not None]
def _excl(t):
    others = " || ".join(e for e in _trig_exprs if e != t)
    return f"({t}) && !({others})"

EXCL_TRIG_LIST = [
    ("NoTrigger",         None),
    ("DST_JetHT_excl",    _excl(DST_JetHT_expr)),
    ("DST_Muon_excl",     _excl(DST_MU_expr)),
    ("DST_Electron_excl", _excl(DST_EL_expr)),
    ("PARKING_HH_excl",   _excl(PARKING_HH_expr)),
]

CUTFLOW_STEPS = [
    ("Trigger",                       None),
    ("≥4 jets",                       "nScoutingPFJetRecluster2 >= 4"),
    ("4 jets pT > 20",                "ak4_pt0 > 20 && ak4_pt1 > 20 && ak4_pt2 > 20 && ak4_pt3 > 20"),
    # TODO: add UParT-based b-tag and tau-tag cuts after gen-matched score study
]

PLOT_VARS = [
    # (branch,     xlabel,                    nbins, xmin, xmax)
    # -- Jet pT --
    ("ak4_pt0",   r"Leading jet $p_T$ [GeV]",        50,    0,  250),
    ("ak4_pt1",   r"Sub-leading jet $p_T$ [GeV]",    50,    0,  250),
    ("ak4_pt2",   r"3rd jet $p_T$ [GeV]",            50,    0,   250),
    ("ak4_pt3",   r"4th jet $p_T$ [GeV]",            50,    0,   250),
    # -- Jet eta --
    ("ak4_eta0",  r"Leading jet $\eta$",              30,   -5,     5),
    ("ak4_eta1",  r"Sub-leading jet $\eta$",          30,   -5,     5),
    ("ak4_eta2",  r"3rd jet $\eta$",                  30,   -5,     5),
    ("ak4_eta3",  r"4th jet $\eta$",                  30,   -5,     5),
    # -- Jet mass --
    ("ak4_mass0", r"Leading jet mass [GeV]",          30,    0,   100),
    ("ak4_mass1", r"Sub-leading jet mass [GeV]",      30,    0,   100),
    # -- Multiplicity --
    ("nJets",     r"Number of AK4 jets",               6,    0,     6),
    ("nLeptons",  r"Number of leptons ($\mu + e$)",   10,    0,    10),
    ("nMuons",    r"Number of muons",                  6,    0,     6),
    ("nElectrons", r"Number of electrons",             6,    0,     6),
    # -- Global event --
    ("HT",        r"$H_T$ [GeV]",                    200,    0,  2000),
    ("MHT",       r"$\slash{H}_T$ [GeV]",             30,    0,  1500),
    ("centrality", r"Centrality ($H_T / \sum E$)",    25,    0,     1),
    # -- Dijet (leading pair) --
    ("mjj_01",    r"$m_{jj}$ (leading dijet) [GeV]", 40,    0,  2000),
    ("dR_01",     r"$\Delta R(j_0, j_1)$",           30,    0,     6),
    ("dEta_01",   r"$|\Delta\eta(j_0, j_1)|$",       25,    0,     5),
    # -- 4-jet --
    ("m4j",       r"$m_{4j}$ (leading 4 jets) [GeV]", 32,   0,   800),
    # -- b-tagged jets --
    ("b0_pt",     r"$b_0$ jet $p_T$ [GeV]",          36,    0,  1000),
    ("b1_pt",     r"$b_1$ jet $p_T$ [GeV]",          36,    0,   600),
    ("b0_score",  r"$b_0$ UParT BvsAll",                25,  0,     1),
    ("b1_score",  r"$b_1$ UParT BvsAll",                25,  0,     1),
    ("b0_raw",    r"$b_0$ UParT raw prob\_b",           25,    0,     1),
    ("b1_raw",    r"$b_1$ UParT raw prob\_b",           25,    0,     1),
    ("mbb",       r"$m_{bb}$ [GeV]",                  30,    0,   300),
    ("dR_bb",     r"$\Delta R(b_0, b_1)$",           30,    0,     6),
    ("ptbb",      r"$p_T^{bb}$ [GeV]",               30,    0,  1000),
    # -- Higgs candidate dijet masses --
    ("mbb_cand",     r"$m_{jj}^{H \to bb}$ candidate [GeV]",    25,   50,   200),
    ("mtautau_cand", r"$m_{jj}^{H \to \tau\tau}$ candidate [GeV]", 25, 0, 200),
    # -- AK8 fat jets (pT > 150 GeV; entries with -1 are filtered at booking) --
    ("nFatJets",  r"Number of AK8 fat jets",  8, 0, 8),
] + [
    entry
    for i, (lbl, pt_max) in enumerate([("Leading", 1000), ("Sub-leading", 800), ("Third", 600)])
    for entry in [
        (f"ak8_pt{i}",        rf"{lbl} AK8 $p_T$ [GeV]",                          50, 150, pt_max),
        (f"ak8_eta{i}",       rf"{lbl} AK8 $\eta$",                                30,  -5,     5),
        (f"ak8_mass{i}",      rf"{lbl} AK8 mass [GeV]",                            40,   0,   400),
        (f"ak8_msd{i}",       rf"{lbl} AK8 $m_{{SD}}$ [GeV]",                      40,   0,   300),
        (f"ak8_Xbb{i}",       rf"{lbl} AK8 $X_{{bb}}$",                            25,   0,     1),
        (f"ak8_Xtt{i}",       rf"{lbl} AK8 $X_{{\tau_h\tau_h}}$",                  25,   0,     1),
        (f"ak8_Xtm{i}",       rf"{lbl} AK8 $X_{{\tau_\mu\tau_h}}$",                25,   0,     1),
        (f"ak8_Xte{i}",       rf"{lbl} AK8 $X_{{\tau_e\tau_h}}$",                  25,   0,     1),
        (f"ak8_QCD{i}",       rf"{lbl} AK8 QCD",                                   25,   0,     1),
        (f"ak8_XbbVsQCD{i}",  rf"{lbl} AK8 $X_{{bb}}$ vs QCD",                     25,   0,     1),
        (f"ak8_XttVsQCD{i}",  rf"{lbl} AK8 $X_{{\tau_h\tau_h}}$ vs QCD",           25,   0,     1),
        (f"ak8_XtmVsQCD{i}",  rf"{lbl} AK8 $X_{{\tau_\mu\tau_h}}$ vs QCD",         25,   0,     1),
        (f"ak8_XteVsQCD{i}",  rf"{lbl} AK8 $X_{{\tau_e\tau_h}}$ vs QCD",           25,   0,     1),
        (f"ak8_massCorr{i}",  rf"{lbl} AK8 regressed mass [GeV]",                  40,   0,   400),
        (f"ak8_massRes{i}",   rf"{lbl} AK8 resonance mass [GeV]",                  40,   0,   400),
    ]
] + [
    entry
    for cand, clbl in [("Hbb", r"$H \to bb$"), ("Htt", r"$H \to \tau_h\tau_h$"),
                        ("Htm", r"$H \to \tau_\mu\tau_h$"), ("Hte", r"$H \to \tau_e\tau_h$")]
    for entry in [
        (f"ak8_{cand}_pt",    rf"{clbl} cand AK8 $p_T$ [GeV]",       50, 150, 1000),
        (f"ak8_{cand}_eta",   rf"{clbl} cand AK8 $\eta$",             30,  -5,    5),
        (f"ak8_{cand}_mass",  rf"{clbl} cand AK8 mass [GeV]",         40,   0,  400),
        (f"ak8_{cand}_msd",   rf"{clbl} cand AK8 $m_{{SD}}$ [GeV]",   40,   0,  300),
        (f"ak8_{cand}_score", rf"{clbl} cand AK8 score",               25,   0,    1),
    ]
] + [
    ("ak8_dR_Hbb_Htt", r"$\Delta R(H_{bb}, H_{\tau\tau})$ cand AK8", 30, 0, 6),
    # -- Gen-matched AK8 tagger scores (signal only, bkg = empty) --
    ("ak8_genHbb_match_dR",       r"$\Delta R$(gen $H_{bb}$, matched AK8)",                    25, 0, 1),
    ("ak8_genHbb_match_Xbb",      r"Gen-matched $H_{bb}$ AK8 raw $X_{bb}$",                   25, 0, 1),
    ("ak8_genHbb_match_XbbVsAll", r"Gen-matched $H_{bb}$ AK8 $X_{bb}$ vs All",                25, 0, 1),
    ("ak8_genHtt_match_dR",       r"$\Delta R$(gen $H_{\tau\tau}$, matched AK8)",              25, 0, 1),
    ("ak8_genHtt_match_Xtt",      r"Gen-matched $H_{\tau\tau}$ AK8 raw $X_{\tau_h\tau_h}$",   25, 0, 1),
    ("ak8_genHtt_match_XttVsAll", r"Gen-matched $H_{\tau\tau}$ AK8 $X_{\tau_h\tau_h}$ vs All", 25, 0, 1),
    # -- Gen-matched AK4 UParT scores (signal only, bkg = empty) --
    ("ak4_genb1_dR",         r"$\Delta R$(gen $b_1$, AK4)",                  25, 0, 0.5),
    ("ak4_genb1_BvsAll",     r"Gen-matched $b_1$ AK4 UParT BvsAll",         25, 0, 1),
    ("ak4_genb2_dR",         r"$\Delta R$(gen $b_2$, AK4)",                  25, 0, 0.5),
    ("ak4_genb2_BvsAll",     r"Gen-matched $b_2$ AK4 UParT BvsAll",         25, 0, 1),
    ("ak4_gentau1_dR",       r"$\Delta R$(gen $\tau_1$, AK4)",              25, 0, 0.5),
    ("ak4_gentau1_TauVsAll",  r"Gen-matched $\tau_1$ AK4 UParT TauVsAll",          25, 0, 1),
    ("ak4_gentau1_TaumVsAll", r"Gen-matched $\tau_1$ AK4 UParT $\tau_\mu$ VsAll", 25, 0, 1),
    ("ak4_gentau1_TaupVsAll", r"Gen-matched $\tau_1$ AK4 UParT $\tau_h$ VsAll",   25, 0, 1),
    ("ak4_gentau2_dR",        r"$\Delta R$(gen $\tau_2$, AK4)",                    25, 0, 0.5),
    ("ak4_gentau2_TauVsAll",  r"Gen-matched $\tau_2$ AK4 UParT TauVsAll",          25, 0, 1),
    ("ak4_gentau2_TaumVsAll", r"Gen-matched $\tau_2$ AK4 UParT $\tau_\mu$ VsAll", 25, 0, 1),
    ("ak4_gentau2_TaupVsAll", r"Gen-matched $\tau_2$ AK4 UParT $\tau_h$ VsAll",   25, 0, 1),
    # -- Gen-matched raw UParT probs --
    ("ak4_genb1_raw_b",       r"Gen-matched $b_1$ AK4 raw prob\_b",               25, 0, 1),
    ("ak4_genb2_raw_b",       r"Gen-matched $b_2$ AK4 raw prob\_b",               25, 0, 1),
    ("ak4_gentau1_raw_b",     r"Gen-matched $\tau_1$ AK4 raw prob\_b",             25, 0, 1),
    ("ak4_gentau1_raw_taup",  r"Gen-matched $\tau_1$ AK4 raw prob $\tau_h^+$",    25, 0, 1),
    ("ak4_gentau1_raw_taum",  r"Gen-matched $\tau_1$ AK4 raw prob $\tau_h^-$",    25, 0, 1),
    ("ak4_gentau1_raw_tau",   r"Gen-matched $\tau_1$ AK4 raw prob $\tau$ (p+m)",  25, 0, 1),
    ("ak4_gentau2_raw_b",     r"Gen-matched $\tau_2$ AK4 raw prob\_b",             25, 0, 1),
    ("ak4_gentau2_raw_taup",  r"Gen-matched $\tau_2$ AK4 raw prob $\tau_h^+$",    25, 0, 1),
    ("ak4_gentau2_raw_taum",  r"Gen-matched $\tau_2$ AK4 raw prob $\tau_h^-$",    25, 0, 1),
    ("ak4_gentau2_raw_tau",   r"Gen-matched $\tau_2$ AK4 raw prob $\tau$ (p+m)",  25, 0, 1),
    ("ak4_genb1_raw_taup",    r"Gen-matched $b_1$ AK4 raw prob $\tau_h^+$",       25, 0, 1),
    ("ak4_genb1_raw_taum",    r"Gen-matched $b_1$ AK4 raw prob $\tau_h^-$",       25, 0, 1),
    ("ak4_genb1_raw_tau",     r"Gen-matched $b_1$ AK4 raw prob $\tau$ (p+m)",     25, 0, 1),
    ("ak4_genb2_raw_taup",    r"Gen-matched $b_2$ AK4 raw prob $\tau_h^+$",       25, 0, 1),
    ("ak4_genb2_raw_taum",    r"Gen-matched $b_2$ AK4 raw prob $\tau_h^-$",       25, 0, 1),
    ("ak4_genb2_raw_tau",     r"Gen-matched $b_2$ AK4 raw prob $\tau$ (p+m)",     25, 0, 1),
    # -- Gen-matched raw UParT probs (c, g, uds, lepb) --
    ("ak4_genb1_raw_c",       r"Gen-matched $b_1$ AK4 raw prob\_c",               25, 0, 1),
    ("ak4_genb1_raw_g",       r"Gen-matched $b_1$ AK4 raw prob\_g",               25, 0, 1),
    ("ak4_genb1_raw_uds",     r"Gen-matched $b_1$ AK4 raw prob\_uds",             25, 0, 1),
    ("ak4_genb1_raw_lepb",    r"Gen-matched $b_1$ AK4 raw prob\_lepb",            25, 0, 1),
    ("ak4_genb2_raw_c",       r"Gen-matched $b_2$ AK4 raw prob\_c",               25, 0, 1),
    ("ak4_genb2_raw_g",       r"Gen-matched $b_2$ AK4 raw prob\_g",               25, 0, 1),
    ("ak4_genb2_raw_uds",     r"Gen-matched $b_2$ AK4 raw prob\_uds",             25, 0, 1),
    ("ak4_genb2_raw_lepb",    r"Gen-matched $b_2$ AK4 raw prob\_lepb",            25, 0, 1),
    ("ak4_gentau1_raw_c",     r"Gen-matched $\tau_1$ AK4 raw prob\_c",            25, 0, 1),
    ("ak4_gentau1_raw_g",     r"Gen-matched $\tau_1$ AK4 raw prob\_g",            25, 0, 1),
    ("ak4_gentau1_raw_uds",   r"Gen-matched $\tau_1$ AK4 raw prob\_uds",          25, 0, 1),
    ("ak4_gentau1_raw_lepb",  r"Gen-matched $\tau_1$ AK4 raw prob\_lepb",         25, 0, 1),
    ("ak4_gentau2_raw_c",     r"Gen-matched $\tau_2$ AK4 raw prob\_c",            25, 0, 1),
    ("ak4_gentau2_raw_g",     r"Gen-matched $\tau_2$ AK4 raw prob\_g",            25, 0, 1),
    ("ak4_gentau2_raw_uds",   r"Gen-matched $\tau_2$ AK4 raw prob\_uds",          25, 0, 1),
    ("ak4_gentau2_raw_lepb",  r"Gen-matched $\tau_2$ AK4 raw prob\_lepb",         25, 0, 1),
]

# Gen-level Higgs plots (signal only, no trigger dependence at gen level)
# (branch, xlabel, nbins, xmin, xmax)
GEN_PLOT_VARS = [
    # -- Per-Higgs kinematics --
    ("gen_pt_Hbb",       r"Gen $H \to bb$ $p_T$ [GeV]",          50,    0,   500),
    ("gen_pt_Htautau",   r"Gen $H \to \tau\tau$ $p_T$ [GeV]",    50,    0,   500),
    ("gen_eta_Hbb",      r"Gen $H \to bb$ $\eta$",               30,   -5,     5),
    ("gen_eta_Htautau",  r"Gen $H \to \tau\tau$ $\eta$",         30,   -5,     5),
    ("gen_phi_Hbb",      r"Gen $H \to bb$ $\phi$",               30, -3.2,   3.2),
    ("gen_phi_Htautau",  r"Gen $H \to \tau\tau$ $\phi$",         30, -3.2,   3.2),
    ("gen_E_Hbb",        r"Gen $H \to bb$ Energy [GeV]",         50,    0,  1000),
    ("gen_E_Htautau",    r"Gen $H \to \tau\tau$ Energy [GeV]",   50,    0,  1000),
    ("gen_mass_Hbb",     r"Gen $H \to bb$ mass [GeV]",           50,  100,   150),
    ("gen_mass_Htautau", r"Gen $H \to \tau\tau$ mass [GeV]",     50,    0,   150),
    # -- Di-Higgs system --
    ("gen_dR_H1H2",      r"Gen $\Delta R(H_{bb}, H_{\tau\tau})$", 30,   0,     6),
    ("gen_pt_HH",        r"Gen $p_T^{HH}$ [GeV]",               50,    0,   500),
    ("gen_eta_HH",       r"Gen $\eta^{HH}$",                     30,   -5,     5),
    ("gen_phi_HH",       r"Gen $\phi^{HH}$",                     30, -3.2,   3.2),
    ("gen_E_HH",         r"Gen $E^{HH}$ [GeV]",                  50,    0,  2000),
    ("gen_mHH",          r"Gen $m_{HH}$ [GeV]",                  50,  200,   800),
    # -- τ charge misidentification (τ₁+τ₂ combined, keyed by true gen charge) --
    # TauVsAll = TaumVsAll + TaupVsAll (charge-agnostic), range [0, 2]
    ("gen_taum_TaumVsAll", r"Gen $\tau^-$ AK4 $\tau_h^-$ vs All (correct sign)",  25, 0, 1),
    ("gen_taum_TaupVsAll", r"Gen $\tau^-$ AK4 $\tau_h^+$ vs All (wrong sign)",    25, 0, 1),
    ("gen_taum_TauVsAll",  r"Gen $\tau^-$ AK4 $\tau_h$ vs All (combined)",         25, 0, 1),
    ("gen_taup_TaupVsAll", r"Gen $\tau^+$ AK4 $\tau_h^+$ vs All (correct sign)",  25, 0, 1),
    ("gen_taup_TaumVsAll", r"Gen $\tau^+$ AK4 $\tau_h^-$ vs All (wrong sign)",    25, 0, 1),
    ("gen_taup_TauVsAll",  r"Gen $\tau^+$ AK4 $\tau_h$ vs All (combined)",         25, 0, 1),
]

# 2D histograms (signal only, gen-matched AK8 tagger studies)
# (name, xvar, yvar, xlabel, ylabel, nx, x0, x1, ny, y0, y1)
PLOT_VARS_2D = [
    ("gen2d_Hbb_dR_vs_Xbb",
     "ak8_genHbb_match_dR", "ak8_genHbb_match_XbbVsAll",
     r"$\Delta R$(gen $H_{bb}$, AK8)", r"$X_{bb}$ vs All",
     25, 0, 1.0, 25, 0, 1.0),
    ("gen2d_Htt_dR_vs_Xtt",
     "ak8_genHtt_match_dR", "ak8_genHtt_match_XttVsAll",
     r"$\Delta R$(gen $H_{\tau\tau}$, AK8)", r"$X_{\tau\tau}$ vs All",
     25, 0, 1.0, 25, 0, 1.0),
    ("gen2d_Hbb_Xbb_vs_Htt_Xtt",
     "ak8_genHbb_match_XbbVsAll", "ak8_genHtt_match_XttVsAll",
     r"$X_{bb}$ vs All (H$\to$bb cand)", r"$X_{\tau\tau}$ vs All (H$\to\tau\tau$ cand)",
     25, 0, 1.0, 25, 0, 1.0),
    ("gen2d_Hbb_pt_vs_Xbb",
     "gen_pt_Hbb", "ak8_genHbb_match_XbbVsAll",
     r"Gen $H \to bb$ $p_T$ [GeV]", r"$X_{bb}$ vs All",
     25, 150, 600, 25, 0, 1.0),
]

# ── Filter variables if --plot-vars given ──
if ARGS.plot_vars:
    _all_var_names = {v[0] for v in PLOT_VARS} | {v[0] for v in GEN_PLOT_VARS}
    PLOT_VARS = [v for v in PLOT_VARS
                 if any(fnmatch.fnmatch(v[0], p) for p in ARGS.plot_vars)]
    GEN_PLOT_VARS = [v for v in GEN_PLOT_VARS
                     if any(fnmatch.fnmatch(v[0], p) for p in ARGS.plot_vars)]
    _unknown = set()
    for p in ARGS.plot_vars:
        if not any(fnmatch.fnmatch(v, p) for v in _all_var_names):
            _unknown.add(p)
    if _unknown:
        print(f"WARNING: no variables matched patterns: {_unknown}")
    if not PLOT_VARS and not GEN_PLOT_VARS:
        print("ERROR: --plot-vars matched no variables")
        sys.exit(1)
    print(f"[--plot-vars] Plotting {len(PLOT_VARS)} main vars, "
          f"{len(GEN_PLOT_VARS)} gen vars: "
          f"{', '.join(v[0] for v in PLOT_VARS + GEN_PLOT_VARS)}")


def asimov_significance(S, B):
    """Asimov significance: Z_A = sqrt(2 * [(S+B)*ln(1 + S/B) - S])."""
    if B <= 0 or S <= 0:
        return float("nan")
    return math.sqrt(2 * ((S + B) * math.log(1 + S / B) - S))


do_stacked = ARGS.plot_type in ("stacked", "both")
do_shape   = ARGS.plot_type in ("shape", "both")


# ── Try loading from histogram cache ──
_cache_path = os.path.join(PLOT_DIR, ".hist_cache.root")
_plot_var_names = [v[0] for v in PLOT_VARS]
_trig_names = [tn for tn, _ in TRIG_LIST]
_expected_meta = _cache_meta(
    max_files=ARGS.max_mc_files,
    plot_vars_names=_plot_var_names,
)

_cache_loaded = False
if not ARGS.recache:
    _cached = load_hist_cache(_cache_path, _expected_meta, _trig_names)
    if _cached is not None:
        (mc_hists_by_trig, mc_denom_hists, mc_items_by_trig, LUMI,
         gen_hists_by_trig, mc_hists_excl_by_trig) = _cached
        # Build trig_selections with just the keys Phase 3 needs
        trig_selections = {}
        for trig_name, trig in TRIG_LIST:
            if trig is not None and data_acc is not None:
                _data_sel = data_acc.Filter(trig)
            else:
                _data_sel = data_acc
            trig_selections[trig_name] = {
                "trig": trig,
                "data_sel": _data_sel, "mc_sel": None,
                "hh": None, "hm": None, "he": None,
                "sig_mHH_ptrs": [],
            }
        cutflow_tables = {}
        _cache_loaded = True
        print(f"Loaded histogram cache: {_cache_path}")
        print(f"  LUMI = {LUMI:.4f} fb^-1, "
              f"{len(mc_hists_by_trig)} triggers, "
              f"{len(next(iter(mc_hists_by_trig.values())))} variables")

# ══════════════════════════════════════════════════════════════════════════════
# Skim — write slim ROOT files if requested (then exit)
# ══════════════════════════════════════════════════════════════════════════════
if ARGS.skim or ARGS.reslim:
    # Collect all expression strings from this script for branch auto-detection.
    # Rather than wrapping every .Define() call, we read the script source and
    # extract all string literals that were passed to Define/Filter. This works
    # because all expressions are Python string literals or f-strings that are
    # fully expanded by the time the script runs.
    import ast as _ast
    with open(__file__, "r") as _f:
        _source = _f.read()
    _tree = _ast.parse(_source)
    for _node in _ast.walk(_tree):
        if isinstance(_node, _ast.Call):
            func = _node.func
            # Match .Define("col", "expr") and .Filter("expr")
            if isinstance(func, _ast.Attribute) and func.attr in ("Define", "Filter"):
                for arg in _node.args:
                    if isinstance(arg, _ast.Constant) and isinstance(arg.value, str):
                        rdf_exprs.append(arg.value)

    # Also add trigger expressions (constructed at runtime, not in source literals)
    rdf_exprs.append(DST_JetHT_expr)
    rdf_exprs.append(DST_MU_expr)
    rdf_exprs.append(DST_EL_expr)
    rdf_exprs.append(PARKING_HH_expr)
    for _, excl_expr in EXCL_TRIG_LIST:
        if excl_expr is not None:
            rdf_exprs.append(excl_expr)

    # Also add f-string expressions that reference AK8 branches (expanded at runtime)
    _AK8_skim = "ScoutingFatPFJetRecluster"
    _SGP_skim = f"{_AK8_skim}_scoutGlobalParT"
    for _i in range(3):
        rdf_exprs.append(f"{_AK8_skim}_pt[{_i}]")
        rdf_exprs.append(f"{_AK8_skim}_eta[{_i}]")
        rdf_exprs.append(f"{_AK8_skim}_mass[{_i}]")
        rdf_exprs.append(f"{_AK8_skim}_msoftdrop[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_prob_Xbb[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_prob_Xtauhtauh[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_prob_Xtauhtaum[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_prob_Xtauhtaue[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_prob_QCD[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_massCorrGeneric[{_i}]")
        rdf_exprs.append(f"{_SGP_skim}_massCorrResonance[{_i}]")
    rdf_exprs.append(f"n{_AK8_skim}")
    # Gen branches
    rdf_exprs.append("GenPart_pdgId GenPart_genPartIdxMother GenPart_statusFlags")
    rdf_exprs.append("GenPart_pt GenPart_eta GenPart_phi GenPart_mass")

    first_name = next(iter(mc_files_map))
    branches = detect_used_branches(rdf_exprs, mc_files_map[first_name][0])
    mc_branches = [b for b in branches if b not in ("run", "luminosityBlock")]
    print(f"[skim] Auto-detected {len(mc_branches)} branches from {len(rdf_exprs)} expressions")
    for b in mc_branches:
        print(f"  {b}")

    run_skim(mc, mc_branches, mc_files_map=mc_files_map,
             force=getattr(ARGS, 'force', False))
    print("Skim complete. Re-run without --skim to use slim files.")
    sys.exit(0)

# ══════════════════════════════════════════════════════════════════════════════
# Phases 1+2 — Book ALL actions, then run ONE event loop
# ══════════════════════════════════════════════════════════════════════════════
#
# Previously Phase 1 (cutflow) and Phase 2 (histograms) each ran separate
# RunGraphs calls per trigger (7 event loops total).  Now we book everything
# first and call RunGraphs once — same result, ~4-5× faster on cold runs.
# Skipped entirely when histogram cache is valid.

if not _cache_loaded:
    cutflow_tables = {}
    trig_selections = {}   # store per-trigger objects for the plotting phase
    unified_ptrs = []      # accumulate ALL lazy actions here
    
    # Book lumi Take actions (lazy — evaluated inside the unified RunGraphs)
    if data_df is not None:
        _lumi_run_take, _lumi_ls_take = book_lumi_actions(data_df)
        unified_ptrs.extend([_lumi_run_take, _lumi_ls_take])
    else:
        _lumi_run_take, _lumi_ls_take = None, None
    
    # ── Phase 1: Book cutflow & yield actions (skip with --skip-cutflow) ──
    _phase1_data = {}  # trig_name -> {cutflow_ptrs, data_n_ptr, sum_w, sum_w2, ...}
    
    for trig_name, trig in TRIG_LIST:
        print(f"Booking Phase 1+2 actions for {trig_name} ...")
    
        # ── Final selection (needed for both cutflow yields AND Phase 2 histograms) ──
        if trig is not None:
            data_sel = data_acc.Filter(trig) if data_acc is not None else None
            mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}
        else:
            data_sel = data_acc
            mc_sel = dict(mc_acc)
    
        sig_mHH_ptrs = [mc_sel[s].Histo1D(
            (f"h_mHH_{s}_{trig_name}",
             "m_{HH} gen-level;m_{HH} [GeV];Events", 32, 0, 800),
            "gen_mHH", "w") for s in sig_samples]
        unified_ptrs.extend(sig_mHH_ptrs)

        # Store selections for the plotting phase
        trig_selections[trig_name] = {
            "trig": trig,
            "data_sel": data_sel, "mc_sel": mc_sel,
            "sig_mHH_ptrs": sig_mHH_ptrs,
        }
    
        if not ARGS.skip_cutflow:
            # ── Cutflow: progressive cuts ──
            if trig is not None:
                data_cf = data_df.Filter(trig) if data_df is not None else None
                mc_cf = {name: df.Filter(trig) for name, df in mc.items()}
            else:
                data_cf = data_df
                mc_cf = dict(mc)
    
            cutflow_ptrs = []
            for step_name, cut_expr in CUTFLOW_STEPS:
                if cut_expr is not None:
                    if data_cf is not None:
                        data_cf = data_cf.Filter(cut_expr)
                    mc_cf = {s: d.Filter(cut_expr) for s, d in mc_cf.items()}
                d_ptr = data_cf.Count() if data_cf is not None else None
                w_ptrs = {s: mc_cf[s].Sum("w") for s in mc_cf}
                cutflow_ptrs.append((step_name, d_ptr, w_ptrs))
                if d_ptr is not None:
                    unified_ptrs.append(d_ptr)
                unified_ptrs.extend(w_ptrs.values())
    
            data_n_ptr = data_sel.Count() if data_sel is not None else None
            sum_w  = {s: mc_sel[s].Sum("w") for s in mc_sel}
            sum_w2 = {s: mc_sel[s].Define("w2", "w*w").Sum("w2") for s in mc_sel}

            # Per-mode yield sums (for cutflow table)
            mode_w = {}  # mode_int -> [Sum ptrs]
            _proc_samples = {"DY": dy_samples, "TT": tt_samples,
                             "HHbbtt": sig_samples, "QCD": qcd_samples}
            for mode in ACTIVE_MODES:
                info = DECAY_MODES[mode]
                samples_for_mode = _proc_samples.get(info["process"], [])
                if not samples_for_mode:
                    mode_w[mode] = []
                    continue
                mode_w[mode] = [mc_sel[s].Filter(f"decayMode=={mode}").Sum("w")
                                for s in samples_for_mode]
                unified_ptrs.extend(mode_w[mode])

            if data_n_ptr is not None:
                unified_ptrs.append(data_n_ptr)
            unified_ptrs.extend(sum_w.values())
            unified_ptrs.extend(sum_w2.values())

            _phase1_data[trig_name] = {
                "cutflow_ptrs": cutflow_ptrs,
                "data_n_ptr": data_n_ptr,
                "sum_w": sum_w, "sum_w2": sum_w2,
                "mode_w": mode_w,
            }
    
    # ── Phase 2: Book MC denominator histograms (no trigger) ──
    # Build mc_denom_groups from global decay mode LUT
    _proc_to_samples = {"DY": dy_samples, "TT": tt_samples,
                         "HHbbtt": sig_samples, "QCD": qcd_samples}

    # Which modes actually have samples? (determines group ordering)
    _BUILT_MODES = [m for m in ACTIVE_MODES
                    if _proc_to_samples.get(DECAY_MODES[m]["process"], [])]

    def _build_groups(sample_dfs):
        """Build mc_groups list from ACTIVE_MODES using decayMode filter."""
        groups = []
        for mode in _BUILT_MODES:
            info = DECAY_MODES[mode]
            proc = info["process"]
            slist = _proc_to_samples[proc]
            if proc == "QCD":
                dfs = [sample_dfs[s] for s in slist]  # no decayMode filter
            else:
                dfs = [sample_dfs[s].Filter(f"decayMode=={mode}") for s in slist]
            groups.append((dfs, info["label"], info["color"]))
        return groups

    mc_denom_groups = _build_groups(mc_acc)
    
    # AK8 variables use -1 sentinel when no jet passes pT cut; filter those out
    _AK8_VARS = {v[0] for v in PLOT_VARS
                 if v[0].startswith("ak8_") or v[0].startswith("ak4_gen")}

    denom_book = {}
    for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
        denom_book[var_name] = []
        for gi, (dfs, _lbl, _col) in enumerate(mc_denom_groups):
            group_ptrs = []
            for si, df in enumerate(dfs):
                uid = f"denom_{var_name}_{gi}_{si}"
                df_v = df.Filter(f"{var_name} > -0.5f") if var_name in _AK8_VARS else df
                ptr = df_v.Histo1D(
                    (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                    var_name, "w")
                group_ptrs.append(ptr)
                unified_ptrs.append(ptr)
            denom_book[var_name].append(group_ptrs)
    
    # ── Phase 2: Book per-trigger MC histograms ──
    mc_items_by_trig = {}
    _histo_books = {}  # trig_name -> {var_name -> [[ptrs per group]]}
    
    # Derive sig/bkg indices from _BUILT_MODES (only modes with samples)
    sig_indices = [i for i, m in enumerate(_BUILT_MODES) if m in SIG_MODES]
    bkg_indices = [i for i, m in enumerate(_BUILT_MODES) if m in BKG_MODES]

    for trig_name in trig_selections:
        sel = trig_selections[trig_name]
        mc_sel = sel["mc_sel"]

        mc_groups = _build_groups(mc_sel)
        mc_items = [{"label": lbl, "color": col} for _, lbl, col in mc_groups]
        mc_items_by_trig[trig_name] = mc_items
    
        histo_book = {}
        for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
            histo_book[var_name] = []
            for gi, (dfs, _lbl, _col) in enumerate(mc_groups):
                group_ptrs = []
                for si, df in enumerate(dfs):
                    uid = f"{trig_name}_{var_name}_{gi}_{si}"
                    df_v = df.Filter(f"{var_name} > -0.5f") if var_name in _AK8_VARS else df
                    ptr = df_v.Histo1D(
                        (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                        var_name, "w")
                    group_ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                histo_book[var_name].append(group_ptrs)
        _histo_books[trig_name] = histo_book
    
    # ── Phase 2b: Book 2D histograms (signal only, gen-matched) ──
    h2d_book = {}
    for pname, xvar, yvar, _xl, _yl, nx, x0, x1, ny, y0, y1 in PLOT_VARS_2D:
        h2d_book[pname] = []
        for trig_name in trig_selections:
            mc_sel_t = trig_selections[trig_name]["mc_sel"]
            for s in sig_samples:
                df_v = mc_sel_t[s].Filter(f"{xvar} > -0.5f && {yvar} > -0.5f")
                uid = f"h2d_{pname}_{s}_{trig_name}"
                ptr = df_v.Histo2D(
                    ROOT.RDF.TH2DModel(uid, f";{xvar};{yvar}",
                                       nx, x0, x1, ny, y0, y1),
                    xvar, yvar, "w")
                h2d_book[pname].append((trig_name, s, ptr))
                unified_ptrs.append(ptr)

    # ── Phase 2c: Book gen-level histograms (signal only, split by decay channel) ──
    gen_histo_book = {}  # {trig_name: {var_name: {mode: [ptrs per sig sample]}}}
    for trig_name in trig_selections:
        mc_sel_t = trig_selections[trig_name]["mc_sel"]
        gen_histo_book[trig_name] = {}
        for var_name, _xlabel, nbins, vmin, vmax in GEN_PLOT_VARS:
            gen_histo_book[trig_name][var_name] = {}
            for mode in SIG_MODES:
                ptrs = []
                for s in sig_samples:
                    uid = f"gen_{trig_name}_{var_name}_{s}_m{mode}"
                    ptr = (mc_sel_t[s]
                           .Filter(f"decayMode=={mode}")
                           .Histo1D((uid, f";{var_name};Events", nbins, vmin, vmax),
                                    var_name, "w"))
                    ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                gen_histo_book[trig_name][var_name][mode] = ptrs

    # ── Phase 2d: Book exclusive-trigger MC histograms ──
    _excl_histo_books = {}
    for trig_name, trig in EXCL_TRIG_LIST:
        if trig is not None:
            mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}
        else:
            mc_sel = dict(mc_acc)
        mc_groups = _build_groups(mc_sel)
        histo_book = {}
        for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
            histo_book[var_name] = []
            for gi, (dfs, _lbl, _col) in enumerate(mc_groups):
                group_ptrs = []
                for si, df in enumerate(dfs):
                    uid = f"excl_{trig_name}_{var_name}_{gi}_{si}"
                    df_v = df.Filter(f"{var_name} > -0.5f") if var_name in _AK8_VARS else df
                    ptr = df_v.Histo1D(
                        (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                        var_name, "w")
                    group_ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                histo_book[var_name].append(group_ptrs)
        _excl_histo_books[trig_name] = histo_book

    # ══════════════════════════════════════════════════════════════════════════════
    # Run ONE unified event loop for all Phase 1 + Phase 2 actions
    # ══════════════════════════════════════════════════════════════════════════════

    print(f"\nBooked {len(unified_ptrs)} total actions (Phase 1+2 combined)")
    print(f"Running unified event loop ...", flush=True)
    sys.stdout.flush()
    t0 = time.time()
    ROOT.RDF.RunGraphs(unified_ptrs)
    dt = time.time() - t0
    print(f"Unified event loop done in {dt:.1f}s")
    
    # ── Extract LUMI from Take vectors ──
    if LUMI is None and _lumi_run_take is not None:
        print("  Extracting lumi from Take vectors ...", flush=True)
        t0_lumi = time.time()
        LUMI, data_runs, n_missing = extract_lumi(
            _lumi_run_take, _lumi_ls_take, _brilcalc_default)
        dt_lumi = time.time() - t0_lumi
        print(f"  Lumi: {len(data_runs)} runs, {n_missing} not in brilcalc, "
              f"extracted in {dt_lumi:.1f}s")
        print(f"  Data lumi (LS-matched) = {LUMI:.4f} fb^-1")
    
    # ── Extract Phase 1 results (cutflow + yields) ──
    if not ARGS.skip_cutflow:
        for trig_name, trig in TRIG_LIST:
            p1 = _phase1_data[trig_name]
            print("=" * 63)
            print(f"  {trig_name}")
            print("=" * 63)
    
            cutflow_rows = []
            for step_name, d_ptr, w_ptrs in p1["cutflow_ptrs"]:
                data_n_cf = d_ptr.GetValue() if d_ptr is not None else -1
                dy_cf  = sum(w_ptrs[s].GetValue() for s in dy_samples) * LUMI
                tt_cf  = sum(w_ptrs[s].GetValue() for s in tt_samples) * LUMI
                qcd_cf = sum(w_ptrs[s].GetValue() for s in qcd_samples) * LUMI if qcd_samples else 0.0
                sig_cf = sum(w_ptrs[s].GetValue() for s in sig_samples) * LUMI
                B_cf = dy_cf + tt_cf + qcd_cf
                s_sqrtb = sig_cf / math.sqrt(B_cf) if B_cf > 0 else float("nan")
                z_a = asimov_significance(sig_cf, B_cf)
                cutflow_rows.append(
                    (step_name, data_n_cf, dy_cf, tt_cf, qcd_cf, sig_cf, s_sqrtb, z_a))
            cutflow_tables[trig_name] = cutflow_rows
    
            print(f"\n{'Cut':<30s} {'Data':>12s} {'DY':>14s} {'TT':>14s} "
                  f"{'QCD':>14s} {'Signal':>14s} {'S/√B':>10s} {'Z_A':>10s}")
            print("-" * 122)
            for row in cutflow_rows:
                step, dn, dy, tt, qcd, sig, sb, za = row
                dn_str = f"{dn:>12d}" if dn >= 0 else f"{'(no data)':>12s}"
                print(f"{step:<30s} {dn_str} {dy:>14.4g} {tt:>14.4g} "
                      f"{qcd:>14.4g} {sig:>14.4g} {sb:>10.4g} {za:>10.4g}")
    
            data_n_ptr = p1["data_n_ptr"]
            sum_w, sum_w2 = p1["sum_w"], p1["sum_w2"]
            mode_w = p1["mode_w"]
            data_n = data_n_ptr.GetValue() if data_n_ptr is not None else -1

            def group_yield(samples):
                y = sum(sum_w[s].GetValue() for s in samples) * LUMI
                yerr = math.sqrt(sum(sum_w2[s].GetValue() for s in samples)) * LUMI
                return y, yerr

            dy_y, dy_yerr = group_yield(dy_samples)
            tt_y, tt_yerr = group_yield(tt_samples)
            qcd_y, qcd_yerr = group_yield(qcd_samples) if qcd_samples else (0.0, 0.0)
            sig_y, sig_yerr = group_yield(sig_samples)
            B = dy_y + tt_y + qcd_y
            S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")

            print(f"\nData events: {data_n if data_n >= 0 else '(no data)'}")
            print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
            print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
            print(f"QCD yield:   {qcd_y:.6g} +/- {qcd_yerr:.3g}")
            print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")
            print(f"S/sqrt(B):   {S_over_sqrtB:.6g}")
            # Per-mode yields
            for mode in ACTIVE_MODES:
                y_mode = sum(d.GetValue() for d in mode_w[mode]) * LUMI
                print(f"  mode {mode:>2d} ({DECAY_MODES[mode]['label']:>25s}): {y_mode:.6g}")
            print()
    
        # ── Write cutflow markdown ──
        cutflow_path = os.path.join(PLOT_DIR, "cutflow.md")
        with open(cutflow_path, "w") as f:
            f.write("# Cutflow Tables\n\n")
            f.write(f"Luminosity: {LUMI} fb$^{{-1}}$\n\n")
            for trig_name, rows in cutflow_tables.items():
                f.write(f"## {trig_name}\n\n")
                f.write("| Cut | Data | DY | TT | QCD | Signal "
                        "| S/sqrt(B) | Z_A (Asimov) |\n")
                f.write("|-----|-----:|---:|---:|----:|-------:"
                        "|----------:|-------------:|\n")
                for step, dn, dy, tt, qcd, sig, sb, za in rows:
                    f.write(f"| {step} | {dn:,d} | {dy:.4g} | {tt:.4g} "
                            f"| {qcd:.4g} | {sig:.4g} | {sb:.4g} | {za:.4g} |\n")
                f.write("\n")
        print(f"Cutflow table saved: {cutflow_path}")
    else:
        print("[--skip-cutflow] Skipping cutflow extraction")
    
    # ── Materialize Phase 2 histograms ──
    
    # Denominator histograms
    mc_denom_hists = {}
    for var_name in denom_book:
        h_list = []
        for gi in range(len(mc_denom_groups)):
            ptrs = denom_book[var_name][gi]
            combined = ptrs[0].GetValue().Clone()
            for ptr in ptrs[1:]:
                combined.Add(ptr.GetValue())
            combined.Scale(LUMI)
            h_list.append(combined)
        mc_denom_hists[var_name] = h_list
    
    # Per-trigger histograms
    mc_hists_by_trig = {}
    for trig_name in trig_selections:
        histo_book = _histo_books[trig_name]
        mc_groups_info = mc_items_by_trig[trig_name]
        h_by_var = {}
        for var_name in histo_book:
            h_mc_list = []
            for gi in range(len(mc_groups_info)):
                ptrs = histo_book[var_name][gi]
                combined = ptrs[0].GetValue().Clone()
                for ptr in ptrs[1:]:
                    combined.Add(ptr.GetValue())
                combined.Scale(LUMI)
                h_mc_list.append(combined)
            h_by_var[var_name] = h_mc_list
        mc_hists_by_trig[trig_name] = h_by_var

    # Gen-level histograms (signal only, split by decay channel)
    gen_hists_by_trig = {}  # {trig_name: {var_name: {mode: TH1D}}}
    for trig_name in trig_selections:
        h_by_var = {}
        for var_name in gen_histo_book[trig_name]:
            h_by_mode = {}
            for mode in SIG_MODES:
                ptrs = gen_histo_book[trig_name][var_name][mode]
                combined = ptrs[0].GetValue().Clone(f"gen_{trig_name}_{var_name}_m{mode}")
                for ptr in ptrs[1:]:
                    combined.Add(ptr.GetValue())
                combined.Scale(LUMI)
                h_by_mode[mode] = combined
            h_by_var[var_name] = h_by_mode
        gen_hists_by_trig[trig_name] = h_by_var

    # Exclusive-trigger histograms
    mc_hists_excl_by_trig = {}
    for trig_name, histo_book in _excl_histo_books.items():
        h_by_var = {}
        for var_name, group_list in histo_book.items():
            h_mc_list = []
            for gi, ptrs in enumerate(group_list):
                combined = ptrs[0].GetValue().Clone()
                for ptr in ptrs[1:]:
                    combined.Add(ptr.GetValue())
                combined.Scale(LUMI)
                h_mc_list.append(combined)
            h_by_var[var_name] = h_mc_list
        mc_hists_excl_by_trig[trig_name] = h_by_var

    # ── Save histogram cache ──
    save_hist_cache(_cache_path, mc_hists_by_trig, mc_denom_hists,
                    mc_items_by_trig, LUMI, _expected_meta,
                    gen_hists=gen_hists_by_trig,
                    excl_hists=mc_hists_excl_by_trig)

    # Drop lazy action pointers — they hold references to internal ROOT objects
    # that can cause double-free after RunGraphs with ImplicitMT.
    del unified_ptrs, _histo_books, _excl_histo_books, denom_book, gen_histo_book
    if '_phase1_data' in dir():
        del _phase1_data
    import gc; gc.collect()


# ══════════════════════════════════════════════════════════════════════════════
# Phase 3 — Draw MC-only plots  (no event loop, just drawing)
# ══════════════════════════════════════════════════════════════════════════════

for trig_name in trig_selections:
    sel = trig_selections[trig_name]
    h_by_var = mc_hists_by_trig[trig_name]
    mc_items = mc_items_by_trig[trig_name]

    # Comparison: gen mHH vs signal m4j (MC only, data added in Phase 5)
    mHH_path = _plot_path("mc", "shape", trig_name, "mHH_m4j")
    if (sel["sig_mHH_ptrs"] and "m4j" in h_by_var and sig_indices
            and (not os.path.exists(mHH_path) or ARGS.overwrite)):
        fig, ax = plt.subplots(figsize=(8, 6))

        # 1) gen mHH (combine signal samples)
        h_gen = sel["sig_mHH_ptrs"][0].GetPtr().Clone("h_gen_mHH_combined")
        for h_ptr in sel["sig_mHH_ptrs"][1:]:
            h_gen.Add(h_ptr.GetPtr())
        edges, vals, _ = th1_to_np(h_gen)
        area = float(np.sum(vals * np.diff(edges)))
        if area > 0:
            vals = vals / area
        ax.step(edges, np.r_[vals, vals[-1]], where="post",
                linewidth=2, color="tab:blue", label=r"$m_{HH}$ gen-level (signal)")

        # 2) signal m4j (reco) — combine all signal channels
        h_sig_m4j = h_by_var["m4j"][sig_indices[0]].Clone("h_sig_m4j_combined")
        for si in sig_indices[1:]:
            h_sig_m4j.Add(h_by_var["m4j"][si])
        edges2, vals2, _ = th1_to_np(h_sig_m4j)
        area2 = float(np.sum(vals2 * np.diff(edges2)))
        if area2 > 0:
            vals2 = vals2 / area2
        ax.step(edges2, np.r_[vals2, vals2[-1]], where="post",
                linewidth=2, color="tab:red", label=r"$m_{4j}$ reco (signal)")

        ax.set_xlabel(r"Mass [GeV]")
        ax.set_ylabel("Normalised to unity")
        ax.legend()
        ax.grid(True, alpha=0.3)
        cms_label(ax, lumi=LUMI)
        fig.tight_layout()
        fig.savefig(mHH_path, dpi=150)
        print(f"Saved: {mHH_path}")
        plt.close(fig)
    elif os.path.exists(mHH_path) and not ARGS.overwrite:
        print(f"Skipping (exists): {mHH_path}")

    # MC-only stacked & shape plots
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        h_mc_list = h_by_var[var_name]

        if do_stacked:
            outpath = _plot_path("mc", "stacked", trig_name, var_name)
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax = plot_stacked_from_hists(
                    h_mc_list, mc_items, h_data=None,
                    logy=True, title=None)
                ax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

        if do_shape:
            outpath = _plot_path("mc", "shape", trig_name, var_name)
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax = plot_shape_from_hists(
                    h_mc_list, mc_items, title=None)
                ax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

        # Stacked + trigger efficiency panel (MC only) — skip for NoTrigger
        if do_stacked and sel["trig"] is not None:
            outpath = _plot_path("mc", "eff", trig_name, var_name)
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                h_denom_list = mc_denom_hists.get(var_name)
                fig, ax, rax = plot_stacked_with_efficiency(
                    h_mc_list, mc_items,
                    h_mc_denom_list=h_denom_list,
                    logy=True, title=None)
                rax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.tight_layout()
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

        # Stacked + S/sqrt(B) panel (MC only)
        if do_stacked:
            outpath = _plot_path("mc", "sig", trig_name, var_name)
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax, rax = plot_stacked_with_significance(
                    h_mc_list, mc_items,
                    sig_indices=sig_indices, bkg_indices=bkg_indices,
                    logy=True, title=None)
                rax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.tight_layout()
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

    print(f"[{trig_name}] Phase 3: MC-only plots done")

    # ── Gen-level plots (signal only, split by decay channel) ──
    if trig_name in gen_hists_by_trig:
        for var_name, xlabel, nbins, vmin, vmax in GEN_PLOT_VARS:
            gen_dir = os.path.join(PLOT_DIR, "mc", "shape", trig_name,
                                   "gen", "HBBTT_SM")
            os.makedirs(gen_dir, exist_ok=True)
            outpath = os.path.join(gen_dir, f"{var_name}.png")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
                continue
            h_by_mode = gen_hists_by_trig[trig_name][var_name]
            fig, ax = plt.subplots(figsize=(8, 6))
            has_entries = False
            for mode in SIG_MODES:
                h = h_by_mode.get(mode)
                if h is None or h.GetEntries() == 0:
                    continue
                has_entries = True
                edges, vals, errs = th1_to_np(h)
                # Normalise to unit area (same logic as plot_trigger_shape_overlay)
                area = float(np.sum(vals * np.diff(edges)))
                if area > 0:
                    vals = vals / area
                    errs = errs / area
                color = DECAY_MODES[mode]["color"]
                label = DECAY_MODES[mode]["label"]
                ax.step(edges, np.r_[vals, vals[-1]], where="post",
                        linewidth=2, color=color, label=label)
                ax.fill_between(edges,
                                np.r_[vals - errs, (vals - errs)[-1]],
                                np.r_[vals + errs, (vals + errs)[-1]],
                                step="post", alpha=0.2, color=color)
            if not has_entries:
                plt.close(fig)
                continue
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Normalised to unity")
            ax.set_xlim(vmin, vmax)
            ax.legend()
            ax.grid(True, alpha=0.3)
            cms_label(ax, lumi=LUMI)
            fig.tight_layout()
            fig.savefig(outpath, dpi=150)
            print(f"Saved: {outpath}")
            plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Phase 3.5 — Trigger shape overlay  (one plot per variable, all triggers)
# ══════════════════════════════════════════════════════════════════════════════

if do_shape:
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        outpath = _plot_path("mc", "overlay", var_name=var_name)
        if os.path.exists(outpath) and not ARGS.overwrite:
            print(f"Skipping (exists): {outpath}")
            continue

        # Sum all MC groups into one total TH1 per trigger
        h_total_by_trig = {}
        for trig_name in trig_selections:
            h_list = mc_hists_by_trig[trig_name][var_name]
            h_total = h_list[0].Clone(f"_overlay_{trig_name}_{var_name}")
            for h in h_list[1:]:
                h_total.Add(h)
            h_total_by_trig[trig_name] = h_total

        fig, ax = plot_trigger_shape_overlay(h_total_by_trig)
        ax.set_xlabel(xlabel)
        cms_label(ax, lumi=LUMI)
        fig.tight_layout()
        fig.savefig(outpath, dpi=150)
        print(f"Saved: {outpath}")
        plt.close(fig)

    # Per-signal-channel trigger overlays
    SIG_CHANNELS = []
    if len(sig_indices) >= 3:
        SIG_CHANNELS = [
            (sig_indices[0], "hh", r"$HH\to bb\tau_h\tau_h$"),
            (sig_indices[1], "hm", r"$HH\to bb\tau_\mu\tau_h$"),
            (sig_indices[2], "he", r"$HH\to bb\tau_e\tau_h$"),
        ]
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        for gi, ch_key, ch_label in SIG_CHANNELS:
            outpath = _plot_path("mc", "overlay", var_name=f"{ch_key}_{var_name}")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
                continue

            h_by_trig = {}
            for trig_name in trig_selections:
                h_by_trig[trig_name] = mc_hists_by_trig[trig_name][var_name][gi]

            fig, ax = plot_trigger_shape_overlay(h_by_trig, title=ch_label)
            ax.set_xlabel(xlabel)
            cms_label(ax, lumi=LUMI)
            fig.tight_layout()
            fig.savefig(outpath, dpi=150)
            print(f"Saved: {outpath}")
            plt.close(fig)

    print("Phase 3.5: Trigger overlay plots done")

    # ── Gen-level overlay: combined → overlay/gen/{var}.png ──
    for var_name, xlabel, nbins, vmin, vmax in GEN_PLOT_VARS:
        if not gen_hists_by_trig or var_name not in next(iter(gen_hists_by_trig.values()), {}):
            continue
        outpath = _plot_path("mc", "overlay", trig_name="gen", var_name=var_name)
        if os.path.exists(outpath) and not ARGS.overwrite:
            print(f"Skipping (exists): {outpath}")
            continue
        # Sum all decay modes into one TH1 per trigger (mirrors AK4 group-summing)
        h_total_by_trig = {}
        for trig_name in trig_selections:
            h_by_mode = gen_hists_by_trig.get(trig_name, {}).get(var_name, {})
            vals = list(h_by_mode.values())
            if not vals:
                continue
            h_total = vals[0].Clone(f"_gen_ov_{trig_name}_{var_name}")
            for h in vals[1:]:
                h_total.Add(h)
            h_total_by_trig[trig_name] = h_total
        if h_total_by_trig:
            fig, ax = plot_trigger_shape_overlay(h_total_by_trig)
            ax.set_xlabel(xlabel)
            cms_label(ax, lumi=LUMI)
            fig.tight_layout()
            fig.savefig(outpath, dpi=150)
            print(f"Saved: {outpath}")
            plt.close(fig)

    # ── Gen-level overlay: per-channel → overlay/gen/{ch}_{var}.png ──
    for var_name, xlabel, nbins, vmin, vmax in GEN_PLOT_VARS:
        if not gen_hists_by_trig or var_name not in next(iter(gen_hists_by_trig.values()), {}):
            continue
        for mode in SIG_MODES:
            ch_key = {20: "hh", 21: "hm", 22: "he"}.get(mode, str(mode))
            outpath = _plot_path("mc", "overlay", trig_name="gen",
                                 var_name=f"{ch_key}_{var_name}")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
                continue
            h_by_trig = {}
            for trig_name in trig_selections:
                h = gen_hists_by_trig.get(trig_name, {}).get(var_name, {}).get(mode)
                if h and h.GetEntries() > 0:
                    h_by_trig[trig_name] = h
            if h_by_trig:
                fig, ax = plot_trigger_shape_overlay(
                    h_by_trig, title=DECAY_MODES[mode]["label"])
                ax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.tight_layout()
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)
    print("Phase 3.5: Gen-level overlay plots done")

    # ── Exclusive-trigger overlay: overlay_exclusive/{var}.png ──
    if mc_hists_excl_by_trig:
        for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
            outpath = _plot_path("mc", "overlay_exclusive", var_name=var_name)
            if not os.path.exists(outpath) or ARGS.overwrite:
                h_by_trig = {}
                for trig_name, h_by_var in mc_hists_excl_by_trig.items():
                    h_list = h_by_var.get(var_name, [])
                    if not h_list:
                        continue
                    h_total = h_list[0].Clone(f"_exclov_{trig_name}_{var_name}")
                    for h in h_list[1:]:
                        h_total.Add(h)
                    h_by_trig[trig_name] = h_total
                if h_by_trig:
                    fig, ax = plot_trigger_shape_overlay(h_by_trig)
                    ax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.tight_layout()
                    fig.savefig(outpath, dpi=150)
                    plt.close(fig)
        # Per-channel exclusive overlay
        for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
            for gi, ch_key, ch_label in SIG_CHANNELS:
                outpath = _plot_path("mc", "overlay_exclusive",
                                     var_name=f"{ch_key}_{var_name}")
                if os.path.exists(outpath) and not ARGS.overwrite:
                    continue
                h_by_trig = {}
                for trig_name, h_by_var in mc_hists_excl_by_trig.items():
                    h_list = h_by_var.get(var_name, [])
                    if gi < len(h_list):
                        h_by_trig[trig_name] = h_list[gi]
                if h_by_trig:
                    fig, ax = plot_trigger_shape_overlay(h_by_trig, title=ch_label)
                    ax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.tight_layout()
                    fig.savefig(outpath, dpi=150)
                    plt.close(fig)
        print("Phase 3.5: Exclusive-trigger overlay plots done")

    # ── Phase 3.6: 2D gen-matched plots (signal only) ──
    for pname, _xvar, _yvar, xlabel, ylabel, *_ in PLOT_VARS_2D:
        # Group pointers by trigger, sum over signal samples
        trigs_seen = {}
        for trig_name, s, ptr in h2d_book.get(pname, []):
            h2 = ptr.GetPtr().Clone()
            h2.Scale(LUMI)
            if trig_name not in trigs_seen:
                trigs_seen[trig_name] = h2
            else:
                trigs_seen[trig_name].Add(h2)

        for trig_name, h2_combined in trigs_seen.items():
            outpath = _plot_path("mc", "gen2d", trig_name, var_name=f"{pname}")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
                continue
            fig = plot_2d_hist(h2_combined, xlabel=xlabel, ylabel=ylabel)
            fig.savefig(outpath, dpi=150)
            plt.close(fig)
            print(f"Saved: {outpath}")

    print("Phase 3.6: 2D gen-matched plots done")


# ══════════════════════════════════════════════════════════════════════════════
# Phase 4 — Book data histograms  (RunGraphs on XCache — slow)
# ══════════════════════════════════════════════════════════════════════════════

data_hists_by_trig = {}  # trig -> {var_name: TH1}
data_denom_hists = {}    # var_name -> TH1 (no trigger, acceptance only)

if not ARGS.no_data:
    # ── Book data denominator histograms (no trigger, acceptance only) ──
    data_denom_ptrs = []
    data_denom_book = {}
    for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
        if var_name == "gen_mHH":
            continue
        uid_d = f"data_denom_{var_name}"
        d_ptr = data_acc.Histo1D(
            (f"h_{uid_d}", f";{var_name};Events", nbins, vmin, vmax),
            var_name)
        data_denom_book[var_name] = d_ptr
        data_denom_ptrs.append(d_ptr)

    for trig_name in trig_selections:
        sel = trig_selections[trig_name]
        data_sel = sel["data_sel"]

        data_histo_ptrs = []
        data_histo_book = {}

        for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
            if var_name == "gen_mHH":
                continue  # gen-only, no data
            uid_d = f"{trig_name}_{var_name}_data"
            d_ptr = data_sel.Histo1D(
                (f"h_{uid_d}", f";{var_name};Events", nbins, vmin, vmax),
                var_name)
            data_histo_book[var_name] = d_ptr
            data_histo_ptrs.append(d_ptr)

        # Include denominator histos in first trigger's RunGraphs
        if not data_denom_hists:
            data_histo_ptrs.extend(data_denom_ptrs)

        print(f"\n[{trig_name}] Phase 4: Booked {len(data_histo_ptrs)} data histograms")
        print(f"[{trig_name}] Running data event loop (XCache) ...", flush=True)
        t0 = time.time()
        ROOT.RDF.RunGraphs(data_histo_ptrs)
        dt = time.time() - t0
        print(f"[{trig_name}] Data event loop done in {dt:.1f}s")

        # Materialize denominator histos (once)
        if not data_denom_hists:
            for var_name, d_ptr in data_denom_book.items():
                data_denom_hists[var_name] = d_ptr.GetValue()

        # Materialize trigger-filtered data histos
        h_data_var = {}
        for var_name, d_ptr in data_histo_book.items():
            h_data_var[var_name] = d_ptr.GetValue()
        data_hists_by_trig[trig_name] = h_data_var
else:
    print("\n[--no-data] Skipping Phase 4 & 5 (data histograms and MC+Data plots)")


# ══════════════════════════════════════════════════════════════════════════════
# Phase 5 — Draw MC+Data plots  (reuse Phase 2 MC + Phase 4 data, no event loop)
# ══════════════════════════════════════════════════════════════════════════════

if not ARGS.no_data:
    for trig_name in trig_selections:
        h_by_var = mc_hists_by_trig[trig_name]
        h_data_var = data_hists_by_trig[trig_name]
        mc_items = mc_items_by_trig[trig_name]

        for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
            h_mc_list = h_by_var[var_name]
            h_data = h_data_var.get(var_name, None)
            if h_data is None:
                continue

            if do_stacked:
                outpath = _plot_path("data", "stacked", trig_name, var_name)
                if os.path.exists(outpath) and not ARGS.overwrite:
                    print(f"Skipping (exists): {outpath}")
                else:
                    fig, ax = plot_stacked_from_hists(
                        h_mc_list, mc_items, h_data=h_data,
                        logy=True, title=None)
                    ax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)

            if do_shape:
                outpath = _plot_path("data", "shape", trig_name, var_name)
                if os.path.exists(outpath) and not ARGS.overwrite:
                    print(f"Skipping (exists): {outpath}")
                else:
                    fig, ax = plot_shape_from_hists(
                        h_mc_list, mc_items, h_data=h_data,
                        title=None)
                    ax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)

            # Stacked + trigger efficiency panel (MC + Data) — skip for NoTrigger
            if do_stacked and trig_selections[trig_name]["trig"] is not None:
                outpath = _plot_path("data", "eff", trig_name, var_name)
                if os.path.exists(outpath) and not ARGS.overwrite:
                    print(f"Skipping (exists): {outpath}")
                else:
                    h_denom_list = mc_denom_hists.get(var_name)
                    h_d_denom = data_denom_hists.get(var_name)
                    fig, ax, rax = plot_stacked_with_efficiency(
                        h_mc_list, mc_items,
                        h_data=h_data,
                        h_mc_denom_list=h_denom_list,
                        h_data_denom=h_d_denom,
                        logy=True, title=None)
                    rax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.tight_layout()
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)

            # Stacked + S/sqrt(B) panel (MC + Data)
            if do_stacked:
                outpath = _plot_path("data", "sig", trig_name, var_name)
                if os.path.exists(outpath) and not ARGS.overwrite:
                    print(f"Skipping (exists): {outpath}")
                else:
                    fig, ax, rax = plot_stacked_with_significance(
                        h_mc_list, mc_items, h_data=h_data,
                        sig_indices=sig_indices, bkg_indices=bkg_indices,
                        logy=True, title=None)
                    rax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.tight_layout()
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)

        # Comparison: gen mHH vs signal m4j vs data m4j
        sel = trig_selections[trig_name]
        cmp_path = _plot_path("data", "shape", trig_name, "mHH_m4j")
        if (sel["sig_mHH_ptrs"] and "m4j" in h_by_var
                and (not os.path.exists(cmp_path) or ARGS.overwrite)):
            fig, ax = plt.subplots(figsize=(8, 6))

            # 1) gen mHH (signal)
            h_gen = sel["sig_mHH_ptrs"][0].GetPtr().Clone("h_gen_mHH_cmp")
            for h_ptr in sel["sig_mHH_ptrs"][1:]:
                h_gen.Add(h_ptr.GetPtr())
            edges, vals, _ = th1_to_np(h_gen)
            area = float(np.sum(vals * np.diff(edges)))
            if area > 0:
                vals = vals / area
            ax.step(edges, np.r_[vals, vals[-1]], where="post",
                    linewidth=2, color="tab:blue", label=r"$m_{HH}$ gen-level (signal)")

            # 2) signal m4j (reco)
            h_sig_m4j = h_by_var["m4j"][sig_indices[0]].Clone("h_sig_m4j_cmp")
            for si in sig_indices[1:]:
                h_sig_m4j.Add(h_by_var["m4j"][si])
            edges2, vals2, _ = th1_to_np(h_sig_m4j)
            area2 = float(np.sum(vals2 * np.diff(edges2)))
            if area2 > 0:
                vals2 = vals2 / area2
            ax.step(edges2, np.r_[vals2, vals2[-1]], where="post",
                    linewidth=2, color="tab:red", label=r"$m_{4j}$ reco (signal)")

            # 3) data m4j
            h_data_m4j = h_data_var.get("m4j", None)
            if h_data_m4j is not None:
                edges3, vals3, _ = th1_to_np(h_data_m4j)
                area3 = float(np.sum(vals3 * np.diff(edges3)))
                if area3 > 0:
                    vals3 = vals3 / area3
                ax.step(edges3, np.r_[vals3, vals3[-1]], where="post",
                        linewidth=2, color="black", linestyle="--",
                        label=r"$m_{4j}$ reco (data)")

            ax.set_xlabel(r"Mass [GeV]")
            ax.set_ylabel("Normalised to unity")
            ax.legend()
            ax.grid(True, alpha=0.3)
            cms_label(ax, lumi=LUMI)
            fig.tight_layout()
            fig.savefig(cmp_path, dpi=150)
            print(f"Saved: {cmp_path}")
            plt.close(fig)

        print(f"[{trig_name}] Phase 5: MC+Data plots done")
