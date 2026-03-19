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
    python cutflow_TrigEff.py --max-files 2                # quick test with 2 MC files per sample
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
import math
import time
import resource

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing plots instead of skipping them")
parser.add_argument("--theme", choices=["light", "dark"], default="light",
                    help="Plot colour theme (default: light)")
parser.add_argument("--max-files", type=int, default=0,
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
    PARKING_HH_expr,
)
from utils.plotting import (
    setup_style,
    cms_label,
    plot_stacked_from_hists,
    plot_shape_from_hists,
    plot_stacked_with_efficiency,
    plot_stacked_with_significance,
    plot_trigger_shape_overlay,
    th1_to_np,
)

import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# Histogram cache helpers
# ──────────────────────────────────────────────────────────────────────────────

_CACHE_VERSION = 1   # bump to invalidate all caches

_INVALIDATION_FILES = [
    "cutflow_TrigEff.py",
    "config/samples.yaml",
    "HHbbtautauAnaElements.C",
    "Particle.h",
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
                    mc_items_by_trig, lumi, meta):
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

    tf.Close()
    return mc_hists_by_trig, mc_denom_hists, mc_items_by_trig, lumi


# ──────────────────────────────────────────────────────────────────────────────
# ROOT setup
# ──────────────────────────────────────────────────────────────────────────────

ANA_DIR = os.path.expanduser(
    "/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana"
)
os.chdir(ANA_DIR)
sys.path.insert(0, ANA_DIR)

ROOT.gInterpreter.AddIncludePath(ANA_DIR)
BUILD_DIR = os.path.join(ANA_DIR, "build")
os.makedirs(BUILD_DIR, exist_ok=True)
ROOT.gSystem.SetBuildDir(BUILD_DIR)
ROOT.gROOT.LoadMacro("Particle.h+")
ROOT.gROOT.LoadMacro("HHbbtautauAnaElements.C+")

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.ROOT.EnableImplicitMT(ARGS.nthreads)

print("CWD =", os.getcwd())
print("sys.path[0] =", sys.path[0])


BTAG_WP = 0.1  # loose b-tag working point for cutflow

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
    """Filter missing files and apply --max-files / MAX_EVENTS file cap."""
    good_files = [f for f in files if os.path.isfile(f)]
    if len(good_files) < len(files):
        print(f"  WARNING: {len(files) - len(good_files)} missing file(s) "
              f"in {sample_name}, using {len(good_files)}/{len(files)}")
    if not good_files:
        raise FileNotFoundError(f"No valid files for {sample_name}")
    if ARGS.max_files > 0:
        good_files = good_files[:ARGS.max_files]
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

for group_name, samples_dict in group_files_by_sample.items():
    print(f"\n[{group_name}]")
    for sample_name, files in samples_dict.items():
        good_files = _limit_files(files, sample_name)
        df = ROOT.RDataFrame("Events", good_files)
        sum_genw_ptrs[sample_name] = df.Sum("genWeight")
        mc_dfs[sample_name] = df
        print(f"  {sample_name[:60]:60s}  ({len(good_files)} files)")

print("\nRunning sum(genWeight) for all samples in parallel...")
ROOT.RDF.RunGraphs(list(sum_genw_ptrs.values()))

for sample_name, df in mc_dfs.items():
    sum_genw = sum_genw_ptrs[sample_name].GetValue()
    xsec_pb = XSEC[sample_name]
    scale_nolumi = xsec_pb * 1000.0 / sum_genw  # LUMI applied post-hoc
    print(f"  {sample_name[:60]:60s}  xsec={xsec_pb:.4g} pb  "
          f"sum_genw={sum_genw:.4g}  scale(no L)={scale_nolumi:.4e}")
    mc[sample_name] = df.Define("w", f"genWeight * {scale_nolumi:.10e}")

# Convenience: lists of sample names per group
dy_samples  = list(group_files_by_sample["DY"].keys())
tt_samples  = list(group_files_by_sample["TT"].keys())
sig_samples = list(group_files_by_sample["HHbbtt"].keys())




# ──────────────────────────────────────────────────────────────────────────────
# Gen-level columns
# ──────────────────────────────────────────────────────────────────────────────

gen_decay_expr = (
    "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)"
    .format("GenPart")
)

for name in mc:
    mc[name] = mc[name].Define("GenDecay", gen_decay_expr)

for name in sig_samples:
    mc[name] = (mc[name]
        .Define("genHbb_idx",       "(int)GenDecay.Htob")
        .Define("genHtautau_idx",   "(int)GenDecay.Htotau")
        .Define("genHbb_p4",        "Ana::getP4(genHbb_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("genHtautau_p4",    "Ana::getP4(genHtautau_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
        .Define("gen_mHH",          "(genHbb_p4 + genHtautau_p4).M()")
        .Define("gen_pt_Hbb",       "genHbb_p4.Pt()")
        .Define("gen_pt_Htautau",   "genHtautau_p4.Pt()")
        .Define("gen_eta_Hbb",      "genHbb_p4.Eta()")
        .Define("gen_eta_Htautau",  "genHtautau_p4.Eta()")
    )


# ──────────────────────────────────────────────────────────────────────────────
# Kinematic definitions
# ──────────────────────────────────────────────────────────────────────────────

def define_kinematics(df):
    """Define derived kinematic columns for jets.

    Intermediate columns are avoided to keep the branch-proxy chain shallow
    and prevent stack overflows with large TChains.
    """
    df = (df
        .Define("ak4_pt0", "ScoutingPFJetRecluster_pt[0]") # Check this (Could be a bug in the processing)
        .Define("ak4_pt1", "ScoutingPFJetRecluster_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster_pt[3]")
        .Define("HT",  "Sum(ScoutingPFJetRecluster_pt)")
        .Define("nJets", "nScoutingPFJetRecluster")
        .Define("nLeptons", "nScoutingMuonVtx + nScoutingElectron")
    )
    df = (df.Define("mjj_01",
                "(float)(ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[0],ScoutingPFJetRecluster_eta[0],"
                "ScoutingPFJetRecluster_phi[0],ScoutingPFJetRecluster_mass[0])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[1],ScoutingPFJetRecluster_eta[1],"
                "ScoutingPFJetRecluster_phi[1],ScoutingPFJetRecluster_mass[1])).M()")
        .Define("dR_01",
                "ROOT::VecOps::DeltaR("
                "ScoutingPFJetRecluster_eta[0],ScoutingPFJetRecluster_eta[1],"
                "ScoutingPFJetRecluster_phi[0],ScoutingPFJetRecluster_phi[1])")
        .Define("MHT",
                "(float)sqrt("
                "pow(Sum(ScoutingPFJetRecluster_pt*cos(ScoutingPFJetRecluster_phi)),2) + "
                "pow(Sum(ScoutingPFJetRecluster_pt*sin(ScoutingPFJetRecluster_phi)),2))")
        .Define("m4j",
                "(float)(ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[0],ScoutingPFJetRecluster_eta[0],"
                "ScoutingPFJetRecluster_phi[0],ScoutingPFJetRecluster_mass[0])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[1],ScoutingPFJetRecluster_eta[1],"
                "ScoutingPFJetRecluster_phi[1],ScoutingPFJetRecluster_mass[1])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[2],ScoutingPFJetRecluster_eta[2],"
                "ScoutingPFJetRecluster_phi[2],ScoutingPFJetRecluster_mass[2])"
                " + ROOT::Math::PtEtaPhiMVector("
                "ScoutingPFJetRecluster_pt[3],ScoutingPFJetRecluster_eta[3],"
                "ScoutingPFJetRecluster_phi[3],ScoutingPFJetRecluster_mass[3])).M()")
    )
    # b-jet selection: ParticleNet discriminators (prob_bb is for AK8 fat jets, not AK4)
    # pnet_b_raw  = raw prob_b score per jet
    # pnet_BvsAll = prob_b / (prob_b + prob_c + prob_cc + prob_g + prob_uds + prob_undef)
    # Sort jets by BvsAll in descending order; [0] = highest, [1] = second-highest
    df = (df
        .Define("pnet_b_raw", "ScoutingPFJetRecluster_particleNet_prob_b")
        .Define("pnet_BvsAll",
                "ScoutingPFJetRecluster_particleNet_prob_b"
                " / (ScoutingPFJetRecluster_particleNet_prob_b"
                " + ScoutingPFJetRecluster_particleNet_prob_c + ScoutingPFJetRecluster_particleNet_prob_cc"
                " + ScoutingPFJetRecluster_particleNet_prob_g + ScoutingPFJetRecluster_particleNet_prob_uds"
                " + ScoutingPFJetRecluster_particleNet_prob_undef)")
        .Define("bsort_idx",
                "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(pnet_BvsAll))")
        .Define("b0_idx", "(int)bsort_idx[0]")
        .Define("b1_idx", "(int)bsort_idx[1]")
        .Define("b0_pt",    "ScoutingPFJetRecluster_pt[b0_idx]")
        .Define("b0_eta",   "ScoutingPFJetRecluster_eta[b0_idx]")
        .Define("b0_phi",   "ScoutingPFJetRecluster_phi[b0_idx]")
        .Define("b0_mass",  "ScoutingPFJetRecluster_mass[b0_idx]")
        .Define("b0_score", "pnet_BvsAll[b0_idx]")
        .Define("b1_pt",    "ScoutingPFJetRecluster_pt[b1_idx]")
        .Define("b1_eta",   "ScoutingPFJetRecluster_eta[b1_idx]")
        .Define("b1_phi",   "ScoutingPFJetRecluster_phi[b1_idx]")
        .Define("b1_mass",  "ScoutingPFJetRecluster_mass[b1_idx]")
        .Define("b1_score", "pnet_BvsAll[b1_idx]")
        .Define("b0_raw",   "pnet_b_raw[b0_idx]")
        .Define("b1_raw",   "pnet_b_raw[b1_idx]")
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
        .Define("ak4_eta0", "ScoutingPFJetRecluster_eta[0]")
        .Define("ak4_eta1", "ScoutingPFJetRecluster_eta[1]")
        .Define("ak4_eta2", "ScoutingPFJetRecluster_eta[2]")
        .Define("ak4_eta3", "ScoutingPFJetRecluster_eta[3]")
        .Define("ak4_mass0", "ScoutingPFJetRecluster_mass[0]")
        .Define("ak4_mass1", "ScoutingPFJetRecluster_mass[1]")
        .Define("nMuons", "nScoutingMuonVtx")
        .Define("nElectrons", "nScoutingElectron")
        .Define("centrality",
                "(float)(Sum(ScoutingPFJetRecluster_pt) / "
                "Sum(ScoutingPFJetRecluster_pt * cosh(ScoutingPFJetRecluster_eta)))")
        .Define("dEta_01",
                "(float)abs(ScoutingPFJetRecluster_eta[0] - ScoutingPFJetRecluster_eta[1])")
    )
    # H→bb candidate: best dijet pair with m_jj ∈ [100, 150] GeV (closest to 125)
    df = (df
        .Define("hbb_pair",
                "Ana::findDijetInWindow("
                "ScoutingPFJetRecluster_pt, ScoutingPFJetRecluster_eta,"
                "ScoutingPFJetRecluster_phi, ScoutingPFJetRecluster_mass,"
                "100.f, 150.f, 125.f)")
        .Define("has_hbb",      "hbb_pair.i1 >= 0")
        .Define("mbb_cand",     "hbb_pair.mass")
        # H→ττ candidate: best dijet pair from remaining jets, m_jj ∈ [40, 150] GeV
        .Define("htautau_pair",
                "Ana::findDijetInWindow("
                "ScoutingPFJetRecluster_pt, ScoutingPFJetRecluster_eta,"
                "ScoutingPFJetRecluster_phi, ScoutingPFJetRecluster_mass,"
                "40.f, 150.f, 80.f, {hbb_pair.i1, hbb_pair.i2})")
        .Define("has_htautau",  "htautau_pair.i1 >= 0")
        .Define("mtautau_cand", "htautau_pair.mass")
    )
    return df

if data_df is not None:
    data_df = define_kinematics(data_df)
for name in mc:
    mc[name] = define_kinematics(mc[name])


# ──────────────────────────────────────────────────────────────────────────────
# Selection & style
# ──────────────────────────────────────────────────────────────────────────────

base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_acc = data_df.Filter(base_cut) if data_df is not None else None
mc_acc = {name: df.Filter(base_cut) for name, df in mc.items()}

setup_style(dark=(ARGS.theme == "dark"))

PLOT_DIR = os.path.join(ANA_DIR, "plots", ARGS.theme)
os.makedirs(PLOT_DIR, exist_ok=True)
for _subdir in ["stacked", "shape", "eff_stacked", "sig_stacked", "overlay",
                "data_stacked", "data_shape", "data_eff_stacked", "data_sig_stacked"]:
    os.makedirs(os.path.join(PLOT_DIR, _subdir), exist_ok=True)

TRIG_LIST = [
    ("NoTrigger",    None),
    ("DST_JetHT",    DST_JetHT_expr),
    ("PARKING_HH",   PARKING_HH_expr),
]

CUTFLOW_STEPS = [
    ("Trigger",                       None),
    ("≥4 jets",                       "nScoutingPFJetRecluster >= 4"),
    ("4 jets pT > 20",                "ak4_pt0 > 20 && ak4_pt1 > 20 && ak4_pt2 > 20 && ak4_pt3 > 20"),
    (f"2 b-tag (BvsAll > {BTAG_WP})", f"b0_score > {BTAG_WP} && b1_score > {BTAG_WP}"),
    ("H→bb cand (m_jj ∈ [100,150])",  "has_hbb"),
    ("H→ττ cand (m_jj ∈ [40,150])",   "has_hbb && has_htautau"),
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
    ("b0_score",  r"$b_0$ PNet $b/(b+c+cc+g+uds+undef)$", 25,  0,     1),
    ("b1_score",  r"$b_1$ PNet $b/(b+c+cc+g+uds+undef)$", 25,  0,     1),
    ("b0_raw",    r"$b_0$ PNet raw prob\_b",           25,    0,     1),
    ("b1_raw",    r"$b_1$ PNet raw prob\_b",           25,    0,     1),
    ("mbb",       r"$m_{bb}$ [GeV]",                  30,    0,   300),
    ("dR_bb",     r"$\Delta R(b_0, b_1)$",           30,    0,     6),
    ("ptbb",      r"$p_T^{bb}$ [GeV]",               30,    0,  1000),
    # -- Higgs candidate dijet masses --
    ("mbb_cand",     r"$m_{jj}^{H \to bb}$ candidate [GeV]",    25,   50,   200),
    ("mtautau_cand", r"$m_{jj}^{H \to \tau\tau}$ candidate [GeV]", 25, 0, 200),
]

# ── Filter variables if --plot-vars given ──
if ARGS.plot_vars:
    _all_var_names = {v[0] for v in PLOT_VARS}
    PLOT_VARS = [v for v in PLOT_VARS
                 if any(fnmatch.fnmatch(v[0], p) for p in ARGS.plot_vars)]
    _matched = {v[0] for v in PLOT_VARS}
    _unknown = set()
    for p in ARGS.plot_vars:
        if not any(fnmatch.fnmatch(v, p) for v in _all_var_names):
            _unknown.add(p)
    if _unknown:
        print(f"WARNING: no variables matched patterns: {_unknown}")
    if not PLOT_VARS:
        print("ERROR: --plot-vars matched no variables")
        sys.exit(1)
    print(f"[--plot-vars] Plotting {len(PLOT_VARS)} variables: "
          f"{', '.join(v[0] for v in PLOT_VARS)}")


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
    max_files=ARGS.max_files,
    plot_vars_names=_plot_var_names,
)

_cache_loaded = False
if not ARGS.recache:
    _cached = load_hist_cache(_cache_path, _expected_meta, _trig_names)
    if _cached is not None:
        mc_hists_by_trig, mc_denom_hists, mc_items_by_trig, LUMI = _cached
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
    
        # Signal decay channels
        hh = [mc_sel[s].Filter("GenDecay.decayType==1") for s in sig_samples]
        hm = [mc_sel[s].Filter("GenDecay.decayType==2") for s in sig_samples]
        he = [mc_sel[s].Filter("GenDecay.decayType==3") for s in sig_samples]
    
        sig_mHH_ptrs = [mc_sel[s].Histo1D(
            (f"h_mHH_{s}_{trig_name}",
             "m_{HH} gen-level;m_{HH} [GeV];Events", 32, 0, 800),
            "gen_mHH", "w") for s in sig_samples]
        unified_ptrs.extend(sig_mHH_ptrs)
    
        # Store selections for the plotting phase
        trig_selections[trig_name] = {
            "trig": trig,
            "data_sel": data_sel, "mc_sel": mc_sel,
            "hh": hh, "hm": hm, "he": he,
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
            hh_w = [d.Sum("w") for d in hh]
            hm_w = [d.Sum("w") for d in hm]
            he_w = [d.Sum("w") for d in he]
    
            if data_n_ptr is not None:
                unified_ptrs.append(data_n_ptr)
            unified_ptrs.extend(sum_w.values())
            unified_ptrs.extend(sum_w2.values())
            unified_ptrs.extend(hh_w + hm_w + he_w)
    
            _phase1_data[trig_name] = {
                "cutflow_ptrs": cutflow_ptrs,
                "data_n_ptr": data_n_ptr,
                "sum_w": sum_w, "sum_w2": sum_w2,
                "hh_w": hh_w, "hm_w": hm_w, "he_w": he_w,
            }
    
    # ── Phase 2: Book MC denominator histograms (no trigger) ──
    mc_denom_groups = [
        ([mc_acc[s] for s in dy_samples],  "DY",  "tab:orange"),
        ([mc_acc[s] for s in tt_samples],  "TT",  "tab:green"),
        ([mc_acc[s].Filter("GenDecay.decayType==1") for s in sig_samples],
         r"$HH\to bb\tau_h\tau_h$", "tab:red"),
        ([mc_acc[s].Filter("GenDecay.decayType==2") for s in sig_samples],
         r"$HH\to bb\tau_\mu\tau_h$", "tab:pink"),
        ([mc_acc[s].Filter("GenDecay.decayType==3") for s in sig_samples],
         r"$HH\to bb\tau_e\tau_h$", "tab:purple"),
    ]
    
    denom_book = {}
    for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
        denom_book[var_name] = []
        for gi, (dfs, _lbl, _col) in enumerate(mc_denom_groups):
            group_ptrs = []
            for si, df in enumerate(dfs):
                uid = f"denom_{var_name}_{gi}_{si}"
                ptr = df.Histo1D(
                    (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                    var_name, "w")
                group_ptrs.append(ptr)
                unified_ptrs.append(ptr)
            denom_book[var_name].append(group_ptrs)
    
    # ── Phase 2: Book per-trigger MC histograms ──
    mc_items_by_trig = {}
    _histo_books = {}  # trig_name -> {var_name -> [[ptrs per group]]}
    
    for trig_name in trig_selections:
        sel = trig_selections[trig_name]
        mc_sel = sel["mc_sel"]
        hh, hm, he = sel["hh"], sel["hm"], sel["he"]
    
        mc_groups = [
            ([mc_sel[s] for s in dy_samples], "DY",                              "tab:orange"),
            ([mc_sel[s] for s in tt_samples], "TT",                              "tab:green"),
            (hh,                              r"$HH\to bb\tau_h\tau_h$",         "tab:red"),
            (hm,                              r"$HH\to bb\tau_\mu\tau_h$",       "tab:pink"),
            (he,                              r"$HH\to bb\tau_e\tau_h$",         "tab:purple"),
        ]
        mc_items = [{"label": lbl, "color": col} for _, lbl, col in mc_groups]
        mc_items_by_trig[trig_name] = mc_items
    
        histo_book = {}
        for var_name, _xlabel, nbins, vmin, vmax in PLOT_VARS:
            histo_book[var_name] = []
            for gi, (dfs, _lbl, _col) in enumerate(mc_groups):
                group_ptrs = []
                for si, df in enumerate(dfs):
                    uid = f"{trig_name}_{var_name}_{gi}_{si}"
                    ptr = df.Histo1D(
                        (f"h_{uid}", f";{var_name};Events", nbins, vmin, vmax),
                        var_name, "w")
                    group_ptrs.append(ptr)
                    unified_ptrs.append(ptr)
                histo_book[var_name].append(group_ptrs)
        _histo_books[trig_name] = histo_book
    
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
                sig_cf = sum(w_ptrs[s].GetValue() for s in sig_samples) * LUMI
                B_cf = dy_cf + tt_cf
                s_sqrtb = sig_cf / math.sqrt(B_cf) if B_cf > 0 else float("nan")
                z_a = asimov_significance(sig_cf, B_cf)
                cutflow_rows.append(
                    (step_name, data_n_cf, dy_cf, tt_cf, sig_cf, s_sqrtb, z_a))
            cutflow_tables[trig_name] = cutflow_rows
    
            print(f"\n{'Cut':<30s} {'Data':>12s} {'DY':>14s} {'TT':>14s} "
                  f"{'Signal':>14s} {'S/√B':>10s} {'Z_A':>10s}")
            print("-" * 108)
            for row in cutflow_rows:
                step, dn, dy, tt, sig, sb, za = row
                dn_str = f"{dn:>12d}" if dn >= 0 else f"{'(no data)':>12s}"
                print(f"{step:<30s} {dn_str} {dy:>14.4g} {tt:>14.4g} "
                      f"{sig:>14.4g} {sb:>10.4g} {za:>10.4g}")
    
            data_n_ptr = p1["data_n_ptr"]
            sum_w, sum_w2 = p1["sum_w"], p1["sum_w2"]
            hh_w, hm_w, he_w = p1["hh_w"], p1["hm_w"], p1["he_w"]
            data_n = data_n_ptr.GetValue() if data_n_ptr is not None else -1
    
            def group_yield(samples):
                y = sum(sum_w[s].GetValue() for s in samples) * LUMI
                yerr = math.sqrt(sum(sum_w2[s].GetValue() for s in samples)) * LUMI
                return y, yerr
    
            dy_y, dy_yerr = group_yield(dy_samples)
            tt_y, tt_yerr = group_yield(tt_samples)
            sig_y, sig_yerr = group_yield(sig_samples)
            y_hh = sum(d.GetValue() for d in hh_w) * LUMI
            y_hm = sum(d.GetValue() for d in hm_w) * LUMI
            y_he = sum(d.GetValue() for d in he_w) * LUMI
            B = dy_y + tt_y
            S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")
    
            print(f"\nData events: {data_n if data_n >= 0 else '(no data)'}")
            print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
            print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
            print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")
            print(f"S/sqrt(B):   {S_over_sqrtB:.6g}")
            print(f"y_hh:y_hm:y_he: {y_hh:.6g}:{y_hm:.6g}:{y_he:.6g}\n")
    
        # ── Write cutflow markdown ──
        cutflow_path = os.path.join(PLOT_DIR, "cutflow.md")
        with open(cutflow_path, "w") as f:
            f.write("# Cutflow Tables\n\n")
            f.write(f"Luminosity: {LUMI} fb$^{{-1}}$\n\n")
            for trig_name, rows in cutflow_tables.items():
                f.write(f"## {trig_name}\n\n")
                f.write("| Cut | Data | DY | TT | Signal "
                        "| S/sqrt(B) | Z_A (Asimov) |\n")
                f.write("|-----|-----:|---:|---:|-------:"
                        "|----------:|-------------:|\n")
                for step, dn, dy, tt, sig, sb, za in rows:
                    f.write(f"| {step} | {dn:,d} | {dy:.4g} | {tt:.4g} "
                            f"| {sig:.4g} | {sb:.4g} | {za:.4g} |\n")
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

    # ── Save histogram cache ──
    save_hist_cache(_cache_path, mc_hists_by_trig, mc_denom_hists,
                    mc_items_by_trig, LUMI, _expected_meta)


# ══════════════════════════════════════════════════════════════════════════════
# Phase 3 — Draw MC-only plots  (no event loop, just drawing)
# ══════════════════════════════════════════════════════════════════════════════

for trig_name in trig_selections:
    sel = trig_selections[trig_name]
    h_by_var = mc_hists_by_trig[trig_name]
    mc_items = mc_items_by_trig[trig_name]

    # Comparison: gen mHH vs signal m4j (MC only, data added in Phase 5)
    mHH_path = os.path.join(PLOT_DIR, "shape", f"{trig_name}_mHH_m4j.png")
    if not sel["sig_mHH_ptrs"]:
        pass  # skip when loaded from cache (no lazy ptrs available)
    elif os.path.exists(mHH_path) and not ARGS.overwrite:
        print(f"Skipping (exists): {mHH_path}")
    else:
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

        # 2) signal m4j (reco) — combine hh + hm + he (groups 2,3,4)
        h_sig_m4j = h_by_var["m4j"][2].Clone("h_sig_m4j_combined")
        h_sig_m4j.Add(h_by_var["m4j"][3])
        h_sig_m4j.Add(h_by_var["m4j"][4])
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

    # MC-only stacked & shape plots
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        h_mc_list = h_by_var[var_name]

        if do_stacked:
            outpath = os.path.join(PLOT_DIR, "stacked", f"{trig_name}_{var_name}.png")
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
            outpath = os.path.join(PLOT_DIR, "shape", f"{trig_name}_{var_name}.png")
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
            outpath = os.path.join(PLOT_DIR, "eff_stacked", f"{trig_name}_{var_name}.png")
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
            outpath = os.path.join(PLOT_DIR, "sig_stacked", f"{trig_name}_{var_name}.png")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax, rax = plot_stacked_with_significance(
                    h_mc_list, mc_items,
                    sig_indices=[2, 3, 4], bkg_indices=[0, 1],
                    logy=True, title=None)
                rax.set_xlabel(xlabel)
                cms_label(ax, lumi=LUMI)
                fig.tight_layout()
                fig.savefig(outpath, dpi=150)
                print(f"Saved: {outpath}")
                plt.close(fig)

    print(f"[{trig_name}] Phase 3: MC-only plots done")


# ══════════════════════════════════════════════════════════════════════════════
# Phase 3.5 — Trigger shape overlay  (one plot per variable, all triggers)
# ══════════════════════════════════════════════════════════════════════════════

if do_shape:
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        outpath = os.path.join(PLOT_DIR, "overlay", f"{var_name}.png")
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
    SIG_CHANNELS = [
        (2, "hh", r"$HH\to bb\tau_h\tau_h$"),
        (3, "hm", r"$HH\to bb\tau_\mu\tau_h$"),
        (4, "he", r"$HH\to bb\tau_e\tau_h$"),
    ]
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        for gi, ch_key, ch_label in SIG_CHANNELS:
            outpath = os.path.join(PLOT_DIR, "overlay", f"{ch_key}_{var_name}.png")
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
                outpath = os.path.join(PLOT_DIR, "data_stacked", f"{trig_name}_{var_name}.png")
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
                outpath = os.path.join(PLOT_DIR, "data_shape", f"{trig_name}_{var_name}.png")
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
                outpath = os.path.join(PLOT_DIR, "data_eff_stacked", f"{trig_name}_{var_name}.png")
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
                outpath = os.path.join(PLOT_DIR, "data_sig_stacked", f"{trig_name}_{var_name}.png")
                if os.path.exists(outpath) and not ARGS.overwrite:
                    print(f"Skipping (exists): {outpath}")
                else:
                    fig, ax, rax = plot_stacked_with_significance(
                        h_mc_list, mc_items, h_data=h_data,
                        sig_indices=[2, 3, 4], bkg_indices=[0, 1],
                        logy=True, title=None)
                    rax.set_xlabel(xlabel)
                    cms_label(ax, lumi=LUMI)
                    fig.tight_layout()
                    fig.savefig(outpath, dpi=150)
                    print(f"Saved: {outpath}")
                    plt.close(fig)

        # Comparison: gen mHH vs signal m4j vs data m4j
        sel = trig_selections[trig_name]
        cmp_path = os.path.join(PLOT_DIR, f"{trig_name}_data_shape_mHH_m4j.png")
        if not sel["sig_mHH_ptrs"]:
            pass  # skip when loaded from cache (no lazy ptrs available)
        elif os.path.exists(cmp_path) and not ARGS.overwrite:
            print(f"Skipping (exists): {cmp_path}")
        else:
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
            h_sig_m4j = h_by_var["m4j"][2].Clone("h_sig_m4j_cmp")
            h_sig_m4j.Add(h_by_var["m4j"][3])
            h_sig_m4j.Add(h_by_var["m4j"][4])
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
