#!/usr/bin/env python3
"""
Cutflow & Trigger Efficiency analysis script.

Output:
    plots/light/  or  plots/dark/   (depending on --theme)

Usage:
    python cutflow_TrigEff.py                              # skip existing, light theme -> plots/light/
    python cutflow_TrigEff.py --overwrite                  # regenerate all plots
    python cutflow_TrigEff.py --theme dark                 # dark theme -> plots/dark/
    python cutflow_TrigEff.py --plot-type stacked          # only stacked plots
    python cutflow_TrigEff.py --plot-type shape            # only shape overlay plots
    python cutflow_TrigEff.py --max-files 2                # quick test with 2 files per sample
    python cutflow_TrigEff.py --all-events                 # process all events (ignore MAX_EVENTS)
    python cutflow_TrigEff.py --overwrite --theme dark --plot-type shape  # combine flags
"""

import os
import sys
import argparse
import math
import resource

import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--overwrite", action="store_true",
                    help="Overwrite existing plots instead of skipping them")
parser.add_argument("--theme", choices=["light", "dark"], default="light",
                    help="Plot colour theme (default: light)")
parser.add_argument("--max-files", type=int, default=0,
                    help="Max files per sample (0 = use MAX_EVENTS limit)")
parser.add_argument("--all-events", action="store_true",
                    help="Process all events (ignore MAX_EVENTS limit)")
parser.add_argument("--plot-type", choices=["stacked", "shape", "both"],
                    default="both",
                    help="Which plot types to produce (default: both)")
ARGS = parser.parse_args()

# Increase stack size to 64 MB to prevent stack overflow from deep
# RDataFrame Define/Filter chains with large TChains.
STACK_SIZE = 64 * 1024 * 1024  # 64 MB
resource.setrlimit(resource.RLIMIT_STACK, (STACK_SIZE, resource.RLIM_INFINITY))

import ROOT
import anaConfig

from scouting_utils.data import (
    ls_nanoaod_files_groups,
    load_scouting_data,
    BASE, GROUPS, XSEC, MAX_EVENTS,
)
from scouting_utils.triggers import (
    DST_JetHT_expr,
    PARKING_HH_expr,
)
from scouting_utils.plotting import (
    setup_style,
    plot_stacked_all_mc,
    plot_shape_overlay,
    th1_to_np,
)

import matplotlib.pyplot as plt


# ──────────────────────────────────────────────────────────────────────────────
# ROOT setup
# ──────────────────────────────────────────────────────────────────────────────

ANA_DIR = os.path.expanduser(
    "/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana"
)
os.chdir(ANA_DIR)
sys.path.insert(0, ANA_DIR)

ROOT.gInterpreter.AddIncludePath(ANA_DIR)
ROOT.gROOT.LoadMacro("Particle.h+")
ROOT.gROOT.LoadMacro("HHbbtautauAnaElements.C+")

ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.ROOT.EnableImplicitMT(32)  # use 32 threads for RDataFrame parallelism

print("CWD =", os.getcwd())
print("sys.path[0] =", sys.path[0])


LUMI = 1.0  # fb^-1
BTAG_WP = 0.1  # loose b-tag working point for cutflow


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
# MC normalisation  (3-phase RunGraphs pattern)
# ──────────────────────────────────────────────────────────────────────────────

print(f"\nNormalising MC to L = {LUMI} fb^-1  (MAX_EVENTS = {MAX_EVENTS})")
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
    scale = xsec_pb * LUMI * 1000.0 / sum_genw
    print(f"  {sample_name[:60]:60s}  xsec={xsec_pb:.4g} pb  "
          f"sum_genw={sum_genw:.4g}  scale={scale:.4e}")
    mc[sample_name] = df.Define("w", f"genWeight * {scale:.10e}")

# Convenience: lists of sample names per group
dy_samples  = list(group_files_by_sample["DY"].keys())
tt_samples  = list(group_files_by_sample["TT"].keys())
sig_samples = list(group_files_by_sample["HHbbtt"].keys())


# ──────────────────────────────────────────────────────────────────────────────
# Data loading
# ──────────────────────────────────────────────────────────────────────────────

YEARS = ["2024"]
RUNS = None  # set e.g. ["Run2024C"] for a single run

data_files, data_summary = load_scouting_data(years=YEARS, runs=RUNS)
if ARGS.max_files > 0:
    data_files = data_files[:ARGS.max_files]
elif not ARGS.all_events and MAX_EVENTS > 0:
    n_files = max(1, MAX_EVENTS // 10_000)
    data_files = data_files[:n_files]
data_df = ROOT.RDataFrame("Events", data_files)
print(f"\nLoaded {len(data_files)} data files into RDataFrame")


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
        .Define("ak4_pt0", "ScoutingPFJetRecluster_pt[0]")
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
    # b-jet selection: ParticleNet BvsAll discriminator = (prob_b + prob_bb) / sum(all probs)
    # Sort jets by BvsAll in descending order; [0] = highest, [1] = second-highest
    df = (df
        .Define("pnet_BvsAll",
                "(ScoutingPFJetRecluster_particleNet_prob_b + ScoutingPFJetRecluster_particleNet_prob_bb)"
                " / (ScoutingPFJetRecluster_particleNet_prob_b + ScoutingPFJetRecluster_particleNet_prob_bb"
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
    return df

data_df = define_kinematics(data_df)
for name in mc:
    mc[name] = define_kinematics(mc[name])


# ──────────────────────────────────────────────────────────────────────────────
# Selection & style
# ──────────────────────────────────────────────────────────────────────────────

base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_acc = data_df.Filter(base_cut)
mc_acc = {name: df.Filter(base_cut) for name, df in mc.items()}

setup_style(dark=(ARGS.theme == "dark"))

PLOT_DIR = os.path.join(ANA_DIR, "plots", ARGS.theme)
os.makedirs(PLOT_DIR, exist_ok=True)

TRIG_LIST = [
    ("DST_JetHT",    DST_JetHT_expr),
    ("PARKING_HH",   PARKING_HH_expr),
]

CUTFLOW_STEPS = [
    ("Trigger",                       None),
    ("≥4 jets",                       "nScoutingPFJetRecluster >= 4"),
    ("4 jets pT > 20",                "ak4_pt0 > 20 && ak4_pt1 > 20 && ak4_pt2 > 20 && ak4_pt3 > 20"),
    (f"2 b-tag (BvsAll > {BTAG_WP})", f"b0_score > {BTAG_WP} && b1_score > {BTAG_WP}"),
]

PLOT_VARS = [
    # (branch,     xlabel,                    nbins, xmin, xmax)
    # -- Jet pT --
    ("ak4_pt0",   r"Leading jet $p_T$ [GeV]",        36,    0,  1800),
    ("ak4_pt1",   r"Sub-leading jet $p_T$ [GeV]",    36,    0,  1000),
    ("ak4_pt2",   r"3rd jet $p_T$ [GeV]",            36,    0,   600),
    ("ak4_pt3",   r"4th jet $p_T$ [GeV]",            36,    0,   400),
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
    ("HT",        r"$H_T$ [GeV]",                    40,    0,  3000),
    ("MHT",       r"$\slash{H}_T$ [GeV]",            30,    0,  1500),
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
    ("b0_score",  r"$b_0$ ParticleNet b-score",       25,    0,     1),
    ("b1_score",  r"$b_1$ ParticleNet b-score",       25,    0,     1),
    ("mbb",       r"$m_{bb}$ [GeV]",                  30,    0,   300),
    ("dR_bb",     r"$\Delta R(b_0, b_1)$",           30,    0,     6),
    ("ptbb",      r"$p_T^{bb}$ [GeV]",               30,    0,  1000),
]


def asimov_significance(S, B):
    """Asimov significance: Z_A = sqrt(2 * [(S+B)*ln(1 + S/B) - S])."""
    if B <= 0 or S <= 0:
        return float("nan")
    return math.sqrt(2 * ((S + B) * math.log(1 + S / B) - S))


do_stacked = ARGS.plot_type in ("stacked", "both")
do_shape   = ARGS.plot_type in ("shape", "both")


# ══════════════════════════════════════════════════════════════════════════════
# Phase 1 — Cutflow & yields  (no plotting here)
# ══════════════════════════════════════════════════════════════════════════════

cutflow_tables = {}
trig_selections = {}   # store per-trigger objects for the plotting phase

for trig_name, trig in TRIG_LIST:
    print("=" * 63)
    print(f"  {trig_name}")
    print("=" * 63)

    # ── Cutflow: progressive cuts ──
    data_cf = data_df.Filter(trig)
    mc_cf = {name: df.Filter(trig) for name, df in mc.items()}

    cutflow_ptrs = []
    for step_name, cut_expr in CUTFLOW_STEPS:
        if cut_expr is not None:
            data_cf = data_cf.Filter(cut_expr)
            mc_cf = {s: d.Filter(cut_expr) for s, d in mc_cf.items()}
        d_ptr = data_cf.Count()
        w_ptrs = {s: mc_cf[s].Sum("w") for s in mc_cf}
        cutflow_ptrs.append((step_name, d_ptr, w_ptrs))

    # ── Final selection (for yields now, plots later) ──
    data_sel = data_acc.Filter(trig)
    mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}

    # Signal decay channels
    hh = [mc_sel[s].Filter("GenDecay.decayType==1") for s in sig_samples]
    hm = [mc_sel[s].Filter("GenDecay.decayType==2") for s in sig_samples]
    he = [mc_sel[s].Filter("GenDecay.decayType==3") for s in sig_samples]

    data_n_ptr = data_sel.Count()
    sum_w  = {s: mc_sel[s].Sum("w") for s in mc_sel}
    sum_w2 = {s: mc_sel[s].Define("w2", "w*w").Sum("w2") for s in mc_sel}
    hh_w = [d.Sum("w") for d in hh]
    hm_w = [d.Sum("w") for d in hm]
    he_w = [d.Sum("w") for d in he]

    sig_mHH_ptrs = [mc_sel[s].Histo1D(
        (f"h_mHH_{s}_{trig_name}",
         "m_{HH} gen-level;m_{HH} [GeV];Events", 40, 200, 1200),
        "gen_mHH", "w") for s in sig_samples]

    # ── Run ALL actions in one parallel pass ──
    cf_all = []
    for _, d_ptr, w_ptrs in cutflow_ptrs:
        cf_all.append(d_ptr)
        cf_all.extend(w_ptrs.values())

    all_ptrs = (cf_all
                + [data_n_ptr] + list(sum_w.values())
                + list(sum_w2.values())
                + hh_w + hm_w + he_w
                + sig_mHH_ptrs)
    ROOT.RDF.RunGraphs(all_ptrs)

    # ── Extract cutflow results ──
    cutflow_rows = []
    for step_name, d_ptr, w_ptrs in cutflow_ptrs:
        data_n_cf = d_ptr.GetValue()
        dy_cf  = sum(w_ptrs[s].GetValue() for s in dy_samples)
        tt_cf  = sum(w_ptrs[s].GetValue() for s in tt_samples)
        sig_cf = sum(w_ptrs[s].GetValue() for s in sig_samples)
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
        print(f"{step:<30s} {dn:>12d} {dy:>14.4g} {tt:>14.4g} "
              f"{sig:>14.4g} {sb:>10.4g} {za:>10.4g}")

    # ── Extract final yields ──
    data_n = data_n_ptr.GetValue()

    def group_yield(samples):
        y = sum(sum_w[s].GetValue() for s in samples)
        yerr = math.sqrt(sum(sum_w2[s].GetValue() for s in samples))
        return y, yerr

    dy_y, dy_yerr = group_yield(dy_samples)
    tt_y, tt_yerr = group_yield(tt_samples)
    sig_y, sig_yerr = group_yield(sig_samples)
    y_hh = sum(d.GetValue() for d in hh_w)
    y_hm = sum(d.GetValue() for d in hm_w)
    y_he = sum(d.GetValue() for d in he_w)
    B = dy_y + tt_y
    S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")

    print(f"\nData events: {data_n}")
    print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
    print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
    print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")
    print(f"S/sqrt(B):   {S_over_sqrtB:.6g}")
    print(f"y_hh:y_hm:y_he: {y_hh:.6g}:{y_hm:.6g}:{y_he:.6g}\n")

    # Store selections for the plotting phase
    trig_selections[trig_name] = {
        "data_sel": data_sel, "mc_sel": mc_sel,
        "hh": hh, "hm": hm, "he": he,
        "sig_mHH_ptrs": sig_mHH_ptrs,
    }


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


# ══════════════════════════════════════════════════════════════════════════════
# Phase 2 — Plotting  (no cutflow/yield computation here)
# ══════════════════════════════════════════════════════════════════════════════

for trig_name in trig_selections:
    sel = trig_selections[trig_name]
    mc_sel = sel["mc_sel"]
    hh, hm, he = sel["hh"], sel["hm"], sel["he"]

    mc_items = [
        {"dfs": [mc_sel[s] for s in dy_samples],
         "label": "DY",   "color": "tab:orange"},
        {"dfs": [mc_sel[s] for s in tt_samples],
         "label": "TT",   "color": "tab:green"},
        {"dfs": hh,
         "label": r"$HH\to bb\tau_h\tau_h$",  "color": "tab:red"},
        {"dfs": hm,
         "label": r"$HH\to bb\tau_e\tau_h$",  "color": "tab:pink"},
        {"dfs": he,
         "label": r"$HH\to bb\tau_\mu\tau_h$", "color": "tab:purple"},
    ]

    # -- Signal-only gen mHH shape --
    mHH_path = os.path.join(PLOT_DIR, f"{trig_name}_shape_gen_mHH.png")
    if os.path.exists(mHH_path) and not ARGS.overwrite:
        print(f"Skipping (exists): {mHH_path}")
    else:
        fig, ax = plt.subplots(figsize=(8, 6))
        for h_ptr in sel["sig_mHH_ptrs"]:
            edges, vals, _ = th1_to_np(h_ptr.GetPtr())
            area = float(np.sum(vals * np.diff(edges)))
            if area > 0:
                vals = vals / area
            ax.step(edges, np.r_[vals, vals[-1]], where="post",
                    linewidth=2, label="HH signal (gen)")
        ax.set_xlabel(r"$m_{HH}$ (gen-level) [GeV]")
        ax.set_ylabel("Normalised")
        ax.set_title(f"{trig_name} — gen-level $m_{{HH}}$")
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        fig.savefig(mHH_path, dpi=150)
        print(f"Saved: {mHH_path}")
        plt.close(fig)

    # -- Stacked & shape plots --
    for var_name, xlabel, nbins, vmin, vmax in PLOT_VARS:
        if do_stacked:
            outpath = os.path.join(PLOT_DIR, f"{trig_name}_stacked_{var_name}.png")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax = plot_stacked_all_mc(
                    mc_items=mc_items, var=var_name, weight="w",
                    nbins=nbins, xmin=vmin, xmax=vmax, logy=True,
                    title=trig_name)
                ax.set_xlabel(xlabel)
                fig.savefig(outpath)
                print(f"Saved: {outpath}")
                plt.close(fig)

        if do_shape:
            outpath = os.path.join(PLOT_DIR, f"{trig_name}_shape_{var_name}.png")
            if os.path.exists(outpath) and not ARGS.overwrite:
                print(f"Skipping (exists): {outpath}")
            else:
                fig, ax = plot_shape_overlay(
                    mc_items=mc_items, var=var_name, weight="w",
                    nbins=nbins, xmin=vmin, xmax=vmax,
                    title=f"{trig_name} (shape)")
                ax.set_xlabel(xlabel)
                fig.savefig(outpath)
                print(f"Saved: {outpath}")
                plt.close(fig)
