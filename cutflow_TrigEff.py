#!/usr/bin/env python3
"""
Cutflow & Trigger Efficiency analysis script.

Usage:
    python cutflow_TrigEff.py
"""

import os
import sys
import math

import ROOT
import anaConfig

from scouting_utils.data import (
    ls_nanoaod_files_groups,
    load_scouting_data,
    BASE, GROUPS, FILE_SIZE,
)
from scouting_utils.triggers import (
    DST_JetHT_expr,
    PARKING_HH_expr,
)
from scouting_utils.plotting import (
    setup_style,
    plot_stacked_all_mc,
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
# ROOT.ROOT.DisableImplicitMT()

print("CWD =", os.getcwd())
print("sys.path[0] =", sys.path[0])


# ──────────────────────────────────────────────────────────────────────────────
# Load MC samples
# ──────────────────────────────────────────────────────────────────────────────

group_files, _ = ls_nanoaod_files_groups(
    base_dir=BASE,
    groups=GROUPS,
    verbose=0,
)

dy_df  = ROOT.RDataFrame("Events", group_files["DY"]).Range(FILE_SIZE)
tt_df  = ROOT.RDataFrame("Events", group_files["TT"]).Range(FILE_SIZE)
sig_df = ROOT.RDataFrame("Events", group_files["HHbbtt"]).Range(FILE_SIZE)


# ──────────────────────────────────────────────────────────────────────────────
# Load real data via XCache
# ──────────────────────────────────────────────────────────────────────────────

YEARS = ["2024"]
RUNS = None  # set e.g. ["Run2024C"] for a single run

data_files, data_summary = load_scouting_data(years=YEARS, runs=RUNS)
data_jetmet_df = ROOT.RDataFrame("Events", data_files).Range(FILE_SIZE)
print(f"\nLoaded {len(data_files)} data files into RDataFrame")


# ──────────────────────────────────────────────────────────────────────────────
# Gen-level definitions (MC only)
# ──────────────────────────────────────────────────────────────────────────────

gen_decay_expr = (
    "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)"
    .format("GenPart")
)

dy_df  = dy_df.Define("GenDecay", gen_decay_expr)
tt_df  = tt_df.Define("GenDecay", gen_decay_expr)
sig_df = sig_df.Define("GenDecay", gen_decay_expr)

sig_df = (sig_df
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

dy_df  = dy_df.Define("w", "genWeight")
tt_df  = tt_df.Define("w", "genWeight")
sig_df = sig_df.Define("w", "genWeight")


# ──────────────────────────────────────────────────────────────────────────────
# AK4 jet pT columns
# ──────────────────────────────────────────────────────────────────────────────

for df_name in ["data_jetmet_df", "dy_df", "tt_df", "sig_df"]:
    df = locals()[df_name]
    df = (df
        .Define("ak4_pt0", "ScoutingPFJetRecluster_pt[0]")
        .Define("ak4_pt1", "ScoutingPFJetRecluster_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster_pt[3]")
    )
    locals()[df_name] = df


# ──────────────────────────────────────────────────────────────────────────────
# Acceptance cut
# ──────────────────────────────────────────────────────────────────────────────

base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_df_acc = data_jetmet_df.Filter(base_cut)
dy_df_acc   = dy_df.Filter(base_cut)
tt_df_acc   = tt_df.Filter(base_cut)
sig_df_acc  = sig_df.Filter(base_cut)

sig_df_hh = sig_df_acc.Filter("GenDecay.decayType==1")
sig_df_hm = sig_df_acc.Filter("GenDecay.decayType==2")
sig_df_he = sig_df_acc.Filter("GenDecay.decayType==3")


# ──────────────────────────────────────────────────────────────────────────────
# Cutflow + stacked plots per trigger
# ──────────────────────────────────────────────────────────────────────────────

setup_style()

PLOT_DIR = os.path.join(ANA_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

TRIG_LIST = [
    ("DST_JetHT",    DST_JetHT_expr),
    ("PARKING_HH",   PARKING_HH_expr),
]

for trig_name, trig in TRIG_LIST:
    print("=" * 63)
    print(trig)
    print("=" * 63)

    data_sel = data_df_acc.Filter(trig)
    dy_sel   = dy_df_acc.Filter(trig)
    tt_sel   = tt_df_acc.Filter(trig)
    sig_sel  = sig_df_acc.Filter(trig)

    # Data: unweighted event count
    data_n = data_sel.Count().GetValue()

    # MC weighted yields
    dy_y  = dy_sel.Sum("w").GetValue()
    tt_y  = tt_sel.Sum("w").GetValue()
    sig_y = sig_sel.Sum("w").GetValue()

    dy_yerr  = math.sqrt(dy_sel.Define("w2", "w*w").Sum("w2").GetValue())
    tt_yerr  = math.sqrt(tt_sel.Define("w2", "w*w").Sum("w2").GetValue())
    sig_yerr = math.sqrt(sig_sel.Define("w2", "w*w").Sum("w2").GetValue())

    print(f"Data events: {data_n}")
    print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
    print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
    print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")

    hh_sel = sig_sel.Filter("GenDecay.decayType==1")
    hm_sel = sig_sel.Filter("GenDecay.decayType==2")
    he_sel = sig_sel.Filter("GenDecay.decayType==3")

    y_hh = hh_sel.Sum("w").GetValue()
    y_hm = hm_sel.Sum("w").GetValue()
    y_he = he_sel.Sum("w").GetValue()

    if abs(sig_y - (y_hh + y_hm + y_he)) > 1e-6 * max(1.0, abs(sig_y)):
        print("WARNING: yield closure mismatch (float/weights?)")

    B = dy_y + tt_y
    S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")

    print(f"S/sqrt(B) (yields): {S_over_sqrtB:.6g}")
    print(f"y_hh:y_hm:y_he: {y_hh:.6g}:{y_hm:.6g}:{y_he:.6g}")
    print()

    # Stacked plot
    mc_items = [
        {"df": dy_sel,  "label": "DY",                       "color": "tab:orange"},
        {"df": tt_sel,  "label": "TT",                       "color": "tab:green"},
        {"df": hh_sel,  "label": r"$HH\to bb\tau_h\tau_h$",  "color": "tab:red"},
        {"df": hm_sel,  "label": r"$HH\to bb\tau_e\tau_h$",  "color": "tab:pink"},
        {"df": he_sel,  "label": r"$HH\to bb\tau_\mu\tau_h$", "color": "tab:purple"},
    ]

    fig, ax = plot_stacked_all_mc(
        data_df=data_sel,
        mc_items=mc_items,
        var="ak4_pt0",
        weight="w",
        nbins=20, xmin=0, xmax=1800,
        logy=True,
    )

    outpath = os.path.join(PLOT_DIR, f"stacked_ak4_pt0_{trig_name}.pdf")
    fig.savefig(outpath)
    print(f"Saved: {outpath}")
    plt.close(fig)
