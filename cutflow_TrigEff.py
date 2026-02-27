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
    BASE, GROUPS, XSEC, FILE_SIZE,
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


LUMI = 1.0  # fb^-1


_, group_files_by_sample = ls_nanoaod_files_groups(
    base_dir=BASE,
    groups=GROUPS,
    verbose=0,
)

# COmment: Add teh ht turn on curve


def make_mc_df(files, sample_name):
    """Create an RDataFrame with weight 'w' normalised to LUMI.

    w_i = genWeight * (xsec_pb * LUMI_fb * 1000) / sum(genWeight)
    The factor 1000 converts fb^-1 → pb^-1.
    """
    df = ROOT.RDataFrame("Events", files).Range(FILE_SIZE)
    sum_genw = df.Sum("genWeight").GetValue()
    xsec_pb = XSEC[sample_name]
    scale = xsec_pb * LUMI * 1000.0 / sum_genw
    print(f"  {sample_name[:60]:60s}  xsec={xsec_pb:.4g} pb  "
          f"sum_genw={sum_genw:.4g}  scale={scale:.4e}")
    return df.Define("w", f"genWeight * {scale:.10e}")


print(f"\nNormalising MC to L = {LUMI} fb^-1  (FILE_SIZE = {FILE_SIZE})")
mc = {}   # sample_name -> RDataFrame (with column "w")
for group_name, samples_dict in group_files_by_sample.items():
    print(f"\n[{group_name}]")
    for sample_name, files in samples_dict.items():
        mc[sample_name] = make_mc_df(files, sample_name)

# Convenience: lists of sample names per group
dy_samples  = list(group_files_by_sample["DY"].keys())
tt_samples  = list(group_files_by_sample["TT"].keys())
sig_samples = list(group_files_by_sample["HHbbtt"].keys())



YEARS = ["2024"]
RUNS = None  # set e.g. ["Run2024C"] for a single run

data_files, data_summary = load_scouting_data(years=YEARS, runs=RUNS)
data_df = ROOT.RDataFrame("Events", data_files).Range(FILE_SIZE)
print(f"\nLoaded {len(data_files)} data files into RDataFrame")


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

def define_ak4_pt(df):
    return (df
        .Define("ak4_pt0", "ScoutingPFJetRecluster_pt[0]")
        .Define("ak4_pt1", "ScoutingPFJetRecluster_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster_pt[3]")
    )

data_df = define_ak4_pt(data_df)
for name in mc:
    mc[name] = define_ak4_pt(mc[name])


base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_acc = data_df.Filter(base_cut)
mc_acc = {name: df.Filter(base_cut) for name, df in mc.items()}

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

    # Apply trigger
    data_sel = data_acc.Filter(trig)
    mc_sel = {name: df.Filter(trig) for name, df in mc_acc.items()}

    # Data: unweighted event count
    data_n = data_sel.Count().GetValue()

    # MC weighted yields per group
    def group_yield(samples):
        y = sum(mc_sel[s].Sum("w").GetValue() for s in samples)
        yerr = math.sqrt(sum(
            mc_sel[s].Define("w2", "w*w").Sum("w2").GetValue()
            for s in samples
        ))
        return y, yerr

    dy_y,  dy_yerr  = group_yield(dy_samples)
    tt_y,  tt_yerr  = group_yield(tt_samples)
    sig_y, sig_yerr = group_yield(sig_samples)

    print(f"Data events: {data_n}")
    print(f"DY yield:    {dy_y:.6g} +/- {dy_yerr:.3g}")
    print(f"TT yield:    {tt_y:.6g} +/- {tt_yerr:.3g}")
    print(f"SIG yield:   {sig_y:.6g} +/- {sig_yerr:.3g}")

    # Signal decay channels
    hh = [mc_sel[s].Filter("GenDecay.decayType==1") for s in sig_samples]
    hm = [mc_sel[s].Filter("GenDecay.decayType==2") for s in sig_samples]
    he = [mc_sel[s].Filter("GenDecay.decayType==3") for s in sig_samples]

    y_hh = sum(d.Sum("w").GetValue() for d in hh)
    y_hm = sum(d.Sum("w").GetValue() for d in hm)
    y_he = sum(d.Sum("w").GetValue() for d in he)

    if abs(sig_y - (y_hh + y_hm + y_he)) > 1e-6 * max(1.0, abs(sig_y)):
        print("WARNING: yield closure mismatch (float/weights?)")

    B = dy_y + tt_y
    S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")

    print(f"S/sqrt(B) (yields): {S_over_sqrtB:.6g}")
    print(f"y_hh:y_hm:y_he: {y_hh:.6g}:{y_hm:.6g}:{y_he:.6g}")
    print()

    # Stacked plot  (MC only — no data overlay)
    mc_items = [
        {"dfs": [mc_sel[s] for s in dy_samples],
         "label": "DY",                        "color": "tab:orange"},
        {"dfs": [mc_sel[s] for s in tt_samples],
         "label": "TT",                        "color": "tab:green"},
        {"dfs": hh, "label": r"$HH\to bb\tau_h\tau_h$",   "color": "tab:red"},
        {"dfs": hm, "label": r"$HH\to bb\tau_e\tau_h$",   "color": "tab:pink"},
        {"dfs": he, "label": r"$HH\to bb\tau_\mu\tau_h$",  "color": "tab:purple"},
    ]

    fig, ax = plot_stacked_all_mc(
        mc_items=mc_items,
        var="ak4_pt0",
        weight="w",
        nbins=36, xmin=0, xmax=1800,
        logy=True,
    )

    outpath = os.path.join(PLOT_DIR, f"stacked_ak4_pt0_{trig_name}.pdf")
    fig.savefig(outpath)
    print(f"Saved: {outpath}")
    plt.close(fig)
