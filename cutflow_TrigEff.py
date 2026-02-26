#!/usr/bin/env python3
"""
Cutflow & trigger efficiency study for HH -> bb tautau (Scouting).
Converted from cutflow_TrigEff.ipynb.

Usage:
    python cutflow_TrigEff.py
"""

import os, sys, glob, math, subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
from matplotlib.lines import Line2D

# ── ROOT + Ana framework setup ───────────────────────────────────────────────
ANA_DIR = os.path.expanduser("/home/das214/HHtobbtautau/Scouting/CMSSW_15_0_15/src/Ana")
os.chdir(ANA_DIR)
sys.path.insert(0, ANA_DIR)

import ROOT
ROOT.gInterpreter.AddIncludePath(ANA_DIR)
ROOT.gROOT.LoadMacro("Particle.h+")
ROOT.gROOT.LoadMacro("HHbbtautauAnaElements.C+")
ROOT.gErrorIgnoreLevel = ROOT.kInfo
ROOT.ROOT.DisableImplicitMT()

import anaConfig
from ROOT import Ana

PLOT_DIR = os.path.join(ANA_DIR, "plots")
os.makedirs(PLOT_DIR, exist_ok=True)

FILE_SIZE = 10_000  # limit events for quick testing


# ── Helper: collect local NanoAOD files (for MC) ─────────────────────────────

def ls_nanoaod_files_from_fs(
    base_dir, sample, dataset_tag=None, production_id=None,
    pattern="*.root", recursive=True, verbose=0,
):
    base_dir = os.path.abspath(base_dir)
    sample_dir = os.path.join(base_dir, sample)
    if not os.path.isdir(sample_dir):
        avail = sorted(d for d in os.listdir(base_dir) if os.path.isdir(os.path.join(base_dir, d)))
        raise FileNotFoundError(
            f"Sample dir not found:\n  {sample_dir}\n"
            f"Available under {base_dir}:\n  - " + "\n  - ".join(avail)
        )

    dataset_tags = sorted(d for d in os.listdir(sample_dir) if os.path.isdir(os.path.join(sample_dir, d)))
    if not dataset_tags:
        raise FileNotFoundError(f"No dataset_tag dirs under:\n  {sample_dir}")
    if dataset_tag is None:
        if len(dataset_tags) == 1:
            dataset_tag = dataset_tags[0]
        else:
            raise ValueError(f"Multiple dataset_tag dirs under:\n  {sample_dir}\n  - " + "\n  - ".join(dataset_tags))
    elif dataset_tag not in dataset_tags:
        raise ValueError(f"dataset_tag '{dataset_tag}' not found under:\n  {sample_dir}")

    dataset_dir = os.path.join(sample_dir, dataset_tag)
    prod_ids = sorted(d for d in os.listdir(dataset_dir) if os.path.isdir(os.path.join(dataset_dir, d)))
    if not prod_ids:
        raise FileNotFoundError(f"No production_id dirs under:\n  {dataset_dir}")
    if production_id is None:
        production_id = prod_ids[-1]
    elif production_id not in prod_ids:
        raise ValueError(f"production_id '{production_id}' not found under:\n  {dataset_dir}")

    prod_dir = os.path.join(dataset_dir, production_id)
    if recursive:
        files = sorted(glob.glob(os.path.join(prod_dir, "**", pattern), recursive=True))
    else:
        files = sorted(glob.glob(os.path.join(prod_dir, pattern)))
    files = [f for f in files if os.path.isfile(f)]

    if verbose:
        print(f"  {sample}: {len(files)} files (tag={dataset_tag}, prod={production_id})")
    if not files:
        raise FileNotFoundError(f"No ROOT files matched in {prod_dir}")
    return files


def ls_nanoaod_files_groups(base_dir, groups, dataset_tag=None, production_id=None,
                            pattern="*.root", recursive=True, verbose=0, allow_missing=False):
    group_files = {}
    for gname, samples in groups.items():
        all_files = []
        for s in samples:
            try:
                files = ls_nanoaod_files_from_fs(
                    base_dir, s, dataset_tag=dataset_tag, production_id=production_id,
                    pattern=pattern, recursive=recursive, verbose=verbose,
                )
                all_files.extend(files)
            except Exception as e:
                if allow_missing:
                    print(f"[WARN] skipping '{s}': {e}")
                    continue
                raise
        group_files[gname] = all_files
    return group_files


# ── Helper: query DAS for data files via XCache ──────────────────────────────

XCACHE_PREFIX = "root://xcache.cms.rcac.purdue.edu/"
DASGOCLIENT = "/cvmfs/cms.cern.ch/common/dasgoclient"


def das_files(dataset, prefix=XCACHE_PREFIX):
    """Query DAS for file list and prepend XCache prefix for remote access."""
    result = subprocess.run(
        [DASGOCLIENT, "-query", f"file dataset={dataset}"],
        capture_output=True, text=True, check=True,
    )
    files = [prefix + f.strip() for f in result.stdout.strip().split("\n") if f.strip()]
    return sorted(files)


# ── Plotting helpers ─────────────────────────────────────────────────────────

def th1_to_np(h):
    nb = h.GetNbinsX()
    edges = np.array([h.GetXaxis().GetBinLowEdge(1 + i) for i in range(nb)] + [h.GetXaxis().GetBinUpEdge(nb)])
    vals = np.array([h.GetBinContent(1 + i) for i in range(nb)])
    errs = np.array([h.GetBinError(1 + i) for i in range(nb)])
    return edges, vals, errs


def teff_to_np(teff):
    htot = teff.GetTotalHistogram()
    nb = htot.GetNbinsX()
    x, y, xerr, yerr_lo, yerr_hi = [], [], [], [], []
    for i in range(1, nb + 1):
        if htot.GetBinContent(i) <= 0:
            continue
        xc = htot.GetXaxis().GetBinCenter(i)
        hw = 0.5 * htot.GetXaxis().GetBinWidth(i)
        eff = float(np.clip(teff.GetEfficiency(i), 0.0, 1.0))
        elo = float(np.clip(teff.GetEfficiencyErrorLow(i), 0.0, 1.0))
        ehi = float(np.clip(teff.GetEfficiencyErrorUp(i), 0.0, 1.0))
        x.append(xc); xerr.append(hw)
        y.append(eff); yerr_lo.append(elo); yerr_hi.append(ehi)
    return np.array(x), np.array(y), np.array(xerr), np.array(yerr_lo), np.array(yerr_hi)


def plot_stacked_all_mc(data_df, mc_items, *, var="ak4_pt0", weight="w",
                        nbins=20, xmin=0, xmax=1800, logy=True, sort_mc_by_yield=True):
    uid = ROOT.TUUID().AsString().replace("-", "_")
    h_data_ptr = data_df.Histo1D((f"h_data_{uid}", f";{var};Events", nbins, xmin, xmax), var)
    h_mc_ptrs = [
        item["df"].Histo1D((f"h_mc_{i}_{uid}", f";{var};Events", nbins, xmin, xmax), var, weight)
        for i, item in enumerate(mc_items)
    ]
    ROOT.RDF.RunGraphs([h_data_ptr] + h_mc_ptrs)

    edges, data_vals, data_errs = th1_to_np(h_data_ptr.GetValue())
    centers = 0.5 * (edges[:-1] + edges[1:])
    widths = np.diff(edges)

    mc_components = []
    for ptr, item in zip(h_mc_ptrs, mc_items):
        _, vals_raw, _ = th1_to_np(ptr.GetValue())
        if np.any(vals_raw < 0):
            print(f"  -> WARNING: '{item['label']}' has {np.sum(vals_raw < 0)} negative bin(s), clipping to 0.")
        vals = np.maximum(vals_raw, 0.0)
        mc_components.append({"item": item, "vals": vals, "yield": float(np.sum(vals))})

    if sort_mc_by_yield:
        mc_components.sort(key=lambda x: x["yield"])

    fig, ax = plt.subplots(figsize=(8, 6))
    bottom = np.zeros_like(centers, dtype=float)
    for comp in mc_components:
        item = comp["item"]
        ax.bar(centers, comp["vals"], width=widths, bottom=bottom, align="center",
               label=item["label"], color=item["color"], edgecolor="black", linewidth=0.2)
        bottom += comp["vals"]

    ax.errorbar(centers, data_vals, yerr=data_errs, fmt="o", color="black",
                label="Data", ms=4, capsize=2, linewidth=1, zorder=10)
    ax.set_xlabel(var); ax.set_ylabel("Events")
    ax.grid(True, axis="y", alpha=0.25); ax.legend(ncol=2)
    if logy:
        ax.set_yscale("log")
        ymax = max(float(np.max(bottom)) if len(bottom) else 1.0,
                   float(np.max(data_vals)) if len(data_vals) else 1.0)
        ax.set_ylim(0.5, max(10.0, 5.0 * ymax))
    return fig, ax


def plot_yield_overlay(dfs, labels, *, var=None, weight="w", pre_expr=None,
                       nbins=20, xmin=200.0, xmax=1800.0,
                       cms_year="2024", cms_com="13.6", cms_label_text="Simulation",
                       figsize=(10, 7), xlabel=None, ylabel="Weighted yields",
                       density=False, logy=False):
    h_ptrs = []
    for i, df in enumerate(dfs):
        df_sel = df if (pre_expr is None or str(pre_expr).strip() in ["", "1"]) else df.Filter(pre_expr)
        h_ptrs.append(df_sel.Histo1D(
            (f"h_{i}", f";{var};Yields", int(nbins), float(xmin), float(xmax)), var, weight
        ))
    try:
        ROOT.RDF.RunGraphs(h_ptrs)
    except Exception:
        pass

    fig, ax = plt.subplots(figsize=figsize)
    line_objs = []
    for hptr in h_ptrs:
        edges, vals, errs = th1_to_np(hptr.GetValue())
        if density:
            bw = np.diff(edges); area = float(np.sum(vals * bw))
            if area > 0: vals = vals / area; errs = errs / area
        (ln,) = ax.step(edges, np.r_[vals, vals[-1]], where="post", linewidth=2.5)
        line_objs.append(ln)

    ax.set_xlabel(xlabel if xlabel else f"{var} [arb]")
    ax.set_ylabel("Density" if density else ylabel)
    if logy: ax.set_yscale("log")
    hep.cms.label(cms_label_text, data=False, ax=ax, year=str(cms_year), com=str(cms_com))
    ax.legend(line_objs, labels, loc="best", frameon=True)
    return fig, ax


def plot_trigger_eff_overlay_channels(
    channels, triggers, *, var="gen_mHH", weight="w", pre_expr=None,
    nbins=20, xmin=200.0, xmax=1800.0, cms_year="2024", cms_com="13.6",
    cms_label_text="", figsize=(10, 10), pre_label="Preselection",
    trigger_colors=None, pre_key="pre", xlabel=None,
):
    if trigger_colors is None:
        trigger_colors = {}
    color_map = {pre_key: trigger_colors.get(pre_key, "tab:blue")}
    for t in triggers:
        k = t["key"]
        color_map[k] = trigger_colors.get(k, {"nom": "tab:orange", "park": "tab:green",
                                                "all": "tab:red"}.get(k, "black"))

    h_pre_ptrs = []
    h_num_ptrs = {t["key"]: [] for t in triggers}
    for i, ch in enumerate(channels):
        df = ch["df"]
        df_pre = df if (pre_expr is None or str(pre_expr).strip() in ["", "1"]) else df.Filter(pre_expr)
        h_pre_ptrs.append(df_pre.Histo1D(
            (f"h_pre_{i}", f";{var} [GeV];Events", nbins, xmin, xmax), var, weight))
        for t in triggers:
            k = t["key"]
            h_num_ptrs[k].append(df_pre.Filter(t["expr"]).Histo1D(
                (f"h_{k}_{i}", f";{var} [GeV];Events", nbins, xmin, xmax), var, weight))

    all_ptrs = list(h_pre_ptrs)
    for k in h_num_ptrs:
        all_ptrs.extend(h_num_ptrs[k])
    try:
        ROOT.RDF.RunGraphs(all_ptrs)
    except Exception:
        pass

    per_ch = []
    for i, ch in enumerate(channels):
        Hpre = h_pre_ptrs[i].GetValue()
        Hnums = {t["key"]: h_num_ptrs[t["key"]][i].GetValue() for t in triggers}
        effs = {}
        for t in triggers:
            te = ROOT.TEfficiency(Hnums[t["key"]], Hpre)
            te.SetUseWeightedEvents(True)
            effs[t["key"]] = te
        edges, pre_vals, pre_errs = th1_to_np(Hpre)
        per_ch.append({
            "key": ch.get("key", f"ch{i}"), "label": ch.get("label", f"ch{i}"),
            "ls": ch.get("ls", "-"), "marker": ch.get("marker", "o"),
            "edges": edges, "pre": (pre_vals, pre_errs),
            "num": {k: th1_to_np(Hnums[k])[1:] for k in Hnums},
            "eff_np": {k: teff_to_np(effs[k]) for k in effs},
        })

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax = fig.add_subplot(gs[0])
    rax = fig.add_subplot(gs[1], sharex=ax)

    for ch in per_ch:
        edges = ch["edges"]
        pre_vals, _ = ch["pre"]
        ax.step(edges, np.r_[pre_vals, pre_vals[-1]], where="post", linewidth=2.5,
                linestyle=ch["ls"], color=color_map[pre_key])
        for t in triggers:
            vals, _ = ch["num"][t["key"]]
            ax.step(edges, np.r_[vals, vals[-1]], where="post",
                    linewidth=float(t.get("lw", 2.5)), linestyle=ch["ls"], color=color_map[t["key"]])

    ax.set_ylabel("Yields"); ax.tick_params(labelbottom=False)
    hep.cms.label(cms_label_text, data=False, ax=ax, year=str(cms_year), com=str(cms_com))

    for ch in per_ch:
        for t in triggers:
            x, y, xerr, ylo, yhi = ch["eff_np"][t["key"]]
            rax.errorbar(x, y, xerr=xerr, yerr=np.vstack([ylo, yhi]),
                         fmt=ch["marker"], markersize=float(t.get("markersize", 5.0)),
                         capsize=2, linewidth=1.2, linestyle="none", color=color_map[t["key"]])

    rax.set_xlabel(xlabel if xlabel else f"{var} [GeV]")
    rax.set_ylabel("Trigger Efficiency"); rax.set_ylim(0.0, 1.1)
    rax.grid(True, axis="y", alpha=0.3)

    ch_handles = [Line2D([0], [0], color="black", linestyle=ch["ls"], marker=ch["marker"],
                         markersize=6, linewidth=2, label=ch["label"]) for ch in channels]
    trig_handles = [Line2D([0], [0], color=color_map[pre_key], linestyle="-", linewidth=2.5, label=pre_label)]
    trig_handles += [Line2D([0], [0], color=color_map[t["key"]], linestyle="-",
                            linewidth=float(t.get("lw", 2.5)), label=t.get("label", t["key"])) for t in triggers]
    leg1 = ax.legend(handles=ch_handles, title="Decay channels", loc="upper left", frameon=True)
    ax.add_artist(leg1)
    ax.legend(handles=trig_handles, loc="upper right", frameon=True)

    return fig, ax, rax, {"per_channel": per_ch, "color_map": color_map}


# ── Matplotlib style ─────────────────────────────────────────────────────────

plt.style.use(hep.style.CMS)
mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.labelsize": 13, "axes.titlesize": 13,
    "xtick.labelsize": 13, "ytick.labelsize": 13,
    "legend.fontsize": 13, "legend.title_fontsize": 13,
    "axes.linewidth": 1.2,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.major.size": 6, "ytick.major.size": 6,
    "xtick.minor.size": 3, "ytick.minor.size": 3,
})


# ── Trigger bit definitions ──────────────────────────────────────────────────

def or_expr(bits):
    return "(" + " || ".join(bits) + ")"


DST_MU_BIT = ["DST_PFScouting_DatasetMuon", "DST_PFScouting_DoubleMuon", "DST_PFScouting_SingleMuon"]
DST_EL_BITS = ["DST_PFScouting_DoubleEG", "DST_PFScouting_SinglePhotonEB"]
DST_JetHT_BITS = ["DST_PFScouting_JetHT"]

DST_MU_expr = or_expr(DST_MU_BIT)
DST_EL_expr = or_expr(DST_EL_BITS)
DST_JetHT_expr = or_expr(DST_JetHT_BITS)
DST_ALL_expr = f"(({DST_MU_expr}) || ({DST_EL_expr}) || ({DST_JetHT_expr}))"

PARKING_HH_BITS = [
    "HLT_PFHT280_QuadPFJet30", "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55",
    "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p60", "HLT_PFHT280_QuadPFJet35_PNet2BTagMean0p60",
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40",
    "HLT_PFHT330PT30_QuadPFJet_75_60_45_40_TriplePFBTagDeepJet_4p5",
    "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70",
    "HLT_PFHT400_SixPFJet32", "HLT_PFHT400_SixPFJet32_PNet2BTagMean0p50",
    "HLT_PFHT450_SixPFJet36", "HLT_PFHT450_SixPFJet36_PNetBTag0p35",
]

PARKING_MUON_BITS_2023 = [
    "HLT_Dimuon0_Jpsi3p5_Muon2", "HLT_Dimuon0_Jpsi_NoVertexing",
    "HLT_Dimuon0_Jpsi_NoVertexing_NoOS", "HLT_Dimuon0_LowMass",
    "HLT_Dimuon0_Upsilon_L1_4p5er2p0", "HLT_Dimuon0_Upsilon_NoVertexing",
    "HLT_Dimuon10_Upsilon_y1p4", "HLT_Dimuon12_Upsilon_y1p4",
    "HLT_Dimuon14_Phi_Barrel_Seagulls", "HLT_Dimuon14_PsiPrime",
    "HLT_Dimuon18_PsiPrime", "HLT_Dimuon18_PsiPrime_noCorrL1",
    "HLT_Dimuon24_Phi_noCorrL1", "HLT_Dimuon24_Upsilon_noCorrL1",
    "HLT_Dimuon25_Jpsi", "HLT_Dimuon25_Jpsi_noCorrL1",
    "HLT_DoubleMu2_Jpsi_DoubleTrk1_Phi1p05",
    "HLT_DoubleMu3_DoubleEle7p5_CaloIdL_TrackIdL_Upsilon",
    "HLT_DoubleMu3_TkMu_DsTau3Mu", "HLT_DoubleMu3_Trk_Tau3mu",
    "HLT_DoubleMu3_Trk_Tau3mu_NoL1Mass", "HLT_DoubleMu4_3_Bs",
    "HLT_DoubleMu4_3_Displaced_Photon4_BsToMMG", "HLT_DoubleMu4_3_Jpsi",
    "HLT_DoubleMu4_3_LowMass", "HLT_DoubleMu4_3_Photon4_BsToMMG",
    "HLT_DoubleMu4_JpsiTrkTrk_Displaced", "HLT_DoubleMu4_JpsiTrk_Bc",
    "HLT_DoubleMu4_Jpsi_Displaced", "HLT_DoubleMu4_Jpsi_NoVertexing",
    "HLT_DoubleMu4_LowMass_Displaced", "HLT_DoubleMu4_MuMuTrk_Displaced",
    "HLT_DoubleMu5_Upsilon_DoubleEle3_CaloIdL_TrackIdL",
    "HLT_Mu25_TkMu0_Phi", "HLT_Mu30_TkMu0_Psi", "HLT_Mu30_TkMu0_Upsilon",
    "HLT_Mu7p5_L2Mu2_Jpsi", "HLT_Mu7p5_L2Mu2_Upsilon",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15", "HLT_Tau3Mu_Mu7_Mu1_TkMu1_IsoTau15_Charge1",
    "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15", "HLT_Tau3Mu_Mu7_Mu1_TkMu1_Tau15_Charge1",
    "HLT_Trimuon5_3p5_2_Upsilon_Muon", "HLT_TrimuonOpen_5_3p5_2_Upsilon_Muon",
]

PARKING_ELECTRON_BITS_2023 = [
    "HLT_DoubleEle10_eta1p22_mMax6", "HLT_DoubleEle6p5_eta1p22_mMax6",
    "HLT_DoubleEle8_eta1p22_mMax6", "HLT_SingleEle8", "HLT_SingleEle8_SingleEGL1",
]

PARKING_HH_expr = or_expr(PARKING_HH_BITS)
PARKING_MU_expr = or_expr(PARKING_MUON_BITS_2023)
PARKING_EG_expr = or_expr(PARKING_ELECTRON_BITS_2023)
PARKING_LEPT_expr = f"({PARKING_MU_expr} || {PARKING_EG_expr})"


# =============================================================================
# LOAD DATA (via XCache / DAS — no local download needed)
# =============================================================================

print("=" * 63)
print("Loading DATA from DAS via XCache...")
print("=" * 63)

DAS_DATASET_2024C = "/ScoutingPFRun3/Run2024C-ScoutNano-v1/NANOAOD"
run2024c_files = das_files(DAS_DATASET_2024C)
print(f"Run2024C files from DAS: {len(run2024c_files)}")
print(f"Example: {run2024c_files[0]}")

data_jetmet_df = ROOT.RDataFrame("Events", run2024c_files).Range(FILE_SIZE)


# =============================================================================
# LOAD MC (local EOS)
# =============================================================================

print("\n" + "=" * 63)
print("Loading MC from local EOS...")
print("=" * 63)

MC_BASE = os.path.join(
    "/eos/purdue/store/user/arghyara", "production", "Scouting",
    "NanoAODv15Scouting24", "mcscouting_2024",
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
SIG = ["GluGluHHto2B2Tau_Par-c2-0p00-kl-0p00-kt-1p00_TuneCP5_13p6TeV_powheg-pythia8"]

group_files = ls_nanoaod_files_groups(
    MC_BASE, {"DY": DY, "TT": TT, "HHbbtt": SIG}, verbose=1,
)

dy_df = ROOT.RDataFrame("Events", group_files["DY"]).Range(FILE_SIZE)
tt_df = ROOT.RDataFrame("Events", group_files["TT"]).Range(FILE_SIZE)
sig_df = ROOT.RDataFrame("Events", group_files["HHbbtt"]).Range(FILE_SIZE)


# =============================================================================
# DEFINE COLUMNS
# =============================================================================

# GenDecay for MC
dy_df = dy_df.Define("GenDecay", "Ana::DecayGenMatching(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
tt_df = tt_df.Define("GenDecay", "Ana::DecayGenMatching(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")
sig_df = sig_df.Define("GenDecay", "Ana::DecayGenMatching(GenPart_pdgId, GenPart_genPartIdxMother, GenPart_statusFlags)")

sig_df = (sig_df
    .Define("genHbb_idx",      "(int)GenDecay.Htob")
    .Define("genHtautau_idx",  "(int)GenDecay.Htotau")
    .Define("genHbb_p4",       "Ana::getP4(genHbb_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("genHtautau_p4",   "Ana::getP4(genHtautau_idx, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    .Define("gen_mHH",         "(genHbb_p4 + genHtautau_p4).M()")
    .Define("gen_pt_Hbb",      "genHbb_p4.Pt()")
    .Define("gen_pt_Htautau",  "genHtautau_p4.Pt()")
    .Define("gen_eta_Hbb",     "genHbb_p4.Eta()")
    .Define("gen_eta_Htautau", "genHtautau_p4.Eta()")
)

dy_df = dy_df.Define("w", "genWeight")
tt_df = tt_df.Define("w", "genWeight")
sig_df = sig_df.Define("w", "genWeight")

# AK4 jet pT aliases — data uses Jet_pt, MC uses ScoutingPFJetRecluster_pt
data_jetmet_df = (data_jetmet_df
    .Define("ak4_pt0", "Jet_pt[0]").Define("ak4_pt1", "Jet_pt[1]")
    .Define("ak4_pt2", "Jet_pt[2]").Define("ak4_pt3", "Jet_pt[3]")
)
for _df_name in ["dy_df", "tt_df", "sig_df"]:
    _df = locals()[_df_name]
    _df = (_df
        .Define("ak4_pt0", "ScoutingPFJetRecluster_pt[0]")
        .Define("ak4_pt1", "ScoutingPFJetRecluster_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster_pt[3]")
    )
    locals()[_df_name] = _df


# =============================================================================
# CUTFLOW: acceptance + trigger yields
# =============================================================================

base_cut = "ak4_pt0 > 20.0 && ak4_pt1 > 20.0 && ak4_pt2 > 20.0 && ak4_pt3 > 20.0"

data_df_acc = data_jetmet_df.Filter(base_cut)
dy_df_acc = dy_df.Filter(base_cut)
tt_df_acc = tt_df.Filter(base_cut)
sig_df_acc = sig_df.Filter(base_cut)

sig_df_hh = sig_df_acc.Filter("GenDecay.decayType==1")
sig_df_hm = sig_df_acc.Filter("GenDecay.decayType==2")
sig_df_he = sig_df_acc.Filter("GenDecay.decayType==3")

TRIG_LIST = [DST_JetHT_expr, PARKING_HH_expr]

plot_idx = 0
for trig in TRIG_LIST:
    print("\n" + "=" * 63)
    print(trig)
    print("=" * 63)

    data_sel = data_df_acc.Filter(trig)
    dy_sel = dy_df_acc.Filter(trig)
    tt_sel = tt_df_acc.Filter(trig)
    sig_sel = sig_df_acc.Filter(trig)

    data_n = data_sel.Count().GetValue()
    dy_y = dy_sel.Sum("w").GetValue()
    tt_y = tt_sel.Sum("w").GetValue()
    sig_y = sig_sel.Sum("w").GetValue()

    dy_yerr = math.sqrt(dy_sel.Define("w2", "w*w").Sum("w2").GetValue())
    tt_yerr = math.sqrt(tt_sel.Define("w2", "w*w").Sum("w2").GetValue())
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

    B = dy_y + tt_y
    S_over_sqrtB = sig_y / math.sqrt(B) if B > 0 else float("nan")
    print(f"S/sqrt(B):   {S_over_sqrtB:.6g}")
    print(f"y_hh:y_hm:y_he = {y_hh:.6g} : {y_hm:.6g} : {y_he:.6g}")

    # ── Stacked data+MC plot ──
    mc_items = [
        {"df": dy_sel, "label": "DY", "color": "tab:orange"},
        {"df": tt_sel, "label": "TT", "color": "tab:green"},
        {"df": hh_sel, "label": r"HH$\to$bb$\tau_h\tau_h$", "color": "tab:red"},
        {"df": hm_sel, "label": r"HH$\to$bb$\tau_e\tau_h$", "color": "tab:pink"},
        {"df": he_sel, "label": r"HH$\to$bb$\tau_\mu\tau_h$", "color": "tab:purple"},
    ]
    fig, ax = plot_stacked_all_mc(data_sel, mc_items, var="ak4_pt0",
                                  weight="w", nbins=20, xmin=0, xmax=1800, logy=True)
    fig.savefig(os.path.join(PLOT_DIR, f"stacked_ak4pt0_trig{plot_idx}.pdf"))
    plt.close(fig)
    print(f"  -> saved plots/stacked_ak4pt0_trig{plot_idx}.pdf")

    # ── Yield overlay (data + MC) ──
    data_plot = data_sel.Define("w_plot", "1.0")
    dy_plot = dy_sel.Define("w_plot", "w")
    tt_plot = tt_sel.Define("w_plot", "w")
    sig_plot = sig_sel.Define("w_plot", "w")
    fig, ax = plot_yield_overlay(
        dfs=[data_plot, dy_plot, tt_plot, sig_plot],
        labels=["Data", "DY", "TT", "HH signal"],
        var="ak4_pt0", weight="w_plot",
        nbins=20, xmin=0, xmax=1800, logy=True,
        cms_label_text="Data + Simulation",
    )
    fig.savefig(os.path.join(PLOT_DIR, f"overlay_ak4pt0_trig{plot_idx}.pdf"))
    plt.close(fig)
    print(f"  -> saved plots/overlay_ak4pt0_trig{plot_idx}.pdf")

    plot_idx += 1


# =============================================================================
# TRIGGER EFFICIENCY vs gen_mHH
# =============================================================================

print("\n" + "=" * 63)
print("Trigger efficiency overlay (gen_mHH)")
print("=" * 63)

sig_df_eff = (sig_df
    .Filter("fabs(gen_eta_Hbb) < 2.5 && fabs(gen_eta_Htautau) < 2.5")
    .Filter("gen_pt_Hbb > 20.0 && gen_pt_Htautau > 40.0")
)

channels = [
    {"key": "hh", "label": r"HH$\rightarrow$bb$\tau_h\tau_h$",
     "df": sig_df_eff.Filter("GenDecay.decayType==1"), "ls": "-", "marker": "o"},
    {"key": "hm", "label": r"HH$\rightarrow$bb$\tau_e\tau_h$",
     "df": sig_df_eff.Filter("GenDecay.decayType==2"), "ls": "--", "marker": "s"},
    {"key": "he", "label": r"HH$\rightarrow$bb$\tau_\mu\tau_h$",
     "df": sig_df_eff.Filter("GenDecay.decayType==3"), "ls": ":", "marker": "^"},
]

triggers = [
    {"key": "dst_all", "expr": DST_ALL_expr, "label": "Data Scouting"},
    {"key": "park", "expr": PARKING_HH_expr, "label": "Data Parking"},
]

trigger_colors = {"pre": "tab:blue", "nom": "tab:orange", "park": "tab:green", "dst_all": "tab:red"}

fig, ax, rax, out = plot_trigger_eff_overlay_channels(
    channels, triggers, var="gen_mHH",
    trigger_colors=trigger_colors,
    nbins=50, xmin=0, xmax=1200, figsize=(10, 10),
)
fig.savefig(os.path.join(PLOT_DIR, "trigger_eff_gen_mHH.pdf"))
plt.close(fig)
print("  -> saved plots/trigger_eff_gen_mHH.pdf")

print("\nDone! All plots saved to:", PLOT_DIR)
