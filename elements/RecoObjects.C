// elements/RecoObjects.C
//
// Reco-level object selection for all processes in the HH->bb tautau
// scouting analysis.
//
// Functions provided (all in namespace Ana):
//
//   SelectMuon()        — DY/TT: select best ScoutingMuonVtx
//   SelectElectron()    — DY/TT: select best ScoutingElectron
//   SelectDiMuon()      — DY:    OS di-muon pair closest to Z mass
//   SelectDiElectron()  — DY:    OS di-electron pair closest to Z mass
//   WTransverseMass()   — TT:    W transverse mass from lepton + MET
//   BvsAll()            — Signal/all: compute BvsAll b-tag discriminant per jet
//   BvsAllSortIdx()     — Signal/all: sort jet indices by BvsAll (descending)
//
// ALL cut values are explicit function parameters — no magic numbers.
// Values come from config/objects.yaml and config/acceptance.yaml,
// interpolated by cutflow_TrigEff.py into the Define() call strings.
//
// Load in cutflow_TrigEff.py after elements/GenMatching.C:
//   ROOT.gROOT.LoadMacro("elements/RecoObjects.C+")
//

#include "common.h"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Ana {


// ═════════════════════════════════════════════════════════════════════════════
//  Muon selection  (ScoutingMuonVtx collection)
// ═════════════════════════════════════════════════════════════════════════════
//
// Selects the highest-pT ScoutingMuonVtx passing all quality cuts.
// Returns collection index or -1 if none pass.
//
// Cut values — all from config/objects.yaml muon section:
//   pt_min    : muon.pt_min
//   eta_max   : muon.eta_max
//   iso_max   : muon.iso_max   (trackIso / pT)
//   dxy_max   : muon.dxy_max   [cm]
//   dz_max    : muon.dz_max    [cm]
//   chi2_max  : muon.chi2_max  (normalised track chi2)
//
// Usage:
//   mu = CUTS_OBJ["muon"]
//   df.Define("mu_idx",
//       f"Ana::SelectMuon(ScoutingMuonVtx_pt, ScoutingMuonVtx_eta, "
//       f"ScoutingMuonVtx_trackIso, ScoutingMuonVtx_trk_dxy, "
//       f"ScoutingMuonVtx_trk_dz, ScoutingMuonVtx_normchi2, "
//       f"{mu['pt_min']}f, {mu['eta_max']}f, {mu['iso_max']}f, "
//       f"{mu['dxy_max']}f, {mu['dz_max']}f, {mu['chi2_max']}f)")
//
int SelectMuon(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& trackIso,
    const ROOT::RVec<float>& trk_dxy,
    const ROOT::RVec<float>& trk_dz,
    const ROOT::RVec<float>& normchi2,
    float pt_min,
    float eta_max,
    float iso_max,
    float dxy_max,
    float dz_max,
    float chi2_max)
{
    int best = -1;
    float best_pt = -1.f;
    for (int i = 0; i < (int)pt.size(); ++i) {
        if (pt[i]                  < pt_min)  continue;
        if (std::abs(eta[i])       > eta_max) continue;
        if (pt[i] > 0 && trackIso[i] / pt[i] > iso_max)  continue;
        if (std::abs(trk_dxy[i])  > dxy_max) continue;
        if (std::abs(trk_dz[i])   > dz_max)  continue;
        if (normchi2[i]            > chi2_max) continue;
        if (pt[i] > best_pt) { best_pt = pt[i]; best = i; }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Electron selection  (ScoutingElectron collection)
// ═════════════════════════════════════════════════════════════════════════════
//
// Selects the highest-pT ScoutingElectron passing loose scouting-tuned cuts.
// ECAL crack (1.4442 < |eta| < 1.566) is excluded.
// Barrel and endcap have separate shower-shape thresholds.
//
// Cut values — all from config/objects.yaml electron section.
//
// Usage:
//   el = CUTS_OBJ["electron"]
//   df.Define("el_idx",
//       f"Ana::SelectElectron(ScoutingElectron_pt, ScoutingElectron_eta, "
//       f"ScoutingElectron_hOverE, ScoutingElectron_sigmaIetaIeta, "
//       f"ScoutingElectron_dEtaIn, ScoutingElectron_dPhiIn, "
//       f"ScoutingElectron_trackIso, "
//       f"ScoutingElectron_bestTrack_d0, ScoutingElectron_bestTrack_dz, "
//       f"{el['pt_min']}f, {el['eta_max']}f, {el['iso_max']}f, "
//       f"{el['d0_max']}f, {el['dz_max']}f, {el['hoe_max']}f, "
//       f"{el['sieta_barrel']}f, {el['sieta_endcap']}f, "
//       f"{el['deta_barrel']}f, {el['deta_endcap']}f, "
//       f"{el['dphi_barrel']}f, {el['dphi_endcap']}f)")
//
int SelectElectron(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& hOverE,
    const ROOT::RVec<float>& sigmaIetaIeta,
    const ROOT::RVec<float>& dEtaIn,
    const ROOT::RVec<float>& dPhiIn,
    const ROOT::RVec<float>& trackIso,
    const ROOT::RVec<float>& bestTrack_d0,
    const ROOT::RVec<float>& bestTrack_dz,
    float pt_min,
    float eta_max,
    float iso_max,
    float d0_max,
    float dz_max,
    float hoe_max,
    float sieta_barrel,
    float sieta_endcap,
    float deta_barrel,
    float deta_endcap,
    float dphi_barrel,
    float dphi_endcap)
{
    constexpr float CRACK_LO = 1.4442f;
    constexpr float CRACK_HI = 1.566f;

    int best = -1;
    float best_pt = -1.f;
    for (int i = 0; i < (int)pt.size(); ++i) {
        float aeta = std::abs(eta[i]);
        if (pt[i]  < pt_min)  continue;
        if (aeta   > eta_max) continue;
        if (aeta   > CRACK_LO && aeta < CRACK_HI) continue;  // ECAL crack

        bool isBarrel = (aeta <= CRACK_LO);
        float sieta_cut = isBarrel ? sieta_barrel : sieta_endcap;
        float deta_cut  = isBarrel ? deta_barrel  : deta_endcap;
        float dphi_cut  = isBarrel ? dphi_barrel  : dphi_endcap;

        if (hOverE[i]                      > hoe_max)   continue;
        if (sigmaIetaIeta[i]               > sieta_cut) continue;
        if (std::abs(dEtaIn[i])            > deta_cut)  continue;
        if (std::abs(dPhiIn[i])            > dphi_cut)  continue;
        if (pt[i] > 0 && trackIso[i]/pt[i] > iso_max)  continue;
        if (std::abs(bestTrack_d0[i])      > d0_max)    continue;
        if (std::abs(bestTrack_dz[i])      > dz_max)    continue;

        if (pt[i] > best_pt) { best_pt = pt[i]; best = i; }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Di-muon Z candidate  (opposite-sign, closest to Z mass)
// ═════════════════════════════════════════════════════════════════════════════
//
// Considers all OS muon pairs, returns the one with invariant mass
// closest to Z_MASS (91.19 GeV). Both muons must pass pt_min and eta_max.
//
// Usage:
//   mu = CUTS_OBJ["muon"]
//   df.Define("dimu",
//       f"Ana::SelectDiMuon(ScoutingMuonVtx_pt, ScoutingMuonVtx_eta, "
//       f"ScoutingMuonVtx_phi, ScoutingMuonVtx_m, ScoutingMuonVtx_charge, "
//       f"{mu['pt_min']}f, {mu['eta_max']}f)")
//   df.Define("mll_mumu", "dimu.mass")
//
DileptonPair SelectDiMuon(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    const ROOT::RVec<float>& mass,
    const ROOT::RVec<int>&   charge,
    float pt_min,
    float eta_max)
{
    DileptonPair best;
    float best_dm = std::numeric_limits<float>::max();

    for (int i = 0; i < (int)pt.size(); ++i) {
        if (pt[i] < pt_min || std::abs(eta[i]) > eta_max) continue;
        for (int j = i + 1; j < (int)pt.size(); ++j) {
            if (pt[j] < pt_min || std::abs(eta[j]) > eta_max) continue;
            if (charge[i] + charge[j] != 0) continue;  // require OS

            TLorentzVector v;
            v.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
            TLorentzVector v2;
            v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
            float m = (float)(v + v2).M();
            float dm = std::abs(m - Z_MASS);
            if (dm < best_dm) {
                best_dm = dm;
                best.i1 = i; best.i2 = j;
                best.mass = m;
                best.charge_sum = charge[i] + charge[j];
            }
        }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Di-electron Z candidate  (opposite-sign, closest to Z mass)
// ═════════════════════════════════════════════════════════════════════════════
//
// Usage:
//   el = CUTS_OBJ["electron"]
//   df.Define("diel",
//       f"Ana::SelectDiElectron(ScoutingElectron_pt, ScoutingElectron_eta, "
//       f"ScoutingElectron_phi, ScoutingElectron_bestTrack_charge, "
//       f"{el['pt_min']}f, {el['eta_max']}f)")
//   df.Define("mll_ee", "diel.mass")
//
// Note: ScoutingElectron mass is very small; we use electron mass = 0.000511 GeV.
//
DileptonPair SelectDiElectron(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    const ROOT::RVec<int>&   charge,
    float pt_min,
    float eta_max)
{
    constexpr float ELECTRON_MASS = 0.000511f;
    constexpr float CRACK_LO = 1.4442f;
    constexpr float CRACK_HI = 1.566f;

    DileptonPair best;
    float best_dm = std::numeric_limits<float>::max();

    for (int i = 0; i < (int)pt.size(); ++i) {
        float aetai = std::abs(eta[i]);
        if (pt[i] < pt_min || aetai > eta_max) continue;
        if (aetai > CRACK_LO && aetai < CRACK_HI) continue;

        for (int j = i + 1; j < (int)pt.size(); ++j) {
            float aetaj = std::abs(eta[j]);
            if (pt[j] < pt_min || aetaj > eta_max) continue;
            if (aetaj > CRACK_LO && aetaj < CRACK_HI) continue;
            if (charge[i] + charge[j] != 0) continue;

            TLorentzVector v, v2;
            v.SetPtEtaPhiM(pt[i], eta[i], phi[i], ELECTRON_MASS);
            v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], ELECTRON_MASS);
            float m = (float)(v + v2).M();
            float dm = std::abs(m - Z_MASS);
            if (dm < best_dm) {
                best_dm = dm;
                best.i1 = i; best.i2 = j;
                best.mass = m;
                best.charge_sum = charge[i] + charge[j];
            }
        }
    }
    return best;
}


// ═════════════════════════════════════════════════════════════════════════════
//  W transverse mass  (TT semi-leptonic)
// ═════════════════════════════════════════════════════════════════════════════
//
// MT = sqrt(2 * pT_lep * MET * (1 - cos(Δφ)))
// Returns -1.f if lepton_idx < 0 (no selected lepton).
//
// Usage (after SelectMuon define):
//   df.Define("mt_mu",
//       "Ana::WTransverseMass(mu_idx, ScoutingMuonVtx_pt, ScoutingMuonVtx_phi,"
//       " ScoutingMET_pt, ScoutingMET_phi)")
//
// Note: confirm MET branch names against your NanoAOD schema
//   (ScoutingMET_pt / ScoutingMET_phi  or  ScoutingPuppiMET_*).
//
float WTransverseMass(
    int   lepton_idx,
    const ROOT::RVec<float>& lep_pt,
    const ROOT::RVec<float>& lep_phi,
    float met_pt,
    float met_phi)
{
    if (lepton_idx < 0 || lepton_idx >= (int)lep_pt.size()) return -1.f;
    double dp = (double)lep_phi[lepton_idx] - (double)met_phi;
    while (dp >  M_PI) dp -= 2.0 * M_PI;
    while (dp < -M_PI) dp += 2.0 * M_PI;
    float dphi = (float)dp;
    return std::sqrt(2.f * lep_pt[lepton_idx] * met_pt * (1.f - std::cos(dphi)));
}


// ═════════════════════════════════════════════════════════════════════════════
//  BvsAll discriminant  (Signal + all processes)
// ═════════════════════════════════════════════════════════════════════════════
//
// Computes the BvsAll b-tag discriminant for each jet:
//   BvsAll = prob_b / (prob_b + prob_c + prob_cc + prob_g + prob_uds + prob_undef)
//
// Note: prob_bb is intentionally excluded — it is an AK8 fat-jet variable
// and is not meaningful for AK4 scouting jets.
//
// This mirrors the Define() expressions currently inlined in cutflow_TrigEff.py,
// moved here for reuse and clarity.
//
// Usage:
//   df.Define("pnet_BvsAll", "Ana::BvsAll("
//       "ScoutingPFJetRecluster_particleNet_prob_b, "
//       "ScoutingPFJetRecluster_particleNet_prob_c, "
//       "ScoutingPFJetRecluster_particleNet_prob_cc, "
//       "ScoutingPFJetRecluster_particleNet_prob_g, "
//       "ScoutingPFJetRecluster_particleNet_prob_uds, "
//       "ScoutingPFJetRecluster_particleNet_prob_undef)")
//
ROOT::RVec<float> BvsAll(
    const ROOT::RVec<float>& prob_b,
    const ROOT::RVec<float>& prob_c,
    const ROOT::RVec<float>& prob_cc,
    const ROOT::RVec<float>& prob_g,
    const ROOT::RVec<float>& prob_uds,
    const ROOT::RVec<float>& prob_undef)
{
    ROOT::RVec<float> result(prob_b.size(), 0.f);
    for (int i = 0; i < (int)prob_b.size(); ++i) {
        float denom = prob_b[i] + prob_c[i] + prob_cc[i]
                    + prob_g[i] + prob_uds[i] + prob_undef[i];
        result[i] = (denom > 0.f) ? prob_b[i] / denom : 0.f;
    }
    return result;
}


// ═════════════════════════════════════════════════════════════════════════════
//  BvsAll sort indices  (Signal + all processes)
// ═════════════════════════════════════════════════════════════════════════════
//
// Returns indices that sort jets by BvsAll in descending order
// (highest b-score first). Use index [0] for leading b-jet candidate, [1] for
// sub-leading, etc.
//
// Usage (after defining pnet_BvsAll):
//   df.Define("bsort_idx", "Ana::BvsAllSortIdx("
//       "ScoutingPFJetRecluster_particleNet_prob_b, "
//       "ScoutingPFJetRecluster_particleNet_prob_c, "
//       "ScoutingPFJetRecluster_particleNet_prob_cc, "
//       "ScoutingPFJetRecluster_particleNet_prob_g, "
//       "ScoutingPFJetRecluster_particleNet_prob_uds, "
//       "ScoutingPFJetRecluster_particleNet_prob_undef)")
//   df.Define("b0_idx", "(int)bsort_idx[0]")
//   df.Define("b1_idx", "(int)bsort_idx[1]")
//
ROOT::RVec<int> BvsAllSortIdx(
    const ROOT::RVec<float>& prob_b,
    const ROOT::RVec<float>& prob_c,
    const ROOT::RVec<float>& prob_cc,
    const ROOT::RVec<float>& prob_g,
    const ROOT::RVec<float>& prob_uds,
    const ROOT::RVec<float>& prob_undef)
{
    auto scores = BvsAll(prob_b, prob_c, prob_cc, prob_g, prob_uds, prob_undef);

    ROOT::RVec<int> idx(scores.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(),
              [&scores](int a, int b) { return scores[a] > scores[b]; });
    return idx;
}


// ═════════════════════════════════════════════════════════════════════════════
//  Dijet mass-window pairing
// ═════════════════════════════════════════════════════════════════════════════
//
// Finds the jet pair with invariant mass inside [mass_lo, mass_hi] that is
// closest to `target`. Optionally excludes jet indices (e.g. jets already
// used by the H→bb candidate when searching for H→ττ).
//
// Used for H→bb (mass window [100,150], target 125) and H→ττ (mass window
// [40,150], target 80).
//
// Cut values come from config/regions.yaml via Python interpolation:
//   reg = CUTS_REG["hbb"]
//   df.Define("hbb_pair",
//       f"Ana::findDijetInWindow(ScoutingPFJetRecluster_pt, "
//       f"ScoutingPFJetRecluster_eta, ScoutingPFJetRecluster_phi, "
//       f"ScoutingPFJetRecluster_mass, "
//       f"{reg['mass_lo']}f, {reg['mass_hi']}f, {reg['target']}f)")
//
DijetPair findDijetInWindow(
    const ROOT::RVec<float>& pt,
    const ROOT::RVec<float>& eta,
    const ROOT::RVec<float>& phi,
    const ROOT::RVec<float>& mass,
    float mass_lo, float mass_hi, float target,
    const std::vector<int>& exclude = {})
{
    DijetPair best;
    float bestDist = 1e9f;
    int n = (int)pt.size();

    for (int i = 0; i < n; ++i) {
        bool skip_i = false;
        for (int ex : exclude) { if (i == ex) { skip_i = true; break; } }
        if (skip_i) continue;

        for (int j = i + 1; j < n; ++j) {
            bool skip_j = false;
            for (int ex : exclude) { if (j == ex) { skip_j = true; break; } }
            if (skip_j) continue;

            TLorentzVector v1, v2;
            v1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]);
            v2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]);
            float mjj = (float)(v1 + v2).M();

            if (mjj >= mass_lo && mjj <= mass_hi) {
                float dist = std::abs(mjj - target);
                if (dist < bestDist) {
                    bestDist = dist;
                    best.i1 = i;
                    best.i2 = j;
                    best.mass = mjj;
                }
            }
        }
    }
    return best;
}

} // namespace Ana
