// elements/GenMatching.C
//
// Gen-level particle matching for all MC sample types used in the
// HH->bb tautau scouting analysis.
//
// Functions provided (all in namespace Ana):
//
//   ZGenMatching()   — DY:     find gen Z boson, classify Z->ll/tautau decay
//   TopGenMatching() — TT:     find gen t+tbar, classify W decay topology
//   HHGenMatching()  — Signal: find gen HH->bbtautau, classify tau decay channels
//
// All functions are safe to call on any MC sample:
//   - return sentinel values (-1 / decayMode=0) if the target particle is absent
//   - no hard-coded thresholds; cut values are passed as parameters from Python
//     (read from config/objects.yaml and config/regions.yaml via cutflow_TrigEff.py)
//
// This file does NOT include HHbbtautauAnaElements.C. It is a standalone
// reimplementation that coexists with the legacy code when both .so files
// are loaded via ROOT::gROOT->LoadMacro().
//
// ACLiC note: LoadMacro cd-s into elements/ before compiling, so
// #include "common.h" resolves correctly (not "elements/common.h").
//
// Load order in cutflow_TrigEff.py (after Particle.h and HHbbtautauAnaElements.C):
//   ROOT.gROOT.LoadMacro("elements/GenMatching.C+")
//

#include "common.h"
#include <algorithm>
#include <vector>

namespace Ana {

// ═════════════════════════════════════════════════════════════════════════════
//  Internal helpers (file-local, not exposed in the public API)
// ═════════════════════════════════════════════════════════════════════════════

namespace {

// Local copies of angular/status-flag helpers to avoid symbol conflicts
// with HHbbtautauAnaElements.C when both .so files are loaded.
double deltaPhi_(double phi1, double phi2) {
    double dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2.0 * M_PI;
    while (dphi < -M_PI) dphi += 2.0 * M_PI;
    return dphi;
}
double deltaR_(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = deltaPhi_(phi1, phi2);
    return std::sqrt(deta * deta + dphi * dphi);
}
bool isLastCopy_(int flags)      { return (flags >> 13) & 1; } // Check from the NanoAOD version
bool isHardProcess_(int flags)   { return (flags >>  7) & 1; }
bool fromHardProcess_(int flags) { return (flags >>  8) & 1; }

// Find all gen particles whose mother is `parent_idx` (direct daughters only).
std::vector<int> directDaughters(int parent_idx,
                                 const ROOT::RVec<int>& mothers) {
    std::vector<int> daughters;
    for (int i = 0; i < (int)mothers.size(); ++i) {
        if (mothers[i] == parent_idx) daughters.push_back(i);
    }
    return daughters;
}

// Find all descendants of `ancestor` with the given |pdgId|, traversing
// the full mother chain. Optionally require a status flag (pass 0 to skip).
std::vector<int> findDescendantsByPdgId(
    int ancestor_idx,
    int target_abs_pdgid,
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& mothers,
    const ROOT::RVec<int>& statusFlags,
    bool require_last_copy = false)
{
    std::vector<int> result;
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) != target_abs_pdgid) continue;
        if (require_last_copy && !isLastCopy_(statusFlags[i])) continue;
        if (isDescendantOf(i, ancestor_idx, mothers)) result.push_back(i);
    }
    return result;
}

// Find first particle (last copy) with |pdgId| == target_abs_pdgid.
// Returns -1 if not found.
int findFirstLastCopy(int target_abs_pdgid,
                      const ROOT::RVec<int>& pdgId,
                      const ROOT::RVec<int>& statusFlags) {
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == target_abs_pdgid && isLastCopy_(statusFlags[i]))
            return i;
    }
    return -1;
}

// Find all particles (last copy) with |pdgId| == target_abs_pdgid.
std::vector<int> findAllLastCopies(int target_abs_pdgid,
                                   const ROOT::RVec<int>& pdgId,
                                   const ROOT::RVec<int>& statusFlags) {
    std::vector<int> result;
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (std::abs(pdgId[i]) == target_abs_pdgid && isLastCopy_(statusFlags[i]))
            result.push_back(i);
    }
    return result;
}

// Classify tau decay: returns 0=hadronic, 1=muonic, 2=electronic
// by checking descendants of the tau for muon (|pdgId|=13) or electron (|pdgId|=11).
int classifyTauDecay(int tau_idx,
                     const ROOT::RVec<int>& pdgId,
                     const ROOT::RVec<int>& mothers) {
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (!isDescendantOf(i, tau_idx, mothers)) continue;
        int aid = std::abs(pdgId[i]);
        if (aid == 13) return 1;  // muonic
        if (aid == 11) return 2;  // electronic
    }
    return 0;  // hadronic
}

} // anonymous namespace


// ═════════════════════════════════════════════════════════════════════════════
//  ZGenMatching — Drell-Yan Z boson gen matching
// ═════════════════════════════════════════════════════════════════════════════
//
// Finds the gen-level Z boson (pdgId=23, isLastCopy) and classifies its decay:
//
//   decayMode = 1 : Z -> e+ e-
//   decayMode = 2 : Z -> mu+ mu-
//   decayMode = 3 : Z -> tau tau (both hadronic)
//   decayMode = 4 : Z -> tau tau (one tau -> mu)
//   decayMode = 5 : Z -> tau tau (one tau -> e)
//
// Safe to call on TT / signal events — returns GenZResult with z_idx=-1 if
// no Z boson is found (DY->ll will have a Z; TT/signal will not).
//
// Usage in cutflow_TrigEff.py Define() call:
//   df.Define("gen_Z", "Ana::ZGenMatching(GenPart_pdgId, "
//                      "GenPart_genPartIdxMother, GenPart_statusFlags)")
//
GenZResult ZGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& mothers,
    const ROOT::RVec<int>& statusFlags)
{
    GenZResult result;

    // Find Z boson (last copy)
    int z_idx = findFirstLastCopy(23, pdgId, statusFlags);
    if (z_idx < 0) return result;
    result.z_idx = z_idx;

    // Collect direct daughters of Z
    auto daughters = directDaughters(z_idx, mothers);

    // Classify by |pdgId| of daughters
    int n_e   = 0, n_mu = 0, n_tau = 0;
    int lep1  = -1, lep2 = -1;
    for (int d : daughters) {
        int aid = std::abs(pdgId[d]);
        if      (aid == 11) { ++n_e;   if (lep1 < 0) lep1 = d; else lep2 = d; }
        else if (aid == 13) { ++n_mu;  if (lep1 < 0) lep1 = d; else lep2 = d; }
        else if (aid == 15) { ++n_tau; if (lep1 < 0) lep1 = d; else lep2 = d; }
    }

    result.lep1 = lep1;
    result.lep2 = lep2;

    if (n_e == 2) {
        result.decayMode = 1;  // Z -> ee
    } else if (n_mu == 2) {
        result.decayMode = 2;  // Z -> mumu
    } else if (n_tau == 2 && lep1 >= 0 && lep2 >= 0) {
        // Sub-classify tau decays
        int decay1 = classifyTauDecay(lep1, pdgId, mothers);
        int decay2 = classifyTauDecay(lep2, pdgId, mothers);
        int n_lep_tau = (decay1 > 0 ? 1 : 0) + (decay2 > 0 ? 1 : 0);
        bool has_mu = (decay1 == 1 || decay2 == 1);
        bool has_e  = (decay1 == 2 || decay2 == 2);

        if (n_lep_tau == 0)       result.decayMode = 3;  // Z -> tautau (hh)
        else if (has_mu)          result.decayMode = 4;  // Z -> tautau (muh)
        else if (has_e)           result.decayMode = 5;  // Z -> tautau (eh)
    }
    // else: Z -> bb or other hadronic; decayMode stays 0

    return result;
}


// ═════════════════════════════════════════════════════════════════════════════
//  TopGenMatching — TTbar top quark gen matching
// ═════════════════════════════════════════════════════════════════════════════
//
// Finds gen top (pdgId=6) and antitop (pdgId=-6) with isLastCopy, then
// classifies the overall ttbar decay topology:
//
//   decayMode = 1 : fully hadronic    (both W -> qq)
//   decayMode = 2 : semi-leptonic     (one W -> lnu, one W -> qq)
//   decayMode = 3 : fully leptonic    (both W -> lnu)
//
// Leptonic includes W -> tau nu (tau counts as leptonic for topology).
//
// Also identifies the b quarks from t -> W b and tbar -> W bbar.
//
// Safe to call on DY / signal events — returns GenTopResult with top_idx=-1.
//
// Usage:
//   df.Define("gen_top", "Ana::TopGenMatching(GenPart_pdgId, "
//             "GenPart_genPartIdxMother, GenPart_statusFlags)")
//
GenTopResult TopGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& mothers,
    const ROOT::RVec<int>& statusFlags)
{
    GenTopResult result;

    // Find top and antitop (last copies)
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (!isLastCopy_(statusFlags[i])) continue;
        if      (pdgId[i] ==  6 && result.top_idx     < 0) result.top_idx     = i;
        else if (pdgId[i] == -6 && result.antitop_idx < 0) result.antitop_idx = i;
    }
    if (result.top_idx < 0 || result.antitop_idx < 0) return result;

    // b quarks: last-copy b (|pdgId|=5) that is a descendant of top/antitop
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (!isLastCopy_(statusFlags[i])) continue;
        if (std::abs(pdgId[i]) != 5) continue;
        if      (result.b_from_top     < 0 && isDescendantOf(i, result.top_idx,     mothers))
            result.b_from_top = i;
        else if (result.b_from_antitop < 0 && isDescendantOf(i, result.antitop_idx, mothers))
            result.b_from_antitop = i;
    }

    // W bosons from top/antitop (last copy, |pdgId|=24)
    int w_from_top     = -1;
    int w_from_antitop = -1;
    for (int i = 0; i < (int)pdgId.size(); ++i) {
        if (!isLastCopy_(statusFlags[i])) continue;
        if (std::abs(pdgId[i]) != 24) continue;
        if      (w_from_top     < 0 && isDescendantOf(i, result.top_idx,     mothers))
            w_from_top = i;
        else if (w_from_antitop < 0 && isDescendantOf(i, result.antitop_idx, mothers))
            w_from_antitop = i;
    }

    // Classify each W: leptonic (e,mu,tau) or hadronic (u,d,s,c)
    auto wIsLeptonic = [&](int w_idx) -> bool {
        if (w_idx < 0) return false;
        auto daus = directDaughters(w_idx, mothers);
        for (int d : daus) {
            int aid = std::abs(pdgId[d]);
            if (aid == 11 || aid == 13 || aid == 15) return true;
        }
        return false;
    };

    bool top_lep     = wIsLeptonic(w_from_top);
    bool antitop_lep = wIsLeptonic(w_from_antitop);

    if (!top_lep && !antitop_lep)      result.decayMode = 1;  // allhad
    else if (top_lep ^ antitop_lep)    result.decayMode = 2;  // semilep
    else                               result.decayMode = 3;  // dilep

    return result;
}


// ═════════════════════════════════════════════════════════════════════════════
//  HHGenMatching — HH -> bb tautau signal gen matching
// ═════════════════════════════════════════════════════════════════════════════
//
// Finds the two gen Higgs bosons and classifies the H->tautau decay channel:
//
//   decayType = 1 : tau_h tau_h  (both taus hadronic)
//   decayType = 2 : tau_mu tau_h (one tau -> mu + neutrinos)
//   decayType = 3 : tau_e  tau_h (one tau -> e + neutrinos)
//
// Mirrors the logic of Ana::DecayGenMatching() in HHbbtautauAnaElements.C
// but returns the new GenHHResult struct and is a standalone implementation.
//
// dr_threshold: max deltaR for gen-reco matching (from config/regions.yaml).
//   Pass as: f"{CUTS_REG['gen_matching']['dr_threshold']}f"
//
// Safe to call on DY / TT events — returns GenHHResult with Htob=-1.
//
// Usage:
//   dr = CUTS_REG['gen_matching']['dr_threshold']
//   df.Define("gen_HH", f"Ana::HHGenMatching(GenPart_pdgId, "
//             "GenPart_genPartIdxMother, GenPart_statusFlags, {dr}f)")
//
GenHHResult HHGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& mothers,
    const ROOT::RVec<int>& statusFlags,
    float /*dr_threshold*/ = 0.05f)   // reserved for future reco-matching use
{
    GenHHResult result;

    // Find all gen Higgs bosons (last copies)
    auto higgs_indices = findAllLastCopies(25, pdgId, statusFlags);
    if (higgs_indices.size() < 2) return result;

    // For each Higgs, classify as H->bb or H->tautau by examining descendants
    int h_bb  = -1;
    int h_tau = -1;

    for (int h : higgs_indices) {
        bool has_b   = false;
        bool has_tau = false;
        auto daus = directDaughters(h, mothers);
        for (int d : daus) {
            int aid = std::abs(pdgId[d]);
            if (aid == 5)  has_b   = true;
            if (aid == 15) has_tau = true;
        }
        if (has_b   && h_bb  < 0) h_bb  = h;
        if (has_tau && h_tau < 0) h_tau = h;
    }

    if (h_bb < 0 || h_tau < 0) return result;  // not a bb+tautau event
    result.Htob   = h_bb;
    result.Htotau = h_tau;

    // Extract b quarks from H->bb daughters
    {
        auto daus = directDaughters(h_bb, mothers);
        for (int d : daus) {
            if (std::abs(pdgId[d]) == 5) {
                if (result.b1 < 0) result.b1 = d;
                else               result.b2 = d;
            }
        }
    }

    // Extract taus from H->tautau daughters
    {
        auto daus = directDaughters(h_tau, mothers);
        for (int d : daus) {
            if (std::abs(pdgId[d]) == 15) {
                if (result.tau1 < 0) result.tau1 = d;
                else                 result.tau2 = d;
            }
        }
    }

    if (result.tau1 < 0 || result.tau2 < 0) return result;

    // Classify tau decay channels
    int decay1 = classifyTauDecay(result.tau1, pdgId, mothers);
    int decay2 = classifyTauDecay(result.tau2, pdgId, mothers);

    // Need exactly one hadronic tau (this analysis targets semi-hadronic modes)
    if (decay1 == 0 && decay2 == 0) {
        result.decayType = 1;  // tau_h tau_h
    } else if (decay1 == 1 || decay2 == 1) {
        result.decayType = 2;  // tau_mu tau_h
        // Record the muon index
        int mu_tau = (decay1 == 1) ? result.tau1 : result.tau2;
        auto mu_descs = findDescendantsByPdgId(mu_tau, 13, pdgId, mothers, statusFlags, true);
        if (!mu_descs.empty()) result.mu = mu_descs[0];
    } else if (decay1 == 2 || decay2 == 2) {
        result.decayType = 3;  // tau_e tau_h
        // Record the electron index
        int e_tau = (decay1 == 2) ? result.tau1 : result.tau2;
        auto e_descs = findDescendantsByPdgId(e_tau, 11, pdgId, mothers, statusFlags, true);
        if (!e_descs.empty()) result.e = e_descs[0];
    }
    // else: e.g. tau_mu tau_mu or other modes — decayType stays 0

    return result;
}


// ═════════════════════════════════════════════════════════════════════════════
//  GlobalDecayMode — unified decay classification for all MC samples
// ═════════════════════════════════════════════════════════════════════════════
//
// Returns a single integer identifying the full decay chain (process + channel).
// Calls ZGenMatching, TopGenMatching, HHGenMatching internally and maps
// to a global numbering scheme:
//
//    0  = unknown / not classified
//    1  = DY  Z→ee
//    2  = DY  Z→μμ
//    3  = DY  Z→τhτh
//    4  = DY  Z→τμτh
//    5  = DY  Z→τeτh
//   10  = TT  fully hadronic
//   11  = TT  semi-leptonic
//   12  = TT  fully leptonic
//   20  = HH  →bb τhτh
//   21  = HH  →bb τμτh
//   22  = HH  →bb τeτh
//
// Gap-based numbering leaves room for future processes (QCD=30, W+jets=40, ...).
//
// The Python-side LUT (DECAY_MODES dict in cutflow_TrigEff.py) maps these
// integers to labels, colors, and process groups.
//
// Usage:
//   df.Define("decayMode", "Ana::GlobalDecayMode(GenPart_pdgId, "
//             "GenPart_genPartIdxMother, GenPart_statusFlags)")
//
int GlobalDecayMode(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<int>& mothers,
    const ROOT::RVec<int>& statusFlags)
{
    // DY: Z boson decay modes 1-5
    auto z = ZGenMatching(pdgId, mothers, statusFlags);
    if (z.found() && z.decayMode > 0) return z.decayMode;

    // TT: top decay topologies → 10-12
    auto top = TopGenMatching(pdgId, mothers, statusFlags);
    if (top.found() && top.decayMode > 0) return 9 + top.decayMode;

    // Signal: HH→bbττ decay types → 20-22
    auto hh = HHGenMatching(pdgId, mothers, statusFlags);
    if (hh.found() && hh.decayType > 0) return 19 + hh.decayType;

    return 0;
}

// ═════════════════════════════════════════════════════════════════════════════
//  NanoAOD type overloads
// ═════════════════════════════════════════════════════════════════════════════
//
// NanoAOD stores GenPart_genPartIdxMother as Short_t and GenPart_statusFlags
// as UShort_t.  ROOT's JIT can silently convert these to RVec<int> in the
// Define() call, but with ImplicitMT the temporary converted vectors cause
// memory corruption ("double free").  These overloads accept the native types
// and convert explicitly.

namespace { // helper to avoid code duplication
inline void _cvt(const ROOT::RVec<short>& s, const ROOT::RVec<unsigned short>& u,
                  ROOT::RVec<int>& m, ROOT::RVec<int>& f) {
    m.assign(s.begin(), s.end());
    f.assign(u.begin(), u.end());
}}

GenZResult ZGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<short>& mothers_s,
    const ROOT::RVec<unsigned short>& flags_us)
{
    ROOT::RVec<int> m, f; _cvt(mothers_s, flags_us, m, f);
    return ZGenMatching(pdgId, m, f);
}

GenTopResult TopGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<short>& mothers_s,
    const ROOT::RVec<unsigned short>& flags_us)
{
    ROOT::RVec<int> m, f; _cvt(mothers_s, flags_us, m, f);
    return TopGenMatching(pdgId, m, f);
}

GenHHResult HHGenMatching(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<short>& mothers_s,
    const ROOT::RVec<unsigned short>& flags_us,
    float dr_threshold = 0.05f)
{
    ROOT::RVec<int> m, f; _cvt(mothers_s, flags_us, m, f);
    return HHGenMatching(pdgId, m, f, dr_threshold);
}

int GlobalDecayMode(
    const ROOT::RVec<int>& pdgId,
    const ROOT::RVec<short>& mothers_s,
    const ROOT::RVec<unsigned short>& flags_us)
{
    ROOT::RVec<int> m, f; _cvt(mothers_s, flags_us, m, f);
    return GlobalDecayMode(pdgId, m, f);
}

// ═════════════════════════════════════════════════════════════════════════════
//  genMatchToRecoJet — match a gen particle to closest reco jet within dR
// ═════════════════════════════════════════════════════════════════════════════

int genMatchToRecoJet(int genIdx,
                      const ROOT::RVec<float>& genPt,
                      const ROOT::RVec<float>& genEta,
                      const ROOT::RVec<float>& genPhi,
                      const ROOT::RVec<float>& genMass,
                      const ROOT::RVec<float>& recoPt,
                      const ROOT::RVec<float>& recoEta,
                      const ROOT::RVec<float>& recoPhi,
                      const ROOT::RVec<float>& recoMass,
                      float drThreshold = 0.4f) {
    if (genIdx < 0 || genIdx >= (int)genPt.size()) return -1;
    TLorentzVector gen;
    gen.SetPtEtaPhiM(genPt[genIdx], genEta[genIdx], genPhi[genIdx], genMass[genIdx]);
    double bestDR = drThreshold;
    int best = -1;
    for (int i = 0; i < (int)recoPt.size(); i++) {
        TLorentzVector reco;
        reco.SetPtEtaPhiM(recoPt[i], recoEta[i], recoPhi[i], recoMass[i]);
        double dR = gen.DeltaR(reco);
        if (dR < bestDR) {
            bestDR = dR;
            best = i;
        }
    }
    return best;
}

} // namespace Ana
