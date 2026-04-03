"""RDataFrame .Define() chains for kinematic variables.

Split into categorized functions so load_and_run() can auto-detect
which groups are needed based on the user's PLOT_VARS.
"""

from analysis.constants import (
    AK4 as _AK4, AK8 as _AK8, AK8_SGP as _SGP, AK4_UPART,
    AK8_PT_MIN as _AK8_PT_MIN, AK8_SGP_PROB_NAMES as _GM_PROB_NAMES,
    UPART_PROB_NAMES,
)


def _vs_all_at(idx_var, prob_name):
    """Build XvsAll expression at a specific jet index."""
    num = f"{_SGP}_prob_{prob_name}[{idx_var}]"
    denom = " + ".join(f"{_SGP}_prob_{p}[{idx_var}]" for p in _GM_PROB_NAMES)
    return f"(float)({num} / ({denom}))"


# ═══════════════════════════════════════════════════════════════════════════════
#  AK4 jet basics
# ═══════════════════════════════════════════════════════════════════════════════

def define_ak4_jets(df):
    """AK4 jet pT/eta/mass, nJets, HT, MHT, nLeptons, nMuons, nElectrons, centrality."""
    df = (df
        .Define("ak4_pt0",  f"{_AK4}_pt[0]")
        .Define("ak4_pt1",  f"{_AK4}_pt[1]")
        .Define("ak4_pt2",  f"{_AK4}_pt[2]")
        .Define("ak4_pt3",  f"{_AK4}_pt[3]")
        .Define("ak4_eta0", f"{_AK4}_eta[0]")
        .Define("ak4_eta1", f"{_AK4}_eta[1]")
        .Define("ak4_eta2", f"{_AK4}_eta[2]")
        .Define("ak4_eta3", f"{_AK4}_eta[3]")
        .Define("ak4_mass0", f"{_AK4}_mass[0]")
        .Define("ak4_mass1", f"{_AK4}_mass[1]")
        .Define("HT",         f"Sum({_AK4}_pt)")
        .Define("nJets",      f"n{_AK4}")
        .Define("nLeptons",   "nScoutingMuonVtx + nScoutingElectron")
        .Define("nMuons",     "nScoutingMuonVtx")
        .Define("nElectrons", "nScoutingElectron")
        .Define("MHT",
                f"(float)sqrt("
                f"pow(Sum({_AK4}_pt*cos({_AK4}_phi)),2) + "
                f"pow(Sum({_AK4}_pt*sin({_AK4}_phi)),2))")
        .Define("centrality",
                f"(float)(Sum({_AK4}_pt) / "
                f"Sum({_AK4}_pt * cosh({_AK4}_eta)))")
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  AK8 fat jet kinematics + tagger scores + candidates
# ═══════════════════════════════════════════════════════════════════════════════

def define_ak8_jets(df):
    """AK8 per-jet kinematics, ScoutGlobalParT tagger scores, VsQCD, candidates, dR."""
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
    # Tagger-based candidates (best XvsAll score per hypothesis)
    _pt_cut = f"{_AK8}_pt > {_AK8_PT_MIN}f"
    _all_probs = " + ".join(f"{_SGP}_prob_{p}" for p in _GM_PROB_NAMES)
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
    df = df.Define("ak8_dR_Hbb_Htt",
                   "(ak8_Hbb_idx >= 0 && ak8_Htt_idx >= 0) ? "
                   "ROOT::VecOps::DeltaR(ak8_Hbb_eta, ak8_Htt_eta, ak8_Hbb_phi, ak8_Htt_phi) : -1.f")
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  B-tagging (UParT discriminators + b-jet kinematics)
# ═══════════════════════════════════════════════════════════════════════════════

def define_b_candidates(df):
    """UParT raw scores, BvsAll/TauVsAll discriminants, b-jet selection, mbb, dR_bb."""
    _U = f"{_AK4}_scoutUParT"
    _denom = (f"{_U}_probb + {_U}_probc + {_U}_probg"
              f" + {_U}_probuds + {_U}_problepb"
              f" + {_U}_probtaum + {_U}_probtaup")
    df = (df
        .Define("upart_b_raw",    f"{_U}_probb")
        .Define("upart_taup_raw", f"{_U}_probtaup")
        .Define("upart_taum_raw", f"{_U}_probtaum")
        .Define("upart_tau_raw",  f"{_U}_probtaup + {_U}_probtaum")
        .Define("upart_c_raw",    f"{_U}_probc")
        .Define("upart_g_raw",    f"{_U}_probg")
        .Define("upart_uds_raw",  f"{_U}_probuds")
        .Define("upart_lepb_raw", f"{_U}_problepb")
        .Define("ak4_BvsAll",    f"{_U}_probb / ({_denom})")
        .Define("ak4_TaupVsAll", f"{_U}_probtaup / ({_denom})")
        .Define("ak4_TaumVsAll", f"{_U}_probtaum / ({_denom})")
        .Define("ak4_TauVsAll",  f"({_U}_probtaup + {_U}_probtaum) / ({_denom})")
        .Define("bsort_idx",
                "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(ak4_BvsAll))")
        # B COI: jets with highest/2nd-highest BvsAll score
        .Define("b_coi0_idx", "nJets >= 1 ? (int)bsort_idx[0] : -1")
        .Define("b_coi1_idx", "nJets >= 2 ? (int)bsort_idx[1] : -1")
        .Define("b_coi0_pt",    f"b_coi0_idx >= 0 ? {_AK4}_pt[b_coi0_idx] : -1.f")
        .Define("b_coi0_eta",   f"b_coi0_idx >= 0 ? {_AK4}_eta[b_coi0_idx] : -99.f")
        .Define("b_coi0_phi",   f"b_coi0_idx >= 0 ? {_AK4}_phi[b_coi0_idx] : -99.f")
        .Define("b_coi0_mass",  f"b_coi0_idx >= 0 ? {_AK4}_mass[b_coi0_idx] : -1.f")
        .Define("b_coi0_score", "b_coi0_idx >= 0 ? ak4_BvsAll[b_coi0_idx] : -1.f")
        .Define("b_coi0_raw",   "b_coi0_idx >= 0 ? upart_b_raw[b_coi0_idx] : -1.f")
        .Define("b_coi1_pt",    f"b_coi1_idx >= 0 ? {_AK4}_pt[b_coi1_idx] : -1.f")
        .Define("b_coi1_eta",   f"b_coi1_idx >= 0 ? {_AK4}_eta[b_coi1_idx] : -99.f")
        .Define("b_coi1_phi",   f"b_coi1_idx >= 0 ? {_AK4}_phi[b_coi1_idx] : -99.f")
        .Define("b_coi1_mass",  f"b_coi1_idx >= 0 ? {_AK4}_mass[b_coi1_idx] : -1.f")
        .Define("b_coi1_score", "b_coi1_idx >= 0 ? ak4_BvsAll[b_coi1_idx] : -1.f")
        .Define("b_coi1_raw",   "b_coi1_idx >= 0 ? upart_b_raw[b_coi1_idx] : -1.f")
        .Define("mbb_coi",
                "b_coi0_idx >= 0 && b_coi1_idx >= 0 ? "
                "(float)(ROOT::Math::PtEtaPhiMVector(b_coi0_pt,b_coi0_eta,b_coi0_phi,b_coi0_mass)"
                " + ROOT::Math::PtEtaPhiMVector(b_coi1_pt,b_coi1_eta,b_coi1_phi,b_coi1_mass)).M() : -1.f")
        .Define("dR_bb_coi",
                "b_coi0_idx >= 0 && b_coi1_idx >= 0 ? "
                "ROOT::VecOps::DeltaR(b_coi0_eta,b_coi1_eta,b_coi0_phi,b_coi1_phi) : -1.f")
        .Define("ptbb_coi",
                "b_coi0_idx >= 0 && b_coi1_idx >= 0 ? "
                "(float)(ROOT::Math::PtEtaPhiMVector(b_coi0_pt,b_coi0_eta,b_coi0_phi,b_coi0_mass)"
                " + ROOT::Math::PtEtaPhiMVector(b_coi1_pt,b_coi1_eta,b_coi1_phi,b_coi1_mass)).Pt() : -1.f")
        # Backwards-compat aliases
        .Define("b0_idx", "b_coi0_idx").Define("b1_idx", "b_coi1_idx")
        .Define("b0_pt", "b_coi0_pt").Define("b0_eta", "b_coi0_eta")
        .Define("b0_phi", "b_coi0_phi").Define("b0_mass", "b_coi0_mass")
        .Define("b0_score", "b_coi0_score").Define("b0_raw", "b_coi0_raw")
        .Define("b1_pt", "b_coi1_pt").Define("b1_eta", "b_coi1_eta")
        .Define("b1_phi", "b_coi1_phi").Define("b1_mass", "b_coi1_mass")
        .Define("b1_score", "b_coi1_score").Define("b1_raw", "b_coi1_raw")
        .Define("mbb", "mbb_coi").Define("dR_bb", "dR_bb_coi").Define("ptbb", "ptbb_coi")
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Tau COI: jets with highest/2nd-highest TauVsAll score
# ═══════════════════════════════════════════════════════════════════════════════

def define_tau_candidates(df):
    """Tau COI: jets sorted by TauVsAll score (highest = tau_coi0)."""
    df = (df
        .Define("tausort_idx",
                "ROOT::VecOps::Reverse(ROOT::VecOps::Argsort(ak4_TauVsAll))")
        .Define("tau_coi0_idx", "nJets >= 1 ? (int)tausort_idx[0] : -1")
        .Define("tau_coi1_idx", "nJets >= 2 ? (int)tausort_idx[1] : -1")
        .Define("tau_coi0_pt",    f"tau_coi0_idx >= 0 ? {_AK4}_pt[tau_coi0_idx] : -1.f")
        .Define("tau_coi0_eta",   f"tau_coi0_idx >= 0 ? {_AK4}_eta[tau_coi0_idx] : -99.f")
        .Define("tau_coi0_phi",   f"tau_coi0_idx >= 0 ? {_AK4}_phi[tau_coi0_idx] : -99.f")
        .Define("tau_coi0_mass",  f"tau_coi0_idx >= 0 ? {_AK4}_mass[tau_coi0_idx] : -1.f")
        .Define("tau_coi0_score", "tau_coi0_idx >= 0 ? ak4_TauVsAll[tau_coi0_idx] : -1.f")
        .Define("tau_coi0_taup",  "tau_coi0_idx >= 0 ? ak4_TaupVsAll[tau_coi0_idx] : -1.f")
        .Define("tau_coi0_taum",  "tau_coi0_idx >= 0 ? ak4_TaumVsAll[tau_coi0_idx] : -1.f")
        .Define("tau_coi1_pt",    f"tau_coi1_idx >= 0 ? {_AK4}_pt[tau_coi1_idx] : -1.f")
        .Define("tau_coi1_eta",   f"tau_coi1_idx >= 0 ? {_AK4}_eta[tau_coi1_idx] : -99.f")
        .Define("tau_coi1_phi",   f"tau_coi1_idx >= 0 ? {_AK4}_phi[tau_coi1_idx] : -99.f")
        .Define("tau_coi1_mass",  f"tau_coi1_idx >= 0 ? {_AK4}_mass[tau_coi1_idx] : -1.f")
        .Define("tau_coi1_score", "tau_coi1_idx >= 0 ? ak4_TauVsAll[tau_coi1_idx] : -1.f")
        .Define("tau_coi1_taup",  "tau_coi1_idx >= 0 ? ak4_TaupVsAll[tau_coi1_idx] : -1.f")
        .Define("tau_coi1_taum",  "tau_coi1_idx >= 0 ? ak4_TaumVsAll[tau_coi1_idx] : -1.f")
        .Define("mtautau_coi",
                "tau_coi0_idx >= 0 && tau_coi1_idx >= 0 ? "
                "(float)(ROOT::Math::PtEtaPhiMVector(tau_coi0_pt,tau_coi0_eta,tau_coi0_phi,tau_coi0_mass)"
                " + ROOT::Math::PtEtaPhiMVector(tau_coi1_pt,tau_coi1_eta,tau_coi1_phi,tau_coi1_mass)).M() : -1.f")
        .Define("dR_tautau_coi",
                "tau_coi0_idx >= 0 && tau_coi1_idx >= 0 ? "
                "ROOT::VecOps::DeltaR(tau_coi0_eta,tau_coi1_eta,tau_coi0_phi,tau_coi1_phi) : -1.f")
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Dijet / 4-jet / Higgs candidate masses
# ═══════════════════════════════════════════════════════════════════════════════

def define_dijet(df):
    """mjj_01, dR_01, dEta_01, m4j, Higgs candidate dijet masses."""
    _p4 = "ROOT::Math::PtEtaPhiMVector"
    df = (df
        .Define("mjj_01",
                f"(float)({_p4}("
                f"{_AK4}_pt[0],{_AK4}_eta[0],{_AK4}_phi[0],{_AK4}_mass[0])"
                f" + {_p4}("
                f"{_AK4}_pt[1],{_AK4}_eta[1],{_AK4}_phi[1],{_AK4}_mass[1])).M()")
        .Define("dR_01",
                f"ROOT::VecOps::DeltaR("
                f"{_AK4}_eta[0],{_AK4}_eta[1],{_AK4}_phi[0],{_AK4}_phi[1])")
        .Define("dEta_01",
                f"(float)abs({_AK4}_eta[0] - {_AK4}_eta[1])")
        .Define("m4j",
                f"(float)({_p4}("
                f"{_AK4}_pt[0],{_AK4}_eta[0],{_AK4}_phi[0],{_AK4}_mass[0])"
                f" + {_p4}("
                f"{_AK4}_pt[1],{_AK4}_eta[1],{_AK4}_phi[1],{_AK4}_mass[1])"
                f" + {_p4}("
                f"{_AK4}_pt[2],{_AK4}_eta[2],{_AK4}_phi[2],{_AK4}_mass[2])"
                f" + {_p4}("
                f"{_AK4}_pt[3],{_AK4}_eta[3],{_AK4}_phi[3],{_AK4}_mass[3])).M()")
    )
    # H→bb candidate: best dijet pair with m_jj ∈ [100, 150] GeV (closest to 125)
    df = (df
        .Define("hbb_pair",
                f"Ana::findDijetInWindow("
                f"{_AK4}_pt, {_AK4}_eta, {_AK4}_phi, {_AK4}_mass,"
                f"100.f, 150.f, 125.f)")
        .Define("has_hbb",      "hbb_pair.i1 >= 0")
        .Define("mbb_cand",     "hbb_pair.mass")
        .Define("htautau_pair",
                f"Ana::findDijetInWindow("
                f"{_AK4}_pt, {_AK4}_eta, {_AK4}_phi, {_AK4}_mass,"
                f"40.f, 150.f, 80.f, {{hbb_pair.i1, hbb_pair.i2}})")
        .Define("has_htautau",  "htautau_pair.i1 >= 0")
        .Define("mtautau_cand", "htautau_pair.mass")
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Convenience wrapper (calls all sub-functions)
# ═══════════════════════════════════════════════════════════════════════════════

def define_kinematics(df):
    """Define all derived kinematic columns. Calls the categorized functions."""
    df = define_ak4_jets(df)
    df = define_ak8_jets(df)
    df = define_b_candidates(df)
    df = define_tau_candidates(df)
    df = define_dijet(df)
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Lepton selection (muon + electron for τμτh and τeτh channels)
# ═══════════════════════════════════════════════════════════════════════════════

def define_lepton_selection(df, obj_cfg):
    """Select best muon and electron using cuts from config/objects.yaml.

    Defines: mu_idx, mu_pt, mu_eta, mu_phi, mu_relIso, has_good_muon
             el_idx, el_pt, el_eta, el_phi, el_relIso, has_good_electron
             dR_mu_tau, dR_el_tau
    """
    _MU = "ScoutingMuonVtx"
    _EL = "ScoutingElectron"

    mu = obj_cfg.get("muon", {})
    el = obj_cfg.get("electron", {})

    # ── Muon selection ──
    df = df.Define("mu_idx",
        f"Ana::SelectMuon({_MU}_pt, {_MU}_eta, "
        f"{_MU}_trackIso, {_MU}_ecalIso, {_MU}_hcalIso, "
        f"{_MU}_trk_dxy, {_MU}_trk_dz, "
        f"{mu.get('pt_min', 20.0)}f, {mu.get('eta_max', 2.4)}f, "
        f"{mu.get('iso_max', 0.15)}f, "
        f"{mu.get('dxy_max', 0.045)}f, {mu.get('dz_max', 0.2)}f)")

    df = (df
        .Define("has_good_muon", "mu_idx >= 0")
        .Define("mu_pt",  f"mu_idx >= 0 ? {_MU}_pt[mu_idx] : -1.f")
        .Define("mu_eta", f"mu_idx >= 0 ? {_MU}_eta[mu_idx] : -99.f")
        .Define("mu_phi", f"mu_idx >= 0 ? {_MU}_phi[mu_idx] : -99.f")
        .Define("mu_relIso",
            f"mu_idx >= 0 && {_MU}_pt[mu_idx] > 0 ? "
            f"(float)(({_MU}_trackIso[mu_idx] + {_MU}_ecalIso[mu_idx] + {_MU}_hcalIso[mu_idx]) "
            f"/ {_MU}_pt[mu_idx]) : -1.f")
    )

    # ── Electron selection ──
    df = df.Define("el_idx",
        f"Ana::SelectElectron({_EL}_pt, {_EL}_eta, "
        f"{_EL}_trackIso, {_EL}_ecalIso, {_EL}_hcalIso, "
        f"{_EL}_bestTrack_d0, {_EL}_bestTrack_dz, "
        f"{el.get('pt_min', 20.0)}f, {el.get('eta_max', 2.5)}f, "
        f"{el.get('crack_lo', 1.4442)}f, {el.get('crack_hi', 1.566)}f, "
        f"{el.get('iso_max', 0.15)}f, "
        f"{el.get('d0_max', 0.045)}f, {el.get('dz_max', 0.2)}f)")

    df = (df
        .Define("has_good_electron", "el_idx >= 0")
        .Define("el_pt",  f"el_idx >= 0 ? {_EL}_pt[el_idx] : -1.f")
        .Define("el_eta", f"el_idx >= 0 ? {_EL}_eta[el_idx] : -99.f")
        .Define("el_phi", f"el_idx >= 0 ? {_EL}_phi[el_idx] : -99.f")
        .Define("el_relIso",
            f"el_idx >= 0 && {_EL}_pt[el_idx] > 0 ? "
            f"(float)(({_EL}_trackIso[el_idx] + {_EL}_ecalIso[el_idx] + {_EL}_hcalIso[el_idx]) "
            f"/ {_EL}_pt[el_idx]) : -1.f")
    )

    # ── ΔR(lepton, tau COI candidate) ──
    # tau_coi0 is the leading tau COI candidate (defined in define_tau_candidates)
    df = (df
        .Define("dR_mu_tau",
            "mu_idx >= 0 && tau_coi0_pt > 0 ? "
            "(float)sqrt(pow(mu_eta - tau_coi0_eta, 2) + "
            "pow(TVector2::Phi_mpi_pi(mu_phi - tau_coi0_phi), 2)) : -1.f")
        .Define("dR_el_tau",
            "el_idx >= 0 && tau_coi0_pt > 0 ? "
            "(float)sqrt(pow(el_eta - tau_coi0_eta, 2) + "
            "pow(TVector2::Phi_mpi_pi(el_phi - tau_coi0_phi), 2)) : -1.f")
    )

    return df


_DATA_AK4 = "ScoutingPFJetRecluster"


def define_kinematics_data(df):
    """Define kinematic columns for data (v1 branches, no UParT/COI)."""
    _D = _DATA_AK4
    df = (df
        .Define("ak4_pt0",  f"{_D}_pt[0]")
        .Define("ak4_pt1",  f"{_D}_pt[1]")
        .Define("ak4_pt2",  f"{_D}_pt[2]")
        .Define("ak4_pt3",  f"{_D}_pt[3]")
        .Define("ak4_eta0", f"{_D}_eta[0]")
        .Define("ak4_eta1", f"{_D}_eta[1]")
        .Define("ak4_eta2", f"{_D}_eta[2]")
        .Define("ak4_eta3", f"{_D}_eta[3]")
        .Define("ak4_mass0", f"{_D}_mass[0]")
        .Define("ak4_mass1", f"{_D}_mass[1]")
        .Define("HT",         f"Sum({_D}_pt)")
        .Define("nJets",      f"n{_D}")
        .Define("nLeptons",   "nScoutingMuonVtx + nScoutingElectron")
        .Define("nMuons",     "nScoutingMuonVtx")
        .Define("nElectrons", "nScoutingElectron")
        .Define("MHT",
                f"(float)sqrt("
                f"pow(Sum({_D}_pt*cos({_D}_phi)),2) + "
                f"pow(Sum({_D}_pt*sin({_D}_phi)),2))")
        .Define("centrality",
                f"(float)(Sum({_D}_pt) / "
                f"Sum({_D}_pt * cosh({_D}_eta)))")
        # Dijet variables from leading pair
        .Define("mjj_01",
                f"(float)sqrt(2*{_D}_pt[0]*{_D}_pt[1]*"
                f"(cosh({_D}_eta[0]-{_D}_eta[1])-cos({_D}_phi[0]-{_D}_phi[1])))")
        .Define("dR_01",
                f"(float)sqrt(pow({_D}_eta[0]-{_D}_eta[1],2)+"
                f"pow(TVector2::Phi_mpi_pi({_D}_phi[0]-{_D}_phi[1]),2))")
        .Define("dEta_01", f"(float)abs({_D}_eta[0]-{_D}_eta[1])")
    )
    return df


# ═══════════════════════════════════════════════════════════════════════════════
#  Gen-level columns (signal-only gen matching)
# ═══════════════════════════════════════════════════════════════════════════════

def define_gen_columns(mc, sig_samples):
    """Define gen-level columns (decayMode, gen HH matching, gen-match AK8)."""
    _GP = "GenPart"
    decay_mode_expr = f"Ana::GlobalDecayMode({_GP}_pdgId, {_GP}_genPartIdxMother, {_GP}_statusFlags)"
    gen_hh_expr = f"Ana::HHGenMatching({_GP}_pdgId, {_GP}_genPartIdxMother, {_GP}_statusFlags)"

    for name in mc:
        mc[name] = mc[name].Define("decayMode", decay_mode_expr)

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
        # Gen-match AK8 jets to gen H→bb and H→ττ
        mc[name] = (mc[name]
            .Define("ak8_dR_to_genHbb",
                    f"ROOT::VecOps::DeltaR("
                    f"{_AK8}_eta, ROOT::RVecF(n{_AK8}, (float)genHbb_p4.Eta()), "
                    f"{_AK8}_phi, ROOT::RVecF(n{_AK8}, (float)genHbb_p4.Phi()))")
            .Define("ak8_genHbb_match_idx",
                    f"n{_AK8} > 0 && genHbb_p4.Pt() > 0 ? "
                    f"(Min(ak8_dR_to_genHbb) < 0.8f ? (int)ArgMin(ak8_dR_to_genHbb) : -1) : -1")
            .Define("ak8_genHbb_match_dR",
                    "ak8_genHbb_match_idx >= 0 ? (float)ak8_dR_to_genHbb[ak8_genHbb_match_idx] : -1.f")
            .Define("ak8_genHbb_match_Xbb",
                    f"ak8_genHbb_match_idx >= 0 ? {_SGP}_prob_Xbb[ak8_genHbb_match_idx] : -1.f")
            .Define("ak8_genHbb_match_XbbVsAll",
                    f"ak8_genHbb_match_idx >= 0 ? {_vs_all_at('ak8_genHbb_match_idx', 'Xbb')} : -1.f")
            .Define("ak8_dR_to_genHtt",
                    f"ROOT::VecOps::DeltaR("
                    f"{_AK8}_eta, ROOT::RVecF(n{_AK8}, (float)genHtautau_p4.Eta()), "
                    f"{_AK8}_phi, ROOT::RVecF(n{_AK8}, (float)genHtautau_p4.Phi()))")
            .Define("ak8_genHtt_match_idx",
                    f"n{_AK8} > 0 && genHtautau_p4.Pt() > 0 ? "
                    f"(Min(ak8_dR_to_genHtt) < 0.8f ? (int)ArgMin(ak8_dR_to_genHtt) : -1) : -1")
            .Define("ak8_genHtt_match_dR",
                    "ak8_genHtt_match_idx >= 0 ? (float)ak8_dR_to_genHtt[ak8_genHtt_match_idx] : -1.f")
            .Define("ak8_genHtt_match_Xtt",
                    f"ak8_genHtt_match_idx >= 0 ? {_SGP}_prob_Xtauhtauh[ak8_genHtt_match_idx] : -1.f")
            .Define("ak8_genHtt_match_XttVsAll",
                    f"ak8_genHtt_match_idx >= 0 ? {_vs_all_at('ak8_genHtt_match_idx', 'Xtauhtauh')} : -1.f")
        )

    # Dummy gen columns for non-signal samples
    _gen_higgs_cols = [
        "gen_mHH", "gen_pt_Hbb", "gen_pt_Htautau",
        "gen_eta_Hbb", "gen_eta_Htautau",
        "gen_phi_Hbb", "gen_phi_Htautau",
        "gen_E_Hbb", "gen_E_Htautau",
        "gen_mass_Hbb", "gen_mass_Htautau",
        "gen_dR_H1H2", "gen_pt_HH", "gen_eta_HH", "gen_phi_HH", "gen_E_HH",
    ]
    _gen_match_cols_ak8 = ["ak8_genHbb_match_dR", "ak8_genHbb_match_Xbb",
                           "ak8_genHbb_match_XbbVsAll", "ak8_genHtt_match_dR",
                           "ak8_genHtt_match_Xtt", "ak8_genHtt_match_XttVsAll"]
    for name in mc:
        if name not in sig_samples:
            for col in _gen_higgs_cols + _gen_match_cols_ak8:
                mc[name] = mc[name].Define(col, "-1.f")


# ═══════════════════════════════════════════════════════════════════════════════
#  Gen-matched AK4 UParT scores (signal only)
# ═══════════════════════════════════════════════════════════════════════════════

def define_gen_matched_ak4(mc, sig_samples):
    """Define gen-matched AK4 jet columns (signal only, after define_btagging)."""
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
                        f"{_AK4}_eta, ROOT::RVecF(n{_AK4}, (float){_p4}.Eta()), "
                        f"{_AK4}_phi, ROOT::RVecF(n{_AK4}, (float){_p4}.Phi()))")
                .Define(_idx,
                        f"n{_AK4} > 0 && {_p4}.Pt() > 0 ? "
                        f"(int)(Min({_dR_vec}) < 0.4f ? ArgMin({_dR_vec}) : -1) : -1")
                .Define(f"ak4_gen{_gp_label}_dR",
                        f"{_idx} >= 0 ? {_dR_vec}[{_idx}] : -1.f")
                .Define(f"ak4_gen{_gp_label}_{_score_name}",
                        f"{_idx} >= 0 ? ak4_{_score_name}[{_idx}] : -1.f")
            )

    # Additional tau scores using existing gen-matched indices
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

    # Charge misidentification columns
    for name in sig_samples:
        def _neg(tl): return (f"gen_{tl}_charge < 0.f && gen_{tl}_charge > -90.f"
                              f" && ak4_gen{tl}_dR > -0.5f")
        def _pos(tl): return f"gen_{tl}_charge > 0.f && ak4_gen{tl}_dR > -0.5f"
        mc[name] = (mc[name]
            .Define("gen_taum_TaumVsAll",
                    f"({_neg('tau1')}) ? ak4_gentau1_TaumVsAll :"
                    f" ({_neg('tau2')}) ? ak4_gentau2_TaumVsAll : -1.f")
            .Define("gen_taum_TaupVsAll",
                    f"({_neg('tau1')}) ? ak4_gentau1_TaupVsAll :"
                    f" ({_neg('tau2')}) ? ak4_gentau2_TaupVsAll : -1.f")
            .Define("gen_taum_TauVsAll",
                    f"({_neg('tau1')}) ? ak4_gentau1_TauVsAll :"
                    f" ({_neg('tau2')}) ? ak4_gentau2_TauVsAll : -1.f")
            .Define("gen_taup_TaupVsAll",
                    f"({_pos('tau1')}) ? ak4_gentau1_TaupVsAll :"
                    f" ({_pos('tau2')}) ? ak4_gentau2_TaupVsAll : -1.f")
            .Define("gen_taup_TaumVsAll",
                    f"({_pos('tau1')}) ? ak4_gentau1_TaumVsAll :"
                    f" ({_pos('tau2')}) ? ak4_gentau2_TaumVsAll : -1.f")
            .Define("gen_taup_TauVsAll",
                    f"({_pos('tau1')}) ? ak4_gentau1_TauVsAll :"
                    f" ({_pos('tau2')}) ? ak4_gentau2_TauVsAll : -1.f")
        )

    # Charge-ID prediction: 1 = tagger correct, 0 = tagger wrong
    for name in sig_samples:
        mc[name] = (mc[name]
            .Define("gen_taum_pred_correct",
                    "gen_taum_TaumVsAll > -0.5f ? (gen_taum_TaumVsAll > gen_taum_TaupVsAll ? 1.f : 0.f) : -1.f")
            .Define("gen_taup_pred_correct",
                    "gen_taup_TaupVsAll > -0.5f ? (gen_taup_TaupVsAll > gen_taup_TaumVsAll ? 1.f : 0.f) : -1.f")
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
    _charge_misid_cols = [
        "gen_taum_TaumVsAll", "gen_taum_TaupVsAll", "gen_taum_TauVsAll",
        "gen_taup_TaupVsAll", "gen_taup_TaumVsAll", "gen_taup_TauVsAll",
        "gen_taum_pred_correct", "gen_taup_pred_correct",
    ]
    _all_ak4_gen_cols = _ak4_gen_cols + _ak4_gen_tau_extra + _ak4_gen_raw_cols + _charge_misid_cols
    for name in mc:
        if name not in sig_samples:
            for col in _all_ak4_gen_cols:
                mc[name] = mc[name].Define(col, "-1.f")


# ═══════════════════════════════════════════════════════════════════════════════
#  COI matching study (signal only — compares COI selection to gen truth)
# ═══════════════════════════════════════════════════════════════════════════════

def define_coi_matching(mc, sig_samples):
    """Define columns comparing COI jets to gen-matched jets (signal only)."""
    _coi_cols = []

    for name in sig_samples:
        mc[name] = (mc[name]
            # B COI matching: did the highest-BvsAll jet match a gen b-quark?
            .Define("b_coi0_is_genb",
                    "b_coi0_idx >= 0 ? "
                    "((ak4_genb1_idx >= 0 && b_coi0_idx == ak4_genb1_idx) || "
                    " (ak4_genb2_idx >= 0 && b_coi0_idx == ak4_genb2_idx) ? 1.f : 0.f) : -1.f")
            .Define("b_coi1_is_genb",
                    "b_coi1_idx >= 0 ? "
                    "((ak4_genb1_idx >= 0 && b_coi1_idx == ak4_genb1_idx) || "
                    " (ak4_genb2_idx >= 0 && b_coi1_idx == ak4_genb2_idx) ? 1.f : 0.f) : -1.f")
            .Define("b_both_coi_are_genb",
                    "b_coi0_is_genb > 0.5f && b_coi1_is_genb > 0.5f ? 1.f : "
                    "(b_coi0_idx >= 0 && b_coi1_idx >= 0 ? 0.f : -1.f)")

            # Tau COI matching: did the highest-TauVsAll jet match a gen tau?
            .Define("tau_coi0_is_gentau",
                    "tau_coi0_idx >= 0 ? "
                    "((ak4_gentau1_idx >= 0 && tau_coi0_idx == ak4_gentau1_idx) || "
                    " (ak4_gentau2_idx >= 0 && tau_coi0_idx == ak4_gentau2_idx) ? 1.f : 0.f) : -1.f")
            .Define("tau_coi1_is_gentau",
                    "tau_coi1_idx >= 0 ? "
                    "((ak4_gentau1_idx >= 0 && tau_coi1_idx == ak4_gentau1_idx) || "
                    " (ak4_gentau2_idx >= 0 && tau_coi1_idx == ak4_gentau2_idx) ? 1.f : 0.f) : -1.f")
            .Define("tau_both_coi_are_gentau",
                    "tau_coi0_is_gentau > 0.5f && tau_coi1_is_gentau > 0.5f ? 1.f : "
                    "(tau_coi0_idx >= 0 && tau_coi1_idx >= 0 ? 0.f : -1.f)")

            # Overlap: is the same jet both best-b and best-tau?
            .Define("b_tau_overlap_0",
                    "b_coi0_idx >= 0 && tau_coi0_idx >= 0 ? "
                    "(b_coi0_idx == tau_coi0_idx ? 1.f : 0.f) : -1.f")
            .Define("any_b_tau_overlap",
                    "b_coi0_idx >= 0 && tau_coi0_idx >= 0 ? "
                    "((b_coi0_idx == tau_coi0_idx || b_coi0_idx == tau_coi1_idx || "
                    "  b_coi1_idx == tau_coi0_idx || b_coi1_idx == tau_coi1_idx) ? 1.f : 0.f) : -1.f")

            # Cross-scores: does the best b-jet look like a tau? Vice versa?
            .Define("b_coi0_TauVsAll",
                    f"b_coi0_idx >= 0 ? ak4_TauVsAll[b_coi0_idx] : -1.f")
            .Define("tau_coi0_BvsAll",
                    f"tau_coi0_idx >= 0 ? ak4_BvsAll[tau_coi0_idx] : -1.f")
        )

    _coi_cols = [
        "b_coi0_is_genb", "b_coi1_is_genb", "b_both_coi_are_genb",
        "tau_coi0_is_gentau", "tau_coi1_is_gentau", "tau_both_coi_are_gentau",
        "b_tau_overlap_0", "any_b_tau_overlap",
        "b_coi0_TauVsAll", "tau_coi0_BvsAll",
    ]
    for name in mc:
        if name not in sig_samples:
            for col in _coi_cols:
                mc[name] = mc[name].Define(col, "-1.f")
