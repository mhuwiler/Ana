"""RDataFrame .Define() chains for kinematic variables."""

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
        # -- Gen-match AK8 jets to gen H→bb and H→ττ --
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

    # Dummy gen-match columns for non-signal samples (AK8 only)
    _gen_match_cols_ak8 = ["ak8_genHbb_match_dR", "ak8_genHbb_match_Xbb",
                           "ak8_genHbb_match_XbbVsAll", "ak8_genHtt_match_dR",
                           "ak8_genHtt_match_Xtt", "ak8_genHtt_match_XttVsAll"]
    for name in mc:
        if name not in sig_samples:
            for col in _gen_match_cols_ak8:
                mc[name] = mc[name].Define(col, "-1.f")


def define_kinematics(df):
    """Define derived kinematic columns for jets.

    Intermediate columns are avoided to keep the branch-proxy chain shallow
    and prevent stack overflows with large TChains.
    """
    df = (df
        .Define("ak4_pt0", "ScoutingPFJetRecluster2_pt[0]")
        .Define("ak4_pt1", "ScoutingPFJetRecluster2_pt[1]")
        .Define("ak4_pt2", "ScoutingPFJetRecluster2_pt[2]")
        .Define("ak4_pt3", "ScoutingPFJetRecluster2_pt[3]")
        .Define("HT",  "Sum(ScoutingPFJetRecluster2_pt)")
        .Define("nJets", "nScoutingPFJetRecluster2")
        .Define("nLeptons", "nScoutingMuonVtx + nScoutingElectron")
    )
    # -- AK8 fat jet kinematics + ScoutGlobalParT tagger + XvsQCD (pT > 150 GeV) --
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
    # b-jet selection: UParT discriminators
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
        .Define("htautau_pair",
                "Ana::findDijetInWindow("
                "ScoutingPFJetRecluster2_pt, ScoutingPFJetRecluster2_eta,"
                "ScoutingPFJetRecluster2_phi, ScoutingPFJetRecluster2_mass,"
                "40.f, 150.f, 80.f, {hbb_pair.i1, hbb_pair.i2})")
        .Define("has_htautau",  "htautau_pair.i1 >= 0")
        .Define("mtautau_cand", "htautau_pair.mass")
    )
    return df


def define_gen_matched_ak4(mc, sig_samples):
    """Define gen-matched AK4 jet columns (signal only, after define_kinematics)."""
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
