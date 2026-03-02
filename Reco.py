#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
ROOT.gROOT.LoadMacro("Particle.h+")
ROOT.gROOT.LoadMacro("HHbbtautauAnaElements.C+")
import os
import math
import collections
import copy
import json
import anaConfig
from ROOT import Ana, TCanvas, TH1D, TPad, TLegend, THStack, RDataFrame
from argparse import ArgumentParser


blacklist = ["Muon_P4", "GenPart_Particle", "GenMuon", "HLT*", "L1*"] # TODO: add autoblacklist


filedict = {"sigggF": "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv15/signalggF.root", }
Ana.filemanager.AddItem("sigggF", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv15/signalggF.root", "Events")
Ana.filemanager.AddItem("official", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv12/Run3Summer22NanoAODv12_1-1.root", "Events")
Ana.filemanager.AddItem("ggfBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/glugluHHto2b2tau/ggf.root", "Events")
Ana.filemanager.AddItem("VBFBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/VBFHHto2b2tau/VBF_SM.root", "Events")
Ana.filemanager.AddItem("QCDBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/QCD/qcd_HT_100-1200.root", "Events")
Ana.filemanager.AddItem("VBFBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/dataset/jetmet.root", "Events")
Ana.filemanager.AddItem("ggfBoostedPrivateLarge", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/glugluHHto2b2tau/ggf_large.root", "Events")


def loadFile(desc): 
	frame = ROOT.RDataFrame(Ana.filemanager.GetItem(desc, True))
	ROOT.SetOwnership(frame, 0)
	return generalise(frame)


def loadFileFull(desc, writemode = False, treename = "Events"): 
	path = filedict[desc]
	flag = "WRITE" if writemode else "READ" 
	file = ROOT.TFile.Open(path, flag)
	tree = copy.deepcopy(file.Get(treename))
	tree.Print()
	frame = ROOT.RDataFrame(tree)
	print(frame)
	ROOT.SetOwnership(frame, 0)
	return frame


def dropBranchNames(frame, filename, exclusionlist = []):
	with open(filename, "w") as file: 
		for name in frame.GetColumnNames(): 
			name = str(name)
			#print("{} {}".format(name, [(excluded in name) for excluded in exclusionlist]))
			if not (any([excluded in name for excluded in exclusionlist])): 
				file.write("{}\n".format(name))

def WriteFile(sample, filename, blacklist, treename = "Events"): 
	sample.Snapshot(treename, filename, Ana.purgeColumns(sample.GetColumnNames(), blacklist))


def generalise(df): 
	return ROOT.ROOT.RDF.AsRNode(df)


def AdditionalVariables(sample): 
	# Adding AK8 jet related varibles

	# Defining n-subjettiness ratios 
	sample = sample.Define("{}_tau21".format("FatJet"), "{0}_tau2/{0}_tau1".format("FatJet"))
	sample = sample.Define("{}_tau32".format("FatJet"), "{0}_tau3/{0}_tau2".format("FatJet"))


	# Constructing normalised discriminators
	sample = sample.Define("{}_{}_QCD".format("FatJet", anaConfig.tauIDvar), "{0}_{1}_QCD0HF+{0}_{1}_QCD1HF+{0}_{1}_QCD2HF".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Top".format("FatJet", anaConfig.tauIDvar), "{0}_{1}_TopbW+{0}_{1}_TopW".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtauh{}".format("FatJet", anaConfig.tauIDvar, "vsQCD"), "{0}_{1}_Xtauhtauh/({0}_{1}_Xtauhtauh+{0}_{1}_QCD)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaum{}".format("FatJet", anaConfig.tauIDvar, "vsQCD"), "{0}_{1}_Xtauhtaum/({0}_{1}_QCD+{0}_{1}_Xtauhtaum)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaue{}".format("FatJet", anaConfig.tauIDvar, "vsQCD"), "{0}_{1}_Xtauhtaue/({0}_{1}_QCD+{0}_{1}_Xtauhtaue)".format("FatJet", anaConfig.tauIDvar))

	sample = sample.Define("{}_{}_Xtauhtauh{}".format("FatJet", anaConfig.tauIDvar, "vsQCDTop"), "{0}_{1}_Xtauhtauh/({0}_{1}_Xtauhtauh+{0}_{1}_QCD+{0}_{1}_Top)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaum{}".format("FatJet", anaConfig.tauIDvar, "vsQCDTop"), "{0}_{1}_Xtauhtaum/({0}_{1}_Xtauhtaum+{0}_{1}_QCD+{0}_{1}_Top)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaue{}".format("FatJet", anaConfig.tauIDvar, "vsQCDTop"), "{0}_{1}_Xtauhtaue/({0}_{1}_Xtauhtaue+{0}_{1}_QCD+{0}_{1}_Top)".format("FatJet", anaConfig.tauIDvar))

	sample = sample.Define("{}_{}_Xbb{}".format("FatJet", anaConfig.tauIDvar, "vsQCD"), "{0}_{1}_Xbb/({0}_{1}_QCD+{0}_{1}_Xbb)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xbb{}".format("FatJet", anaConfig.tauIDvar, "vsQCDTop"), "{0}_{1}_Xbb/({0}_{1}_QCD+{0}_{1}_Xbb+{0}_{1}_Top)".format("FatJet", anaConfig.tauIDvar))

	sample = sample.Define("{}_{}_Xtauhtauh{}".format("FatJet", anaConfig.tauIDvar, "vsbb"), "{0}_{1}_Xtauhtauh/({0}_{1}_Xtauhtauh+{0}_{1}_Xbb)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaum{}".format("FatJet", anaConfig.tauIDvar, "vsbb"), "{0}_{1}_Xtauhtaum/({0}_{1}_Xtauhtaum+{0}_{1}_Xbb)".format("FatJet", anaConfig.tauIDvar))
	sample = sample.Define("{}_{}_Xtauhtaue{}".format("FatJet", anaConfig.tauIDvar, "vsbb"), "{0}_{1}_Xtauhtaue/({0}_{1}_Xtauhtaue+{0}_{1}_Xbb)".format("FatJet", anaConfig.tauIDvar))

	# TODO: add those variables https://github.com/LPC-HH/bbtautau/blob/591d8d1acea32896ef652157d5547ffccd861130/src/bbtautau/processors/objects.py#L75


	return sample


def DefineTauTaggerVarsForJet(theJet, sample): 
	sample = sample.Define("{}_{}_Xtauhtaum".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtaum, {})".format(anaConfig.tauIDvar, theJet)).Define("{}_{}_Xtauhtauh".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtauh, {})".format(anaConfig.tauIDvar, theJet)).Define("{}_{}_Xtauhtaue".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtaue, {})".format(anaConfig.tauIDvar, theJet))
	sample = sample.Define("{}_{}_XtauhtaumvsQCD".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_XtauhtaumvsQCD, {})".format(anaConfig.tauIDvar, theJet)).Define("{}_{}_XtauhtauhvsQCD".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_XtauhtauhvsQCD, {})".format(anaConfig.tauIDvar, theJet)).Define("{}_{}_XtauhtauevsQCD".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_XtauhtauevsQCD, {})".format(anaConfig.tauIDvar, theJet))
	return sample


def DefineBTaggerVarsForJet(theJet, sample): 
	sample = sample.Define("{}_{}_Xbb".format(theJet, anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xbb, {})".format(anaConfig.tauIDvar, theJet)).Define("{}_btagDeepB".format(theJet), "Ana::overflowProtected(FatJet_btagDeepB, {})".format(theJet))
	return sample


def ApplyTriggerSelection(sample): 
	hadronicTriggers = ["HLT_AK8PFJet250_SoftDropMass40_PFAK8ParticleNetBB0p35", "HLT_AK8PFJet230_SoftDropMass40_PFAK8ParticleNetTauTau0p30", "HLT_AK8PFJet230_SoftDropMass40_PNetBB0p06", "HLT_AK8PFJet230_SoftDropMass40_PNetTauTau0p03", "HLT_AK8PFJet420_MassSD30", "HLT_AK8PFJet425_SoftDropMass40",]
	resolvedTriggers = ["HLT_QuadPFJet70_50_40_35_PFBTagParticleNet_2BTagSum0p65", "HLT_QuadPFJet70_50_40_35_PNet2BTagMean0p65", "HLT_QuadPFJet103_88_75_15_PFBTagDeepJet_1p3_VBF2", "HLT_QuadPFJet103_88_75_15_DoublePFBTagDeepJet_1p3_7p7_VBF1",]
	hadronicTriggers += resolvedTriggers
	tauTriggers = ["HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1", "HLT_DoubleMediumDeepTauPFTauHPS35_L2NN_eta2p1", "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet60", "HLT_DoubleMediumDeepTauPFTauHPS30_L2NN_eta2p1_PFJet75",]
	muonTriggers = ["HLT_IsoMu24", "HLT_Mu50", "HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1",]
	electronTriggers = ["HLT_Ele30_WPTight_Gsf", "HLT_Ele115_CaloIdVT_GsfTrkIdT", "HLT_Ele50_CaloIdVT_GsfTrkIdT_PFJet165", "HLT_Photon200", "HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1",]
	METTriggers = ["HLT_PFMET120_PFMHT120_IDTight", "HLT_PFHT280_QuadPFJet30_PNet2BTagMean0p55", "HLT_PFHT340_QuadPFJet70_50_40_40_PNet2BTagMean0p70",]

	# define which triggers to be selected on
	triggers = hadronicTriggers + tauTriggers + muonTriggers + electronTriggers + METTriggers

	
	branches = [str(name) for name in sample.GetColumnNames()] # Get the column names
			

	activetriggers = []

	for path in triggers: # check which triggers exist in the sample
		if ((path in branches)): 
			print(path)
			activetriggers.append(path)


	triggerselection = "||".join(activetriggers)
	print(triggerselection)


	sample = sample.Filter(triggerselection)

	return sample





if __name__ == "__main__":

	parser = ArgumentParser(description="Selection") 
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v7", help="Which version (cycle) of files to run on")
	parser.add_argument("--prefix", dest="xrdpfx", action="store", type=str, default="root://cms-xrd-global.cern.ch//", help="XRootD prefix to be used to access files")
	parser.add_argument("-o", "--out", dest="outputpath", action="store", type=str, default="./temp", help="Local output path")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("--test", dest="test", action="store_true", default=False, help="Process a reduced number of files for testing purposes")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")

	

	options = parser.parse_args()

	
	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 


	Ana.Init()


	#ROOT.gSystem.Load("HHbbtautauAnaELements.so")
	#ROOT.gSystem.Load("MyDict.so")
	

	sig = loadFile("ggfBoostedPrivate")

	if (options.test): 
		sig = generalise(sig.Range(0, 500000))

	ROOT.gStyle.SetOptStat(0)

	print(sig)

	dropBranchNames(sig, "branchnamesNanoAODv12.txt", ["L1", "HLT", "DST"])

	#sig = Ana.GetP4(sig, "Muon") #sig = Ana.GetP4["float"](sig, "Muon")
	#sig = Ana.GetGenParticles(sig, "Muon")
	#sig = Ana.GetGenParticles(sig, "Electron")

	genprefix = "GenPart"

	cutflow = ROOT.CutFlow("cutflow", "Selection cutflow")

	cutflow.Increment("Initial", sig.Count().GetValue())

	sig = sig.Filter("nFatJet>=2").Filter("FatJet_pt[0]>250&&FatJet_pt[1]>200")

	cutflow.Add("jet selection", sig.Count().GetValue())

	sig = ApplyTriggerSelection(sig)

	#sig = sig.Define("GenDecay", "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)".format("GenPart"))

	#sig = sig.Define("TheGenMuon_pt", "Ana::overflowProtected(GenPart_pt, GenDecay.mu)").Define("TheGenMuon_eta", "GenPart_eta[GenDecay.mu]").Define("TheGenMuon_phi", "GenPart_phi[GenDecay.mu]")


	sig = AdditionalVariables(sig)
	

	# TODO: try out taking the max from tauhtauh,tauhtaumu, tauhtaue
	differentialTagger = "vsbb" #vsQCDTop
	sig = sig.Define("{}_{}_tautau{}".format("FatJet", anaConfig.tauIDvar, differentialTagger), "{0}_{1}_Xtauhtaum{2}+{0}_{1}_Xtauhtauh{2}+{0}_{1}_Xtauhtaue{2}".format("FatJet", anaConfig.tauIDvar, differentialTagger))
	sig = sig.Define("TheTauFatJet", "Ana::RecoTauJet( {0}_pt, {0}_eta, {0}_phi, {0}_{1})".format("FatJet", "{}_tautau{}".format(anaConfig.tauIDvar, differentialTagger)))
	#sig = sig.Define("TheTauFatJet", "1")
	#sig = sig.Define("ThebFatJet", "0")

	sig = sig.Define("ThebFatJet", "Ana::RecoBJet( {0}_pt, {0}_eta, {0}_phi, {0}_{1}, TheTauFatJet)".format("FatJet", "{}_Xbb{}".format(anaConfig.tauIDvar, "vsQCD")))
	#sig = sig.Define("TheTauFatJet", "1-ThebFatJet")

	sig = sig.Define("TheTauFatJet_eta", "Ana::overflowProtected(FatJet_eta, TheTauFatJet)").Define("TheTauFatJet_phi", "Ana::overflowProtected(FatJet_phi, TheTauFatJet)")
	sig = sig.Define("TheTauFatJet_{}_Xtauhtaum".format(anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtaum, TheTauFatJet)".format(anaConfig.tauIDvar)).Define("TheTauFatJet_{}_Xtauhtauh".format(anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtauh, TheTauFatJet)".format(anaConfig.tauIDvar)).Define("TheTauFatJet_{}_Xtauhtaue".format(anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xtauhtaue, TheTauFatJet)".format(anaConfig.tauIDvar)).Define("TheBFatJet_{}_Xbb".format(anaConfig.tauIDvar), "Ana::overflowProtected(FatJet_{}_Xbb, ThebFatJet)".format(anaConfig.tauIDvar))
	#sig = sig.Define("TheTauFatJet_eta", "Ana::overflowProtected(FatJet_eta, TheTauFatJet)").Define("TheTauFatJet_phi", "Ana::overflowProtected(FatJet_phi, TheTauFatJet)")
	sig = sig.Define("TheMuon", "Ana::RecoMuon( {0}_pt, {0}_eta, {0}_phi, {0}_tightId, {0}_dz, {0}_dxy, {1}_eta, {1}_phi)".format("Muon", "TheTauFatJet"))


	

	sig = sig.Filter("TheTauFatJet>=0&&ThebFatJet>=0") #&&TheMuon>=0

	cutflow.Add("topology reco", sig.Count().GetValue())

	sig = sig.Filter("(TheTauFatJet!=ThebFatJet)") # removing overlap between bb and tautau jets

	cutflow.Add("no overlap b tau", sig.Count().GetValue())


	gensig = sig.Define("GenDecay", "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)".format("GenPart"))

	gensig = gensig.Define("GenMatchedTauFatJet", "Ana::closestMatch(GenDecay.Htotau, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))
	gensig = gensig.Define("GenMatchedbFatJet", "Ana::closestMatch(GenDecay.Htob, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))
	gensig = gensig.Define("OverlapbTauJet", "GenMatchedTauFatJet==GenMatchedbFatJet")
	gensig = gensig.Define("dR_gen_reco_bb", "Ana::deltaR(GenDecay.Htob, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedbFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Define("dR_gen_reco_tautau", "Ana::deltaR(GenDecay.Htotau, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedTauFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Define("dR_gen_reco_b1", "Ana::deltaR(GenDecay.b1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedbFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Define("dR_gen_reco_b2", "Ana::deltaR(GenDecay.b2, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedbFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Define("dR_gen_reco_tau1", "Ana::deltaR(GenDecay.tau1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedTauFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Define("dR_gen_reco_tau2", "Ana::deltaR(GenDecay.tau2, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedTauFatJet)".format("GenPart", "FatJet"))
	gensig = gensig.Filter("dR_gen_reco_bb<0.8&&dR_gen_reco_tautau<0.8&&dR_gen_reco_b1<0.8&&dR_gen_reco_b2<0.8")
	#gensig = gensig.Filter("dR_gen_reco_tau1<0.8&&dR_gen_reco_tau2<0.8")
	gensig = gensig.Define("ClosestFatJettoGenMatchedbFatJet", "Ana::closestMatch(GenMatchedbFatJet, {0}_pt, {0}_eta, {0}_phi, {0}_mass)".format("FatJet"))
	gensig = gensig.Define("dR_reco_bb_closest", "Ana::deltaR(ClosestFatJettoGenMatchedbFatJet, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenMatchedbFatJet)".format("FatJet"))
	gensig = gensig.Define("ClosestJettoGenMatchedbFatJet", "Ana::closestMatchAboveThreshold(GenMatchedbFatJet, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("FatJet", "Jet"))
	gensig = gensig.Define("dR_reco_bb_closest_jet", "Ana::deltaR(ClosestJettoGenMatchedbFatJet, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, GenMatchedbFatJet)".format("FatJet", "Jet"))


	gensig = DefineTauTaggerVarsForJet("GenMatchedTauFatJet", gensig)
	gensig = DefineBTaggerVarsForJet("GenMatchedbFatJet", gensig)
	gensig = DefineTauTaggerVarsForJet("GenMatchedbFatJet", gensig)
	gensig = DefineBTaggerVarsForJet("GenMatchedTauFatJet", gensig)
	gensig = gensig.Define("SubSubLeadingFatJet", "3")
	gensig = DefineTauTaggerVarsForJet("SubSubLeadingFatJet", gensig)
	gensig = DefineBTaggerVarsForJet("SubSubLeadingFatJet", gensig)



	gensig = gensig.Define("dR_gen_reco_mu", "Ana::deltaR(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, TheMuon)".format("GenPart", "Muon"))
	gensig = gensig.Define("TheGenMuon_pt", "Ana::overflowProtected(GenPart_pt, GenDecay.mu)").Define("TheGenMuon_eta", "GenPart_eta[GenDecay.mu]").Define("TheGenMuon_phi", "GenPart_phi[GenDecay.mu]")
	gensig = gensig.Define("TheRecoMuon_pt", "Ana::overflowProtected(Muon_pt, TheMuon)").Define("TheRecoMuon_eta", "Muon_eta[TheMuon]").Define("TheRecoMuon_phi", "Muon_phi[TheMuon]")

	gensig = gensig.Define("dR_HH", "Ana::deltaR(GenDecay.Htotau, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.Htob)".format("GenPart"))

	gensig = gensig.Define("dR_tautau", "Ana::deltaR(GenDecay.tau1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.tau2)".format("GenPart"))

	gensig = gensig.Define("dR_mu_gen", "Ana::deltaR(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, {1}_pdgId, \"{1}\")".format("GenPart", "Muon"))

	gensig = gensig.Define("closest_mu_gen", "Ana::closestMatchBelowThreshold(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, {1}_pdgId, \"{1}\")".format("GenPart", "Muon"))

	gensig = gensig.Define("closest_mu_FatJet", "Ana::closestMatch(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))

	gensig = gensig.Define("dR_mu_FatJet", "Ana::deltaR(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))

	
	hh = gensig.Filter("GenDecay.decayType==1")

	hm = gensig.Filter("GenDecay.decayType==2")

	cutflow.Add("tauh taumu", hm.Count().GetValue())

	#hi = hm.Filter("dR_mu_FatJet<1.6&&dR_mu_FatJet>0")

	#cutflow.Add("gen mu in jet", hi.Count().GetValue())

	#hi = hi.Define("TheRecoMuon_pt", "Muon_pt[closest_mu_gen]")

	#hi = hi.Filter("Muon_tightId[closest_mu_gen]")

	#cutflow.Add("reco mu matched", hi.Count().GetValue())

	#hi = hi.Filter("Muon_pt[closest_mu_gen]>20&&abs(Muon_eta[closest_mu_gen])<2.4&&Muon_dz[closest_mu_gen]<0.2&&Muon_dxy[closest_mu_gen]<0.05")

	#cutflow.Add("muon sel", hi.Count().GetValue())

	#hh = hh.Define("dR_tautau", "Ana::deltaR(GenDecay.tau1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.tau2)".format("GenPart"))


	
	WriteFile(sig, "./SigAll.root", blacklist)

	WriteFile(gensig, "./SigAllwithGen.root", blacklist)

	gensig = gensig.Filter("{0}_{1}_Xbb[{2}]<0.4&&{0}_particleNetLegacy_Xbb[{2}]>0.6".format("FatJet", anaConfig.tauIDvar, "GenMatchedbFatJet"))
	WriteFile(gensig, "./SigAllwithGen_misclassified.root", blacklist)

	#hh.Snapshot("Events", "./Sighh.root", Ana.purgeColumns(hh.GetColumnNames(), blacklist))

	#hm.Snapshot("Events", "./Sighm.root", Ana.purgeColumns(hm.GetColumnNames(), blacklist))

	#hi.Snapshot("Events", "./Sigwithproxy.root", Ana.purgeColumns(hi.GetColumnNames(), blacklist))

	#print("Initial: {}, 2 FatJets: {}, tauhtaumu: {}, gen mu within jet: {}, reco mu within jet: {}, other selection requirements: {}".format(n0, n1, n2, n3, n4, n5))

	# TODO: implement this into the cutflow class 
	cutflow.Print()


	canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
	histo = cutflow.GenerateHistogram()
	histo.Draw("HIST")
	histo.SetLineColor(ROOT.kBlue)
	histo.SetLineWidth(2)
	canvas.Draw()
	canvas.SaveAs("cutflow.png")

	efficiencies = ROOT.TCanvas("efficiencies", "efficiencies", 800, 600)
	eff = cutflow.MakeEffHistogram()
	eff.Draw("HIST")
	eff.SetLineColor(ROOT.kBlue)
	eff.SetLineWidth(2)
	efficiencies.Draw()
	efficiencies.SaveAs("efficiencies.png")

	Ana.filemanager.CloseAll()

	print("Hello")


	
	




	

