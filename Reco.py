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


filedict = {"sigggF": "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv15/signalggF.root", }
Ana.filemanager.AddItem("sigggF", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv15/signalggF.root", "Events")
Ana.filemanager.AddItem("official", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv12/Run3Summer22NanoAODv12_1-1.root", "Events")
Ana.filemanager.AddItem("ggfBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/glugluHHto2b2tau/ggf.root", "Events")
Ana.filemanager.AddItem("VBFBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/VBFHHto2b2tau/VBF_SM.root", "Events")
Ana.filemanager.AddItem("QCDBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/QCD/qcd_HT_100-1200.root", "Events")
Ana.filemanager.AddItem("VBFBoostedPrivate", "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/dataset/jetmet.root", "Events")


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


def generalise(df): 
	return ROOT.ROOT.RDF.AsRNode(df)



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
		sig = generalise(sig.Range(0, 5000))

	ROOT.gStyle.SetOptStat(0)

	print(sig)

	dropBranchNames(sig, "branchnames.txt", ["L1", "HLT", "DST"])

	#sig = Ana.GetP4(sig, "Muon") #sig = Ana.GetP4["float"](sig, "Muon")
	#sig = Ana.GetGenParticles(sig, "Muon")
	#sig = Ana.GetGenParticles(sig, "Electron")

	genprefix = "GenPart"

	cutflow = ROOT.CutFlow("cutflow", "Selection cutflow")

	n0 = sig.Count().GetValue()
	cutflow.Increment("Initial", n0)

	sig = sig.Filter("nFatJet>=2").Filter("FatJet_pt[0]>250&&FatJet_pt[1]>200")

	n1 = sig.Count().GetValue()
	cutflow.Add("jet selection", n1)

	#sig = sig.Define("GenDecay", "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)".format("GenPart"))

	#sig = sig.Define("TheGenMuon_pt", "Ana::overflowProtected(GenPart_pt, GenDecay.mu)").Define("TheGenMuon_eta", "GenPart_eta[GenDecay.mu]").Define("TheGenMuon_phi", "GenPart_phi[GenDecay.mu]")

	sig = sig.Define("TheTauFatJet", "Ana::RecoTauJet( {0}_pt, {0}_eta, {0}_phi, {0}_{1})".format("FatJet", "{}_Xtauhtaum".format(anaConfig.tauIDvar)))

	sig = sig.Define("TheRecoMuon", "Ana::RecoMuon( {0}_pt, {0}_eta, {0}_tightId, {0}_dz, {0}_dxy)".format("Muon"))


	gensig = sig.Define("GenDecay", "Ana::DecayGenMatching({0}_pdgId, {0}_genPartIdxMother, {0}_statusFlags)".format("GenPart"))

	gensig = gensig.Define("dR_HH", "Ana::deltaR(GenDecay.Htotau, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.Htob)".format("GenPart"))

	gensig = gensig.Define("dR_tautau", "Ana::deltaR(GenDecay.tau1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.tau2)".format("GenPart"))

	gensig = gensig.Define("dR_mu_gen", "Ana::deltaR(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, {1}_pdgId, \"{1}\")".format("GenPart", "Muon"))

	gensig = gensig.Define("closest_mu_gen", "Ana::closestMatch(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass, {1}_pdgId, \"{1}\")".format("GenPart", "Muon"))

	gensig = gensig.Define("closest_mu_FatJet", "Ana::closestMatch(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))

	gensig = gensig.Define("dR_mu_FatJet", "Ana::deltaR(GenDecay.mu, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {1}_pt, {1}_eta, {1}_phi, {1}_mass)".format("GenPart", "FatJet"))

	
	hh = gensig.Filter("GenDecay.decayType==1")

	hm = gensig.Filter("GenDecay.decayType==2")

	n2 = hm.Count().GetValue()
	cutflow.Add("tauh taumu", n2)

	hi = hm.Filter("dR_mu_FatJet<1.6&&dR_mu_FatJet>0")

	n3 = hi.Count().GetValue()
	cutflow.Add("gen mu in jet", n3)

	hi = hi.Define("TheRecoMuon_pt", "Muon_pt[closest_mu_gen]")

	hi = hi.Filter("Muon_tightId[closest_mu_gen]")

	n4 = hi.Count().GetValue()
	cutflow.Add("reco mu matched", n4)

	hi = hi.Filter("Muon_pt[closest_mu_gen]>20&&abs(Muon_eta[closest_mu_gen])<2.4&&Muon_dz[closest_mu_gen]<0.2&&Muon_dxy[closest_mu_gen]<0.05")

	n5 = hi.Count().GetValue()
	cutflow.Add("muon sel", n5)

	#hh = hh.Define("dR_tautau", "Ana::deltaR(GenDecay.tau1, {0}_pt, {0}_eta, {0}_phi, {0}_mass, {0}_pt, {0}_eta, {0}_phi, {0}_mass, GenDecay.tau2)".format("GenPart"))


	blacklist = ["Muon_P4", "GenPart_Particle", "GenMuon", "HLT*", "L1*"] # TODO: add autoblacklist

	sig.Snapshot("Events", "./SigAll.root", Ana.purgeColumns(sig.GetColumnNames(), blacklist))

	#hh.Snapshot("Events", "./Sighh.root", Ana.purgeColumns(hh.GetColumnNames(), blacklist))

	#hm.Snapshot("Events", "./Sighm.root", Ana.purgeColumns(hm.GetColumnNames(), blacklist))

	#hi.Snapshot("Events", "./Sigwithproxy.root", Ana.purgeColumns(hi.GetColumnNames(), blacklist))

	print("Initial: {}, 2 FatJets: {}, tauhtaumu: {}, gen mu within jet: {}, reco mu within jet: {}, other selection requirements: {}".format(n0, n1, n2, n3, n4, n5))

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


	
	




	

