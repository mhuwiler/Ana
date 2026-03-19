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


def loadFile(desc): 
	frame = ROOT.RDataFrame(Ana.filemanager.GetItem(desc, True))
	ROOT.SetOwnership(frame, 0)
	return generalise(frame)


def loadFileFull(desc, writemode = False, treename = "Events"): 
	path = filedict[desc]
	return loadFileBase(path, writemode, treename)


def loadFileBase(path, writemode = False, treename = "Events"): 
	flag = "WRITE" if writemode else "READ" 
	globals()["file"] = ROOT.TFile.Open(path, flag) # need to make this global to preserve pointer outside function 
	#tree = copy.deepcopy(file.Get(treename))
	tree = file.Get(treename)
	#tree.Print()
	frame = generalise(ROOT.RDataFrame(tree))
	#print(frame)
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
	parser.add_argument("--file", dest="file", action="store", type=str, default="", help="File name")
	parser.add_argument("--tree", dest="tree", action="store", type=str, default="Events", help="Path of tree within file")
	parser.add_argument("--variables", dest="variables", action="store", nargs="+", default="", help="List of variables to plot")
	parser.add_argument("-o", "--out", dest="outputpath", action="store", type=str, default="./temp", help="Local output path")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("--folder", dest="folder", action="store", type=str, default="./", help="Folder to store the output plots")
	parser.add_argument("--selection", dest="selection", action="store", type=str, default="1", help="Selection to be applied to events")
	parser.add_argument("--test", dest="test", action="store_true", default=False, help="Process a reduced number of files for testing purposes")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")

	

	options = parser.parse_args()

	
	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 


	Ana.Init()


	#ROOT.gSystem.Load("HHbbtautauAnaELements.so")
	#ROOT.gSystem.Load("MyDict.so")
	

	#sample = loadFile("sigggF")
	sample = loadFileBase(options.file)

	if (options.test): 
		sample = generalise(sample.Range(0, 5000))

	print(sample)

	nBins = 50
	hist = ("", "#mu p_{T};#mu p_{T};", nBins, 0., 30.)

	from anaPlotting import PlotSimple

	if options.selection != "1": 
		sample = generalise(sample.Filter(options.selection))

	for var in options.variables: 
		if "[" in var: 
			varname = "autoVar_{}".format(var.replace("[", "_").replace("]", "_"))
			#print(var)
			sampleToPlot = generalise(sample.Define(varname, var)) #sample.Define("TheRecoMuon_pt", "Muon_pt[TheRecoMuon]") 
			PlotSimple(sampleToPlot, varname, options.folder)
		else: 
			PlotSimple(sample, var, options.folder)

	
	Ana.filemanager.CloseAll()

	print("Hello")


	
	




	

