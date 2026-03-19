#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import collections
import copy
import json
from argparse import ArgumentParser
from ROOT import TCanvas, TTree, TFile, TH1D, TPad, TLegend, THStack, RDataFrame
from uncertainties import ufloat
from uncertainties.umath import *



class CacheFile: 
	@classmethod
	def SetFile(self, somefile): 
		self.file = somefile
		self.open = True

	file = ""
	folder = "./cache/"
	open = False

def PrintEvents(tree, options): 
	cache = CacheFile
	if (options.cut != ""): 
		if (options.cache): 
			ROOT.gSystem.Exec("mkdir -p {}".format(cache.folder))
			cachefile = TFile.Open(cache.folder+"/TreeForDecayString.root", "RECREATE")
			cachefile.cd()
			cache.SetFile(cachefile)
			#ache.file = cachefile
		else:
			ROOT.gROOT.cd(); # Making sure the tree with the cut is memory resident 
		tree = tree.CopyTree(options.cut)

	eventLoopFactory = ROOT.EventLoopFactory(tree, options.full, options.max)
	eventLoopFactory.PrintEvents()

	if (cache.open): 
		cache.file.Close(); 
		ROOT.gSystem.Exec("rm -r {}".format(cache.folder))



if __name__ == "__main__":

	parser = ArgumentParser(description="EventPlotter")
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("--out", dest="out", action="store", type=str, default="BackgroundEstimateUpdate/", help="Directory where the plots shuld go")
	parser.add_argument("--file", dest="file", action="store", type=str, default="test", help="Turn on debug output")
	parser.add_argument("--sample", dest="sample", action="store", type=str, default="", help="Turn on debug output")
	parser.add_argument("--full", dest="full", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("--tree", dest="tree", action="store", type=str, default="ntuplizer/tree", help="Turn on debug output")
	parser.add_argument("--cut", dest="cut", action="store", type=str, default="", help="Cut to apply to the tree")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v5", help="Which version (cycle) of files to run on")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("-f", "--forcepath", dest="forcepath", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
	parser.add_argument("--stats", dest="stats", action="store_true", default=False, help="Show stats box in ROOT")
	parser.add_argument("--max", dest="max", action="store", default=-1, help="How many events to plot")
	#parser.add_argument("--cache", dest="cache", action="store", type=str, default="./cache/", help="Cached tree for cuts")
	parser.add_argument("--cache", dest="cache", action="store_true", default=False, help="Cache tree for cuts")

	options = parser.parse_args()

	
	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 

	if (options.stats): 
		ROOT.gStyle.SetOptStat(1111111)
	else: 
		# Don't plot stats box
		ROOT.gStyle.SetOptStat(0) 


	ROOT.gInterpreter.Declare("""
	int PrintEventsEvtLoop(std::string decay)
	{
		std::cout << decay << std::endl; 
		return 1; 
	}
	""")

	#if (options.cut != ""): 
	#	ROOT.gSystem.Exec("mkdir -p {}".format(options.cache))
	#	cachefile = TFile.Open(options.cache+"TreeForDecayString.root", "WRITE")


	ROOT.gROOT.LoadMacro("FileFlow.h")
	ROOT.gROOT.LoadMacro("EventLoopFactory.C")
	from ROOT import Ana

	if (options.sample == ""): 
		file = TFile.Open(options.file, "READ")
		tree = file.Get(options.tree)
		
		PrintEvents(tree, options)

		file.Close()

	#frame = RDataFrame(tree)

	else: 
		Ana.Init(options.version)

		Ana.filemanager.OpenItem(options.sample)

		#filtered = frame.Filter("pttau_tau_m>1.5").Define("numPrinted", ROOT.PrintDecayString, ("genstring"))
		tree = Ana.filemanager.GetItem(options.sample)
		
		PrintEvents(tree, options)

		#ROOT.PrintDecayString(tree, options.full) #+"_ntuple"

		Ana.filemanager.CloseAll()

	print("Done")
	#file.Close()

