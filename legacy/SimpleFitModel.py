#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import collections
import copy
import json
from argparse import ArgumentParser
from ROOT import TCanvas, TH1D, TPad, TLegend, THStack, RDataFrame
from uncertainties import ufloat
from uncertainties.umath import * 


#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/libFunctions.C+")#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/ExperimentSpecificLayer.C")
#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/FileManager/CFileManager.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/CMS/tdrstyle.C")
#ROOT.gROOT.LoadMacro("FileFlow.h")
ROOT.gROOT.LoadMacro("Tau.h")
#ROOT.setTDRStyle()
#import CMS_lumi
import anaConfig
from ROOT import Ana



webpublication =False



if __name__ == "__main__":

	parser = ArgumentParser(description="SaveRegions")
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v1", help="Which version (cycle) of files to run on")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
	
	options = parser.parse_args()


	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 

	# Global initialisations
	Ana.Init(options.version)

	nBins = 6
	rangeMin = 0.2 #0.37
	rangeMax = 1.5 #1.43


	# Starting script 
	samples = {}
	frames = collections.defaultdict(dict)
	baseline = collections.defaultdict(dict)
	anasamples = anaConfig.samples

	#anasamples = {"Sig", "data", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "WS"}

	# for item in anasamples:  
	# 	Ana.filemanager.OpenItem(item)
	# 	samples[item] = ROOT.RDataFrame(Ana.filemanager.GetItem(item))
	# 	frames[item]["baseline"] =  samples[item].Filter((Ana.cut["base"]+Ana.samples.at(item).cut).GetTitle()) #Ana.cut["base"].GetTitle() "1."
	# 	frames[item]["all"] = samples[item].Filter("1.")
	# 	#baseline[item]
	# 	for region in anaConfig.regions: 
	# 		cut = (Ana.cut[region]+Ana.samples.at(item).cut).GetTitle()
	# 		#if "WS" in item: 
	# 		#	cut = Ana.cutstandalone[region].GetTitle()
	# 		if (options.debug): print("Using following cut string (from TCut): {}".format(cut))
	# 		frames[item][region] = samples[item].Filter(cut)
	# 		ROOT.SetOwnership(frames[item][region], 0)
	# 		# For histogram legacy compatibility

	from anaPrepareRegions import PrepareSamples
	frames, effs = PrepareSamples(anasamples) #samples = Ana.GetSamples(anasamples)

	if options.debug: print(frames)


	ROOT.gInterpreter.Declare('''
		template<typename T>
		void fixStringVariables(T &dataframe)
		{
			#include "stringbranches.gcf"

			for (auto branch : stringbranches) // Hack to fix string branche 
			{
				dataframe = dataframe.Redefine(branch, [](const ROOT::RVec<std::string> &v) {return std::vector<std::string>(v.begin(), v.end());}, {branch}); 
			}
		}
		'''
	)

	from libFitting import WriteWorkspace, WriteDatacardSimple
	from libEfficiencies import ReadEffs2D

	yields = ReadEffs2D("./data/etc/RegionEffs.json")

	regions = ["baseline"] #anaConfig.regions

	variables = ["b_B_m"] #anaConfig.variables

	fitvariable = "b_B_m"

	WriteWorkspace(frames, effs, variables, regions, anaConfig.data, "workspaceFromExportLatest.root")

	WriteDatacardSimple(frames, effs, variables, regions, fitvariable, anaConfig.data, "datacardGeneratedSimpleLatest.txt", "workspaceFromExportLatest.root")

	Ana.filemanager.CloseAll()



