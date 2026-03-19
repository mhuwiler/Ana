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


	# Global initialisations
	Ana.Init(options.version)

	nBins = 6
	rangeMin = 0.2 #0.37
	rangeMax = 1.5 #1.43


	# Starting script 
	from anaPrepareRegions import PrepareRegionsSimple, PrepareSamples
	#frames = PrepareRegionsSimple()
	frames, effs = PrepareSamples(anaConfig.samples)
	

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

	from libFitting import WriteWorkspace, WriteWorkspaceDataset, WriteWorkspaceRooType
	from libEfficiencies import ReadEffs2D, MultiplyFinalEffs, ReadEffs

	# yields = ReadEffs2D("./data/etc/RegionEffs.json")
	selectioneffs = ReadEffs(Ana.folder+"/Expectedyields.json")
	yields = MultiplyFinalEffs(selectioneffs, effs)

	WriteWorkspace(frames, effs, ["b_B_q2"], ["baseline"], anaConfig.data, "workspaceCombinedWODirq2.root")
	#WriteWorkspaceDataset(frames, effs, anaConfig.variables, ["baseline"], anaConfig.data, "workspaceFromExportHistsUpdate.root")

	Ana.filemanager.CloseAll()



