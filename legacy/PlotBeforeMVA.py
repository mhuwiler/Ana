#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import collections
import copy
import json
import anaConfig
ROOT.gROOT.LoadMacro("FileFlow.h+")
from ROOT import Ana, TCanvas, TH1D, TPad, TLegend, THStack, RDataFrame
from collections import defaultdict
from argparse import ArgumentParser
from anaPrepareRegions import PrepareRegions, PrepareRegionsSimple, PrepareSamples, UnrollHist, GetEventList


if __name__ == "__main__":

	parser = ArgumentParser(description="PlotBeforeMVA") 
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v7", help="Which version (cycle) of files to run on")
	parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="TauSelectionv7", help="Which version (cycle) of files to run on")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
	parser.add_argument('-d', "--denom", dest="denom", action="store_true", default=False, help="Denominator analysis")
	parser.add_argument('-f', "--full", dest="allvars", action="store_true", default=False, help="Plot all variables")
	parser.add_argument("-g", "--cut", dest="cut", action="store", type=str, default="1.", help="Custom cut to be included added")
	parser.add_argument("--veto", dest="vetofile", action="store", type=str, default="", help="Name of JSON file containing ID of selected events")
	parser.add_argument("--comp", dest="comp", action="store_true", default=False, help="Plot comparisons")
	parser.add_argument("--stack", dest="stack", action="store_true", default=False, help="Plot stacked distributions")

	

	options = parser.parse_args()

	
	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 

	outputfolder = "./plots/"+options.out+"/"
	#if (options.forcepath): 
	#	print("WARNING: You have used option '-f' or '--forcepath'. Files will be written to: {}".format(options.out))
	#	outputfolder = options.out+"/"

	os.system("mkdir -p "+outputfolder)


	postfix = "" #"_DNN_m"
	if (options.denom): 
		anaConfig.Denominator()
		postfix = "_DNN"




	# Global initialisations
	Ana.Init(options.version, options.denom)

	variables = anaConfig.variables
	if options.allvars: 
		variables = [item.first for item in Ana.binning]
	print(variables)

	samples = anaConfig.samples #[anaConfig.data, anaConfig.Sig] #anaConfig.samples
	sample = {}
	for item in samples:
		sample[item] = item+postfix
	

	#for file in sample.values(): 
	#	Ana.filemanager.OpenItem(file)



	frames, effs = PrepareSamples(list(sample.values()), options.cut)
	#if (not options.denom): print("WS yield {}".format(frames[anaConfig.dataWS]["baseline"].Count().GetValue()))

	if options.debug: print(frames)

	if ((options.vetofile != "") and (not options.denom)): 
		events = GetEventList(frames[anaConfig.data]["baseline"])

		vetofloder = Ana.folderbase + options.vetofile + "/numveto.json"
		with open(vetofloder, "w") as vetofile: 
			json.dump(events, vetofile, ensure_ascii=False, sort_keys=False) #encoding="utf8", 
	#import pickle
	#with open(options.vetofile, "w") as vetofile: 
	#	pickle.dump(events, vetofile)


	from anaPrepareRegions import SaveDataframe
	#SaveDataframe(frames[sample["data"]]["baseline"].Filter("b_B_fsig>2."), "dataFullWithCuts.root")

	# from here on starts teting

	from anaPlotting import PlotComparison, PlotStack

	if (options.comp): 
		PlotComparison(frames, sample[anaConfig.data], sample[anaConfig.Sig], ["baseline"], variables, outputfolder+"/comp/", False)

	if (options.stack):
		from libEfficiencies import MultiplyFinalEffs, ReadEffs, PrintEfficiencies2D
		from uncertainties import ufloat

		selectioneffs = ReadEffs(Ana.folder+"/Expectedyields.json")
		selectioneffs["dataB2WS"] = ufloat(1.39e-5*-12.0*7*8, 0.)
		selectioneffs["dataDWS"] = ufloat(1.39e-5*-12.0*11*3.3*1.2, 0.)
		selectioneffs["dataD1WS"] = ufloat(1.39e-5*-12.0*11*3.3*1.2*20, 0.)
		selectioneffs["WS"] = ufloat(310000/3.*1.1,0.)
		if (options.denom): 
			selectioneffs["WS"] = ufloat(35000.,0.)

		PrintEfficiencies2D(effs)

		regioneffs = MultiplyFinalEffs(selectioneffs, effs)

		PrintEfficiencies2D(regioneffs)

		PlotStack(frames, sample[anaConfig.data], list(sample.values()), ["baseline"], variables, regioneffs, outputfolder+"/stack/", False)
		from anaPlotting import PlotOverlay
		#PlotOverlay(frames, sample[anaConfig.data], list(sample.values()), ["baseline"], variables, regioneffs, outputfolder+"/comp/", False)


	Ana.filemanager.CloseAll()


	from libUtils import GetDate
	datestring = GetDate()

	from webInterface import PublishToWeb
	PublishToWeb(outputfolder, "Variables_{}_{}".format(options.out, datestring)) #Variables_23_8_14_beforeBDT_Sigvsdata



