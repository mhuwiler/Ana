#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import copy
import json
from argparse import ArgumentParser
from uncertainties import ufloat
from uncertainties.umath import * 
import os # TODO: remove 


ROOT.gROOT.LoadMacro("FileFlow.h")

from ROOT import Ana, TFile

from libEfficiencies import getEffFromInfo, getGenmatchingEff, DumpEffs, ReadEffs, FormatLatex, getOfflineEff




# Efficiency calculations 



if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser(description="GetEfficiency")
	#parser.add_argument("filename", action="store", type=str, default="", help="Name of file")
	#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v7", help="Which version (cycle) of files to run on")
	parser.add_argument("-e", "--object", dest="object", action="store", type=str, default="ntuplizer/EffCalc", help="Efficiency info object within file")
	parser.add_argument("-g", "--cut", dest="cut", action="store", type=str, default="1", help="Custom cut to be included in eff calculation")

	options = parser.parse_args()


	filtereffs = ReadEffs("./data/etc/FilterEfficiencies.json")

	constants = ReadEffs("./data/etc/Constants.json")

	forcedbr = ReadEffs("./data/etc/ForcedBranchingFractions.json")




	samples = ["B0toDstar3pi", "B0toDstarrho0pi", "B0toDstara1"] 

	
	Ana.Init(options.version)


	effs = {}

	Eff = ufloat(0., 0.)

	for sample in samples: 
		print(sample)
		# Getting selection efficiency from file 
		item = sample+"_ntuple"
		Ana.filemanager.OpenItem(item)
		file = ROOT.TFile.Open(Ana.filemanager.GetFile(item), "READ")
		if file.IsZombie():
			print("ERROR: File {} of item {} does not exist... Skipping.".format(Ana.filemanager.GetFile(item), item))
			continue
		info = file.Get(options.object)

		if not info: 
			raise ValueError("ERROR: No efficiency info found in file. Are you sure this file should contain efficiency information at {} ?".format(options.object))

		eff = getEffFromInfo(info)

		offlinesample = sample+"_DNN"
		if not os.path.isfile(Ana.filemanager.GetFile(offlinesample)):
			print("ERROR: File {} of item {} does not exist... Skipping.".format(Ana.filemanager.GetFile(offlinesample), offlinesample))
			continue
		Ana.filemanager.OpenItem(offlinesample)

		if not Ana.filemanager.GetItem(offlinesample):
			print("ERROR: Object {} of item {} does not exist... Skipping.".format(Ana.filemanager.GetFile(offlinesample), offlinesample))
			continue

		genmatcheff = getGenmatchingEff(Ana.filemanager.GetItem(offlinesample), Ana.samples.at(sample).cut.GetTitle())

		customeff = getGenmatchingEff(Ana.filemanager.GetItem(item), options.cut)

		file = TFile.Open(Ana.filemanager.GetFile(offlinesample))
		if file.IsZombie():
			print("ERROR: Object {} of item {} does not exist... Skipping.".format(Ana.filemanager.GetFile(offlinesample), offlinesample))
			continue
		vec = file.Get("ntuplizer/EffUpdateTau")
		offlineeff = getOfflineEff(vec)
		#print(offlineeff)

		effs[sample] = filtereffs[sample]*eff*offlineeff 

		print("Efficiency: {}".format(effs[sample]))


	Br3pi = forcedbr["B0toDstar3pi"]
	Brrho0pi = forcedbr["B0toDstarrho0pi"]
	Bra1 = (forcedbr["B0toDstar3pitotal"]-forcedbr["B0toDstar3pi"]-forcedbr["B0toDstarrho0pi"])

	Eff = Br3pi*effs["B0toDstar3pi"] + Brrho0pi*effs["B0toDstarrho0pi"]+Bra1*effs["B0toDstara1"]
			

	print("Denominator signal efficiency: {}".format(Eff))

	print("Denominator 3pi efficiency: {}".format(forcedbr["B0toDstar3pi"]*effs["B0toDstar3pi"]))
	print("Denominator rho0 pi efficiency: {}".format(forcedbr["B0toDstarrho0pi"]*effs["B0toDstarrho0pi"]))
	print("Denominator a1 efficiency: {}".format(forcedbr["B0toDstara1"]*effs["B0toDstara1"]))

	Ana.filemanager.CloseAll()


