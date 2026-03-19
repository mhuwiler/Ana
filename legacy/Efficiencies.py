#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import copy
import json
from argparse import ArgumentParser
from uncertainties import ufloat
from uncertainties.umath import * 
import os # TODO: remove 


#ROOT.gROOT.LoadMacro("FileFlow.h")
import anaConfig

from ROOT import Ana, TFile

from libEfficiencies import getEffFromInfo, getGenmatchingEff, DumpEffs, ReadEffs, FormatLatex, getOfflineEff




# Efficiency calculations 



if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser(description="GetEfficiency")
	#parser.add_argument("filename", action="store", type=str, default="", help="Name of file")
	#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v7", help="Which version (cycle) of files to run on")
	parser.add_argument("-l", "--lumi", dest="lumi", action="store", type=float, default=41.5, help="Luminostiy processed")
	parser.add_argument("-e", "--object", dest="object", action="store", type=str, default="ntuplizer/EffCalc", help="Efficiency info object within file")
	parser.add_argument("-g", "--cut", dest="cut", action="store", type=str, default="1", help="Custom cut to be included in eff calculation")
	parser.add_argument("-t", "--table", dest="table", action="store", type=str, default="/Users/mhuwiler/cernbox/DoctoralThesis/Analysis/Presentations/PresentationVFS_24_3_12/efftable.tex", help="Latex fragment with summary table")
	parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="Expectedyields.json", help="Path for json with yield info")
	parser.add_argument("-n", "--target", dest="target", action="store", type=float, default=10000., help="Target number of events after selection")
	parser.add_argument('-m', "--denom", dest="denom", action="store_true", default=False, help="Denominator analysis")

	options = parser.parse_args()


	filtereffs = ReadEffs("./data/etc/FilterEfficiencies.json")

	constants = ReadEffs("./data/etc/Constants.json")

	forcedbr = ReadEffs("./data/etc/ForcedBranchingFractions.json")

	finalbr = ReadEffs("./data/etc/FinalStateBranchingFractionsGenerated.json")

	if (options.denom): anaConfig.Denominator()

	#print(constants)

	#print(filtereffs)

	#N = constants["sigmabb"]*constants["fB0"]*forcedbr["Sig"]

	#print("N expected: {}".format(N))

	expected = {}


	samples = anaConfig.samples #["Sig", "B0toDstarD0K", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarDs1", "B0toDstarD0Kstar", "B0toDstarDs0star", "B0toDstara1", "B0toDstarD", "B0toDstar3pi", "B0toDstar3pipi0", "B0toDstar5pi", "ButoDstarDK", "B0toDstarDsX", "ButoDstarXc", "ButoDstarDK", "BstoDD", "B0toDstarrho0pi", "B0toDstarKKstar", "B0toDstarKpipi"] #["Sig", "BkgDstarDs", "BkgDstarDsstar", "BkgB0DD", "BkgBuDXc"] , "BkgDstara1Part"
	samples.remove(anaConfig.data)
	#samples.remove(anaConfig.dataWS)

	Ana.Init(options.version, options.denom)

	outputfile = Ana.folder+"/"+options.out 

	template = "{} & ${}$ & ${}$ & ${}$ & ${}$ & ${:fL}$ \\\\\n" #"{} & ${}$ & ${}$ & ${}$ & ${}$ & ${:fL}$ & ${}$ \\\\\n" #{:.1e} "{} & ${:.{precision}eL}$ & ${:.{precision}eL}$ & ${:.{precision}eL}$ & ${:fL}$ & ${}$ \\\\\n"

	numberafterselection = options.target

	Nexpected = {}

	Sum = 0

	with open(options.table, "w") as outfile: 
		outfile.write("\\begin{tabular}{lcccccr}\n")
		outfile.write("sample & $\\epsilon_{filter}$ & $\\epsilon_{ana}$ & $\\epsilon_{match}$ & Br & $N_{exp.}$ & N requested \\\\\n\\hline\n")

		print("Expected yields")
		for sample in samples: 
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

			finaleff = ufloat(1., 0.)
			if (options.denom): getGenmatchingEff(Ana.filemanager.GetItem(offlinesample), Ana.cut.at("base").GetTitle())

			customeff = getGenmatchingEff(Ana.filemanager.GetItem(item), options.cut)

			file = TFile.Open(Ana.filemanager.GetFile(offlinesample))
			if file.IsZombie():
				print("ERROR: Object {} of item {} does not exist... Skipping.".format(Ana.filemanager.GetFile(offlinesample), offlinesample))
				continue
			vec = file.Get("ntuplizer/EffUpdateTau")
			offlineeff = getOfflineEff(vec)
			print(offlineeff)

			N = options.lumi*constants["sigmabb"]*constants["fB0"]*2*forcedbr[sample]*constants["BrDstar2D0pi"]*constants["BrD02Kpi"]*1000*filtereffs[sample]*eff*offlineeff #*genmatcheff*customeff
			print("Eff for fit: {}".format(options.lumi*constants["sigmabb"]*constants["fB0"]*2*constants["BrDstar2D0pi"]*constants["BrD02Kpi"]*1000*filtereffs[sample]*eff*offlineeff*finaleff))
			# Adding Br error ad hoc. 
			try: 
				finalstatebr = finalbr[sample]
			except: 
				finalstatebr = ufloat(1., 0.1)
				print("Warning: No final state branching fraction found for sample {}. Setting uncertainty to 0.1".format(sample))
			BRerror = ufloat(1., finalstatebr.s/finalstatebr.n)
			N = N*BRerror
			Nexpected[sample] = N
			Sum += N
			#print(10000./eff)
			print("\tN expected for {}: {} (filter eff: {}, ana eff: {}, genmatch eff: {}, number requested: {})".format(sample, N, filtereffs[sample], eff, genmatcheff, numberafterselection/eff))
			if (genmatcheff==0): genmatcheff=ufloat(1., 0) # Hack to avoid division by zero
			numrequested = numberafterselection/(eff*offlineeff*genmatcheff)
			n = round(numrequested.n, -3)
			samplename = sample # TODO; implement handlig of casw where sample name not found
			outfile.write(template.format(Ana.samples.at(sample).latex, FormatLatex(filtereffs[sample].n), FormatLatex(eff.n), FormatLatex(genmatcheff.n), FormatLatex(forcedbr[sample].n), N, FormatLatex(numrequested.n), precision=2).replace("\\times", "\\cdot"))

		outfile.write("\\end{tabular}\n")

		# Add here computation for WS sample 
		# e.g. data - sum 

	DumpEffs(Nexpected, outputfile)

	Ana.filemanager.CloseAll()


