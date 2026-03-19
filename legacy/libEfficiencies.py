#!/usr/bin/env python

from __future__ import division, print_function

import os
from ROOT import RDataFrame, TFile, TVectorD, gInterpreter
from uncertainties import ufloat
from uncertainties.umath import * 
import json
from collections import OrderedDict, defaultdict



def getEff(n, N): 
	if (N==0): 
		print("(getEff) WARNING: N is 0, returning efficiency = 0")
		return ufloat(0., 0.)
	eff = float(n)/float(N)
	#print eff
	err = sqrt(eff*(1.-eff)/float(N))
	#print err
	#return eff, err
	return ufloat(eff, err)

def getEffFromInfo(tree): 
	frame = RDataFrame(tree)

	n = frame.Sum("numSelected").GetValue()

	N = frame.Sum("numTotal").GetValue()

	eff = getEff(n, N)
	return eff

def getGenmatchingEff(tree, cut):
	frame = RDataFrame(tree)

	initialCount = frame.Count().GetValue()

	genmatchedframe = frame.Filter(cut)
	finalCount = genmatchedframe.Count().GetValue()

	return getEff(finalCount, initialCount)

gInterpreter.Declare("TVectorD castVectorD(TObject* obj) { return *static_cast<TVectorD*>(obj); };") #"TH1F * convertHisto(TH1D *histo) { return static_cast<TH1F*>(histo); } "

def getOfflineEff(vector):
	from ROOT import castVectorD
	vec = TVectorD(castVectorD(vector))
	initialCount = vec(0)
	finalCount = vec(1)

	return getEff(finalCount, initialCount)


def DumpEffs(effs, path): 
	effsForWrite = OrderedDict()
	for item, content in effs.items(): 
			eff = effs[item]
			effsForWrite[item] = (eff.n, eff.s)
	with open(path, "w") as file: 
		#effsForWrote = effsForWrite.encode("utf-8")
		json.dump(effsForWrite, file, ensure_ascii=False, encoding="utf8", sort_keys=False) #indent=4, 

def ReadEffs(path): 
	effs = OrderedDict()
	with open(path, "r") as file: 
		effsFromFile = json.load(file, object_pairs_hook=OrderedDict) #json.load(file, encoding="utf8", object_pairs_hook=OrderedDict)
		for item, content in effsFromFile.items(): 
				eff = effsFromFile[item]
				assert(len(eff)==2)
				effs[item] = ufloat(eff[0], eff[1])
	return effs

def DumpEffs2D(effs, path): 
	effsForWrite = defaultdict(dict)
	for item, content in effs.items(): 
		for key, value in content.items(): 
			eff = effs[item][key]
			effsForWrite[item][key] = (eff.n, eff.s)
	with open(path, "w") as file: 
		json.dump(effsForWrite, file, ensure_ascii=False, sort_keys=False) #encoding="utf8", 


def ReadEffs2D(path): 
	effs = defaultdict(dict)
	with open(path, "r") as file: 
		effsFromFile = json.load(file) #, encoding="utf8"
		for item, content in effsFromFile.items(): 
			for key, value in content.items(): 
				eff = effsFromFile[item][key]
				assert(len(eff)==2)
				effs[item][key] = ufloat(eff[0], eff[1])
	return effs

def PrintEfficiencies2D(effs): 
	for item, content in effs.items(): 
		print("{}:".format(item))
		for key, value in content.items(): 
			print(("\t{}: {}".format(key, value)))

def FormatLatex(number, precision = 2): 
	numstring = "{:.{precision}e}".format(number, precision=precision)
	collection = numstring.split("e")
	assert(len(collection)==2)
	num = float(collection[0])
	err = int(collection[1])
	#err = err.replace("+", "")
	latexstring = "{:.{precision}} \\times 10^{{{}}}".format(num, err, precision=precision)
	return latexstring

# For Ana regions
def GetEfficiencies(frames): 
	keys = [] #TODO: generate a json
	#for item in frames: 
		#localkeys = frames[item].keys()
		#print localkeys
	
	#Credit to: https://stackoverflow.com/questions/35491223/inverting-a-dictionary-with-list-values
	inv_frames = {}
	for k,v in list(frames.items()):
		for x in v:
			inv_frames.setdefault(x,[]).append(k)
	print(inv_frames)

	for key, value in list(inv_frames.items()): 
		print(key)
		for item in value: 
			n = frames[item][key].Count().GetValue()
			N = frames[item]["all"].Count().GetValue()
			eff = n/N
			print(("\t{}: {} ({}/{})".format(item, eff, n, N)))
		print("\n")

def ComputeEfficiencies(frames): 
	efficiencies = defaultdict(dict)
	for item, content in frames.items(): 
		#print("{}:".format(item))
		for key, value in content.items(): 
			n = frames[item][key].Count().GetValue()
			N = frames[item]["all"].Count().GetValue()
			# see: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
			ne = ufloat(n, sqrt(n))
			NE = ufloat(N, sqrt(N))
			efferr = ne/NE
			#print(efferr)
			#print(("\t{}: {} ({}/{})".format(key, eff, n, N)))
			efficiencies[item][key] = getEff(n, N) #ufloat(eff, err)
	return efficiencies

def InitialEffs(lumi): 
	# Computing the initial efficiencies 
	#effs = collections.defaultdict(dict)
	# Efficiencies relative to the 
	effs = { "Sig": {"br":ufloat(1.84e-2, 2.2e-3), "geneff":ufloat(3.72e-4, 0.), "eff":ufloat(1.4e-3, 0.)}, # TODO: group these with the others into another file 
		"BkgDstarDs": {"br":ufloat(8e-3, 1.1e-3), "geneff":ufloat(1.458e-3, 0.), "eff":ufloat(1.4e-3, 0.)},
		"BkgDstarDsstar": {"br":ufloat(1.77e-2, 1.4e-3), "geneff":ufloat(5.38e-4, 0.), "eff":ufloat(1.4e-3, 0.)}, 
		"BkgDstar3pi": {"br":ufloat(7.21e-3, 2.9e-4), "geneff":ufloat(2.e-5, 0.), "eff":ufloat(1.4e-3, 0.)}, # TODO: obtain ana eff from other script
		"SigPart": {"br":ufloat(1.84e-2, 2.2e-3), "geneff":ufloat(3.72e-4, 0.), "eff":ufloat(1.4e-3, 0.)},
		"BkgDstara1": {"br":ufloat(1.3e-2, 2.7e-3), "geneff":ufloat(3.800e-04, 0.), "eff":ufloat(1.4e-3, 0.)},
		"BkgDstara1Part": {"br":ufloat(1.3e-2, 2.7e-3), "geneff":ufloat(3.800e-04, 0.), "eff":ufloat(1.4e-3, 0.)},
	}
	bbxsec = ufloat(4.72e8, 0.)
	fB0 = fB = ufloat(0.404, 0.006)
	Br_Dstar_D0pi = ufloat(6.77e-1, 0.)
	Br_D0_KPI = ufloat(3.88e-2, 0.)

	expected = {}
	for key, eff in effs.items(): 
		key.replace("Part", "")
		expected[key]= lumi*bbxsec*fB0*2.*Br_Dstar_D0pi*Br_D0_KPI*1000.*eff["br"]*eff["geneff"] #*eff["eff"]

	return expected

def CompleteEffsFromFile(effs, version, filemanager): 
	anaeffs = {"Sig":ufloat(1.4e-3, 0.), "BkgDstarDs":ufloat(1.44e-3, 0.), "BkgDstarDsstar":ufloat(2.26e-3, 0.), "BkgDstar3pi":ufloat(5.3e-4, 0.), "SigPart":ufloat(1.27e-3, 0.), "dataD2WS":ufloat(-0.392, 0.), "dataD21TauWS": ufloat(-0.53, 0.)}
	for key, eff in effs.items(): 
		print(key)
		
		efficiency = 1.
		try: 
			file = TFile.Open(filemanager.GetFile(key+"_ntuple"), "READ")
			efftree = file.Get("ntuplizer/EffCalc")
			print(efftree)
			efficiency = getEffFromInfo(efftree)
			file.Close()
		except: 
			efficiency = anaeffs[key]
		effs[key] = eff*efficiency
		#if anaeffs[key].nominal_value < 0.: 
			#effs[key] = anaeffs[key]
	return effs

def ReadEffsFromFile(item, version, filemanager): 
	anaeffs = {"Sig":ufloat(1.4e-3, 0.), "BkgDstarDs":ufloat(1.44e-3, 0.), "BkgDstarDsstar":ufloat(2.26e-3, 0.), "BkgDstar3pi":ufloat(5.3e-4, 0.), "SigPart":ufloat(1.27e-3, 0.), "dataD2WS":ufloat(-0.392, 0.), "dataD21TauWS": ufloat(-0.53, 0.)}
	
	print(item)
		
	efficiency = 1.
	try: 
		file = TFile.Open(filemanager.GetFile(item+"_ntuple"), "READ")
		efftree = file.Get("ntuplizer/EffCalc")
		print(efftree)
		efficiency = getEffFromInfo(efftree)
		file.Close()
	except: 
		efficiency = anaeffs[item]
		#if anaeffs[key].nominal_value < 0.: 
			#effs[key] = anaeffs[key]
	return efficiency

def MultiplyEffs(effs, regioneffs, region):
	from anaPrepareRegions import GetBaseName
	result = {}
	for item, content in regioneffs.items(): 
		print(item)
		strippeditem = GetBaseName(item)
		try:
			result[strippeditem] = effs[strippeditem]*regioneffs[item][region]
			if (effs[strippeditem].nominal_value < 0.): 
				result[strippeditem]=effs[strippeditem]
		except:
			result[strippeditem] = ufloat(-1., 0.)
	return result

def MultiplyFinalEffs(effs, regioneffs):
	from anaPrepareRegions import GetBaseName
	for item, content in regioneffs.items(): 
		#item = GetBaseName(item)
		print("{}:".format(item))
		for key, value in content.items(): 
			try:
				regioneffs[item][key] = effs[GetBaseName(item)]*regioneffs[item][key]
				if (effs[GetBaseName(item)].nominal_value < 0.): 
					regioneffs[item][key]=effs[GetBaseName(item)]
			except:
				regioneffs[item][key] = ufloat(-1., 0.)
	return regioneffs


if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser(description="GetEfficiency")
	parser.add_argument("filename", action="store", type=str, default="", help="Name of file")
	#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
	#parser.add_argument("-i", "--input", dest="input", action="store", type=str, default="file.root", help="Number of points to be tested")
	parser.add_argument("-t", "--total", dest="total", action="store_true", default=False, help="Print total efficiency, including genmatching")
	parser.add_argument("-c", "--cut", dest="cut", action="store", type=str, default="Dstar_match&&pttau_tau_match", help="PCut to enable genmatching")
	parser.add_argument("-o", "--object", dest="object", action="store", type=str, default="ntuplizer/EffCalc", help="Efficiency info object within file")
	parser.add_argument("-e", "--tree", dest="tree", action="store", type=str, default="ntuplizer/tree", help="Event ntuple object within the file")

	options = parser.parse_args()


	file = TFile.Open(options.filename, "READ")

	info = file.Get(options.object)

	if not info: 
		raise ValueError("ERROR: No efficiency info found in file. Are you sure this file should contain efficiency information at {} ?".format(options.object))

	eff = getEffFromInfo(info)

	#print("Number of selected events: {}/{}".format(n, N))

	print("Efficiency: {}".format(eff))

	if options.total: 
		eventstree = file.Get(options.tree)
		if not eventstree: 
			raise ValueError("ERROR: No events tree '{}' found in file, please provide the correct events tree. ".format(options.tree))
		genmatcheff = getGenmatchingEff(eventstree, options.cut)
		totaleff = eff*genmatcheff
		print("Genmatching efficiency: {}\n".format(genmatcheff))
		print("Total efficiency: {}\n".format(totaleff))

