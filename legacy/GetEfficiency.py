#!/usr/bin/env python

from __future__ import division, print_function

import os
from ROOT import RDataFrame, TFile
from uncertainties import ufloat
from uncertainties.umath import * 


def getEff(n, N): 
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

