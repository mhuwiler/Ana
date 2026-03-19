#!/usr/bin/env python
from __future__ import division, print_function

import os
import math
import numpy as np
import copy
from ROOT import RDataFrame, TGraph, TH1D



debugmode = False


def FindVariableRange(sample, variable, threshold=0.): 
	#hist = sample.Histo1D(variable)
	#maxvar = hist.GetXaxis().GetBinCenter(hist.FindFirstBinAbove(threshold))
	#minvar = hist.GetXaxis().GetBinCenter(hist.FindLastBinAbove(threshold))
	npvariable = sample.AsNumpy([variable])
	print(len(npvariable[variable]))
	if (len(npvariable[variable])):
		maxvar = max(npvariable[variable])
		minvar = min(npvariable[variable])
	else: 
		maxvar = 0. 
		minvar = 0.
	return minvar, maxvar

def GetROCgeneral(sig, bkg, variable, direction=True, bkgInSample=1., sigInSample=1.): 
	minvar, maxvar = FindVariableRange(sig, variable)
	if (maxvar == minvar): 
		graph = TGraph(100, np.zeros(100), np.ones(100))
		FOM = TH1D("FOM{}".format(variable), "", 100, 0., 100.)
		return -9999., -1., -999., graph, FOM
	from sklearn.metrics import roc_curve, auc
	columns = [variable]
	signal = sig.AsNumpy(columns)
	background = bkg.AsNumpy(columns)

	siglabels = signal[variable]
	if (maxvar == minvar): 
		dx=0.5
		maxvar = maxvar+dx
		minvar = minvar-dx
	siglabels = (siglabels - minvar)/(maxvar - minvar)
	#print(siglabels)
	#print(min(siglabels))
	#print(max(siglabels))

	#print(siglabels)

	sigtruth = np.ones(len(siglabels))

	#print(sigtruth)

	bkglabels = background[variable]
	bkgtruth = np.zeros(len(bkglabels))
	bkglabels = (bkglabels - minvar)/(maxvar - minvar)

	#print(bkglabels)
	
	#print(bkgtruth)
	if (not direction):
		siglabels = (-1.*siglabels)+1.
		bkglabels = (-1.*bkglabels)+1.

	labels = np.concatenate([siglabels, bkglabels])

	#print(labels)

	truths = np.concatenate([sigtruth, bkgtruth])

	bkgeff, sigeff, _ = roc_curve(truths, labels)

	graph = TGraph(len(bkgeff), np.asarray(bkgeff, "d"), np.asarray(sigeff, "d"))

	area = auc(bkgeff, sigeff)

	#Computing the FOM
	sigeffs = np.sort(sigeff)
	bkgeffs = np.sort(bkgeff)
	assert(len(sigeffs) == len(bkgeffs))
	numPoints = len(sigeffs)

	FOM = TH1D("FOM{}".format(variable), "", numPoints, minvar, maxvar)

	for point in range(0, numPoints): 
		sigEff = sigeffs[point]
		bkgEff = bkgeffs[point]
		
		if (debugmode): print("Sig eff: {}, bkg eff: {}".format(sigEff, bkgEff))

		B = bkgInSample*bkgEff
		S = sigInSample*sigEff

		Sigma = 0 if (B == 0) else S/math.sqrt(B) #Sigma = 0 if (S+B == 0) else S/math.sqrt(S+B)

		factor = 0.1
		factordenom = 0.0000001
		corrB = factor if (bkgInSample*bkgEff < factor) else bkgInSample*bkgEff
		corrS = factor if (sigInSample*sigEff < factor) else sigInSample*sigEff
		corrSigma = factordenom if (Sigma < factordenom) else Sigma

		# Error commputation taken from slide 13 in: /https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
		error = math.sqrt((S+1)*(S+2) - (S+1)*(S+1))/((B+2)*(B+3)-(B+2)*(B+2)) #(math.sqrt(corrB)/corrB)*(math.sqrt(corrS)/corrS)*corrSigma # *bkgInSample*sigInSample

		FOM.SetBinContent(numPoints -1 - point, Sigma)
		FOM.SetBinError(numPoints -1 - point, error)
		#FOM.SetBinError(point, error)

		maxvalue = FOM.GetMaximum()

		maxbin = FOM.GetMaximumBin()
		maxcut = FOM.GetXaxis().GetBinCenter(maxbin)

	return maxvalue, area, maxcut, graph, FOM

def GetROC(sig, bkg, mvavar = "mvaScore", sigtarget=1., bkgtarget=-1.): 
	from sklearn.metrics import auc
	bkgeff, sigeff = ComputeRoc(sig, bkg, mvavar, sigtarget, bkgtarget)

	# creating a TGraph from the efficiency points
	graph = TGraph(len(bkgeff), np.asarray(bkgeff, "d"), np.asarray(sigeff, "d"))

	area = auc(bkgeff, sigeff)

	return graph, area #copy.deepcopy(graph)


def ComputeRoc(sig, bkg, mvavar = "mvaScore", sigtarget=1., bkgtarget=-1.): 
	from sklearn.metrics import roc_curve
	columns = [mvavar]
	signal = sig.AsNumpy(columns)
	background = bkg.AsNumpy(columns)

	siglabels = signal[mvavar]

	#print(siglabels)

	sigtruth = sigtarget*np.ones(len(siglabels))

	#print(sigtruth)

	bkglabels = background[mvavar]
	bkgtruth = bkgtarget*np.ones(len(bkglabels))

	#print(bkglabels)
	
	#print(bkgtruth)

	labels = np.concatenate([siglabels, bkglabels])

	#print(labels)

	truths = np.concatenate([sigtruth, bkgtruth])

	bkgeff, sigeff, _ = roc_curve(truths, labels)

	return bkgeff, sigeff


def GetFom(sig, bkg,  mvavar = "mvaScore", sigInSample=1., bkgInSample=1.,sigtarget=1., bkgtarget=-1.): 
	bkgeffroc, sigeffroc = ComputeRoc(sig, bkg, mvavar, sigtarget, bkgtarget)
	sigeffs = np.sort(sigeffroc)
	bkgeffs = np.sort(bkgeffroc)
	assert(len(sigeffs) == len(bkgeffs))
	numPoints = len(sigeffs)

	if (debugmode): print(sigeffs)

	FOM = TH1D("FOM", "", numPoints, -1., 1.)

	for point in range(0, numPoints): 
		sigEff = sigeffs[point]
		bkgEff = bkgeffs[point]
		
		if (debugmode): print("Sig eff: {}, bkg eff: {}".format(sigEff, bkgEff))

		B = bkgInSample*bkgEff
		S = sigInSample*sigEff

		Sigma = 0 if (B == 0) else S/math.sqrt(B) #Sigma = 0 if (S+B == 0) else S/math.sqrt(S+B)

		factor = 0.1
		factordenom = 0.0000001
		corrB = factor if (bkgInSample*bkgEff < factor) else bkgInSample*bkgEff
		corrS = factor if (sigInSample*sigEff < factor) else sigInSample*sigEff
		corrSigma = factordenom if (Sigma < factordenom) else Sigma

		# Error commputation taken from slide 13 in: /https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
		error = math.sqrt((S+1)*(S+2) - (S+1)*(S+1))/((B+2)*(B+3)-(B+2)*(B+2)) #(math.sqrt(corrB)/corrB)*(math.sqrt(corrS)/corrS)*corrSigma # *bkgInSample*sigInSample

		FOM.SetBinContent(numPoints -1 - point, Sigma)
		FOM.SetBinError(numPoints -1 - point, error)
		#FOM.SetBinError(point, error)

		maxvalue = FOM.GetMaximum()

		maxbin = FOM.GetMaximumBin()
		maxcut = FOM.GetXaxis().GetBinCenter(maxbin)

	return FOM, maxvalue, maxcut


if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser(description="GetEfficiency")
	#parser.add_argument("filename", action="store", type=str, default="", help="Name of file")
	#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v6.8", help="Which version (cycle) of files to run on")
	parser.add_argument("-g", "--cut", dest="cut", action="store", type=str, default="1", help="Custom cut to be included in eff calculation")
	parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="./plots/testroc/", help="Path for saving plots")
	parser.add_argument('-l', "--variable", dest="variable", action="store", default="mvaScore", help="Variable to reprocess")
	parser.add_argument("-n", "--target", dest="target", action="store", type=float, default=10000., help="Target number of events after selection")

	options = parser.parse_args()


	samples = ["Sig", "dataD1"] #["Sig", "BkgDstarDs", "BkgDstarDsstar", "BkgB0DD", "BkgBuDXc"] , "BkgDstara1Part" # "dataB2", "dataD2WS", "dataD2_SB", 

	import ROOT

	ROOT.gROOT.LoadMacro("FileFlow.h")

	from ROOT import Ana
			
	Ana.Init(options.version)

	for file in samples: 
		Ana.filemanager.OpenItem(file)

	os.system("mkdir -p {}".format(options.out))

	

	# Loading signal and background 
	cutsig = (Ana.cut["base"]+Ana.samples.at("Sig").cut).GetTitle()
	cutbkg = (Ana.cut["base"]+Ana.samples.at("dataB2").cut).GetTitle()

	signal = RDataFrame(Ana.filemanager.GetItem("Sig")).Filter(cutsig)
	background = RDataFrame(Ana.filemanager.GetItem("dataD1")).Filter(cutbkg)

	print("Starting to compute ROC curve... ")
	roc, auc = GetROC(signal, background, options.variable)
	print("Computed ROC curve. ")

	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	roc.Draw("AP")
	canv.Draw()
	roc.SetMarkerColor(ROOT.kBlue)
	#roc.SetMarkerSize(1)
	#roc.SetMarkerStyle(8)
	canv.Draw()
	canv.Print(options.out+"ROC.pdf")
	print("Area under curve (A.U.C.): {}".format(auc))

	from libUtils import HoldUntilKeyPress
	HoldUntilKeyPress()

	fom, maxsig, cutvalue = GetFom(signal, background, options.variable)
	canv2 = ROOT.TCanvas("fom", "fom", 800, 600)
	fom.Draw("E")
	canv2.Draw()
	fom.GetXaxis().SetRangeUser(-1., 1.)
	#fom.SetMarkerColor(ROOT.kGreen+4)
	fom.SetMarkerStyle(8)
	fom.SetMarkerSize(0.2)
	fom.SetLineColor(ROOT.kGreen+2)
	#roc.SetMarkerSize(1)
	#roc.SetMarkerStyle(8)
	canv2.Draw()
	canv2.Print(options.out+"FOM.pdf")
	print("Maximum significance of {} with cut at value {}".format(maxsig, cutvalue))

	HoldUntilKeyPress()


	Ana.filemanager.CloseAll()


