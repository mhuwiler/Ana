#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import uproot
import numpy as np
import collections
import copy
from datetime import datetime
from ROOT import RooRealVar, RooArgSet, RooDataHist, RooArgList, RooFormulaVar, RooAddition

uproot.default_library = "np"


#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/libFunctions.C+")#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/ExperimentSpecificLayer.C")
#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/FileManager/CFileManager.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/CMS/tdrstyle.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
ROOT.gROOT.LoadMacro("FileFlow.h")
ROOT.gROOT.LoadMacro("RooParametricHist.h")
#ROOT.setTDRStyle()
#import CMS_lumi

from ROOT.Ana import filemanager


ROOT.Ana.Init("v1") 


bkgInSample = 9400000000*0.0000376
sigInSample = 1130


graphcollection = ROOT.vector('std::pair<TGraph*,TString>')()


plotstats = False

webpublication =False


outputfolder = "./plots/BackgroundModeling/"

os.system("mkdir -p "+outputfolder)

logscale = [""]

#presentationfolder = "../Presentations/Presentation_22_5_10/reveal.js-master/Figures.root"


#ROOT.ExperimentSpecificLayer.SetStyle("CMS")

#layer = ROOT.ExperimentSpecificLayer("CMS", 21, 11, "Preliminary"); 

#factor = 1.0
#CMS_lumi.writeExtraText = True
#CMS_lumi.extraText = "Simulation Preliminary"
#CMS_lumi.lumiTextSize = 0.45*factor
#CMS_lumi.lumiTextOffset = 0.2*factor
#CMS_lumi.cmsTextSize = 0.75*factor
#CMS_lumi.lumi_sqrtS = "14 TeV"
#CMS_lumi.lumi_13TeV = "2021"

def formatHistoCommon(hist): 
	hist.SetLineWidth(2)
	hist.SetLineStyle(1)
	hist.SetMarkerStyle(8) # None 0
	hist.SetMarkerSize(0.4)
	#hist.SetLineColor(4)
	hist.GetYaxis().SetRangeUser(0, hist.GetMaximum()*1.4)
	#hist.GetYaxis().SetTitle("Events")
	hist.GetYaxis().SetTitleOffset(1.3)
	if quantity in logscale:
		hist.GetYaxis().SetRangeUser(0.1, hist.GetMaximum()*1000)
		canvas.SetLogy()
	#hist.GetXaxis().SetTitleSize(0.05)
	hist.GetXaxis().SetTitleOffset(1.1)
	canvas.SetTopMargin(0.1)
	#hist.GetXaxis().SetTitle("Probe jet #tau_{21}")
	#hist.GetXaxis().SetTitleSize(0.05)
	#hist.GetXaxis().SetTitleOffset(1.2)
	#hist.GetYaxis().SetTitleSize(0.05)
	#hist.GetYaxis().SetTitleOffset(1.2)
	#hist.GetXaxis().SetLabelSize(0.)
	#hist.GetXaxis().SetLabelOffset(-999.)

	return #histo


def formatAxisCommon(axis): 
	axis.SetTitleSize(0.05)
	axis.SetTitleOffset(1.2)

	return

def PromptYesNo(answerasbool=False): 
		# Inspired from Fabrice Couderc 
		rep = ''
		while not rep in [ 'yes', 'no' ]:
			rep = input( "(type 'yes' or 'no'): " ).lower()
		if (answerasbool): 
			if (rep == 'yes'): 
				return True
			else: 
				return False
		return rep

def getStatsBox(histogram):  #Works 
	dummycanvas = ROOT.TCanvas("dummycanvas", "dummycanvas", 800, 600)
	#histogram.SetStats(True)
	newhistogram = histogram.Clone(histogram.GetName())
	newhistogram.SetStats(True)
	newhistogram.Draw()
	dummycanvas.Update()
	stats = newhistogram.GetListOfFunctions().FindObject("stats").Clone("stats"+histogram.GetName())
	ROOT.SetOwnership(stats, 0)
	return stats

def exportHistFromTree(tree, variable, cut = "1", nbins = 0, minbin = None, maxbin = None): 
	dummycanvas = ROOT.TCanvas("dummycanvas", "dummycanvas", 800, 600)
	binning = ""; 
	if (nbins): 
		if (minbin == None): minbin = tree.CopyTree(cut).GetMinimum(variable) 
		if (maxbin == None): maxbin = tree.CopyTree(cut).GetMaximum(variable)
		binning = "({}, {}, {})".format(nbins, minbin, maxbin)
	tree.Draw(variable+">>h"+binning, cut);
	histo = tree.GetHistogram(); 
	histo.SetDirectory(0); 
	return histo; 

def UnrollHist(histo2D, inverted=True): 
	nx = histo2D.GetNbinsX()
	ny = histo2D.GetNbinsY()
	if inverted: 
		ny = histo2D.GetNbinsX()
		nx = histo2D.GetNbinsY()

	nTotal = nx*ny

	unrolled = ROOT.TH1D("unrolled", "unrolled", nTotal, 0, 100)

	print("Nunmber of bins: {}, {}".format(nx, ny))

	for i in range(0, nx):
		for j in range(0, ny): # TODO: check overflow is handled properly 
			if inverted: 
				binContent = histo2D.GetBinContent(j, i)
			else: 
				binContent = histo2D.GetBinContent(i, j)
			unrolled.SetBinContent(i+j*nx, binContent)

	return unrolled


def BackgroundShapeUnrolled(reference, estimate, additionalhists): 
	#estimate = histos["data"]["CR"]
	#reference = histos["data"]["SR"]
	#additionalhists = [histos["MC"]["SR"], histos["DstarDs"]["SR"]]
	#estimate = UnrollHist(estimate.Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2"))
	#reference = UnrollHist(reference.Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2"))
	#additionalhists = []
	#for item in additionalframes: 
		#additionalhists.append(UnrollHist(SR[key].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")))

	colors = [ROOT.kGreen+3, ROOT.kBlue+3] #ROOT.kOrange
	legends = ["Signal MC (signal region)",  "B#rightarrowD*D_{s} MC (signal region)"]
	norms = [1., 1.]

	assert(len(colors) >= len(additionalhists))
	assert(len(legends) >= len(additionalhists))
	assert(len(norms) >= len(additionalhists))

	# Plot the histogram
	estimate.Draw("LE")
	reference.Draw("LE SAME")
	for i, hist in enumerate(additionalhists): 
		additionalhists[i].Draw("LE SAME")

	estimate.GetXaxis().SetRangeUser(18, 100)

	reference.Scale(1./reference.Integral())
	estimate.Scale(1./estimate.Integral())
	estimate.Sumw2()
	reference.Sumw2()
	for i, hist in enumerate(additionalhists): 
		hist.Scale(norms[i]/hist.Integral())
		hist.Sumw2()

	maxes = [estimate.GetMaximum(), reference.GetMaximum()]
	for item in additionalhists: 
		maxes.append(item.GetMaximum())
	themax = max(maxes)

	estimate.SetMaximum(1.1*themax)

	estimate.SetLineColor(ROOT.kRed)
	reference.SetLineColor(ROOT.kBlue)
	for i, hist in enumerate(additionalhists): 
		hist.SetLineColor(colors[i]) 

	estimate.SetTitle("Unrolled 2D #rho mass distribution")

	estimateErr = estimate.DrawCopy("HIST SAME")
	referenceErr = reference.DrawCopy("HIST SAME")
	additionalhistsErr = []
	for hist in additionalhists: 
		additionalhistsErr.append(hist.DrawCopy("HIST SAME"))

	legend = ROOT.TLegend(.65,.60,.90,.85)
	legend.SetBorderSize(0)
	legend.SetFillColor(0)
	legend.SetFillStyle(0)
	legend.SetTextFont(42)
	legend.SetTextSize(0.04)
	legend.AddEntry(referenceErr,"data (signal region)","l")
	legend.AddEntry(estimateErr,"data (sideband)","l")
	for i, hist in enumerate(additionalhists): 
		legend.AddEntry(hist,legends[i],"l")
	legend.Draw()

	histoRatio = estimate.Clone()
	histoRatio.Divide(reference.Clone())

	ratiocanvas = ROOT.TCanvas("ratiocanvas", "ratiocanvas", 800, 600)

	histoRatio.Draw("HIST")
	histoRatio.SetLineColor(1)
	histoRatio.DrawCopy("LE SAME")


	canvas.Draw()

	#histoRatio.SetMaximum(0.2)
	#histoRatio.SetMinimum(1.8)
	histoRatio.GetYaxis().SetRangeUser(0.2, 1.8)
	histoRatio.GetYaxis().SetNdivisions(5)
	#canvas.RemoveMiddleAxis()

	canvas.Update()

		#FOM = filemanager.GetItem("cutflowGen")

	fomcanvas = ROOT.TCanvas("fomcanvas", "fomcanvas", 800, 600)

	

	#formatHistoCommon(histgwf)
#	ROC.SetLineColor(ROOT.kBlue+2)
#	ROC.SetLineWidth(2)
#	ROC.GetXaxis().SetLabelSize(textsize)
#	ROC.GetYaxis().SetLabelSize(textsize)

#	FOM.SetLineColor(ROOT.kGreen)
#	#FOM.SetLineWidth(2)
#	FOM.SetMarkerStyle(7)
#	FOM.SetMarkerColor(ROOT.kGreen)
#	FOM.GetXaxis().SetLabelSize(textsize)
#	FOM.GetYaxis().SetLabelSize(textsize)


	splitfraction = 0.8

	sc = 1.05

	# Making the ratio plots 




	#if plotstats: 
	#	scale = 1./1.2
	#else: 
	#	scale = 1.
	#legend = ROOT.TLegend(.62*scale,.60,.90*scale,.85)
	#legend.SetBorderSize(0)
	#legend.SetFillColor(0)
	#legend.SetFillStyle(0)
	#legend.SetTextFont(42)
	#legend.SetTextSize(0.04)
	#legend.AddEntry(histgwf,"GPU workflow","l")
	#legend.Draw()

	#histlwf.Scale(histgwf.Integral()/histlwf.Integral())



   
	canvas.Update()

	canvas.Print(outputfolder+quantity+".pdf")
	if (webpublication): canvas.Print(webfolder+quantity+".png")


	#f = ROOT.TFile.Open(presentationfolder, "Update")
	#fomcanvas.Write()
	#canvas.Write()
	#f.Write()
	#f.Close()

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
			print("\t{}: {}".format(item, frames[item][key].Count().GetValue()))
		print("\n")


def AtomicDraw(histo, name, options = ""): 
	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	histo.DrawCopy(options)
	canv.Draw()
	canv.Print(name)

ROOT.gInterpreter.Declare("TH1F * convertHisto(TH1D *histo) { TH1F *newHisto; newHisto = new TH1F(); histo->Copy(*newHisto); return newHisto; } ") #"TH1F * convertHisto(TH1D *histo) { return static_cast<TH1F*>(histo); } "

def SuperimposeRegions(frames, sample, variable, outputname, model=("model", "", 50, 0., 2.)): # TODO: make histo stack
	count = 0
	shapes = {}
	maxes=[]
	drawn = {}
	for region in regions: 
		shapes[region] = frames[sample][region].Histo1D(model, variable)
		shapes[region].Sumw2()
		shapes[region].SetLineColor(ROOT.kAzure+1+count)
		shapes[region].SetLineWidth(2)
		shapes[region].Scale(1./shapes[region].Integral())
		maxes.append(shapes[region].GetMaximum())
		count+=1

	newcanvas = ROOT.TCanvas("newcanvas", "newcanvas", 800, 600)
	legend = ROOT.TLegend(	canvas.GetLeftMargin()+0.35, 
								1-canvas.GetTopMargin()-.2, 
								canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
								1-canvas.GetTopMargin())
	drawn["SR"] = shapes["SR"].DrawCopy("HIST E")
	drawn["CR"] = shapes["CR"].DrawCopy("HIST E SAME")
	drawn["SB"] = shapes["SB"].DrawCopy("HIST E SAME")
	drawn["SR"].SetMaximum(max(maxes)*1.3)
	legend.AddEntry(drawn["SR"], "SR", "L")
	legend.AddEntry(drawn["CR"], "CR", "L")
	legend.AddEntry(drawn["SB"], "SB", "L")
	legend.Draw()
	newcanvas.Draw()
	newcanvas.Print(outputname) #outputfolder+"DstarDsShapes.pdf"


if plotstats: 
	ROOT.gStyle.SetOptStat(1111111)

# Web publication
if (webpublication): 
	webfolder = "/eos/home-m/mhuwiler/www/Analysis/BackgroundModellingNew/"
	os.system("mkdir -p "+webfolder)
	webenginesource = "/eos/home-m/mhuwiler/software/php-plots/"
	os.system("cp -r "+webenginesource+"res "+webfolder)
	os.system("cp "+webenginesource+"index.php "+webfolder)
	with open(webenginesource+"example/htaccess", "r") as permissionfile: 
		content = permissionfile.read()
		content = content.replace("/<me>/<my-project>/", webfolder)
		file = open(webfolder+".htaccess", "w")
		file.write(content)
		file.close()
	os.system("cp "+webfolder+".htaccess "+webfolder+"htaccess")
	webfolder = webfolder+"plots/"
	os.system("mkdir -p "+webfolder)

textsize = 0.04

numEvents = -1 



#ROOT.gInterpreter.Declare("""
#	double Rhomass2DUnrolled(Float_t rhomass1, Float_t rhomass2)
#	{
#		return int((std::min(rhomass2, 1.3) - 0.2)/0.22) + 6*int((std::min(tau_rhomass1, 1.3) - 0.2)/0.22); 
#	}
#""")

# Don't plot stats box
ROOT.gStyle.SetOptStat(0)

for quantity in ["Rhomass2Dunrolled"]: 
	print("Plotting {}".format(quantity))
	#canvas = ROOT.TCanvas("romassunrolled", "Unrolled 2D distribution of rho mass", 800, 600) #ROOT.RatioCanvas(quantity, quantity, 950, 800) #800, 800
	#canvas.SetMiddleMargin(0.13)
	#canvas.SetPadDelimitation(0.34)
	#canvas.SetRightMargin(0.12)
	ROOT.gROOT.SetBatch(1)

	filesUsed = ["dataD2", "Sig", "SigOld", "SigPart", "BkgDstarDs", "BkgDstar3pi", "BkgDstarDsstar", "dataD2WS"] #, "DstarDsMCfirst", "Data2018BFirst"

	filemap = {"data":"dataD2", "MC":"Sig", "DstarDs":"BkgDstarDs", "DstarDsstar":"BkgDstar3pi", "Dstar3pi":"BkgDstarDsstar"}

	colors = {	"dataD2": ROOT.kBlue, 
				"Sig": ROOT.kRed, 
				"SigPart": ROOT.kRed+2, # TODO: add color scheme to FileFlow.h
				"BkgDstarDs": ROOT.kGreen, 
				"BkgDstar3pi": ROOT.kGreen+3, 
				"BkgDstarDsstar": ROOT.kOrange+2
			}

	legends = {	"dataD2": "data", 
				"Sig": "signal MC (genmatched)", 
				"SigPart": "partially reconstructed signal",
				"BkgDstarDs": "B^{0}#rightarrow D*D_{s} Inclusive", 
				"BkgDstar3pi": "B^{0}#rightarrow D*D_{s}* Inclusive", 
				"BkgDstarDsstar": "B^{0}#rightarrow D*3#pi Nonresonant"
			}

	regions = ["SR", "CR", "SB"]

	for item in filesUsed: 
		print("Opening file: {}".format(item)) 
		filemanager.OpenItem(item); 


	canvas = ROOT.TCanvas("romassunrolled", "Unrolled 2D distribution of rho mass", 800, 600)

	# Histo parameters
	nBins = 6
	rangeMin = 0.37
	rangeMax = 1.43

	samples = {}
	frames = collections.defaultdict(dict)
	histos = collections.defaultdict(dict)
	histosunrolled = collections.defaultdict(dict)
	cut = {}

	# Signal region definitions
	mvaThreshold = 0.9
	mvaLowThreshold = 0.0
	mvaLowerBound = -0.5

	cut["base"] = "(b_Ds_vprob>0.1) && (b_D0_vprob>0.1)"
	cut["SR"] = cut["base"]+" && (mvaScore>={})".format(mvaThreshold) # TODO: Use TCut 
	cut["CR"] = cut["base"]+" && (mvaScore<{})&&(mvaScore>{})".format(mvaThreshold, mvaLowThreshold) #cut["base"]+" && "+"(mvaScore>={})".format(mvaThreshold)
	cut["SB"] = cut["base"]+" && (mvaScore<{})&&(mvaScore>{})".format(mvaLowThreshold, mvaLowerBound)

	
	norm = collections.defaultdict(dict)
	norm["Sig"]["SR"] = 638. # TODO: load from json
	norm["Sig"]["CR"] = 300.
	norm["Sig"]["SB"] = 62.8

	norm["SigPart"]["SR"] = 638.
	norm["SigPart"]["CR"] = 300.
	norm["SigPart"]["SB"] = 62.8

	norm["BkgDstarDs"]["SR"] = 615.
	norm["BkgDstarDs"]["CR"] = 607.
	norm["BkgDstarDs"]["SB"] = 186.

	norm["BkgDstarDsstar"]["SR"] = 1.1
	norm["BkgDstarDsstar"]["CR"] = 1.1
	norm["BkgDstarDsstar"]["SB"] = 1.1

	norm["BkgDstar3pi"]["SR"] = 1430.
	norm["BkgDstar3pi"]["CR"] = 800.
	norm["BkgDstar3pi"]["SB"] = 188.

	norm["dataD2"]["SR"] = 1.
	norm["dataD2"]["CR"] = 1.
	norm["dataD2"]["SB"] = 1.

	norm["dataD2WS"]["SR"] = 549.
	norm["dataD2WS"]["CR"] = 1772.
	norm["dataD2WS"]["SB"] = 1263.

	partFraction = norm["SigPart"]["CR"]/norm["Sig"]["CR"]


	print(cut["SR"])


	for item in filesUsed:  
		samples[item] = ROOT.RDataFrame(filemanager.GetItem(item))
		for region in regions: 
			frames[item][region] = samples[item].Filter(cut[region])
			ROOT.SetOwnership(frames[item][region], 0)
			# For histogram legacy compatibility
			histos[item][region] = frames[item][region].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
			histosunrolled[item][region] = UnrollHist(histos[item][region])

	print(frames)

	#print "Number of events (SR, CR, SB): {}, {}, {}".format(histos["data"]["SR"].GetEntries(), histos["data"]["CR"].GetEntries(), histos["data"]["SB"].GetEntries())

	#BackgroundShapeUnrolled(histos["dataD2"]["SR"], histos["dataD2"]["CR"], [histos["Sig"]["SR"], histos["BkgDstarDs"]["SR"]])

	GetEfficiencies(frames)

	#BackgroundShapeUnrolled(histosunrolled["dataD2"]["SR"], histosunrolled["dataD2"]["CR"], [histosunrolled["Sig"]["SR"], histosunrolled["BkgDstarDs"]["SR"]])



	variable = "b_tau_rhomass1"
	model = ("model", "", 20, 0.2, 1.6)
	datadesc = "dataD2"

	data = frames["dataD2"]["SR"].Histo1D(model, variable)
	bkg = frames["BkgDstarDs"]["SR"].Histo1D(model, variable)

	AtomicDraw(data, outputfolder+"data_before.pdf")

	background = data.Clone("backgroundModel") # Works (does not change initial histo)
	background.Add(bkg.GetPtr(), -1.)

	AtomicDraw(data, outputfolder+"data_after.pdf")

	AtomicDraw(background, outputfolder+"BackgroundModel.pdf")


	MC = copy.deepcopy(filesUsed)
	MC.remove(datadesc)
	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	count = 0
	maxes = []
	refhist = 0
	# Plot the different regions
	"""
	for region in regions: 
		canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
		legend = ROOT.TLegend(	canvas.GetLeftMargin()+0.35, 
								1-canvas.GetTopMargin()-.2, 
								canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
								1-canvas.GetTopMargin())
		legend.SetFillStyle(1)
		legend.SetBorderSize(1)
		legend.SetMargin(0.3)
		legend.SetTextSize(0.04)
		hists = {}
		reference = frames[datadesc][region].Histo1D(model, variable)
		reference.SetLineColor(ROOT.kBlack)
		reference.SetMarkerColor(ROOT.kBlack)
		reference.SetLineWidth(2)
		reference.SetMarkerStyle(8)
		datahist = reference.DrawCopy("E")
		normalisation = reference.Integral()
		legend.AddEntry(datahist,"data","P")
		for item in MC: 
			hist = frames[item][region].Histo1D(model, variable)
			hist.Scale(norm[item][region]/hist.Integral())
			hist.SetLineColor(colors[item])
			hist.SetMarkerColor(colors[item])
			hist.SetLineWidth(2)
			histo = hist.DrawCopy("HIST SAME")
			hists[item] = hist

			legend.AddEntry(histo, legends[item], "L")


		legend.Draw()

		canvas.Draw()
		canvas.Print(outputfolder+"rhomass1"+region+".pdf")


		canv.cd()
		background = reference.Clone("backgroundModel{}".format(count)) # Works (does not change initial histo)
		bkg = frames["BkgDstar3pi"][region].Histo1D(model, variable)
		normfactor = norm["BkgDstar3pi"][region]/bkg.Integral()
		background.Add(bkg.GetPtr(), -1.*normfactor)
		background.SetLineColor(ROOT.kOrange+1+count)
		currenthist = background.DrawCopy("HIST SAME")
		if (count == 0): 
			refhist = currenthist
		maxes.append(background.GetMaximum())
		count+=1

	refhist.SetMaximum(max(maxes)*1.3)
	canv.Draw()
	canv.Print(outputfolder+"BackgroundModel.pdf")
	"""


	SuperimposeRegions(frames, "BkgDstarDs", variable, outputfolder+"DstarDsShapesTest.pdf", model)

	count = 0
	shapes = {}
	maxes=[]
	drawn = {}
	for region in regions: 
		background = frames["dataD2"][region].Histo1D(model, variable) #.Clone("backgroundModel{}".format(count)) # Works (does not change initial histo)
		background.Sumw2()
		bkg = frames["BkgDstarDs"][region].Histo1D(model, variable)
		normfactor = norm["BkgDstarDs"][region]/bkg.Integral()
		background.Add(bkg.GetPtr(), -1.*normfactor)
		background.SetLineColor(ROOT.kOrange+1+count)
		shapes[region] = background
		background.Scale(1./background.Integral())
		background.SetLineWidth(2)
		maxes.append(background.GetMaximum())
		count+=1

	newcanvas = ROOT.TCanvas("newcanvas", "newcanvas", 800, 600)
	legend = ROOT.TLegend(	canvas.GetLeftMargin()+0.35, 
								1-canvas.GetTopMargin()-.2, 
								canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
								1-canvas.GetTopMargin())
	drawn["SR"] = shapes["SR"].DrawCopy("HIST E")
	drawn["CR"] = shapes["CR"].DrawCopy("HIST E SAME")
	drawn["SB"] = shapes["SB"].DrawCopy("HIST E SAME")
	drawn["SR"].SetMaximum(max(maxes)*1.3)
	legend.AddEntry(drawn["SR"], "SR", "L")
	legend.AddEntry(drawn["CR"], "CR", "L")
	legend.AddEntry(drawn["SB"], "SB", "L")
	legend.Draw()
	newcanvas.Draw()
	newcanvas.Print(outputfolder+"BackgroundModel.pdf")


	count = 0
	shapes = {}
	maxes=[]
	drawn = {}
	for region in regions: 
		shapes[region] = frames["BkgDstarDs"][region].Histo1D(model, variable)
		shapes[region].Sumw2()
		shapes[region].SetLineColor(ROOT.kAzure+1+count)
		shapes[region].SetLineWidth(2)
		shapes[region].Scale(1./shapes[region].Integral())
		maxes.append(shapes[region].GetMaximum())
		count+=1

	newcanvas = ROOT.TCanvas("newcanvas", "newcanvas", 800, 600)
	legend = ROOT.TLegend(	canvas.GetLeftMargin()+0.35, 
								1-canvas.GetTopMargin()-.2, 
								canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
								1-canvas.GetTopMargin())
	drawn["SR"] = shapes["SR"].DrawCopy("HIST E")
	drawn["CR"] = shapes["CR"].DrawCopy("HIST E SAME")
	drawn["SB"] = shapes["SB"].DrawCopy("HIST E SAME")
	drawn["SR"].SetMaximum(max(maxes)*1.3)
	legend.AddEntry(drawn["SR"], "SR", "L")
	legend.AddEntry(drawn["CR"], "CR", "L")
	legend.AddEntry(drawn["SB"], "SB", "L")
	legend.Draw()
	newcanvas.Draw()
	newcanvas.Print(outputfolder+"DstarDsShapes.pdf")


	# Making datacards 
	"""
	with open("datacard.txt", "w") as datacard: 
		datacard.write("# Datacard generated automatically with {}{} on {}.\n".format(os.getcwd(), __file__, datetime.today().strftime("%d.%m.%y %H:%M:%S")))
		datacard.write("# Rhomass fit with DstarDs component\n\n")

		innerlengths = []
		for item, element in frames.items(): 
			innerlengths.append(len(element))
		print(innerlengths)
		assert(len(innerlengths)>=1), "ERROR: No region/channel defined."
		assert(all(element == innerlengths[0] for element in innerlengths)), "ERROR: Collection with different numbers of regions provided. "
		imax = innerlengths[0]
		datacard.write("imax {}\n".format(imax))
		datacard.write("jmax {}\n".format(len(filesUsed)-2))
		datacard.write("kmax {}\n".format(0))

		datacard.write("\n#Observed events (data)\n")
		regionstring = "bin "
		for item in regions: 
			regionstring += (item+" ")
		regionstring+="\n"
		datacard.write(regionstring)

		observationstring = "observation "
		for item in frames["dataD2"]: 
			observationstring += ("{} ".format(frames["dataD2"]["CR"].Count().GetValue()))
		datacard.write(observationstring+"\n")

		datacard.write("\n#Expected events (MC/model)\n")
		binstring = "bins "
		labelstring = "process "
		indexstring = "process "
		expectedstring = "rate "
		count = 1
		for region in regions: 
			for item in MC: 
				binstring += "{} ".format(region)
				labelstring += "{} ".format(item)
				indexstring += "{} ".format(count)
				expectedstring += "{} ".format(norm[item][region])
				count += 1
		datacard.write(binstring+"\n")
		datacard.write(labelstring+"\n")
		datacard.write(indexstring+"\n")
		datacard.write(expectedstring+"\n")

		datacard.write("\n#Shapes and RooFit workspace\n")
		workspacefile = "workspace.root"
		workspacename = "workspace"
		file = ROOT.TFile.Open(workspacefile, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		for region in regions: 
			for item in MC: 
				histname = item+"_"+region
				datacard.write("shapes {} {} {} {}\n".format(item, region, workspacefile, histname))
				hist = frames[item][region].Histo1D(model, variable)
				hist.SetName(histname)
				hist.Write()
		file.Write()
		file.Close()
	"""


	# Deriving the background shape in the SB 
	region = "SB" # we work in the sideband for now
	MC = ["SigPart", "dataD2WS"] #"Sig", 
	BKG = ["dataD2WS"]
	background = histosunrolled["dataD2"][region].Clone("backgroundModel") # Works (does not change initial histo)
	background.Sumw2()
	for item in BKG: 
		hist = histosunrolled[item][region]
		hist.Sumw2()
		normfactor = norm[item][region]/hist.Integral()
		background.Add(hist, -1.*normfactor)

	histograms = collections.defaultdict(dict)
	for region in regions: 
		for item in MC+["dataD2"]: 
			hist = histosunrolled[item][region] #ROOT.convertHisto(frames[item][region].Histo1D(model, variable).GetPtr())
			hist.Sumw2()
			print(hist.Integral())
			histograms[item][region] = hist

	for region in regions: 
		histograms["bkg"][region] = background.Clone()

	# Now we use it to fit the data in the CR
	regions = ["CR", "SB"]
	with open("datacard.txt", "w") as datacard: 
		datacard.write("# Datacard generated automatically with {}{} on {}.\n".format(os.getcwd(), __file__, datetime.today().strftime("%d.%m.%y %H:%M:%S")))
		datacard.write("# Rhomass fit in CR with DstarDs component and bacgkround model from SB\n\n")

		datacard.write("imax {}\n".format(len(regions)))
		datacard.write("jmax {}\n".format(len(BKG)))
		datacard.write("kmax {}\n".format(0)) # For now no systematics

		datacard.write("\n"+"-"*50+"\n")
		datacard.write("# Shapes and RooFit workspace\n")
		workspacefile = "workspace.root"
		workspacename = "w"
		file = ROOT.TFile.Open(workspacefile, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		# Creating the variable on which we fit
		var = RooRealVar("rhomass1", "rhomass1", 0., 100.)
		fitspace = RooArgSet(var)
		datacard.write("shapes data_obs {} {} {}\n".format("CR", workspacefile, workspacename+":data_obs_CR"))
		datacard.write("shapes data_obs {} {} {}\n".format("SB", workspacefile, workspacename+":data_obs_SB"))
		histograms["dataD2"]["CR"].SetName("data_obs_CR")
		dataCR = RooDataHist("data_obs_CR", "data_obs_CR", fitspace, histograms["dataD2"]["CR"])
		histograms["dataD2"]["CR"].Write()
		getattr(workspace, "import")(dataCR)
		histograms["dataD2"]["SB"].SetName("data_obs_SB")
		dataSB = RooDataHist("data_obs_SB", "data_obs_SB", fitspace, histograms["dataD2"]["SB"])
		histograms["dataD2"]["SB"].Write()
		getattr(workspace, "import")(dataSB)
		histograms["dataD2"]["SB"].Write()
		print(histograms["dataD2"]["CR"].Integral())
		localnorm = copy.deepcopy(norm)
		localnorm["BkgDstarDs"]["CR"] = 1. #"Ds_norm"
		localnorm["bkg"]["CR"] = 1. #"bkg_norm"
		localnorm["bkg"]["SB"] = 1.
		localnorm["BkgDstarDs"]["SB"] = 1.
		for region in regions: 
			for item in MC: 
				histname = item+"_"+region
				datacard.write("shapes {} {} {} {}\n".format(item, region, workspacefile, workspacename+":"+histname))
				hist = histograms[item][region] #TODO: fix availablility of histos
				hist.Scale(localnorm[item][region]/hist.Integral())
				hist.SetName(histname)
				roohist = ROOT.RooDataHist(histname, histname, fitspace, hist)
				hist.Write()
				getattr(workspace, "import")(roohist)
		bins = RooArgList()
		binsdest = RooArgList()
		variables = []
		variablesdest = []
		maxval = histograms["dataD2"][region].GetMaximum()
		print(background.GetNbinsX())
		"""
		transferfactor = RooRealVar("bkg_transferfactor_SB_CR", "bkg_transferfactor_SB_CR", 0., 10.)
		for i in range(background.GetNbinsX()): 
			# Creating bins for the shape in SB
			name = "bin_SB_{}".format(i)
			print("bin {} {} {}".format(i, min(1., background.GetBinContent(i)), name))
			mybin = RooRealVar(name, name, min(1.,background.GetBinContent(i)), 1., maxval)
			bins.add(mybin)
			variables.append(mybin) # needed as the ArgList contains references only 
			# Creating the destination bins as RooFormulaVar of transferfactor * source bin
			name = name.replace("SB", "CR")
			destbin = RooFormulaVar(name, name, "@0*@1", RooArgList(transferfactor, mybin))
			binsdest.add(destbin)
			variablesdest.append(destbin)
		parametricshape = ROOT.RooParametricHist("bkg_SB", "bkg_SB", var, bins, background)
		getattr(workspace, "import")(parametricshape)
		bkgNormSB = RooAddition("bkg_SB_norm", "bkg_SB_norm", bins)
		getattr(workspace, "import")(bkgNormSB)
		parametricshapefloat =  ROOT.RooParametricHist("bkg_CR", "bkg_CR", var, binsdest, background)
		getattr(workspace, "import")(parametricshapefloat)
		bkgNormCR = RooAddition("bkg_CR_norm", "bkg_CR_norm", binsdest)
		getattr(workspace, "import")(bkgNormCR, ROOT.RooFit.RecycleConflictNodes())
		datacard.write("shapes {} {} {} {}\n".format("bkg", "CR", workspacefile, workspacename+":"+"bkg_CR"))
		datacard.write("shapes {} {} {} {}\n".format("bkg", "SB", workspacefile, workspacename+":"+"bkg_SB"))
		"""
		workspace.Write()
		file.Write()
		file.Close()

		#MC.append("bkg")
		datacard.write("\n"+"-"*50+"\n")
		datacard.write("# Observed events (data)\n")
		regionstring = "bin "
		for item in regions: 
			regionstring += (item+" ")
		regionstring+="\n"
		datacard.write(regionstring)

		observationstring = "observation "
		for item in regions: 
			print(frames["dataD2"][item].Count().GetValue())
			observationstring += ("-1 ") # "{} ".format(frames["dataD2"][item].Count().GetValue()) # TODO: fix
		datacard.write(observationstring+"\n")

		# We want to leave a few components floating 
		datacard.write("\n"+"-"*50+"\n")
		datacard.write("# Expected events (MC/model)\n")
		binstring = "bin "
		labelstring = "process "
		indexstring = "process "
		expectedstring = "rate "
		count = 0
		for region in regions: 
			for item in MC: 
				binstring += "{} ".format(region)
				labelstring += "{} ".format(item)
				factor = 1
				if "Sig" in item: 
					factor = -1 # make signal negative
				indexstring += "{} ".format(factor*count)
				expectedstring += "{} ".format(localnorm[item][region])
				count += 1
		datacard.write(binstring+"\n")
		datacard.write(labelstring+"\n")
		datacard.write(indexstring+"\n")
		datacard.write(expectedstring+"\n")

		#datacard.write("\n"+"-"*50+"\n")
		#datacard.write("lumi     lnN    1.10       1.0 		1.0\n")

		datacard.write("\n"+"-"*50+"\n")
		#datacard.write("BkgDstarDs_CR_norm rateParam CR BkgDstarDs {} [{},{}]\n".format(norm["BkgDstarDs"]["CR"], 0, norm["BkgDstarDs"]["CR"]*5.))
		#datacard.write("BkgDstarDs_SB_norm rateParam SB BkgDstarDs {} [{},{}]\n".format(norm["BkgDstarDs"]["SB"], 0, norm["BkgDstarDs"]["SB"]*5.))
		#datacard.write("bkg_SB_norm rateParam CR bkg {} [{},{}]\n".format(norm["dataD2"]["CR"]/2., 0., norm["dataD2"]["CR"]))
		#datacard.write("bkg_transferfactor_SB_CR rateParam SB bkg 0.5 [0.0,10]\n")
		#datacard.write("bkg_CR_norm rateParam SB bkg (@0*@1) bkg_SB_norm,bkg_transferfactor_SB_CR\n")
		#datacard.write("Sig_CR_norm rateParam CR Sig {} [{},{}]\n".format(norm["Sig"]["CR"], 0, norm["Sig"]["CR"]*5.))
		#datacard.write("SigPart_CR_fraction rateParam CR SigPart {}\n".format(partFraction))
		#datacard.write("SigPart_CR_norm rateParam CR SigPart (@0*@1) r,SigPart_CR_fraction\n")
		datacard.write("dataD2WS_CR_norm rateParam CR dataD2WS {} [{},{}]\n".format(norm["dataD2WS"]["CR"], 0., norm["dataD2WS"]["CR"]*2.))
		datacard.write("dataD2WS_SB_norm rateParam SB dataD2WS {} [{},{}]\n".format(norm["dataD2WS"]["SB"], 0., norm["dataD2WS"]["SB"]*2.))
		datacard.write("\n"+"-"*50+"\n")
		for item in variables: 
			datacard.write("{} flatParam\n".format(item.GetName()))


	

	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	hist = frames["Sig"]["SR"].Histo2D(("rhomass1", "rhomass2", 15, rangeMin, rangeMax, 15, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
	hist.SetTitle("m(#rho_{1}) vs m(#rho_{2})")
	hist.DrawCopy("COLZ")
	canv.Draw()
	canv.Print(outputfolder+"2Drhomass.pdf")









filemanager.CloseAll()



