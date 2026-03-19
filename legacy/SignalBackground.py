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
ROOT.gROOT.LoadMacro("FileFlow.h")
ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
#ROOT.setTDRStyle()
#import CMS_lumi
from ROOT import Ana



webpublication =False



def LoadFiles(filesUsed, cycle="v1"): # TODO: add into common python include 
	Ana.Init(cycle)
	for item in filesUsed: 
		print("Opening file: {}".format(item))
		Ana.filemanager.OpenItem(item); 
	return Ana.filemanager

#presentationfolder = "../Presentations/Presentation_22_5_10/reveal.js-master/Figures.root"



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

	name = "unrolled"+histo2D.GetName()
	unrolled = ROOT.TH1D(name, name, nTotal, 0, 100)

	print("Nunmber of bins: {}, {}".format(nx, ny))

	for i in range(0, nx):
		for j in range(0, ny): # TODO: check overflow is handled properly 
			if inverted: 
				binContent = histo2D.GetBinContent(j, i)
			else: 
				binContent = histo2D.GetBinContent(i, j)
			unrolled.SetBinContent(i+j*nx, binContent)

	return unrolled

def BackgroundShapeUnrolled(reference, estimate, additionalhists, filename): # Works as in PlotDistributions
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

	canvas.Lower()

	histoRatio.Draw("HIST")
	histoRatio.SetLineColor(1)
	histoRatio.DrawCopy("LE SAME")


	canvas.Draw()

	#histoRatio.SetMaximum(0.2)
	#histoRatio.SetMinimum(1.8)
	histoRatio.GetYaxis().SetRangeUser(0.2, 1.8)
	histoRatio.GetYaxis().SetNdivisions(5)
	canvas.RemoveMiddleAxis()

	canvas.Update()

		#FOM = filemanager.GetItem("cutflowGen")

	fomcanvas = ROOT.TCanvas("fomcanvas", "fomcanvas", 800, 600)

	#canv = ROOT.TCanvas("canv", "canv", 800, 600)
	#hist = SR["MC"].Histo2D(("rhomass1", "rhomass2", 15, rangeMin, rangeMax, 15, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
	#hist.SetTitle("m(#rho_{1}) vs m(#rho_{2})")
	#hist.DrawCopy("COLZ")
	#canv.Draw()
	#canv.Print(outputfolder+"2Drhomass.pdf")

	

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

	canvas.Print(filename+".pdf")
	if (webpublication): canvas.Print(filename+".png")


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
			n = frames[item][key].Count().GetValue()
			N = frames[item]["all"].Count().GetValue()
			eff = n/N
			print(("\t{}: {} ({}/{})".format(item, eff, n, N)))
		print("\n")


def getEff(n, N): 
	eff = float(n)/float(N)
	#print eff
	err = sqrt(abs(eff*(1.-eff)/float(N)))
	#print err
	#return eff, err
	return ufloat(eff, err)


def ComputeEfficiencies(frames, baseline="baseline"): 
	efficiencies = collections.defaultdict(dict)
	for item, content in frames.iteritems(): 
		#print("{}:".format(item))
		for key, value in content.iteritems(): 
			n = frames[item][key].Count().GetValue()
			N = frames[item][baseline].Count().GetValue()
			# see: chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://indico.cern.ch/event/66256/contributions/2071577/attachments/1017176/1447814/EfficiencyErrors.pdf
			ne = ufloat(n, sqrt(n))
			NE = ufloat(N, sqrt(N))
			efferr = ne/NE
			#print(efferr)
			#print(("\t{}: {} ({}/{})".format(key, eff, n, N)))
			efficiencies[item][key] = getEff(n, N) #ufloat(eff, err)
	return efficiencies


def PrintEfficiencies(effs): 
	for item, content in effs.iteritems(): 
		print("{}:".format(item))
		for key, value in content.iteritems(): 
			print(("\t{}: {}".format(key, value)))


def DumpEffs(effs, path): 
	effsForWrite = collections.defaultdict(dict)
	for item, content in effs.iteritems(): 
		for key, value in content.iteritems(): 
			eff = effs[item][key]
			effsForWrite[item][key] = (eff.n, eff.s)
	with open(path, "w") as file: 
		json.dump(effsForWrite, file, ensure_ascii=False, encoding="utf8", sort_keys=False)


def ReadEffs(path): 
	effs = collections.defaultdict(dict)
	with open(path, "r") as file: 
		effsFromFile = json.load(file, encoding="utf8")
		for item, content in effsFromFile.iteritems(): 
			for key, value in content.iteritems(): 
				eff = effsFromFile[item][key]
				assert(len(eff)==2)
				effs[item][key] = ufloat(eff[0], eff[1])
	return effs


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
	for key, eff in effs.iteritems(): 
		key.replace("Part", "")
	 	expected[key]= lumi*bbxsec*fB0*2.*Br_Dstar_D0pi*Br_D0_KPI*1000.*eff["br"]*eff["geneff"] #*eff["eff"]

	return expected


def getEffFromInfo(tree): 
	frame = RDataFrame(tree)

	n = frame.Sum("numSelected").GetValue()

	N = frame.Sum("numTotal").GetValue()

	eff = getEff(n, N)
	return eff

def CompleteEffsFromFile(effs, version, filemanager): 
	anaeffs = {"Sig":ufloat(1.4e-3, 0.), "BkgDstarDs":ufloat(1.44e-3, 0.), "BkgDstarDsstar":ufloat(2.26e-3, 0.), "BkgDstar3pi":ufloat(5.3e-4, 0.), "SigPart":ufloat(1.27e-3, 0.), "dataD2WS":ufloat(-0.392, 0.), "dataD21TauWS": ufloat(-0.53, 0.)}
	for key, eff in effs.iteritems(): 
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


def MultiplyFinalEffs(effs, regioneffs):
	resulteffs = collections.defaultdict(dict)
	for item, content in regioneffs.iteritems(): 
		#print("{}:".format(item))
		for key, value in content.iteritems(): 
			try:
				resulteffs[item][key] = effs[item]*regioneffs[item][key] #*regioneffs[item]["baseline"]
				if (effs[item].nominal_value < 0.): 
					resulteffs[item][key]=effs[item]
			except:
				resulteffs[item][key] = ufloat(-1., 0.)
	return resulteffs


def MakeDatacardSimple(frames, yields, name=""): 
	region = "SB"
	control = "CR"

def InitialiseWebfolder(webfolder): 
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


def PublishToWeb(folder, webfolder): 
	webbasepath = "/eos/home-m/mhuwiler/www/Analysis/"
	publicationfolder = "{}{}/".format(webbasepath, webfolder)
	print(publicationfolder)
	InitialiseWebfolder(publicationfolder)
	os.system("cp -r {} {}".format(folder, publicationfolder))
	publicationurl = "https://mhuwiler.web.cern.ch/Analysis/{}/?match=&depth=3".format(webfolder)
	print("Published plots in '{}' to: \n{}".format(folder, publicationurl))
	


# Web publication
if (webpublication): 
	webfolder = "/eos/home-m/mhuwiler/www/Analysis/BackgroundModellingUpdate/"
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

def AtomicDraw(histo, name, options = ""): 
	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	histo.DrawCopy(options)
	canv.Draw()
	canv.Print(name)


def PlotOverlay(frames, dataname, initialcomponents, regions, variables, yields, outfolder, drawlegend=True, normalise=False): 
	# Plotting distributions over each other 
	outfolder+="overlay/"
	os.system("mkdir -p "+outfolder)
	factor = 1.1 # how much overhead to add to the histos 
	components = copy.deepcopy(initialcomponents)
	components.reverse()
	if (dataname in components): components.remove(dataname)
	for region in regions: 
		for variable in variables: 
			name = "{}_{}".format(variable, region)
			canvas = TCanvas(name, "{} in {}".format(variable, region), 800, 600)

			legend = TLegend(canvas.GetLeftMargin()+0.35, 
	                         	1.-canvas.GetTopMargin()-.2, 
	                            canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
	                           	1.-canvas.GetTopMargin() )

			examplehist = Ana.binning[variable]
			data = frames[dataname][region].Histo1D(examplehist, variable)
			data.SetMarkerStyle(8) # Large scalable dot
			data.SetMarkerSize(0.5)
			data.SetLineColor(ROOT.kBlack)
			data.SetTitle("{}_{}".format(variable, region)) # TODO: Delete once the binning is centralised
			#data.SetFillColor(ROOT.kBlack)
			legend.AddEntry(data.GetPtr(), "data", "PE")
			data.Draw("E")

			maxes = [data.GetMaximum()]

			i = 0
			for component in components: 
				histo = frames[component][region].Histo1D(examplehist, variable)
				#ROOT.SetOwnership(histo, 0)
				if (histo.Integral()==0): 
					continue
				print(component)
				histo.SetLineStyle(1) # plain
				histo.SetLineWidth(2)
				color = Ana.samples.at(component.replace("Part", "")).color
				if (not color): 
					color = colors[i]
				histo.SetLineColor(color) #colors[i]Ana.color[component.replace("Part", "")]
				#histo.SetFillStyle(3003)
				histo.SetFillColorAlpha(color, 0.4)
				#histo = frames[component][region].Histo1D(examplehist, variable).GetPtr()
				if normalise: 
					if (yields[component][region] > 0. and histo.Integral() > 0): 
						histo.Scale(yields[component][region].nominal_value/histo.Integral())
					elif (yields[component][region].nominal_value != -1.):
						histo.Scale(-yields[component][region].nominal_value*histo.Integral())
				else:
					if not (histo.Integral()==0): 
						histo.Scale(data.Integral()/histo.Integral())
				#histo.SetMarkerColor(Ana.color[component])
				histo.DrawCopy("HIST SAME")
				maxes.append(histo.GetMaximum())
				legend.AddEntry(histo.GetPtr(), component)
				i+=1

			if (drawlegend): legend.Draw()
			legend.SetBorderSize(1)
			legend.SetMargin(0.3)
			legend.SetTextSize(0.04)

			if (not normalise): data.SetMaximum(factor*max(maxes))
			canvas.Draw()

			canvas.Print(outfolder+name+".png")
			canvas.Print(outfolder+name+".pdf")


def PlotStack(frames, dataname, initialcomponents, regions, variables, yields, outfolder, drawlegend=True): 
	# Plotting distributions over each other 
	outfolder+="stacked/"
	os.system("mkdir -p "+outfolder)
	factor = 1.1 # how much overhead to add to the histos 
	components = copy.deepcopy(initialcomponents)
	if (dataname in components): components.remove(dataname)
	components.reverse()
	notYetDrawn = True
	for region in regions: 
		for variable in variables: 
			name = "{}_{}".format(variable, region)
			canvas = TCanvas(name, "{} in {}".format(variable, region), 800, 600)

			legend = TLegend(canvas.GetLeftMargin()+0.35, 
	                         	1.-canvas.GetTopMargin()-.2, 
	                            canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
	                           	1.-canvas.GetTopMargin() )

			examplehist = Ana.binning[variable]
			data = frames[dataname][region].Histo1D(examplehist, variable)
			data.SetMarkerStyle(8) # Large scalable dot
			data.SetMarkerSize(0.5)
			data.SetLineColor(ROOT.kBlack)
			data.SetTitle("") #data.SetTitle("{}_{}".format(variable, region))
			#data.SetFillColor(ROOT.kBlack)
			legend.AddEntry(data.GetPtr(), "data", "PE")
			data.Draw("E")

			stack = THStack("stack", "Background modelling")
			hists = {}
			for component in components: 
				if "-" in component: 
					comps = component.split("-")
					assert(len(comps)>=2)
					print(component)
					histo = copy.deepcopy(frames[comps[0]][region].Histo1D(examplehist, variable).GetPtr())
					if (yields[comps[0]][region] > 0. ): 
						histo.Scale(yields[comps[0]][region].nominal_value/histo.Integral())
					elif (yields[comps[0]][region].nominal_value != -1.):
						histo.Scale(-yields[comps[0]][region].nominal_value*histo.Integral())
					comps.remove(comps[0])
					for comp in comps: 
						hist = copy.deepcopy(frames[comp][region].Histo1D(examplehist, variable).GetPtr())
						if (yields[comp][region] > 0. ): 
							hist.Scale(yields[comp][region].nominal_value/hist.Integral())
						elif (yields[comp][region].nominal_value != -1.):
							hist.Scale(-yields[comp][region].nominal_value*hist.Integral())
						histo.Add(hist, -1.)
				else: 
					histo = frames[component][region].Histo1D(examplehist, variable).GetPtr()
					if (yields[component][region] > 0. and histo.Integral() > 0): 
						histo.Scale(yields[component][region].nominal_value/histo.Integral())
					elif (yields[component][region].nominal_value != -1.):
						histo.Scale(-yields[component][region].nominal_value*histo.Integral())
				ROOT.SetOwnership(histo, 0)
				histo.SetLineStyle(1) # plain
				histo.SetLineWidth(2)
				color = Ana.samples.at(component.replace("Part", "")).color
				if (not color): 
					color = colors[0]
				histo.SetLineColor(color)
				#histo.SetMarkerColor(Ana.color[component])
				histo.SetFillStyle(1)
				histo.SetFillColor(color)
				hists[component] = histo
				stack.Add(histo)
				#legend.AddEntry(histo.GetPtr(), Ana.legends[component], "F")

			stack.Draw("HIST SAME") #"SAME"
			data.Draw("E SAME") # Plot on top
			if (drawlegend): legend.Draw()
			legend.SetBorderSize(1)
			legend.SetMargin(0.3)
			legend.SetTextSize(0.04)

			#data.GetXaxis().SetTitle(Ana.labels[variable])
			data.GetXaxis().SetTitleSize(0.06)
			data.GetXaxis().SetLabelSize(0.06)
			data.GetYaxis().SetLabelSize(0.06)
			#data.GetYaxis().SetTitle("Counts")
			data.GetYaxis().SetTitleSize(0.06)
			data.GetXaxis().SetTitleOffset(1.2)
			canvas.SetBottomMargin(0.15)
			canvas.SetTopMargin(0.01)
			canvas.SetLeftMargin(0.15)

			maxes = [data.GetMaximum(), stack.GetMaximum()]

			data.SetMaximum(factor*max(maxes))
			canvas.Draw()

			canvas.Print(outfolder+name+".png")
			canvas.Print(outfolder+name+".pdf")

			if ((not drawlegend) and notYetDrawn): 
				#components.reverse()
				canv = TCanvas("legendCanvas", "legenCanvas", 800, 200*len(hists))
				dummy = TCanvas("dummy", "dummy", 800, 600)
				#legend.AddEntry(data.GetPtr(), "data", "PE")
				for component, histo in hists.iteritems(): 
					#histo = frames[component][region].Histo1D(examplehist, variable)
					ROOT.SetOwnership(histo, 0)
					histo.SetLineStyle(1) # plain
					histo.SetLineWidth(2)
					color = Ana.samples.at(component.replace("Part", "")).color
					if "WS" in component: 
						color = Ana.samples.at("WS").color
					histo.SetLineColor(color)
					histo.SetFillStyle(1)
					histo.SetFillColor(color)
					legend.AddEntry(histo, Ana.samples.at(component.replace("Part", "")).legend, "F")
				canv.cd()
				data.SetMarkerSize(4.)
				data.SetLineWidth(4)
				legend.SetX1(0.)
				legend.SetY1(0.)
				legend.SetX2(1.)
				legend.SetY2(1.)
				legend.SetBorderSize(0)
  				legend.SetFillColor(0)
  				legend.SetFillStyle(0)
				legend.SetTextFont(43)
				legend.SetTextSize(canv.GetWh()/(2*stack.GetNhists()))
				legend.Draw()
				canv.Draw()
				canv.Print(outfolder+"legend.png")
				canv.Print(outfolder+"legend.pdf")
				notYetDrawn = False


def PlotComparison(frames, referencename, comparisonname, regions, variables, outfolder, normalise=False): 
	# Plotting distributions over each other 
	outfolder+="comparisons/"
	os.system("mkdir -p "+outfolder)
	factor = 1.3 # how much overhead to add to the histos 
	for region in regions: 
		for variable in variables: 
			name = "{}_{}".format(variable, region)
			canvas = TCanvas(name, "{} in {}".format(variable, region), 800, 600)

			legend = TLegend(canvas.GetLeftMargin()+0.55, 
	                         	1.-canvas.GetTopMargin()-.11, 
	                            canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
	                           	1.-canvas.GetTopMargin() )

			examplehist = Ana.binning[variable]
			reference = frames[referencename][region].Histo1D(examplehist, variable)
			reference.SetTitle("{}_{}".format(variable, region))
			#reference.SetMarkerStyle(8) # Large scalable dot
			#reference.SetMarkerSize(0.5)
			reference.SetLineWidth(2)
			reference.SetLineColor(ROOT.kBlue)
			reference.SetFillStyle(3003)
			reference.SetFillColor(ROOT.kBlack)
			reference.SetTitle("{}_{}".format(variable, region))
			legend.AddEntry(reference.GetPtr(), referencename, "L")
			reference.Draw("HIST E")

			comparison = frames[comparisonname][region].Histo1D(examplehist, variable)
			comparison.SetLineWidth(2)
			comparison.SetLineColor(ROOT.kRed)
			comparison.SetFillStyle(3356)
			comparison.SetFillColor(ROOT.kRed)
			legend.AddEntry(comparison.GetPtr(), comparisonname, "L")
			if normalise: 
				if "data" in comparisonname: 
					comparison.Scale(dirtynorm[comparisonname][region])
				else: 
					comparison.Scale(dirtynorm[comparisonname][region]/comparison.Integral())
			else: 
				comparison.Scale(reference.Integral()/comparison.Integral())

			comparison.Draw("HIST SAME E") #"SAME"
			legend.Draw()
			legend.SetBorderSize(1)
			legend.SetMargin(0.3)
			legend.SetTextSize(0.04)

			maxes = [reference.GetMaximum(), comparison.GetMaximum()]

			reference.SetMaximum(factor*max(maxes))
			canvas.Draw()

			canvas.Print(outfolder+name+".png")
			canvas.Print(outfolder+name+".pdf")


from ROOT import TFile
def LoadFile(filename, treename="ntuplizer/tree"): 
	globals()[filename] = TFile.Open(filename, "READ")
	tree = globals()[filename].Get(treename)
	return tree


def PrepareCustomFiles(dictf, regions): 
	frames = collections.defaultdict(dict)
	histos = collections.defaultdict(dict)
	histosunrolled = collections.defaultdict(dict)
	for key, item in dictf.iteritems(): 
		frame = RDataFrame(item)
		frames[key]["all"] = frame.Filter("(B_mu_alpha > 1.)") #Ana.cut["base"].GetTitle()
		allregions = copy.deepcopy(regions)
		allregions.append("all")
		for region in regions: 
			cut = "(B_mu_alpha > 1.)" #Ana.cut[region].GetTitle()
			if (options.debug): print("Using following cut string (from TCut): {}".format(cut))
			frames[key][region] = frames[key]["all"].Filter(cut)
			ROOT.SetOwnership(frames[key][region], 0)
			# For histogram legacy compatibility
			#histos[key][region] = frames[key][region].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
			#histosunrolled[key][region] = UnrollHist(histos[key][region])
	return frames, histos, histosunrolled

def WriteDatacard(frames, yields, variables, regions, dataname, datacardname="datacard.txt", workspacefile="workspace.root", workspacename="w"): 
	from datetime import datetime
	from ROOT import RooRealVar, RooArgSet, RooWorkspace, RooDataHist
	MC = frames.keys()
	MC.remove(dataname)
	variable = "b_tau_rhomass1"

	with open(datacardname, "w") as datacard: 
		datacard.write("# Datacard generated automatically with {}{} on {}.\n".format(os.getcwd(), __file__, datetime.today().strftime("%d.%m.%y %H:%M:%S")))
		datacard.write("# Simple fit \n\n")
		datacard.write("imax {}\n".format(len(regions)))
		datacard.write("jmax {}\n".format(len(MC)))
		datacard.write("kmax {}\n".format(0)) # For now no systematics
		datacard.write("\n"+"-"*50+"\n")

		datacard.write("# Shapes and RooFit workspace\n")
		file = ROOT.TFile.Open(workspacefile, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		# Creating the variable on which we fit
		var = RooRealVar(variable, variable, 0., 100.)
		fitspace = RooArgSet(var)
		for region in regions: 
			# Writing the data shapes for each region
			name = "data_obs_{}".format(region)
			datacard.write("shapes data_obs {} {} {}\n".format(region, workspacefile, workspacename+":{}".format(name)))
			examplehist = Ana.binning[variable]
			hist = frames[dataname][region].Histo1D(examplehist, variable).GetPtr()
			hist.SetName(name)
			histogram = RooDataHist(name, name, fitspace, hist)
			hist.Write() # Also saving the ROOT hist
			getattr(workspace, "import")(histogram)
			# Writing the MC shapes 
			for item in MC: 
				histname = item+"_"+region
				datacard.write("shapes {} {} {} {}\n".format(item, region, workspacefile, workspacename+":"+histname))
				hist = frames[item][region].Histo1D(examplehist, variable).GetPtr()
				#hist.Scale(yields[item][region].n/hist.Integral())
				hist.SetName(histname)
				roohist = ROOT.RooDataHist(histname, histname, fitspace, hist)
				hist.Write()
				getattr(workspace, "import")(roohist)

		workspace.Write()
		file.Write()
		file.Close()
		datacard.write("\n"+"-"*50+"\n")
		
		datacard.write("# Observed events (data)\n")
		regionstring = "bin "
		for item in regions: 
			regionstring += (item+" ")
		regionstring+="\n"
		datacard.write(regionstring)

		observationstring = "observation "
		for item in regions: 
			print(frames[dataname][item].Count().GetValue())
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
				expectedstring += "{} ".format(yields[item][region])
				count += 1
		datacard.write(binstring+"\n")
		datacard.write(labelstring+"\n")
		datacard.write(indexstring+"\n")
		datacard.write(expectedstring+"\n")
		#datacard.write("\n"+"-"*50+"\n")
		#datacard.write("lumi     lnN    1.10       1.0 		1.0\n")
		datacard.write("\n"+"-"*50+"\n")
		# Writing out the constraints
		for region in regions: 
			for item in MC: 
				datacard.write("{}_{}_norm rateParam {} {} {} [{},{}]\n".format(item, region, region, item, yields[item][region].n, 0, yields[item][region].n*5.))
		datacard.write("\n"+"-"*50+"\n")
		#for item in variables: 
			#datacard.write("{} flatParam\n".format(item.GetName()))


def ReadEffsSimple(path): 
	effs = {}
	with open(path, "r") as file: 
		effsFromFile = json.load(file, encoding="utf8")
		for item, content in effsFromFile.iteritems(): 
				eff = effsFromFile[item]
				assert(len(eff)==2)
				effs[item] = ufloat(eff[0], eff[1])
	return effs


#ROOT.gInterpreter.Declare("""
#	double Rhomass2DUnrolled(Float_t rhomass1, Float_t rhomass2)
#	{
#		return int((std::min(rhomass2, 1.3) - 0.2)/0.22) + 6*int((std::min(tau_rhomass1, 1.3) - 0.2)/0.22); 
#	}
#""")


if __name__ == "__main__":

	parser = ArgumentParser(description="SignalBackground")
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("--out", dest="out", action="store", type=str, default="BackgroundEstimateTauFLNewTrainingWSlatestBinning/", help="Directory where the plots shuld go")
	parser.add_argument("--name", dest="name", action="store", type=str, default="test", help="Turn on debug output")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v1", help="Which version (cycle) of files to run on")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("-f", "--forcepath", dest="forcepath", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
	parser.add_argument("--stats", dest="stats", action="store_true", default=False, help="Show stats box in ROOT")

	options = parser.parse_args()

	
	if (options.batch): 
		ROOT.gROOT.SetBatch(1) 

	if (options.stats): 
		ROOT.gStyle.SetOptStat(1111111)
	else: 
		# Don't plot stats box
		ROOT.gStyle.SetOptStat(0) 

	outputfolder = "./plots/"+options.out+"/"
	if (options.forcepath): 
		print("WARNING: You have used option '-f' or '--forcepath'. Files will be written to: {}".format(options.out))
		outputfolder = options.out+"/"

	os.system("mkdir -p "+outputfolder)


	# Global initialisations
	Ana.Init(options.version)


	data = "dataD2" #"dataB2"
	filesUsed = ["Sig", data, "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "B0toDstar3pi", "dataD1WS"] #"SigPart", "dataD2WS", "dataD2TauWS", , "B0toDstar5pi" ["Sig", "dataD1", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "BkgDstara1", "B0toDstar3pi", "dataD2WS"]

	regions = ["SR", "CR", "SB"]

	variables = ["b_tau_rhomass1", "b_tau_rhomass2", "b_B_q2", "b_B_m", "b_tau_m", "b_B_proper_xi_rho1", "b_B_proper_xi_rho2", "b_tau_proper_alpha_rho1_pi", "b_tau_proper_alpha_rho2_pi", "b_tau_proper_theta_rho1", "b_tau_proper_theta_rho2"]


	nBins = 6
	rangeMin = 0.2 #0.37
	rangeMax = 1.5 #1.43
	

	for file in filesUsed: 
		Ana.filemanager.OpenItem(file)


	# Starting script 
	canvas = ROOT.RatioCanvas("romassunrolled", "Unrolled 2D distribution of rho mass", 800, 600)

	samples = {}
	frames = collections.defaultdict(dict)
	baseline = collections.defaultdict(dict)
	histos = collections.defaultdict(dict)
	histosunrolled = collections.defaultdict(dict)

	for item in filesUsed:  
		samples[item] = ROOT.RDataFrame(Ana.filemanager.GetItem(item))
		frames[item]["baseline"] =  samples[item].Filter((Ana.cut["base"]+Ana.samples.at(item).cut).GetTitle()) #Ana.cut["base"].GetTitle() "1."
		frames[item]["all"] = samples[item].Filter("1.")
		#baseline[item]
		for region in regions: 
			cut = (Ana.cut[region]+Ana.samples.at(item).cut).GetTitle()
			#if "WS" in item: 
			#	cut = Ana.cutstandalone[region].GetTitle()
			if (options.debug): print("Using following cut string (from TCut): {}".format(cut))
			frames[item][region] = samples[item].Filter(cut)
			ROOT.SetOwnership(frames[item][region], 0)
			# For histogram legacy compatibility
			histos[item][region] = frames[item][region].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
			histosunrolled[item][region] = UnrollHist(histos[item][region])

	if options.debug: print(frames)


	#BackgroundShapeUnrolled(histosunrolled["dataD2"]["SR"], histosunrolled["dataD2"]["SB"], [histosunrolled["Sig"]["SR"], histosunrolled["BkgDstarDs"]["SR"]], outputfolder+"UnrolledRhoMass")

	#AtomicDraw(frames["Sig"]["SR"].Histo1D("b_tau_rhomass1"), outputfolder+"/SignalFromNew.png")

	#colors = {"Sig":2, "BkgDstarDs":3, "BkgDstarDsstar":8, "BkgDstara1":4, "dataD2WS":6, "dataD2TauWS":7, "other":9, "yetanother":1} #4, 3, 6, 7, 9 colors = [2, 3, 8, 4, 6, 7, 9, 1]
	colors = [2, 3, 8, 4, 6, 7, 9, 1]

	#GetEfficiencies(frames)

	effs = ComputeEfficiencies(frames, "all")

	#for key, item in effs.iteritems(): 
	#	print("{}: {}".format(key, item["baseline"]))

	#PrintEfficiencies(effs)

	#DumpEffs(effs, "efficienciesTest.json")

	#newEffs = ReadEffs("efficienciesTest.json")

	#PrintEfficiencies(newEffs)

	#value = InitialEffs(7.98)

	#print(value)

	#selectioneffs = CompleteEffsFromFile(value, "v3", Ana.filemanager)

	#print(selectioneffs)

	#selectioneffs["dataD2WS"] = ufloat(-0.0003, 0.)  #29
	#selectioneffs["dataD2TauWS"] = ufloat(-0.0005, 0.) #21.5

	selectioneffs = ReadEffsSimple("./data/etc/Expectedyields.json")

	print(selectioneffs)
	print(effs)

	selectioneffs["dataD2WS"] = ufloat(1.39e-5*-12.0*7*8, 0.)
	selectioneffs["dataD1WS"] = ufloat(1.39e-5*-12.0*11*22*3.3, 0.)

	regioneffs = MultiplyFinalEffs(selectioneffs, effs)

	PrintEfficiencies(regioneffs)

	from libEfficiencies import DumpEffs2D

	DumpEffs2D(regioneffs, "./data/etc/RegionEffs.json")

	# Quick and dirty significance computation
	s = regioneffs["Sig"]["SR"]
	b = frames[data]["SR"].Count().GetValue()
	significance = s/sqrt(b) # b is here s+b since taken from data

	print("Expected significance: {}/sqrt({}) = {}".format(s, b, significance))

	print("Signal over background ratio: {}".format(s/b))


	files = filesUsed #["Sig", "BkgDstarDs", "BkgDstarDsstar", "dataD2WS", "dataD2TauWS"]
	#PlotOverlay(frames, "dataD2", ["Sig", "dataD2", "BkgDstarDs", "BkgDstarDsstar", "BkgDstara1Part"], regions, variables, regioneffs, outputfolder)

	PlotStack(frames, data, files, regions, variables, regioneffs, outputfolder, False)

	#PlotComparison(frames, data, "dataD2WS", regions, variables, outputfolder, False)

	#PublishToWeb(outputfolder, "Modelling_23_9_27")

	#WriteDatacard(frames, regioneffs, variables, regions, "dataB2")

	#files = {"Sig":LoadFile("/Users/mhuwiler/eos/DoctoralThesis/Analysis/data/v3/Sig.root"), "BkgDstara1":LoadFile("/Users/mhuwiler/eos/DoctoralThesis/Analysis/data/v3/BkgDstara1.root") }
	#newframes, histos, histisunrolled = PrepareCustomFiles(files, regions)

	#print(newframes)

	#PlotComparison(frames, "BkgDstara1Part", "Sig", regions, variables, outputfolder, False) #["t_B_mu_alpha", "t_B_m", "t_tau_m"]["tau_rhomass1", "tau_rhomass2", "B_m", "B_q2"]

	
	#with open("normalisationsFromFit.json", "r") as fitefffile: 
	fiteffs = ReadEffs("normalisationsFromFit.json")
	fiteffs.update(regioneffs)
	fiteffs["dataD2"]["SB"] = ufloat(-1., 0.)
	fiteffs["dataD2"]["CR"] = ufloat(-1., 0.)
	PlotStack(frames, "dataD2", ["Sig", "dataD2", "dataD2-BkgDstarDs"], ["SB", "CR"], variables, fiteffs, "./plots/closuretest/", False)


	Ana.filemanager.CloseAll()



