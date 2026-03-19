#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import collections
import copy
from argparse import ArgumentParser
#ROOT.gROOT.LoadMacro("FileFlow.h+")
import anaConfig
from ROOT import Ana, TCanvas, TH1D, TPad, TLegend, THStack, RDataFrame
from libUtils import HoldUntilKeyPress


def PlotStack(session, dataname, initialcomponents, variables, outfolder, drawlegend=False): 
	# Plotting distributions over each other 
	#outfolder+="stacked/"
	os.system("mkdir -p "+outfolder)
	name = "FitResult"
	factor = 1.1 # how much overhead to add to the histos 
	components = copy.deepcopy(initialcomponents)
	if (dataname in components): components.remove(dataname)
	components.reverse()
	notYetDrawn = True

	# Setting up plot 
	canvas = TCanvas("fitresult", "Fit result", 800, 600)

	legend = TLegend(canvas.GetLeftMargin()+0.35, 
	                         	1.-canvas.GetTopMargin()-.2, 
	                            canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
	                           	1.-canvas.GetTopMargin() )


	stack = THStack("stack", "Fit result")

	data = session.Get(dataname)

	data.SetMarkerStyle(8) # Large scalable dot
	data.SetMarkerSize(0.5)
	data.SetLineColor(ROOT.kBlack)
	data.SetTitle("") #data.SetTitle("{}_{}".format(variable, region))
	#data.SetFillColor(ROOT.kBlack)
	legend.AddEntry(data, "data", "PE")
	data.Draw("E")

	numcomponents = 0 #1
	for component in components: 
		print(component)
		#component = component.replace("Dist", "")
		histo = session.Get(component)
		print("Norm {}: {}".format(component, histo.Integral()))
		#histo.Draw()
		#HoldUntilKeyPress()

		histo.SetLineStyle(1) # plain
		histo.SetLineWidth(2)
		color = Ana.samples.at(component.replace("Part", "")).color
		if (not color): 
			color = colors[0]
		histo.SetLineColor(color)
		#histo.SetMarkerColor(Ana.color[component])
		histo.SetFillStyle(1001)
		histo.SetFillColor(color)
		histo.SetTitle("")
		#hists[component] = histo
		legend.AddEntry(histo, Ana.samples.at(component).legend, "F")

		stack.Add(histo)

		#ROOT.SetOwnership(histo, 0)
		numcomponents += 1

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
	canvas.Draw()

	maxes = [data.GetMaximum(), stack.GetMaximum()]

	print(maxes)

	data.SetMaximum(factor*max(maxes))
	canvas.Update()

	canvas.Print(outfolder+name+".png")
	canvas.Print(outfolder+name+".pdf")


	if ((not drawlegend) and notYetDrawn): 
		#components.reverse()
		canv = TCanvas("legendCanvas", "legenCanvas", 800, 200*numcomponents)
		dummy = TCanvas("dummy", "dummy", 800, 600)
		#legend.AddEntry(data, "data", "PE")
		# for component, histo in hists.iteritems(): 
		# 	#histo = frames[component][region].Histo1D(examplehist, variable)
		# 	ROOT.SetOwnership(histo, 0)
		# 	histo.SetLineStyle(1) # plain
		# 	histo.SetLineWidth(2)
		# 	color = Ana.samples.at(component.replace("Part", "")).color
		# 	if "WS" in component: 
		# 		color = Ana.samples.at("WS").color
		# 	histo.SetLineColor(color)
		# 	histo.SetFillStyle(1)
		# 	histo.SetFillColor(color)
		# 	legend.AddEntry(histo, Ana.samples.at(component.replace("Part", "")).legend, "F")
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
	#outfolder+="comparisons/"
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

			print(variable)
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
			legend.AddEntry(reference, referencename, "L")
			reference.Draw("HIST E")

			comparison = frames[comparisonname][region].Histo1D(examplehist, variable)
			comparison.SetLineWidth(2)
			comparison.SetLineColor(ROOT.kRed)
			comparison.SetFillStyle(3356)
			comparison.SetFillColor(ROOT.kRed)
			legend.AddEntry(comparison.GetPtr(), comparisonname, "L")
			if (comparison.Integral() != 0):
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


if __name__ == "__main__":

	parser = ArgumentParser(description="SignalBackground")
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("--out", dest="out", action="store", type=str, default="TestSimpleFitFullStatsnorms/", help="Directory where the plots shuld go")
	parser.add_argument("--name", dest="name", action="store", type=str, default="test", help="Turn on debug output")
	parser.add_argument("--file", dest="file", action="store", type=str, default="postFitPlots.root", help="Directory where the plots shuld go")
	parser.add_argument("--directory", dest="directory", action="store", type=str, default="SR_postfit", help="Directory where the plots shuld go")
	parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v1", help="Which version (cycle) of files to run on")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("-f", "--forcepath", dest="forcepath", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
	parser.add_argument('-l', "--variable", dest="variable", action="store", default="", help="Variable to reprocess")
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


	samples = ["Sig", "data_obs", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "ABCD"] #anaConfig.samples #["Sig", "data_obs", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "dataD2WS"] # "B0toDstar3pi", #"SigPart", "dataD2WS", "dataD2TauWS", , "B0toDstar5pi" ["Sig", "dataD1", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "BkgDstara1", "B0toDstar3pi", "dataD2WS"]
	for i in range(0, len(samples)): 
		if (("data" in samples[i]) and (not "WS" in samples[i])): 
			samples[i] = "data_obs" 

	regions = ["SR", "CR", "SB"]

	#variables = ["b_tau_rhomass1", "b_tau_rhomass2", "b_B_q2", "b_B_m", "b_tau_m"]

	ROOT.gInterpreter.Declare("std::vector<std::string> getKeys(const std::unordered_map<std::string, ROOT::RDF::TH1DModel>& map) { std::vector<std::string> result; result.reserve(map.size()); for (auto item = map.begin(); item != map.end(); item++) { result.push_back(item->first); } return result; };") #"TH1F * convertHisto(TH1D *histo) { return static_cast<TH1F*>(histo); } "

	varvector = ROOT.getKeys(Ana.binning)
	variables = [] #variables = [ var for var ]
	for i in range(varvector.size()): 
		variables.append(varvector.at(i))

	if (options.variable != ""): 
		variables = options.variable.split(",")

	print(variables)


	file = ROOT.TFile.Open(options.file, "READ")

	session = file.Get(options.directory)


	print(session)

	PlotStack(session, "data_obs", samples, variables, outputfolder)
	

	from webInterface import PublishToWeb
	#PublishToWeb(outputfolder, "Variables_23_8_16_beforeBDT_data") #Variables_23_8_14_beforeBDT_Sigvsdata



