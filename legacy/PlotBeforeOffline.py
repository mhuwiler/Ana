#!/usr/bin/env python
from __future__ import division, print_function

import ROOT
import os
import math
import collections
import copy
from argparse import ArgumentParser
ROOT.gROOT.LoadMacro("FileFlow.h+")
from ROOT import Ana, TCanvas, TH1D, TPad, TLegend, THStack, RDataFrame


def PlotOverlay(frames, dataname, initialcomponents, regions, variables, yields, outfolder, drawlegend=True): 
	# Plotting distributions over each other 
	#outfolder+="overlay/"
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
				ROOT.SetOwnership(histo, 0)
				histo.SetLineStyle(1) # plain
				histo.SetLineWidth(2)
				histo.SetLineColor(Ana.colorold[component.replace("Part", "")]) #colors[i]Ana.color[component.replace("Part", "")]
				#histo.SetFillStyle(3003)
				#histo.SetFillColorAlpha(Ana.color[component], 0.4)
				if (yields[component][region].nominal_value > 0.): 
					histo.Scale(yields[component][region].nominal_value/histo.Integral())
				elif (yields[component][region].nominal_value != -1.):
					histo.Scale(-yields[component][region].nominal_value*histo.Integral())
				#histo.SetMarkerColor(Ana.color[component])
				histo.Draw("HIST SAME")
				maxes.append(histo.GetMaximum())
				legend.AddEntry(histo.GetPtr(), component)
				i+=1

			if (drawlegend): legend.Draw()
			legend.SetBorderSize(1)
			legend.SetMargin(0.3)
			legend.SetTextSize(0.04)

			data.SetMaximum(factor*max(maxes))
			canvas.Draw()

			canvas.Print(outfolder+name+".png")
			canvas.Print(outfolder+name+".pdf")


def PlotStack(frames, dataname, initialcomponents, regions, variables, yields, outfolder, drawlegend=True): 
	# Plotting distributions over each other 
	#outfolder+="stacked/"
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
			legend.AddEntry(reference.GetPtr(), referencename, "L")
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


if __name__ == "__main__":

	parser = ArgumentParser(description="SignalBackground")
	#parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
	parser.add_argument("--out", dest="out", action="store", type=str, default="SigCombPreBDT/", help="Directory where the plots shuld go")
	parser.add_argument("--name", dest="name", action="store", type=str, default="test", help="Turn on debug output")
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


	filesUsed = ["Sig", "dataD2WS", "dataB2", "B0toDstarDsPart", "B0toDstarDsDist"] #"SigPart", "dataD2WS", "dataD2TauWS", , "B0toDstar5pi" ["Sig", "dataD1", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "BkgDstara1", "B0toDstar3pi", "dataD2WS"]

	regions = ["SR", "CR", "SB"]

	variables = ["b_tau_rhomass1", "b_tau_rhomass2", "b_B_q2", "b_B_m", "b_tau_m"]

	ROOT.gInterpreter.Declare("std::vector<std::string> getKeys(const std::unordered_map<std::string, ROOT::RDF::TH1DModel>& map) { std::vector<std::string> result; result.reserve(map.size()); for (auto item = map.begin(); item != map.end(); item++) { result.push_back(item->first); } return result; };") #"TH1F * convertHisto(TH1D *histo) { return static_cast<TH1F*>(histo); } "

	varvector = ROOT.getKeys(Ana.binning)
	variables = [] #variables = [ var for var ]
	for i in range(varvector.size()): 
		variables.append(varvector.at(i))

	if (options.variable != ""): 
		variables = options.variable.split(",")

	print(variables)


	nBins = 6
	rangeMin = 0.2 #0.37
	rangeMax = 1.5 #1.43
	

	for file in filesUsed: 
		Ana.filemanager.OpenItem(file)



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

	# from here on starts teting

	PlotComparison(frames, "Sig", "dataB2", ["baseline"], variables, outputfolder, False)

	#variables = ["D0_pt", "D0_eta", "D0_eta", "D0_alpha", "D0_fl", "D0_gen_pt", "D0_gen_eta", "D0_lip", "D0_mu_alpha", "D0_prefit_pt", "D0_prefit_eta", "D0_prefit_kpt", "D0_prefit_keta", "D0_prefit_kphi", "D0_prefit_pipt", "D0_prefit_pieta", "D0_prefit_piphi", "D0_pvip", "Dstar_pt", "Dstar_pt", "Dstar_phi", "Dstar_alpha", "Dstar_gen_pt", "Dstar_gen_eta", "Dstar_lip", "Dstar_mu_alpha", "Dstar_prefit_pt", "Dstar_prefit_eta", "Dstar_prefit_phi", "Dstar_prefit_pispt", "Dstar_prefit_piseta", "Dstar_prefit_pisphi", "Dstar_pvip", "b_tau_pt", "b_tau_eta", "b_tau_phi", "b_tau_alpha", "b_tau_fl", "b_tau_lip", "b_tau_mu_alpha", "b_tau_prefit_pt", "b_tau_prefit_eta", "b_tau_prefit_phi", "b_tau_pvip", "b_B_pt", "b_B_eta", "b_B_phi", "b_B_alpha", "b_B_fl", "b_B_lip", "b_B_mu_alpha", "b_B_prefit_pt", "b_B_prefit_eta", "b_B_prefit_phi", "b_B_pvip"]
	#PlotComparison(frames, "B0toDstarDsPart", "B0toDstarDsDist", ["baseline"], variables, outputfolder, False)

	from webInterface import PublishToWeb
	PublishToWeb(outputfolder, "Variables_23_8_16_beforeBDT_data") #Variables_23_8_14_beforeBDT_Sigvsdata



