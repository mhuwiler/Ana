#!/usr/bin/env python
from __future__ import division, print_function

import os
import copy
from ROOT import TCanvas, TLegend, THStack, gROOT


#gROOT.LoadMacro("FileFlow.h")
import ROOT
from ROOT import Ana


ROOT.gStyle.SetOptStat(0) 

class Colors: 
	def __init__(self, colorlist): 
		self.colorlist = colorlist
		self.idx = 0
		self.maxidx = len(colorlist)
		assert(self.maxidx > 0)
	def next(self): 
		color = self.idx
		self.idx+=1
		if (self.idx>self.maxidx): 
			self.idx = 0
		return self.colorlist[color]

	colorlist = []
	idx = 0
	maxidx = 0


defaultcolors = Colors([2, 3, 8, 4, 6, 7, 9, 1])

normalisebinwidth = False



def binning(variable): 
	#print(variable)
	variables = [item[0] for item in Ana.binning]
	#print(variables)
	if variable in variables: 
		model = Ana.binning[variable]
	else: 
		proxy = ""
		for element in variables: 
			if element in variable: 
				proxy = element
				break
		#print(proxy)
		model = Ana.binning[proxy]
	
	return model


def PlotOverlay(frames, dataname, initialcomponents, regions, variables, yields, outfolder, drawlegend=True, normalise=False): 
	# Plotting distributions over each other 
	from anaPrepareRegions import GetBaseName
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
				#ROOT.SetOwnership(histo, 0)
				if (histo.Integral()==0): 
					continue
				print(component)
				histo.SetLineStyle(1) # plain
				histo.SetLineWidth(2)
				color = Ana.samples.at(GetBaseName(component.replace("Part", ""))).color
				if (not color): 
					color = defaultcolors.next()
				histo.SetLineColor(color) #colors[i]Ana.color[component.replace("Part", "")]
				#histo.SetFillStyle(3003)
				#histo.SetFillColorAlpha(color, 0.4)
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
				legend.AddEntry(histo.GetPtr(), GetBaseName(component))
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
	from anaPrepareRegions import GetBaseName, GetABCDcomponent
	import anaConfig
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
			if (normalisebinwidth): Ana.normaliseBinContent(data.GetPtr())
			legend.AddEntry(data.GetPtr(), "data", "PE")
			data.Draw("E")

			comb = GetABCDcomponent(anaConfig.data, "(b_tau_sumdnn>2.)", "b_B_nmu<1&&b_B_ne<1&&b_B_nh<1", examplehist, variable, [element for element in components if element != anaConfig.Sig])

			stack = THStack("stack", "Background modelling")
			hists = {}
			MCstats = 0.
			WS = comb
			for component in components: 
				if not "WS" in component: 
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
					if (not "WS" in component):
						MCstats += histo.Integral()
					else:
						WS = comb #histo
					ROOT.SetOwnership(histo, 0)
					histo.SetLineStyle(1) # plain
					histo.SetLineWidth(2)
					color = Ana.samples.at(GetBaseName(component.replace("Part", ""))).color
					if (not color): 
						print("No color for sample {}".format(component))
						color = defaultcolors.next()
					histo.SetLineColor(color)
					#histo.SetMarkerColor(Ana.color[component])
					histo.SetFillStyle(1001)
					histo.SetFillColor(color)
					hists[component] = histo
					if (normalisebinwidth): Ana.normaliseBinContent(histo)
					stack.Add(histo)
					#legend.AddEntry(histo.GetPtr(), Ana.legends[component], "F")

			WS.SetLineStyle(1) # plain
			WS.SetLineWidth(2)
			color = Ana.samples.at(GetBaseName("WS")).color
			if (not color): 
				print("No color for sample {}".format(component))
				color = defaultcolors.next()
			WS.SetLineColor(color)
			#WS.SetMarkerColor(Ana.color[component])
			WS.SetFillStyle(1001)
			WS.SetFillColor(color)
			hists["WS"] = WS
			if (normalisebinwidth): Ana.normaliseBinContent(WS)
			stack.Add(WS)
			#ABCDnorm = (data.Integral()-MCstats)/WS.Integral()
			WS.Scale((data.Integral()-MCstats)/WS.Integral())
			with open("ABCDnorm.txt", "w") as file: 
				file.write("{}\n".format(data.Integral()-MCstats))

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
			canvas.SetTopMargin(0.1)
			canvas.SetLeftMargin(0.15)
			OverheadText()

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
					color = Ana.samples.at(GetBaseName(component.replace("Part", ""))).color
					if "WS" in component: 
						color = Ana.samples.at("WS").color
					histo.SetLineColor(color)
					histo.SetFillStyle(1001)
					histo.SetFillColor(color)
					legend.AddEntry(histo, Ana.samples.at(GetBaseName(component.replace("Part", ""))).legend, "F")
				canv.cd()
				data.SetMarkerSize(4.)
				data.SetLineWidth(4)
				legend.SetX1(0.)
				legend.SetY1(0.)
				legend.SetX2(1.)
				legend.SetY2(1.)
				legend.SetBorderSize(0)
  				legend.SetFillColor(0)
  				legend.SetFillStyle(1001)
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

			examplehist = Ana.binning[variable]
			reference = frames[referencename][region].Histo1D(examplehist, variable)
			reference.SetTitle("") #"{}_{}".format(variable, region)
			#reference.SetMarkerStyle(8) # Large scalable dot
			#reference.SetMarkerSize(0.5)
			reference.SetLineWidth(2)
			reference.SetLineColor(ROOT.kBlue)
			reference.SetFillStyle(3003)
			reference.SetFillColor(ROOT.kBlack)
			#reference.SetTitle("{}_{}".format(variable, region))
			legend.AddEntry(reference.GetPtr(), referencename, "F")
			reference.Draw("HIST E")

			comparison = frames[comparisonname][region].Histo1D(examplehist, variable)
			comparison.SetLineWidth(2)
			comparison.SetLineColor(ROOT.kRed)
			comparison.SetFillStyle(3356)
			comparison.SetFillColor(ROOT.kRed)
			legend.AddEntry(comparison.GetPtr(), comparisonname, "F")
			if normalise: 
				if (yields[component][region] > 0. and comparison.Integral() > 0): 
					comparison.Scale(yields[component][region].nominal_value/comparison.Integral())
				elif (yields[component][region].nominal_value != -1.):
					comparison.Scale(-yields[component][region].nominal_value*comparison.Integral())
			else:
				if not (comparison.Integral()==0): 
					comparison.Scale(reference.Integral()/comparison.Integral())

			comparison.Draw("HIST SAME E") #"SAME"
			legend.Draw()
			legend.SetBorderSize(1)
			legend.SetMargin(0.3)
			legend.SetTextSize(0.04)

			maxes = [reference.GetMaximum(), comparison.GetMaximum()]

			reference.GetXaxis().SetTitleSize(0.06)
			reference.GetXaxis().SetLabelSize(0.06)
			reference.GetYaxis().SetLabelSize(0.06)
			reference.GetYaxis().SetTitleSize(0.06)
			reference.GetXaxis().SetTitleOffset(1.2)
			canvas.SetBottomMargin(0.15)
			canvas.SetTopMargin(0.1)
			canvas.SetLeftMargin(0.15)
			#OverheadText()
			reference.SetMaximum(factor*max(maxes))
			canvas.Draw()

			canvas.Print(outfolder+name+".png")
			canvas.Print(outfolder+name+".pdf")


def PlotSimple(frame, variable, outfolder = ".", name = "", histomodel = None): 
	if name == "": 
		name = "/{}".format(variable)
	print(name)
	os.system("mkdir -p "+outfolder)
	drawlegend = False
	if (histomodel == None): # Take the default binning for this variable defined in the analysis
		histomodel = binning(variable)
	canvas = TCanvas(name, variable, 800, 600)

	legend = TLegend(canvas.GetLeftMargin()+0.35, 
		1.-canvas.GetTopMargin()-.15, 
		canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
		1.-canvas.GetTopMargin() )

	color = ROOT.kBlue
	histo = frame.Histo1D(histomodel, variable)
	histo.SetMarkerStyle(8) # Large scalable dot
	histo.SetMarkerSize(0.5)
	histo.SetLineColor(color)
	histo.SetLineWidth(2)
	#histo.SetTitle("{}_{}".format(variable, region)) # TODO: Delete once the binning is centralised
	#histo.SetFillColor(ROOT.kBlack)
	histo.SetMarkerColor(color)
	legend.AddEntry(histo.GetPtr(), variable, "LPE")
	histo.Draw("HIST E")

	
	if (drawlegend): legend.Draw()
	legend.SetBorderSize(1)
	legend.SetMargin(0.3)
	legend.SetTextSize(0.04)

	canvas.SaveAs(outfolder+name+".png")
	canvas.SaveAs(outfolder+name+".pdf")


def PlotFitResult(session, dataname, initialcomponents, variables, outfolder, drawlegend=False): 
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
		#histo.Draw()
		#HoldUntilKeyPress()

		histo.SetLineStyle(1) # plain
		histo.SetLineWidth(2)
		color = Ana.samples.at(component.replace("Part", "")).color
		if (not color): 
			color = defaultcolors.next()
		histo.SetLineColor(color)
		#histo.SetMarkerColor(Ana.color[component])
		histo.SetFillStyle(1)
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


def PlotSingle(frame, variables, outfolder, color = 2, legendtext="", normalise=1., drawoptions="HIST"): 
	# Plotting distributions over each other 
	os.system("mkdir -p "+outfolder)
	factor = 1.1 # how much overhead to add to the histos 
	for variable in variables: 
		name = "{}".format(variable) #, region)
		canvas = TCanvas(name, variable, 800, 600)

		legend = TLegend(canvas.GetLeftMargin()+0.35, 
                         	1.-canvas.GetTopMargin()-.2, 
                            canvas.GetLeftMargin()+(1.-(canvas.GetLeftMargin()+canvas.GetRightMargin())),
                           	1.-canvas.GetTopMargin() )

		variablename = VariableName(variable)
		if (Ana.binning.find(variablename) == Ana.binning.end()): 
			histo = frame.Histo1D(variable)
		else:
			examplehist = Ana.binning[variablename]
			histo = frame.Histo1D(examplehist, variable)
		#ROOT.SetOwnership(histo, 0)
		#print(variable)
		#print(frame.GetColumnNames())
		#print(examplehist)
		ROOT.SetOwnership(histo, 0)
		histo.SetLineStyle(1) # plain
		histo.SetLineWidth(2)
		#color = Ana.samples.at(GetBaseName(component.replace("Part", ""))).color
		histo.SetLineColor(color)
		histo.SetMarkerColor(color)
		#histo.SetFillColor(color)
		histo.Scale(normalise/histo.Integral())
		histo.SetTitle("{}".format(variable)) # TODO: Delete once the binning is centralised
		histo.Draw(drawoptions)

		maxes = [histo.GetMaximum()]

		
		legend.AddEntry(histo.GetPtr(), legendtext) # , "PE"

		if (legendtext): legend.Draw()
		legend.SetBorderSize(1)
		legend.SetMargin(0.3)
		legend.SetTextSize(0.04)

		histo.SetMaximum(factor*max(maxes))
		canvas.Draw()

		canvas.Print(outfolder+name+".png")
		canvas.Print(outfolder+name+".pdf")


def VariableName(name): 
	name = name.replace("v_", "b_")
	name = name.replace("pttau_", "b_")
	return name


def OverheadText(additionaltext = "Private Work", plot = None, position = 1): 
	if (plot == None): 
		plot = ROOT.gPad
	ROOT.gROOT.LoadMacro("cliPlotting.C+")
	ROOT.SetHeader(plot, position, additionaltext)

