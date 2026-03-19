import ROOT
import os
import math
import uproot
import numpy as np
import collections

uproot.default_library = "np"


#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/libFunctions.C+")#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/ExperimentSpecificLayer.C")
#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/FileManager/CFileManager.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/CMS/tdrstyle.C")
ROOT.gROOT.LoadMacro("FileFlow.h")
ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
#ROOT.setTDRStyle()
#import CMS_lumi


ROOT.Ana.Init("v1")


bkgInSample = 9400000000*0.0000376
sigInSample = 1130


graphcollection = ROOT.vector('std::pair<TGraph*,TString>')()

filesUsed = ["dataD2", "Sig", "BkgDstarDs"]


for file in filesUsed: 
	ROOT.Ana.filemanager.OpenItem(file)


plotstats = False

webpublication =False


outputfolder = "./plots/BackgroundEstimateUpdate/"

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
			rep = raw_input( "(type 'yes' or 'no'): " ).lower()
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

	print "Nunmber of bins: {}, {}".format(nx, ny)

	for i in range(0, nx):
		for j in range(0, ny): # TODO: check overflow is handled properly 
			if inverted: 
				binContent = histo2D.GetBinContent(j, i)
			else: 
				binContent = histo2D.GetBinContent(i, j)
			unrolled.SetBinContent(i+j*nx, binContent)

	return unrolled


if plotstats: 
	ROOT.gStyle.SetOptStat(1111111)

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



#ROOT.gInterpreter.Declare("""
#	double Rhomass2DUnrolled(Float_t rhomass1, Float_t rhomass2)
#	{
#		return int((std::min(rhomass2, 1.3) - 0.2)/0.22) + 6*int((std::min(tau_rhomass1, 1.3) - 0.2)/0.22); 
#	}
#""")

# Don't plot stats box
ROOT.gStyle.SetOptStat(0)

for quantity in ["Rhomass2Dunrolled"]: 
	print "Plotting {}".format(quantity)
	#canvas = ROOT.TCanvas("romassunrolled", "Unrolled 2D distribution of rho mass", 800, 600) #ROOT.RatioCanvas(quantity, quantity, 950, 800) #800, 800
	#canvas.SetMiddleMargin(0.13)
	#canvas.SetPadDelimitation(0.34)
	#canvas.SetRightMargin(0.12)

	canvas = ROOT.RatioCanvas("romassunrolled", "Unrolled 2D distribution of rho mass", 800, 600)

	# Signal region definitions
	mvaThreshold = 0.9
	mvaLowThreshold = 0.0
	mvaLowerBound = -0.5

	# Histo parameters
	nBins = 6
	rangeMin = 0.2
	rangeMax = 1.5

	dataframes = {}
	histos = collections.defaultdict(dict)
	dataframes["data"] = ROOT.RDataFrame(ROOT.Ana.filemanager.GetItem("dataD2"))
	#data = ROOT.RDataFrame(filemanager.GetItem("ParkingBPH1Run2018D")) #"ParkingBPH4-6Run2018B"
	dataframes["MC"] = ROOT.RDataFrame(ROOT.Ana.filemanager.GetItem("Sig"))
	dataframes["DstarDs"] = ROOT.RDataFrame(ROOT.Ana.filemanager.GetItem("BkgDstarDs"))

	SR = {}
	CR = {}
	SB = {}

	for key in dataframes.keys(): 
		dataframes[key] = dataframes[key].Filter("mvaScore>-2.")
		dataframes[key] = dataframes[key].Define("rhomass2D", "int((min(b_tau_rhomass2, float(1.3)) - 0.2)/0.22) + 6*int((min(b_tau_rhomass1, float(1.3)) - 0.2)/0.22)")
		SR[key] = dataframes[key].Filter("mvaScore>={}".format(mvaThreshold))
		CR[key] = dataframes[key].Filter("(mvaScore<{})&&(mvaScore>{})".format(mvaThreshold, mvaLowThreshold))
		SB[key] = dataframes[key].Filter("(mvaScore<{})&&(mvaScore>{})".format(mvaLowThreshold, mvaLowerBound))
		#histoDataSB = dataSB.Histo1D(("rhomass2D", "rhomass2D", 30, 0., 30.), "rhomass2D")
		histos[key]["SR"] = UnrollHist(SR[key].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2"))
		histos[key]["CR"] = UnrollHist(CR[key].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2"))
		histos[key]["SB"] = UnrollHist(SB[key].Histo2D(("rhomass1", "rhomass2", nBins, rangeMin, rangeMax, nBins, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2"))

	
	#histo2D = data.Histo2D(("rhomass1", "rhomass2", 10, 0., 3., 10, 0., 3.), "b_tau_rhomass1", "b_tau_rhomass2")

	#hist = UnrollHist(histo2D)

	showBKG = False


	print "Number of events (SR, CR, SB): {}, {}, {}".format(histos["data"]["SR"].GetEntries(), histos["data"]["CR"].GetEntries(), histos["data"]["SB"].GetEntries())

	estimate = histos["data"]["CR"]
	reference = histos["data"]["SR"]
	additionalhists = [histos["MC"]["SR"], histos["DstarDs"]["SR"]]
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


	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	hist = SR["MC"].Histo2D(("rhomass1", "rhomass2", 15, rangeMin, rangeMax, 15, rangeMin, rangeMax), "b_tau_rhomass1", "b_tau_rhomass2")
	hist.SetTitle("m(#rho_{1}) vs m(#rho_{2})")
	hist.DrawCopy("COLZ")
	canv.Draw()
	canv.Print(outputfolder+"2Drhomass.pdf")

	#f = ROOT.TFile.Open(presentationfolder, "Update")
	#fomcanvas.Write()
	#canvas.Write()
	#f.Write()
	#f.Close()

	AtomicDraw(SR["MC"].Histo1D("b_tau_rhomass1"), outputfolder+"/SignalFromOld.png")




ROOT.Ana.filemanager.CloseAll()



