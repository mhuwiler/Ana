import ROOT
import os
import math
import uproot
import numpy as np

uproot.default_library = "np"


#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/libFunctions.C+")#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/ExperimentSpecificLayer.C")
ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/FileManager/CFileManager.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/CMS/tdrstyle.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
#ROOT.setTDRStyle()
#import CMS_lumi


filemanager = ROOT.FileManager()



filemanager.AddItem("cutflowRef", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/TauCutflowReference.root", "ntuplizer/Taucutflow")
filemanager.AddItem("cutflowRefEff", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/TauCutflowReference.root", "ntuplizer/Taucutflow_eff")

filemanager.AddItem("cutflowGen", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/TauCutflowGen.root", "ntuplizer/Taucutflow")
filemanager.AddItem("cutflowGenEff", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/TauCutflowGen.root", "ntuplizer/Taucutflow_eff")

filemanager.AddItem("TMVAROC", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "xgboOptimized/TMVA-like_ROC")
filemanager.AddItem("bkgEff", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "xgboOptimized/bkgEff(sigEff)")
filemanager.AddItem("bkgEffEval", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/ROCsummary.root", "bkgEff(sigEff)")

filemanager.AddItem("SignalOfficialMC50M", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted_mva.root", "tree")
filemanager.AddItem("ParkingBPHAllRun2018B", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPHRun2018B_converted_mva.root", "tree")

filemanager.AddItem("SignalOfficialMC50M", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted_mva.root", "tree")
filemanager.AddItem("ParkingBPH4-6Run2018B", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPH4-6Run2018B_converted_mva.root", "tree")

filemanager.AddItem("signalTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "signalTest")
filemanager.AddItem("backgroundTrain", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "backgroundTrain")
filemanager.AddItem("backgroundTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "backgroundTest")
filemanager.AddItem("DataBackground", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge_converted.root", "tree")

filemanager.AddItem("signal", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted_mvanew.root", "tree")
filemanager.AddItem("background", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPHRun2018B_converted_mvanew.root", "tree")



filemanager.OpenAllItems()

bkgInSample = 9400000000*0.0000376
sigInSample = 1130


graphcollection = ROOT.vector('std::pair<TGraph*,TString>')()


plotstats = False

webpublication =False


outputfolder = "./plots/"

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

for quantity in ["ROC"]: 
	print "Plotting {}".format(quantity)
	canvas = ROOT.TCanvas("figure", "Figure of merit", 800, 600) #ROOT.RatioCanvas(quantity, quantity, 950, 800) #800, 800
	#canvas.SetMiddleMargin(0.13)
	#canvas.SetPadDelimitation(0.34)
	#canvas.SetRightMargin(0.12)

	ROC = filemanager.GetItem("TMVAROC")

	#ROC = filemanager.GetItem("bkgEff")

	ROC.Draw()

	ROC.SetMarkerColor(ROOT.kGreen)
	ROC.SetLineColor(ROOT.kGreen)

	#canvas.SetLogy()
	#canvas.SetLogx()

	effs = filemanager.GetItem("bkgEffEval")

	numPoints = effs.GetN()

	lowestbin = -1.
	highestbin = 1.
	nBins = 1000

	signal = exportHistFromTree(filemanager.GetItem("signal"), "mvaScore", "1", nBins, lowestbin, highestbin)
	background = exportHistFromTree(filemanager.GetItem("background"), "mvaScore", "1", nBins, lowestbin, highestbin)

	#signal = filemanager.GetItem("signalTest")
	#backgroud = filemanager.GetItem("backgroundTest")

	assert(nBins == signal.GetNbinsX())
	assert(nBins == background.GetNbinsX())



	FOM = ROOT.TH1D("FOM", "", nBins, lowestbin, highestbin) ##frac{S}{#sqrt{S+B}}

	for point in range(0, nBins+1): 
		#sigEff = effs.GetPointX(point)
		#bkgEff = effs.GetPointY(point)

		B = background.Integral(point, nBins+1)
		S = signal.Integral(point, nBins+1)

		Sigma = 0 if (B == 0) else S/math.sqrt(B)

		xval = signal.GetBinCenter(point)

		FOM.SetBinContent(point, Sigma)
		#FOM.SetBinError(point, error)




	#FOM = filemanager.GetItem("cutflowGen")

	fomcanvas = ROOT.TCanvas("fomcanvas", "fomcanvas", 800, 600)

	FOM.Draw("E SAME")

	#formatHistoCommon(histgwf)
	ROC.SetLineColor(ROOT.kBlue+2)
	ROC.SetLineWidth(2)
	ROC.GetXaxis().SetLabelSize(textsize)
	ROC.GetYaxis().SetLabelSize(textsize)

	FOM.SetLineColor(ROOT.kGreen)
	#FOM.SetLineWidth(2)
	FOM.SetMarkerStyle(7)
	FOM.SetMarkerColor(ROOT.kGreen)
	FOM.GetXaxis().SetLabelSize(textsize)
	FOM.GetYaxis().SetLabelSize(textsize)


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

	fomcanvas.Draw("NOSTAT")

	canvas.Print(outputfolder+quantity+".pdf")
	if (webpublication): canvas.Print(webfolder+quantity+".png")
	fomcanvas.Print(outputfolder+"/Fomstandalone"+".pdf")
	fomcanvas.SaveAs(outputfolder+"/fomstandalone.root")

	#f = ROOT.TFile.Open(presentationfolder, "Update")
	#fomcanvas.Write()
	#canvas.Write()
	#f.Write()
	#f.Close()




filemanager.CloseAll()



