import ROOT
import os
import math
from sklearn.metrics import roc_curve
import numpy as np


#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/libFunctions.C+")#ROOT.gROOT.LoadMacro("/eos/home-m/mhuwiler/plugins/libFunctions.C")
#ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/ExperimentSpecificLayer.C")
ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C")
ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/CMS/tdrstyle.C")
ROOT.gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
ROOT.setTDRStyle()
#import CMS_lumi





plotstats = False

webpublication =False


outputfolder = "./plots/MVAFiguresTauFLIsoLatest/"

os.system("mkdir -p "+outputfolder)

logscale = [""]

presentationfolder = "../Presentations/PresentationVFS_15_6_22/reveal.js-master/Figures.root"


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


def GetFom(sigeffs, bkgeffs, sigInSample=1., bkgInSample=1): 

	assert(len(sigeffs) == len(bkgeffs))
	numPoints = len(sigeffs)

	FOM = ROOT.TH1D("FOM", "", numPoints, -1., 1.)

	for point in range(0, numPoints): 
		sigEff = sigeffs[point]
		bkgEff = bkgeffs[point]
		
		print "Sig eff: {}, bkg eff: {}".format(sigEff, bkgEff)

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

	return FOM

def ScaleCutToPad(cut, hist, pad):
	leftmargin = pad.GetLeftMargin()
	rightmargin = pad.GetRightMargin()
	print leftmargin
	print rightmargin
	return leftmargin + (cut1 - hist.GetXaxis().GetXmin())/(hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin())*(1-leftmargin-rightmargin)


# Web publication
if (webpublication): 
	webfolder = "/eos/home-m/mhuwiler/www/Analysis/MVAinputVariables/"
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


if __name__ == "__main__": 

	filemanager = ROOT.FileManager()

	#ROOT.gROOT.ForceStyle()
	ROOT.gSystem.Load("Tau.h")


	filemanager.AddItem("signal", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v6.8/Sig_tauDNN_mva.root", "tree")
	#filemanager.AddItem("bkgnoresponse", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPH1-3Run2018B_converted.root", "tree")
	filemanager.AddItem("background", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v6.8/dataB2_tauDNN_mva.root", "tree")
	filemanager.AddItem("oldroc",  "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "xgboOptimized/scikit-like_ROC")
	#filemanager.AddItem("TMVAROC", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "xgboOptimized/TMVA-like_ROC")
	#filemanager.AddItem("bkgEff", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "xgboOptimized/bkgEff(sigEff)")

	filemanager.OpenAllItems()

	#from anaMVA.loadSamples import loadSamples
	#samples = loadSamples(featuresToLoad, numLoad)

	columns = ["mvaScore"]

	scheme = "original"


	#signal = [samples["signal"]]
	#background = [samples["background"]]

	signal = ROOT.RDataFrame(filemanager.GetItem("signal")).Filter("mvaScore>-1").AsNumpy(columns)
	background = ROOT.RDataFrame(filemanager.GetItem("background")).Filter("mvaScore>-1").AsNumpy(columns)

	print signal.keys()

	siglabels = signal["mvaScore"]

	print siglabels

	sigtruth = np.ones(len(siglabels))

	print sigtruth

	bkglabels = background["mvaScore"]
	bkgtruth = -np.ones(len(bkglabels))

	print bkglabels
	
	print bkgtruth

	labels = np.concatenate([siglabels, bkglabels])

	print labels

	truths = np.concatenate([sigtruth, bkgtruth])

	print "Starting to compute ROC curve... "
	bkgeff, sigeff, _ = roc_curve(truths, labels)
	print "Computed ROC curve. "


	canv = ROOT.TCanvas("canv", "canv", 800, 600)
	graph = ROOT.TGraph(len(bkgeff), np.asarray(bkgeff, "d"), np.asarray(sigeff, "d"))
	oldgraph = filemanager.GetItem("oldroc")
	oldgraph.Draw("AP")
	graph.Draw("SAME P")
	canv.Draw()
	graph.SetMarkerColor(ROOT.kBlue)
	oldgraph.SetMarkerColor(ROOT.kGreen)
	#oldgraph.SetLineColor(ROOT.kGreen)
	oldgraph.SetMarkerSize(1)
	oldgraph.SetMarkerStyle(8)

	legend = ROOT.TLegend(0.6, 0.2, 0.9, 0.3)
	legend.AddEntry(oldgraph, "Initial variables")
	legend.AddEntry(graph, "New variables")
	legend.SetBorderSize(0)
  	legend.SetFillColor(0)
  	legend.SetTextSize(0.04)
	legend.Draw()

	oldgraph.GetXaxis().SetTitleSize(0.04)
	oldgraph.GetYaxis().SetTitleSize(0.04)
	canv.Print(outputfolder+"NewRoc.pdf")


	sigeffcorr = np.sort(sigeff)
	bkgeffcorr = np.sort(bkgeff) #np.ones(len(bkgeff)) - np.sort(bkgeff)

	print sigeffcorr 
	print bkgeffcorr

	fom = GetFom(sigeffcorr, bkgeffcorr, 220, 1296-118) #626, 2698 updated values from Oct. 2022

	# Margins 
	leftmargin = 0.1
	rightmargin = 0.1
	bottommargin = 0.15
	topmargin = 0.1

	cut1 = 0.5
	cut2 = -0.2
	cut3 = -0.5

	ROOT.gStyle.SetPadTickY(0) # Disactivate axes on both sides

	if (scheme == "original"): 
		signalColor = ROOT.kBlue+1
		backgroundColor = ROOT.kRed
	elif (scheme == "match"): 
		signalColor = ROOT.kRed
		backgroundColor = ROOT.kBlue+1
	else: 
		print "Unknown schemme"

	fomcanvas = ROOT.TCanvas("fomcanvas", "fomcanvas", 800, 600)
	fomcanvas.SetMargin(leftmargin, rightmargin, bottommargin, topmargin); 
	signalDist = ROOT.RDataFrame(filemanager.GetItem("signal")).Filter("(mvaScore>-1)&&Dstar_match&&b_tau_match").Histo1D(("signalDist", "signalDist", 30, -1., 1.), "mvaScore") # TODO: use cut >= -1.
	bkgDist = ROOT.RDataFrame(filemanager.GetItem("background")).Filter("mvaScore>-1").Histo1D(("bkgDist", "bkgDist", 30, -1., 1.), "mvaScore")
	fomcanvas.cd()
	fompad = ROOT.TPad("fompad", "fompad", 0., 0., 1., 1.)
	fompad.SetMargin(leftmargin, rightmargin, bottommargin, topmargin); 
	fompad.SetBorderSize(0)
	fompad.SetFrameLineWidth(0)
	fompad.SetFrameBorderMode(0)
	fompad.SetLogy()
	fompad.cd()
	fom.Draw("E")
	fom.GetXaxis().SetRangeUser(-1., 1.)
	fomcanvas.cd()
	fompad.Draw()
	fomcanvas.cd()
	pad = ROOT.TPad("pad", "pad", 0., 0., 1., 1.)
	pad.SetMargin(leftmargin, rightmargin, bottommargin, topmargin); 
	pad.SetBorderSize(0)
	pad.SetFrameLineWidth(0)
	pad.SetFrameBorderMode(0)
	pad.SetBorderMode(0)
	pad.cd()
	pad.SetFillColorAlpha(ROOT.kWhite, 0.); 
	signalDist.Draw("HIST Y+") #"SAME"
	signalDist.Scale(5./signalDist.Integral()) #1./signalDist.Integral()
	signalDist.SetMarkerColor(signalColor)
	signalDist.SetLineColor(signalColor)
	signalDist.SetLineWidth(2)
	bkgDist.Draw("HIST SAME Y+")
	bkgDist.Scale(5./bkgDist.Integral())
	bkgDist.SetMarkerColor(backgroundColor)
	bkgDist.SetLineColor(backgroundColor)
	bkgDist.SetLineWidth(2)
	signalDist.GetXaxis().SetRangeUser(-1., 1.)
	fomcanvas.cd()
	pad.SetLogy()
	pad.Draw()
	fom.SetMarkerSize(0.2)
	fom.SetLineColor(ROOT.kGreen+2)
	#fom.SetMarkerColor(ROOT.kGreen+3)
	#fom.SetMarkerColor(ROOT.kGreen)
	#fom.SetMarkerStyle(1)
	#fom.SetMarkerSize(2)
	fomcanvas.Draw()

	fomlegend = ROOT.TLegend(0.17, 0.7, 0.37, 0.85) #0.125, 0.7, 0.305, 0.85
	fomlegend.AddEntry(fom, "significance")
	fomlegend.AddEntry(signalDist.GetPtr(), "signal")
	fomlegend.AddEntry(bkgDist.GetPtr(), "background")
	fomlegend.SetBorderSize(0)
  	fomlegend.SetFillColor(0)
  	fomlegend.SetTextSize(0.03)
	fomlegend.Draw()

	frac = (cut1 - fom.GetXaxis().GetXmin())/(fom.GetXaxis().GetXmax() - fom.GetXaxis().GetXmin())
	factor = (1-leftmargin-rightmargin)
	value1 = leftmargin + (cut1 - fom.GetXaxis().GetXmin())/(fom.GetXaxis().GetXmax() - fom.GetXaxis().GetXmin())*(1-leftmargin-rightmargin)
	value2 = leftmargin + (cut2 - fom.GetXaxis().GetXmin())/(fom.GetXaxis().GetXmax() - fom.GetXaxis().GetXmin())*(1-leftmargin-rightmargin)
	value3 = leftmargin + (cut3 - fom.GetXaxis().GetXmin())/(fom.GetXaxis().GetXmax() - fom.GetXaxis().GetXmin())*(1-leftmargin-rightmargin)
	print fom.GetXaxis().GetXmax()
	print fom.GetXaxis().GetXmin()
	print frac
	print factor
	print frac*factor
	print value1
	line1 = ROOT.TLine(value1, bottommargin, value1, 1-topmargin)
	#line1.SetLineStyle(7)
	line1.SetLineColor(ROOT.kGray+2)
	line1.Draw()

	line2 = ROOT.TLine(value2, bottommargin, value2, 1-topmargin)
	line2.SetLineColor(ROOT.kGray+2)
	line2.Draw()

	line3 = ROOT.TLine(value3, bottommargin, value3, 1-topmargin)
	line3.SetLineColor(ROOT.kGray+2)
	#line3.SetLineStyle(7)
	line3.Draw()

	label1 = ROOT.TPaveText(0.33, 0.4, 0.4, 0.6)
	label1.AddText("control")
	label1.AddText("region")
	label1.AddText("(CR)")
	label1.SetBorderSize(0)
  	label1.SetFillColor(0)
  	label1.SetFillColorAlpha(ROOT.kWhite, 0.); 
  	label1.SetTextSize(0.04)
  	label1.SetTextColor(ROOT.kGray+2)
	label1.Draw()

	label2 = ROOT.TPaveText(0.5, 0.5, 0.6, 0.55)
	label2.AddText("sideband (SB)")
	#label2.AddText("region")
	label2.SetBorderSize(0)
  	label2.SetFillColor(0)
  	label2.SetFillColorAlpha(ROOT.kWhite, 0.); 
  	label2.SetTextSize(0.04)
  	label2.SetTextColor(ROOT.kGray+2)
	label2.Draw()

	label3 = ROOT.TPaveText(0.7, 0.6, 0.86, 0.7)  
	text = label3.AddText("signal")
	label3.AddText("region (SR)")
	#label3.AddText("signal")
	#label3.AddText("region")
	#label3.SetAllWith("Angle", 90)
	#text.SetTextAngle(90)
	label3.SetBorderSize(0)
  	label3.SetFillColor(0)
  	label3.SetFillColorAlpha(ROOT.kWhite, 0.); 
  	label3.SetTextSize(0.04)
  	label3.SetTextColor(ROOT.kGray+2)
	label3.Draw()
	fomlegend.Draw()

	#label3 = ROOT.TPaveText(0.83, 0.7, 0.93, 0.8)  
	#text = label3.AddText("signal region")
	##label3.AddText("signal")
	##label3.AddText("region")
	##label3.SetAllWith("Angle", 90)
	#text.SetTextAngle(90)
	#label3.SetBorderSize(0)
 	#label3.SetFillColor(0)
 	#label3.SetFillColorAlpha(ROOT.kWhite, 0.); 
 	#label3.SetTextSize(0.04)
 	#label3.SetTextColor(ROOT.kGray+2)
	#label3.Draw()
	textsize = 0.05
	fom.GetXaxis().SetTitle("BDT score")
	fom.GetXaxis().SetTitleSize(textsize)
	fom.GetYaxis().SetTitle("significance")
	fom.GetYaxis().SetTitleSize(textsize)
	fom.GetXaxis().SetTitleOffset(1.2)
	signalDist.GetYaxis().SetTitle("N events")
	signalDist.GetYaxis().SetTitleSize(textsize)
	signalDist.GetYaxis().SetTitleOffset(0.8)
	fom.GetYaxis().SetTitleOffset(1.)

	fomcanvas.Update()
	fomcanvas.Print(outputfolder+"Fom.pdf")






