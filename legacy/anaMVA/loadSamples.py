import os
import pandas as pd
try: 
	import uproot3 as uproot #if the latest uproot is installed we want to use a previous version
except: 
	import uproot
import ROOT
ROOT.gROOT.LoadMacro("../FileFlow.h")



def loadSamples(features, version, eventFraction=-1) : 


	path = "/eos/home-m/mhuwiler/data/Analysis/v5/" #"../localdata/" #"/pnfs/psi.ch/cms/trivcat/store/user/mhuwiler/Analysis/submit/productionProductionFastApplyDnnNoDNNdebug/" #"/eos/home-m/mhuwiler/data/Analysis/v5/"


	if(eventFraction<0):
		numEvents=-1
	else : 
		numEvents=eventFraction

	print "Number of events: {}".format(numEvents)

	# Poweg to be used if available

	samples={}

	ROOT.Ana.Init(version)
	from ROOT.Ana import filemanager

	#samples["signal"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents) #Dict-like structure 
	#samples["background"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPH1-3Run2018B_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents)

	#samples["signal"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_tauDNN.root"))["ntuplizer/tree"].pandas.df(branches=features, entrystop=numEvents) #Dict-like structure 
	#samples["background"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/prod2018BFirst_tauDNN.root"))["ntuplizer/tree"].pandas.df(branches=features, entrystop=numEvents)
	#samples["backgroundDs"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/firstDstarDsMultipleTau_tauDNN.root"))["ntuplizer/tree"].pandas.df(branches=features, entrystop=numEvents)

	#samples["signal"] = uproot.open(os.path.expandvars(filemanager.GetFile("SigTest_DNN")))["ntuplizer/tree"].pandas.df(branches=features, entrystop=numEvents) #Dict-like structure #"/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v4/SigTest_tauDNN.root"
	#samples["background"] = uproot.open(os.path.expandvars(filemanager.GetFile("dataD1WS_DNN")))["ntuplizer/tree"].pandas.df(branches=features, entrystop=numEvents) #"/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v4/ParkingBPHULD1_tauDNN.root"

	#samples["signalfloatgendstar"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/PrivateProductionGenDstar_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents)
	#samples["signalfloat"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents)
	#samples["backgroundfloat"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPH1-3Run2018B_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents)
	#samples["backgroundfloatprevious"] = uproot.open(os.path.expandvars("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge_converted.root"))["tree"].pandas.df(branches=features, entrystop=numEvents)

	from ROOT import Ana

	trainingsamples = { "signal": "SigTest_DNN", "background": "dataB_DNN"}

	for key, item in trainingsamples.iteritems(): 
		filemanager.OpenItem(item)
	
	samples["signal"] = pd.DataFrame(ROOT.RDataFrame(filemanager.GetItem(trainingsamples["signal"])).Filter((Ana.samples.at("Sig").cut+Ana.cut["training"]).GetTitle()).Range(int(numEvents)).AsNumpy(features))
	samples["background"] = pd.DataFrame(ROOT.RDataFrame(Ana.GetSample(trainingsamples["background"])).Filter((Ana.samples.at("data").cut+Ana.cut["training"]).GetTitle()).Range(int(numEvents)).AsNumpy(features))



	return samples


