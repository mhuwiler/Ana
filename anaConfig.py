import ROOT


import os
environment=""
try: 
	environment = os.environ["PLATFORM"]
except: 
	pass
if (environment == "lxplus"): 
	ROOT.gROOT.LoadMacro("FileFlow.h+")
else:
	ROOT.gROOT.LoadMacro("FileFlow.h+")
#ROOT.gROOT.LoadMacro("FileFlow.h+")
from ROOT import Ana


data = "data" #"dataB2"

dataWS = "dataDWS"

Sig = "Sig"

samples = [Sig, data, "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "B0toDstar3pi"] ##[Sig, data, "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "B0toDstar3pi", "B0toDstar3pipi0", "B0toDstar5pi", "B0toDstarDs1", "B0toDstarDs0star", "B0toDstarD0Kstar", "BstoDD", "B0toDstara1"] #, dataWS #"SigPart", "dataD2WS", "dataD2TauWS", , "B0toDstar5pi" ["Sig", "dataD1", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "B0toDstarD0K", "BkgDstara1", "B0toDstar3pi", "dataD2WS"]

fitsamples = [Sig, "data_obs", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarD", "ButoDstarDK", "ABCD"] #"B0toDstara1"

tauIDvar = "globalParT3" #globalParT3_Xtauhtaum

regions = ["SR", "CR", "SB"]

variables = ["b_B_q2", "b_B_m", "b_tau_m", "b_B_mm2"] #"b_tau_rhomass1", "b_tau_rhomass2", 

fitvariables = ["b_B_m"]

MVAsamples = ["Sig", "dataD1"] #, "dataA1"

blacklist = ROOT.vector("std::string")(["v_taucandidates", "b_tau"])

from anaPrepareRegions import ABCDcuts

ABCDconfig = ABCDcuts("(b_tau_sumdnn>2.)", "(b_tau_sumdnn<1.5)", "b_B_nmu<1&&b_B_ne<1&&b_B_nh<1", "b_B_nmu>1||b_B_ne>1||b_B_nh>1")

categories = {
	'SR': [(1, 'SR')],
	'CR': [(2, 'CR')],
	'SB': [(3, 'SB')],
	'baseline': [(4, 'baseline')],
	'all': [(5, 'all')],
 }


def OpenFiles(version): 
	Ana.Init(version)
	for file in samples: 
		Ana.filemanager.OpenItem(file)

def CloseFiles(): 
	Ana.filemanager.CloseAll()

def DebugMode(): 
	global samples
	samples = ["Sig", "B0toDstarDs", "B0toDstarDsstar", "dataD1"] # Using a reduced set of files for debugging
	global regions
	regions = ["SR"]
	global variables
	variables = ["b_B_m"]
	global data
	data = "dataD1"
	import libMLTools
	libMLTools.debugmode = True

def Denominator(): 
	global Sig
	Sig = "B0toDstarrho0pi"
	global samples
	global data
	global fitsamples
	data = "dataD"
	samples = [data, Sig, "B0toDstar3pi", "B0toDstara1", "B0toDstar5pi", "B0toDstar3pipi0", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarKpipi", "B0toDstarKpipi", "B0toDstarKKstar", "ButoDstarpipipipi0", "ButoDstarpipipi"]
	fitsamples = [Sig, "data_obs", "B0toDstar3pi", "B0toDstara1", "B0toDstar5pi", "B0toDstar3pipi0", "B0toDstarDs", "B0toDstarDsstar", "B0toDstarKpipi", "B0toDstarKpipi", "B0toDstarKKstar", "ButoDstarpipipipi0", "ButoDstarpipipi"]
	import anaPlotting 
	anaPlotting.normalisebinwidth = True

def UpdateColorsDenom(): 
	colors = Ana.mycolors
	Ana.samplelist["B0toDstara1"].color = mycolors[2]
	Ana.samplelist["B0toDstarrho0pi"].color = mycolors[3]
