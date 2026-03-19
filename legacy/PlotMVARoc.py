#!/usr/bin/env python

from __future__ import division, print_function


import os
import anaConfig
from ROOT import Ana, RDataFrame, TCanvas, kBlue, kGreen
from libMLTools import GetROC, GetFom
from libUtils import HoldUntilKeyPress



from argparse import ArgumentParser

parser = ArgumentParser(description="GetEfficiency")
#parser.add_argument("filename", action="store", type=str, default="", help="Name of file")
#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v6.8", help="Which version (cycle) of files to run on")
parser.add_argument("-g", "--cut", dest="cut", action="store", type=str, default="1", help="Custom cut to be included in eff calculation")
parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="./plots/testroc/", help="Path for saving plots")
parser.add_argument('-l', "--variable", dest="variable", action="store", default="mvaScore", help="Variable to reprocess")
parser.add_argument("-n", "--target", dest="target", action="store", type=float, default=10000., help="Target number of events after selection")

options = parser.parse_args()


samples = anaConfig.MVAsamples

		
Ana.Init(options.version)

for file in samples: 
	Ana.filemanager.OpenItem(file)

os.system("mkdir -p {}".format(options.out))



# Loading signal and background 
cutsig = (Ana.cut["base"]+Ana.samples.at("Sig").cut).GetTitle()
cutbkg = (Ana.cut["base"]+Ana.samples.at("dataB2").cut).GetTitle()

signal = RDataFrame(Ana.filemanager.GetItem("Sig")).Filter(cutsig)
background = RDataFrame(Ana.filemanager.GetItem("dataD1")).Filter(cutbkg)

print("Starting to compute ROC curve... ")
roc, auc = GetROC(signal, background, options.variable)
print("Computed ROC curve. ")

canv = TCanvas("canv", "canv", 800, 600)
roc.Draw("AP")
canv.Draw()
roc.SetMarkerColor(kBlue)
#roc.SetMarkerSize(1)
#roc.SetMarkerStyle(8)
canv.Draw()
canv.Print(options.out+"ROC.pdf")
print("Area under curve (A.U.C.): {}".format(auc))

HoldUntilKeyPress()

fom, maxsig, cutvalue = GetFom(signal, background, options.variable)
canv2 = TCanvas("fom", "fom", 800, 600)
fom.Draw("E")
canv2.Draw()
fom.GetXaxis().SetRangeUser(-1., 1.)
#fom.SetMarkerColor(ROOT.kGreen+4)
fom.SetMarkerStyle(8)
fom.SetMarkerSize(0.2)
fom.SetLineColor(kGreen+2)
#roc.SetMarkerSize(1)
#roc.SetMarkerStyle(8)
canv2.Draw()
canv2.Print(options.out+"FOM.pdf")
print("Maximum significance of {} with cut at value {}".format(maxsig, cutvalue))

HoldUntilKeyPress()


Ana.filemanager.CloseAll()

