import pandas as pd
from anaMVA.loadSamples import loadSamples
import ROOT
from anaMVA.BDTvariables import features


# TODO: put this in an extra file
# features = [ "BsDstarTauNu_D0_pt", "BsDstarTauNu_D0_eta", "BsDstarTauNu_D0_phi", "BsDstarTauNu_D0_vprob", "BsDstarTauNu_D0_fl3d", "BsDstarTauNu_D0_fls3d",
# 	"BsDstarTauNu_Ds_pt", "BsDstarTauNu_Ds_eta", "BsDstarTauNu_Ds_phi", "BsDstarTauNu_Ds_vprob", "BsDstarTauNu_Ds_fl3d", "BsDstarTauNu_Ds_fls3d",
# 	"BsDstarTauNu_D0_lip", "BsDstarTauNu_D0_lips", "BsDstarTauNu_D0_pvip", "BsDstarTauNu_Ds_lip", "BsDstarTauNu_Ds_lips", "BsDstarTauNu_Ds_pvip", 
# 	"BsDstarTauNu_tau_pt", "BsDstarTauNu_tau_eta", "BsDstarTauNu_tau_phi", "BsDstarTauNu_tau_fl3d", "BsDstarTauNu_tau_fls3d", "BsDstarTauNu_tau_vprob", 
# 	"BsDstarTauNu_tau_lip", "BsDstarTauNu_tau_lips", "BsDstarTauNu_tau_pvip", "BsDstarTauNu_tau_pvips", "BsDstarTauNu_tau_alpha", "BsDstarTauNu_tau_max_dr_3prong", 
# 	"BsDstarTauNu_tau_pi1_pt", "BsDstarTauNu_tau_pi1_eta", "BsDstarTauNu_tau_pi1_phi", 
# 	"BsDstarTauNu_tau_pi2_pt", "BsDstarTauNu_tau_pi2_eta", "BsDstarTauNu_tau_pi2_phi", 
# 	"BsDstarTauNu_tau_pi3_pt", "BsDstarTauNu_tau_pi3_eta", "BsDstarTauNu_tau_pi3_phi", 
# 	"BsDstarTauNu_tau_sumofdnn",  
# 	 ]

samples = loadSamples(features, "v6.8", 10000)

ROOT.gStyle.SetOptStat(0)

signal = samples["signal"]
background = samples["background"]

sigcorrelations = signal.corr(method="pearson")

print sigcorrelations

# Now we do havve the correlation matrix, let's plot it
Nx = len(sigcorrelations.columns)
Ny = len(sigcorrelations)

print sigcorrelations.iloc[1, 1]

histo = ROOT.TH2D("signalcorrelations", "Signal correlations", Nx, 0, Nx, Ny, 0, Ny)

labels = list(sigcorrelations.columns.values)
print labels
labelsy = list(sigcorrelations.index.values)
print labelsy

for i in range(0, Nx): 
	for j in range(0, Ny): 
		entry = sigcorrelations.iloc[i, j]
		histo.SetBinContent(i, j, entry)

for i, label in enumerate(labels): 
	label = label.replace("BsDstarTauNu_", "")
	histo.GetXaxis().SetBinLabel(i+1, label)

for i, label in enumerate(labelsy): 
	label = label.replace("BsDstarTauNu_", "")
	histo.GetYaxis().SetBinLabel(i+1, label)

histo.GetXaxis().LabelsOption("v")


canvas = ROOT.TCanvas("canvas", "canvas", 800, 600)
histo.Draw("COLZ")
canvas.SetLeftMargin(0.15)
canvas.SetBottomMargin(0.15)
canvas.Draw()
canvas.Print("plots/CorrelationsMVA.pdf")



