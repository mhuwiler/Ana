#!/usr/bin/env python

import CombineHarvester.CombineTools.ch as ch
#import CombineHarvester.CombinePdfs.morphing as morphing
from ROOT import TFile, TH1F, Double_t
import os
import sys

print sys.argv


# TODO: Add argument parsing 

from collections import namedtuple

DummyOptions = namedtuple('DummyOptions', ['debug'])

options = DummyOptions(True)

shape_file = 'workspaceFromExportTest.root'
HFSys = [] #"ABCD-sys-HFDown","ABCD-sys-HFUp"
#HighNchSys = "highNch_" + str(sys.argv[2])

file = TFile(shape_file)

cb = ch.CombineHarvester()
if (options.debug): cb.SetVerbosity(3)
cb.SetVerbosity(3)

sig_procs = ["Sig"] #['Sig']

bkg_procs = ["B0toDstarDs", "B0toDstarDsstar", "B0toDstarD0K"]

categories = {
    'SR': [(1, 'SR')],
    'CR': [(2, 'CR')],
    'SB': [(3, 'SB')],
    }


channels = ['SB'] # Channels for which to write datacards
prefix = ['fitFirst'] # also called analysis
era = ['2018']


#if int(sys.argv[1]) > 0:
#  bkg_procs = [HFSys[int(sys.argv[1])-1]]

#if int(sys.argv[2]) > 0:
#  bkg_procs = [HighNchSys]


for chn in channels:

    cb.AddObservations(['*'], prefix, era, [chn], categories[chn])

    cb.AddProcesses(['120'], prefix, era, [chn], sig_procs, categories[chn], True)

    cb.AddProcesses(['*'], prefix, era, [chn], bkg_procs, categories[chn], False)

    

cb.cp().process(sig_procs).AddSyst(cb, 'CMS_lumi', 'lnN', ch.SystMap()(1.05))
print "CMS_lumi added."

"""
print '>> Adding systematic uncertainties...'


cb.cp().process(sig_procs).AddSyst(    cb, 'CMS_lumi', 'lnN', ch.SystMap()(1.05))
print "CMS_lumi added."

cb.cp().process(sig_procs).AddSyst(    cb, 'Tau_BR', 'lnN', ch.SystMap()(1.006))
print "Tau_BR added."

cb.cp().process(sig_procs).AddSyst(    cb, 'MC_size_for_eff', 'lnN', ch.SystMap()(1.01062))
print "MC_size_for_eff added."

# Adding pion SF uncertainty
eta_bins = [0, 0.8, 1.5]
pt_bins = [0, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1]
pionSFgroup = []
for ebin in range(3):
  for pbin in range(14):
    cb.cp().process(sig_procs).AddSyst(    cb, 'pionSF_pt'+str(int(10*pt_bins[pbin]))+'_eta'+str(int(10*eta_bins[ebin])), 'shape', ch.SystMap()(1.0))
    pionSFgroup.append('pionSF_pt'+str(int(10*pt_bins[pbin]))+'_eta'+str(int(10*eta_bins[ebin])))
#cb.cp().process(sig_procs).AddSyst(    cb, 'pion-SF', 'lnN', ch.SystMap()(1.01))
print "pion SF uncertainty added."


cb.cp().process(sig_procs).AddSyst(    cb, 'muon-SF', 'shape', ch.SystMap()(1.0))
#cb.cp().process(sig_procs).AddSyst(    cb, 'muon-SF', 'lnN', ch.SystMap()(1.01))
print "muon-SF added."

if int(sys.argv[1]) == 0:
  cb.cp().process(bkg_procs).AddSyst(    cb, 'ABCD-sys-HF', 'shape', ch.SystMap()(1.0))
  print "ABCD-sys-HF added."
if int(sys.argv[2]) == 0:
  cb.cp().process(bkg_procs).AddSyst(    cb, 'ABCD-sys-nch', 'shape', ch.SystMap()(1.0))
  print "ABCD-sys-Nch added."

#cb.cp().process(bkg_procs).AddSyst(
#    cb, 'CMS_bkg', 'lnN', ch.SystMap()(1.1))

#cb.cp().process(bkg_procs).AddSyst(cb,'rbkg','rateParam',ch.SystMap()(14.0))
"""

print '>> Extracting histograms from input root files...'

for chn in channels:
    cb.cp().channel([chn]).ExtractShapes(
        '%s' % (shape_file),
#        '$BIN/$PROCESS', '$BIN/$PROCESS_$SYSTEMATIC')
        'w:$PROCESS_$BIN', 'w:$PROCESS_$BIN_$SYSTEMATIC') #, '$BIN/$SYSTEMATIC' 'w:$PROCESS_$BIN_$SYSTEMATIC'


"""
cb.SetGroup('syst', ['.*'])
#cb.SetGroup('lumi', ['CMS_lumi'])
cb.SetGroup('pionSF', pionSFgroup)
cb.SetGroup('muonSF', ['muon-SF'])
cb.SetGroup('ABCD_HF', ['ABCD-sys-HF'])
cb.SetGroup('ABCD_nch', ['ABCD-sys-nch'])
#cb.RemoveGroup('syst', ['CMS_lumi'])
"""

print '>> Setting standardised bin names...'
ch.SetStandardBinNames(cb)
cb.PrintAll()

writer = ch.CardWriter('$TAG/$ANALYSIS_$CHANNEL_$BINID_$ERA_$MASS.txt',
                       '$TAG/common/$ANALYSIS_$CHANNEL_$BINID_$ERA_$MASS.input.root')

writer.SetVerbosity(1)

outdir = 'output/sm_cards/LIMITS'

for chn in channels:  # plus a subdir per channel
    print 'writing', chn, cb.cp().channel([chn])
    writer.WriteCards(outdir, cb.cp().channel([chn]))


print '>> Done!'


outcard = outdir + '/ggtautau_3prong_1_2015_120.txt'

#os.system(command)

# overwrite extra rateParam
if os.path.isfile(outcard):

    f = open(outcard, 'a')
#    f.write(extraStr)
    f.write('* autoMCStats 0 1\n')
    f.close()

command = 'text2workspace.py ' + outcard + ' -o ' + outdir + '/workspace.root -m 120'
os.system(command)
