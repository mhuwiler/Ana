#!/usr/bin/env python
from __future__ import division, print_function

import CombineHarvester.CombineTools.ch as ch
import anaConfig
from ROOT import TFile, TH1F, Double_t
import os
import sys
import copy as cp
from argparse import ArgumentParser



if __name__ == "__main__":

    parser = ArgumentParser(description="Model") 
    #parser.add_argument("tool", action="store", type=str, help="Which time list you want to analyse")
    parser.add_argument("-c", "--version", dest="version", action="store", type=str, default="v7", help="Which version (cycle) of files to run on")
    parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="fitmodel", help="Which version (cycle) of files to run on")
    parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
    parser.add_argument('-b', "--batch", dest="batch", action="store_true", default=False, help="Run in batch mode")
    parser.add_argument('-d', "--denom", dest="denom", action="store_true", default=False, help="Denominator analysis")
    parser.add_argument("-t", "--tag", dest="tag", action="store", type=str, default="fitFirst", help="Tag fir this version of fit")
    parser.add_argument("-w", "--workspace", dest="workspace", action="store", type=str, default="workspaceFromExportHists.root", help="Name of the workspace file with the distributions (and fit model)")
    

    options = parser.parse_args()

    
    if (options.batch): 
        ROOT.gROOT.SetBatch(1) 


    #Ana.Init(options.version, options.denom)
    print("Starting datacard creation")


    variables = anaConfig.fitvariables


    os.system("mkdir -p {}".format(options.out))


    cb = ch.CombineHarvester()
    cb.SetVerbosity(3)

    sig_procs = [anaConfig.Sig] #['Sig']

    allsamples = cp.deepcopy(anaConfig.samples)
    allsamples.remove(anaConfig.Sig)
    allsamples.remove(anaConfig.data)
    bkg_procs = allsamples

    categories = anaConfig.categories


    channels = ['baseline'] # Channels for which to write datacards
    prefix = 'fitFirst' # also called analysis
    era = '2018'


    for chn in channels:

        cb.AddObservations(['*'], [prefix], [era], [chn], categories[chn])

        cb.AddProcesses(sig_procs, [prefix], [era], [chn], sig_procs, categories[chn], True)

        cb.AddProcesses(['*'], [prefix], [era], [chn], bkg_procs, categories[chn], False)

        
    # Systematics
    cb.cp().process(sig_procs).AddSyst(cb, 'CMS_lumi', 'lnN', ch.SystMap()(1.05))
    print("CMS_lumi added.")



    print('>> Extracting histograms from input root files...')

    # for chn in channels:
    #     cb.cp().channel([chn]).ExtractShapes(
    #         '%s' % (file),
    # #        '$BIN/$PROCESS', '$BIN/$PROCESS_$SYSTEMATIC')
    #         'w:$PROCESS', 'w:$PROCESS_$BIN_$SYSTEMATIC') #, '$BIN/$SYSTEMATIC' 'w:$PROCESS_$BIN_$SYSTEMATIC'

    cb.cp().backgrounds().ExtractShapes(
        options.workspace, 
        "$BIN/$PROCESS", "");
    cb.cp().signals().ExtractShapes(
        options.workspace, 
        "$BIN/$PROCESS", "");



    print('>> Setting standardised bin names...')
    ch.SetStandardBinNames(cb)
    cb.PrintAll()


    writer = ch.CardWriter('{}/$ANALYSIS_$CHANNEL_$BINID.txt'.format(options.out),
                           '{}/$ANALYSIS_$CHANNEL_$BINID.root'.format(options.out))

    writer.SetVerbosity(1)

    #options.out = 'output/sm_cards/LIMITS'

    for chn in channels:  # plus a subdir per channel
        print('writing', chn, cb.cp().channel([chn]))
        writer.WriteCards(options.out, cb.cp().channel([chn]))


    print('>> Done!')


#os.system(command)

# overwrite extra rateParam
# for channel in channels: 
#     print(categories[channel])
#     print(categories[channel][0])
#     print(categories[channel][0][0])
#     outcard = options.out + "/{}_{}_{}.txt".format(prefix, channel, categories[channel][0][0])
#     print(outcard)
#     if os.path.isfile(outcard):

#         f = open(outcard, 'a')
#     #    f.write(extraStr)
#         f.write('* autoMCStats 0 1\n')
#         f.close()

#     command = 'text2workspace.py ' + outcard + ' -o ' + options.out + '/workspace{}.root -m 120'.format(channel)
#     os.system(command)
