import os, copy, ctypes
from ROOT import gStyle, gROOT, TCanvas, TLegend, TH1F, TH1D, gROOT, TFile, TColor, THStack, Double_t
#from common.officialStyle import officialStyle
#from common.DisplayManager_postfit import DisplayManager
#from common.DataMCPlot import *

gROOT.LoadMacro("/Users/mhuwiler/coding/plugins/Drawing/RatioCanvas.h")
from ROOT import RatioCanvas


gROOT.SetBatch(True)
gStyle.SetOptStat(False)
#gROOT.SetBatch(False)
#officialStyle(gStyle)
gStyle.SetOptTitle(0)
textsize = 0.04

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)


def AtomicDraw(histo, name, options = ""): 
    canv = TCanvas("canv", "canv", 800, 600)
    histo.DrawCopy(options)
    canv.Draw()
    canv.Print(name)


def ConvertGraphToHist(graph): 
    h = TH1D(graph.GetName()+"_converted", graph.GetTitle(), graph.GetN(), graph.GetXaxis().GetXmin(), graph.GetXaxis().GetXmax())
    err = TH1D(graph.GetName()+"_err", graph.GetTitle(), graph.GetN(), graph.GetXaxis().GetXmin(), graph.GetXaxis().GetXmax())
    for i in range(0, graph.GetN()): 
        x = ctypes.c_double()
        y = ctypes.c_double()
        graph.GetPoint(i, x, y);
        h.Fill(x,y)
        error = graph.GetErrorY(i) # Returns the average of high and low
        err.SetBinContent(i, error)
    return h , err


#shape_file = '/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/anal/dev/datacard_MUSF_blind/tau_rhomass_unrolled_coarse_new.root'
#shape_file = '/work/ytakahas/work/analysis/CMSSW_10_2_10/src/rJpsi/anal/dev/datacard_MUSF_blind/tau_rhomass_unrolled_new.root'
#file_shape = TFile(shape_file)
#data_sb = file_shape.Get('sb/data_obs')
#nbin = data_sb.GetXaxis().GetNbins()

# Color blind friendly palette
# https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=7 
mycolors = [TColor(215,48,39), TColor(252,141,89), TColor(254,224,144), TColor(255,255,191), TColor(224,243,248), TColor(145,191,219), TColor(69,117,180)]
#[red, oragnge, dark yellow, yellow, pale blue, blue, darker blue] # TODO: make an enum for them
colors = []
for color in mycolors: 
    colors.append(color.GetNumber())



filename = 'FitRegularNoSigCR.root'

file = TFile(filename)

print file

process = {
    'data':{'name':'data'},
    'total_signal':{'name':'SigPart'},
    'total_background':{'name':'total_background'}, 
    'comb_bkg':{'name':'bkg'},
    'part_bkg':{'name':'BkgDstarDs'},
}

#for ibin in range(1, nbin+1):
#    process['bg_bin' + str(ibin)] = {'name':'bg_bin' + str(ibin)}


ensureDir('plots/combine/')

#nbin = 20

for ftype in ['fit_s', 'fit_b']: #'prefit'
#for ftype in ['prefit', 'fit_b']:


    for cr in ['CR', 'SB']:
#    for cr in ['sr', 'sb']:
#    for cr in ['rJpsi_sr_1_2018']:

        hists = {}
    
        ymax = -1

        
        for ii, var in process.iteritems():
    
            hist = file.Get('shapes_' + ftype + '/' + cr + '/' + var['name'])

            if ((ii.find('signal')!=-1) and (cr == "SB")):
                hist = file.Get('shapes_' + ftype + '/' + "CR" + '/' + var['name'])
                hist.Scale(0.)
            
            print hist, 'shapes_' + ftype + '/' + cr + '/' + var['name']

            hist.SetFillStyle(0)
#            print file
            
            nbin = hist.GetXaxis().GetNbins()
#            print nbin
            
#            if var['name'].find('bin')!=-1:
#                hist.SetLineColor(2)
#                hist.SetLineStyle(2)

            if ii.find('signal')!=-1:
                hist.SetLineColor(2)
                #hist.SetLineStyle(2)
            elif ii.find('data')!=-1:
                hist.SetMarkerStyle(20)
                hist.SetMarkerSize(1)

            hist.SetLineWidth(3)
            #hist.Sumw2()
        
            hists[ii] = copy.deepcopy(hist)

            print hist.GetMaximum()

            if ymax < hist.GetMaximum()*1.1: ymax = hist.GetMaximum()*1.1

            
        name = 'canvas_' + ftype + '_' + cr
        canvas = RatioCanvas(name, name, 800, 600)
#        canvas.SetLogy()
        pad = canvas.Upper()
        legend = TLegend(  pad.GetLeftMargin()+0.45, # TODO: Implement GetMargin in RatioCanvas
                                1-pad.GetTopMargin()-.25, 
                                pad.GetLeftMargin()+(1.-(pad.GetLeftMargin()+pad.GetRightMargin())),
                                1-pad.GetTopMargin())

        frame = TH1F('frame_' + ftype + '_' + cr, 'fname_' + ftype + '_' + cr, nbin, 15., 100.)
        frame.GetXaxis().SetTitle('Tau rhomasses unrolled bin ID')
        frame.GetYaxis().SetTitle('Events')
  
#        if ftype.find('prefit')!=-1:
        frame.SetMaximum(ymax*2.)
        frame.SetMinimum(0.)
#        else:
#            frame.SetMaximum(ymax*8)
#            frame.SetMinimum(0.8)

        frame.Draw()
        
        hs = copy.deepcopy(hists['total_background'])
        hs.Add(copy.deepcopy(hists['total_signal']))
        hs.SetFillStyle(1)
        hs.SetFillColor(colors[0])
        hs.SetLineColor(colors[0])
        hs.Draw('hsame')
        legend.AddEntry(hists['data'], "data", "P")
        legend.AddEntry(hs, "Signal", "F")

        hists['total_background'].SetFillStyle(1)
        hists['total_background'].SetFillColor(colors[2])
        hists['total_background'].SetLineColor(colors[2])
        hists['total_background'].Draw('hsame')
        legend.AddEntry(hists['total_background'], "Background from estimate", "F")
        
        hists['part_bkg'].SetFillStyle(1)
        hists['part_bkg'].SetFillColor(colors[6])
        hists['part_bkg'].SetLineColor(colors[6])
        hists['part_bkg'].Draw('hsame')
        legend.AddEntry(hists['part_bkg'], "B#rightarrowD*D_{s}", "F")
        
        #hists['total_signal'].SetLineColor(colors[0])
        hists['total_signal'].Draw('hsame')
        hists['data'].Draw('epzsame')
        legend.Draw()

        #canvas.Upper().RedrawAxis()                # TODO: Implement RedraAxis in RatioCanvas

        canvas.Lower()
        #stack = THStack("stack", "Pulls")
        #stack.Add(hists["signal"])
        #stack.Add(hists["comb_bkg"])
        #stack.Add(hists["part_bkg"])
        #ratio = copy.deepcopy(hists["data"].GetHistogram())
        ratio, err = ConvertGraphToHist(hists["data"])
        ratio.Sumw2()
        histfit = copy.deepcopy(hs)
        histfit.Sumw2()
        ratio.Add(histfit, -1.)
        ratio.Divide(err)
        ratio.SetMarkerStyle(8)
        ratio.SetMarkerSize(0.5)
        ratio.SetMarkerColor(1)
        ratio.SetLineColor(1)
        ratio.GetYaxis().SetRangeUser(-10., 10.)
        ratio.GetYaxis().SetNdivisions(5)
        #AtomicDraw(ratio, "plots/combine/ratio.png")
        canvas.Lower()
        ratio.Draw("PE SAME")
        canvas.RemoveMiddleAxis()

        frame.GetYaxis().SetLabelSize(textsize)
        ratio.GetYaxis().SetLabelSize(textsize)
        frame.GetYaxis().SetTitle("")
        ratio.GetXaxis().SetLabelSize(textsize)
        canvas.Draw()

        #canvas.SaveAs('plots/combine/' + ftype + '_' + cr + '.gif') #TODO: Implement SaveAs in RatioCanvas
        canvas.Print('plots/combine/' + ftype + '_' + cr + '.png')
        canvas.Print('plots/combine/' + ftype + '_' + cr + '.pdf')


    # Checking consistency of shapes 
    newcanvas = TCanvas("newcanvas", "newcanvas", 800, 600)
    legend = TLegend(  newcanvas.GetLeftMargin()+0.55, 
                                1-newcanvas.GetTopMargin()-.15, 
                                newcanvas.GetLeftMargin()+(1.-(newcanvas.GetLeftMargin()+newcanvas.GetRightMargin())),
                                1-newcanvas.GetTopMargin())
    histSB = file.Get("shapes_fit_s/SB/bkg").DrawCopy("HIST E")
    histCR = file.Get("shapes_fit_s/CR/bkg").DrawCopy("HIST E SAME")
    histCR.SetLineColor(3)
    histSB.SetLineColor(4)
    histCR.Scale(1./histCR.Integral())
    histSB.Scale(1./histSB.Integral())
    histCR.SetFillColor(3)
    histSB.SetFillColor(4)
    histCR.SetFillStyle(3003) #3003
    histSB.SetFillStyle(3356)
    legend.AddEntry(histCR, "CR", "L")
    legend.AddEntry(histSB, "SB", "L")
    legend.SetTextSize(0.04)
    legend.Draw()
    newcanvas.Draw()
    newcanvas.Print("plots/combine/ShapeDstarDs.pdf")
