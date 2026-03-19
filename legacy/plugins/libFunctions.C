#ifndef libFunctions_lib
#define libFunctions_lib
#include <cstdlib>
#include <vector>
#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "THStack.h"
#include <sstream>
#include "FileManager/CFileManager.C"
#include "CustomBins/CustomBins.C"
//#include "TMVA/tmvaglob.h"  // For TMVA style, may be removed once an own style is implemented


Double_t getEfficiency(const Double_t n, const Double_t N) {
	Double_t eff = n/N; 
	return eff; 
}

Double_t getEfficiencyError(const Double_t eff, const Double_t N) {
	return sqrt(eff*(1.-eff)/N); 
}

// Did not prove useful for returning the size of a C array, because still here the size needs to be known at compile time (may be deleted)
template<size_t SIZE, class T> inline size_t array_size(T (&arr)[SIZE]) {
    return SIZE;
}


// Can be deleted, now the class CustomBins takes care of this (or may be rethought to integrate CustomBins) 
TH1D* makeTH1D(TString name, TString title, std::vector<Double_t> bins) 
{
	TH1D *hist = new TH1D(name, title, bins.size()-1, bins.data()); 
	return hist; 
}

void setStyle(TH1* histo) {
	Double_t titleSize = 0.05; 
	Double_t labelSize = 0.04; 
	histo->GetXaxis()->SetTitleSize(titleSize);
	histo->GetXaxis()->SetTitleOffset(1.05);
	histo->GetXaxis()->SetLabelSize(labelSize); 
	histo->GetYaxis()->SetTitleSize(titleSize);
	histo->GetYaxis()->SetTitleOffset(1.05);
	histo->GetYaxis()->SetLabelSize(labelSize); 

	histo->SetLineWidth(2); 
	//histo->SetMarkerSize(1.2); 
	//histo->SetMarkerStyle(8); 


	TObject* isCurrentCanvas = gPad->GetPrimitive((histo->GetName())); 	// Making sure the graph is plot in gPad 
	if (isCurrentCanvas) 
	{
		gPad->SetLeftMargin(0.11); 
		gPad->SetBottomMargin(0.11); 
		//gPad->SetRightMargin(0.15); 
		gPad->Modified();  
    	gPad->Update();

	} 
}

void plotHistograms(std::vector<std::pair<TH1*,TString> > histograms, TCanvas *canvas) 
{
	if (histograms.size() > 10) {
		std::cout << "You are trying to draw more than 10 graphs on the same plot. You may need a more elaborate color handling. " << std::endl; 
		return; 
	}

	vector<Int_t> colors = {4, 2, 1, 3, 6, 7, 5, 9, 8, 15};

	THStack *stack = new THStack("stack", ""); 

	TLegend *legend= new TLegend( canvas->GetLeftMargin(), 
                                    1-canvas->GetTopMargin()-.15, 
                                    //canvas->GetLeftMargin()+.4, 
                                    canvas->GetLeftMargin()+(1.-(canvas->GetLeftMargin()+canvas->GetRightMargin())),
                                    1-canvas->GetTopMargin() );

	for (unsigned int i=0; i<histograms.size(); i++) 
	{
		TH1 *histogram = histograms.at(i).first;

		histogram->SetLineColor(colors.at(i)); 
		histogram->SetMarkerColor(colors.at(i)); 

		setStyle(histogram); 

		stack->Add(histogram); 

		legend->AddEntry(histogram, histograms.at(i).second, "p"); 


	}

	canvas->cd(); 

	stack->Draw("nostack"); 


	legend->SetTextSize(0.04);

	legend->Draw();
	canvas->Modified(); 
	canvas->Update(); 


}

THStack* plotTwoHistograms(TH1 *h1, TH1 *h2, TCanvas *canvas, TString legh1, TString legh2, Int_t norm = 0) 
{
	THStack *hs = new THStack("hs","");

  	if (norm == 1) 
  	{
  		double norm = 1.;
  		h1->Scale(norm/h1->Integral()); 
  		h2->Scale(norm/h2->Integral()); 
  	}
  	else if (norm == 2) 
  	{
  		h2->Scale(h1->Integral()/h2->Integral()); 
  	}
  	else if (norm != 0)
  	{
  		std::cerr << "ERROR: The normalisation parameter is odd. Please specify an integer between 0 and 2. The options are following: \n\n 0: Don't normalize\n 1: Normalize both histograms to 1\n 2: Normalize the area under curve of the second histogram to the one of the first\n\n The parameter specified was ignored and the histograms are not normalised (same behaviour as default option 0).\n" << std::endl; 
  	}

  	setStyle(h1); 
  	setStyle(h2); 

  	Float_t sc = 1.2;
  	Double_t maximumRange = TMath::Max(h1->GetMaximum(), h2->GetMaximum()) * sc;
  	h1->SetMaximum(maximumRange); 
  	h2->SetMaximum(maximumRange); 

  	hs->Add(h1); 
  	hs->Add(h2); 

  	h1->SetLineColor(kBlue); 
  	h2->SetLineColor(kRed); 

  	h1->SetLineWidth(2); 
  	h2->SetLineWidth(2); 

  	//TMVA::TMVAGlob::SetSignalAndBackgroundStyle(h1, h2); 


	TLegend *legend= new TLegend( canvas->GetLeftMargin(), 
                                    1.-canvas->GetTopMargin()-.15, 
                                    //canvas->GetLeftMargin()+.4, 
                                    canvas->GetLeftMargin()+ 0.25,
                                    1.-canvas->GetTopMargin() );

	/* 	// Top right corner 
	TLegend *legend= new TLegend( canvas->GetLeftMargin()+0.5, 
                                    1-canvas->GetTopMargin()-.15, 
                                    //canvas->GetLeftMargin()+.4, 
                                    1.-canvas->GetRightMargin(),
                                    1-canvas->GetTopMargin() );
    */ 


	legend->AddEntry(h1, legh1, "F"); 
	legend->AddEntry(h2, legh2, "F"); 


	canvas->cd(); 

	hs->Draw("nostack, e");

	legend->SetTextSize(0.04);

	legend->Draw();
	canvas->Modified(); 
	canvas->Update(); 

	return hs; 


}

// This is crap, should be done with style 
void setStyle(TGraph* graph) {
	Double_t titleSize = 0.05; 
	Double_t labelSize = 0.04; 
	graph->GetXaxis()->SetTitleSize(titleSize);
	graph->GetXaxis()->SetTitleOffset(1.05);
	graph->GetXaxis()->SetLabelSize(labelSize); 
	graph->GetYaxis()->SetTitleSize(titleSize);
	graph->GetYaxis()->SetTitleOffset(1.05);
	graph->GetYaxis()->SetLabelSize(labelSize); 

	graph->SetLineWidth(2); 
	graph->SetMarkerSize(1.2); 
	graph->SetMarkerStyle(8); 


	TObject* isCurrentCanvas = gPad->GetPrimitive((graph->GetName())); 	// Making sure the graph is plot in gPad 
	if (isCurrentCanvas) 
	{
		gPad->SetLeftMargin(0.11); 
		gPad->SetBottomMargin(0.11); 
		//gPad->SetRightMargin(0.15); 
		gPad->Modified();  
    	gPad->Update();

	} 
}
	

TGraph* plotMultipleGraphs(vector<std::pair<TGraph*, TString> > graphs, TCanvas* canvas) {

	//TCanvas *canvas = new TCanvas("canvas", "ROC curves", 800, 600); 
	vector<Int_t> colors = {1, 4, 2, 3, 6, 7, 5, 9, 8, 15}; 

	if (graphs.size() > 10) {
		std::cout << "You are trying to draw more than 10 graphs on the same plot. You may need a more elaborate color handling. " << std::endl; 
		return nullptr; 
	}

	canvas->cd(); 

	TLegend *legend = new TLegend(0.11, 0.11, 0.7, 0.3);

	//TMultiGraph *graphCollection = new TMultiGraph(); 
	TGraph* referenceCurve = nullptr;

	for (unsigned int i=0; i < graphs.size(); i++) 
	{
		TGraph *curve = graphs.at(i).first; 

		curve->SetTitle("");
		curve->SetLineColor(colors.at(i));
		curve->SetMarkerColor(colors.at(i)); 
 
		if (i==0) {
			curve->SetName("reference"); 
			curve->Draw("AP."); 
			//curve->GetYaxis()->SetRangeUser(0.055, 0.2)	// Setting axis range must be done on the first graph to be plotted (with option A);
			referenceCurve = curve; 

		}
		else 
		{
			curve->Draw("P.");
		}

		setStyle(curve);

		//graphCollection->Add(curve); 
		
		legend->AddEntry(curve, graphs.at(i).second, "p");		
	
	
		

	}

	legend->SetTextSize(0.04);


	//graphCollection->Draw("AP."); 

	

	legend->Draw();
	canvas->Modified(); 
	canvas->Update(); 
	//canvas->Draw();

	return referenceCurve; 

	
}

// This method is not very good, since maps are internally ordered, one has no control of the order in which the graphs are plotted. 
void plotMultipleGraphs(std::map<TGraph*, TString> graphs, TCanvas* canvas) 
{
	std::vector<std::pair<TGraph*, TString> > graphCollection; 
	for(std::map<TGraph*, TString>::iterator it=graphs.begin(); it!=graphs.end(); ++it) 
	{
		graphCollection.push_back(std::pair<TGraph*, TString>(it->first, it->second)); 
	}
	plotMultipleGraphs(graphCollection, canvas); 
}

// HTML writer
// stores string
// saves
// setautosave (save at each modification)

#endif
