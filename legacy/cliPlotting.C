#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLatex.h"
#include <iostream>
#include "DrawTMVAHistogram.C"
#include "GetSeparation.C"
#include "FileFlow.h"


using namespace ROOT; 
using namespace Ana; 


std::vector<TString> formats = {".pdf", ".png"};



template <typename T> 
T* PrettyPlot(const TString& name = "plot", const TString& options = "", T *obj = nullptr) 
{
	if (!obj)
	{
		obj = gPad; 
	}

	TCanvas *canv = new TCanvas("canvas", "canvas", 800, 600); 

	obj->DrawCopy(); 

	canv->Draw(options); 

	for (auto format : formats) 
	{
		canv->Print(name+format); 
	}

	canv->SaveAs(name+".root"); 

	return obj; 	
}



TObject* PrettyPlot(const TString& name = "plot", const TString& options = "", TObject *obj = nullptr) 
{
	if (!obj)
	{
		obj = gPad; 
	}

	TCanvas *canv = new TCanvas("canvas", "canvas", 800, 600); 

	obj->Draw(); 

	canv->Draw(options); 

	for (auto format : formats) 
	{
		canv->Print(name+format); 
	}

	canv->SaveAs(name+".root"); 

	return obj; 	
}


void SetHeader( TVirtualPad* pad, int pos, TString extraText = "Private Work")
{            
	int iPosX = pos; //2; 

	bool writeExtraText = true;
	float extraOverCmsTextSize  = 0.76;	

	TString cmsText = "CMS"; 
	float cmsTextFont   = 61;  // default is helvetic-bold
	//TString extraText   = "Private Work";
	float extraTextFont = 52;  // default is helvetica-italics	
	TString lumiText = "33.6 fb^{-1} (13 TeV)";

	// text sizes and text offsets with respect to the top frame
	// in unit of the top margin size
	float lumiTextSize     = 0.5; //0.6
	float lumiTextOffset   = 0.2;
	float cmsTextSize      = 0.75;
	float cmsTextOffset    = 0.1;  // only used in outOfFrame version	

	float relPosX    = 0.045;
	float relPosY    = 0.035;
	float relExtraDY = 1.2;	

	float extraTextSize = extraOverCmsTextSize*cmsTextSize;


  int alignX_=2;
  int alignY_=3;
  	
  if(iPosX==0) 
  {
  	alignX_=1;
  	alignY_=1;
  }

  if( iPosX/10==0 ) alignX_=1;
  
  if( iPosX/10==1 ) alignX_=1;
  if( iPosX/10==2 ) alignX_=2;
  if( iPosX/10==3 ) alignX_=3;
  //if( iPosX == 0  ) relPosX = 0.12;

  int align_ = 10*alignX_ + alignY_;

 	float H = pad->GetWh();
	float W = pad->GetWw();
	float l = pad->GetLeftMargin();
	float t = pad->GetTopMargin();
	float r = pad->GetRightMargin();
	float b = pad->GetBottomMargin();
  //  float e = 0.025;
  

  pad->cd();

  	

  TLatex latex;
  latex.SetNDC();
  latex.SetTextAngle(0);
  latex.SetTextColor(kBlack);    


  latex.SetTextFont(42);
  latex.SetTextAlign(31); 
  latex.SetTextSize(lumiTextSize*t);    
  latex.DrawLatex(1-r,1-t+lumiTextOffset*t,lumiText);


  if( pos == 0 )
  {
  	float posX_ =   l+0.05; //l +  relPosX*(1-l-r);
	  float posY_ =   1-2.*t+lumiTextOffset*t;
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(l,1-t+lumiTextOffset*t,cmsText);

    if (writeExtraText) 
    {
    	latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(align_);
      latex.DrawLatex(posX_, posY_, extraText);     
    }
  }

  if( pos == 1 )
  {
  	relPosX = cmsTextSize*t*0.75*cmsText.Length(); // Factor defining how much of the CMS text size the extratext should be dispèlaced
  	float xpos = l + relPosX*(1-l-r);
  	float ypos = 1-t+lumiTextOffset*t; 
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(l, ypos, cmsText);

    if (writeExtraText) 
    {
    	latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(xpos, ypos, extraText);     
    }
  }

  if( pos == 2 )
  {
  	float xpos = l+0.05;
  	float ypos = 1-2*t+lumiTextOffset*t; 
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(xpos, ypos, cmsText);

    if (writeExtraText) 
    {
    	latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(xpos, ypos - extraTextSize*t*1.1, extraText);     
    }
  }

  if( pos == 3 )
  {
  	relPosX = 0.4; // Factor defining how much of the CMS text size the extratext should be dispèlaced
  	std::cout << cmsTextSize*t << " " << relPosX << " " << W << std::endl; 
  	float xpos = l+relPosX;
  	float ypos = 1-2*t+lumiTextOffset*t; 
    latex.SetTextFont(cmsTextFont);
    latex.SetTextAlign(11); 
    latex.SetTextSize(cmsTextSize*t);    
    latex.DrawLatex(xpos, ypos, cmsText);

    if (writeExtraText) 
    {
    	latex.SetTextFont(extraTextFont);
      latex.SetTextSize(extraTextSize*t);
      latex.SetTextAlign(11);
      latex.DrawLatex(xpos, ypos - extraTextSize*t*1.1, extraText);     
    }
  }
  
  pad->cd();


  return;
}


void xLegend(const std::string& text) 
{
	auto plot = gPad; 
	auto hist = static_cast<TH1*>(plot->GetPrimitive("htemp"));
	hist->GetXaxis()->SetTitle(text.data()); 
	hist->GetXaxis()->SetTitleSize(0.06);
	hist->GetXaxis()->SetLabelSize(0.06);
	hist->GetXaxis()->SetTitleOffset(1.2); 
	plot->SetBottomMargin(0.15);
}


void yLegend(const std::string& text) 
{
	auto plot = gPad; 
	auto hist = static_cast<TH1*>(plot->GetPrimitive("htemp"));
	hist->GetYaxis()->SetTitle(text.data());
	hist->GetYaxis()->SetLabelSize(0.06);
	hist->GetYaxis()->SetTitleSize(0.06);
	hist->GetXaxis()->SetTitleOffset(1.2); 
	plot->SetLeftMargin(0.15); 
}


void Style(int position = 1) 
{
	auto plot = gPad; 
	gStyle->SetOptStat(0);
	plot->SetTitle(""); 
	static_cast<TH1*>(plot->GetPrimitive("htemp"))->SetTitle("");
	plot->Draw(); 
	SetHeader(plot, position); 

}
