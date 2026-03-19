#include "TMVA/TMVAGui.h"
#include "TMVA/tmvaglob.h"
#include "TCanvas.h"
#include "TH1D.h"



Int_t canvasWidth = 550; 
Int_t canvasHeight = 0.65*canvasWidth; 

Float_t textsize = 0.05; 


TCanvas* DrawTMVAHistogram(TH1 *sig, TH1 *bgd, Int_t n) {	

	TCanvas *canv = new TCanvas( Form("canvas%d", n), "Input variable distribution", canvasWidth, canvasHeight );

      Float_t sc = 1.3;
      TH1 *signal = sig->DrawCopy( "hist" );
      signal->SetStats(false); 

      signal->GetYaxis()->SetTitleOffset( 0.9 );
      TH1 *background = bgd->DrawCopy("histsame");
      signal->SetMaximum( TMath::Max( signal->GetMaximum(), background->GetMaximum() )*sc );
      background->SetStats(false); 
      TString ytit = TString("(1/N) ") + signal->GetYaxis()->GetTitle();
      signal->GetYaxis()->SetTitle( ytit ); // histograms are normalised
      signal->GetYaxis()->SetLabelSize(0.05);
      signal->GetYaxis()->SetTitleSize(textsize); 
      signal->GetXaxis()->SetTitleSize(textsize);
      signal->GetXaxis()->SetLabelSize(0.05);
      signal->Draw("sameaxis"); 

      TMVA::TMVAGlob::SetSignalAndBackgroundStyle( signal, background ); 
      

      //canv->Print(directory+"/"+Form("Variable_no_%i.pdf", n)); 

      bool displaylegend = true;    // We could set this for instance to !n in order to display the legend only on the first graph 
                                    // Or to !(n%r) to display it on all r graphs. 

      if (displaylegend) {
      TLegend *legend= new TLegend( canv->GetLeftMargin(), 
                                    1-canv->GetTopMargin()-.15, 
                                    canv->GetLeftMargin()+.4, 
                                    1-canv->GetTopMargin() );
      legend->SetFillStyle(1);
      legend->AddEntry(signal,"Signal","F");
      legend->AddEntry(background,"Background","F");
      legend->SetBorderSize(1);
      legend->SetMargin( 0.3 );
      legend->Draw("same");
     }

      Int_t    nbin = signal->GetNbinsX();
      Double_t dxu  = signal->GetBinWidth(0);
      Double_t dxo  = signal->GetBinWidth(nbin+1);
      TString uoflow = "";
      uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
                       signal->GetBinContent(0)*dxu*100, background->GetBinContent(0)*dxu*100,
                       signal->GetBinContent(nbin+1)*dxo*100, background->GetBinContent(nbin+1)*dxo*100 );

      uoflow = Form( "U/O-flow (S,B): (%.1f, %.1f)%% / (%.1f, %.1f)%%", 
                       signal->GetBinContent(0)*dxu*100, background->GetBinContent(0)*dxu*100,
                       signal->GetBinContent(nbin+1)*dxo*100, background->GetBinContent(nbin+1)*dxo*100 ); 



     return canv; 

  }