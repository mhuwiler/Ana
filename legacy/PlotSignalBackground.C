#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
ClassImp(FileManager)
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include "FileFlow.h"


using namespace ROOT; 

using namespace Ana; 


void PauseUntilAnyKey() 
{
	std::cout << "Press any key to continue... " << std::endl; 
	std::cin.get(); 
}

void PauseUntilEnter() 
{
	std::cout << "Press 'enter' to continue..." << std::endl;
	std::cin.ignore(); 
	//std::cin.ignore(std::numeric_limits<streamsize>::max(),'\n'); // #include <limits>
}

void Pause(Int_t timeInSec) 
{
	// Better way, taken from: https://stackoverflow.com/questions/23609507/pause-program-execution-for-5-seconds-in-c
	#include <chrono>
	#include <thread>
	//std::this_thread::sleep_for(static_cast<std::chrono::seconds>(timeInSec));
	sleep(timeInSec); 
}

TLorentzVector LV(double pt, double eta, double phi, double m) 
{
	TLorentzVector V; 
	V.SetPtEtaPhiM(pt, eta, phi, m); 
	return V; 
}


void PlotSignalBackground(TString campaignName = "PlotsBackgroundWithWS/") 
{
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C+");
	//gSystem->Load("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.so"); 

	Init("v1"); 

	std::vector<TString> filesUsed = {"dataD2", "Sig", "BkgDstarDs", "BkgDstarDsstar", "BkgDstar3pi", "dataD2WS", "dataD2TauWS"}; //"DstarDsMCfirst", "Data2018BFirst"


	for (auto item : filesUsed) 
	{
		std::cout << "Opening file: " << item << std::endl; 
		filemanager.OpenItem(item); 
	}

	gStyle->SetOptStat(0); 


	//auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>("DstarDsMCfirst")); // tree100k

	//auto sampleMC = RDataFrame(*filemanager.GetItem<TTree*>("Data2018BFirst")); 
	

	// Defining the delta
	auto P = [](double pt, double eta, double phi, double m) {
		TLorentzVector V; 
		V.SetPtEtaPhiM(pt, eta, phi, m); 
		return V; 
	}; 

	auto P_v = [](std::vector<float>& pt, std::vector<float>& eta, std::vector<float>& phi, std::vector<float>& m) {
		size_t vecSize = pt.size(); 
		assert(eta.size() == vecSize and phi.size() == vecSize and m.size() == vecSize); 
		std::vector<TLorentzVector> result; 
		for (unsigned int i=0; i<vecSize; i++) 
		{
			TLorentzVector V; 
			V.SetPtEtaPhiM(pt.at(i), eta.at(i), phi.at(i), m.at(i)); 
			result.push_back(V); 
		}
		return result; 
	}; 

	auto invMass = [](TLorentzVector V1, TLorentzVector V2) {
		auto V = V1 + V2; 
		return V.M(); 
	}; 

	auto invMass_v = [](std::vector<TLorentzVector>& V1, std::vector<TLorentzVector>& V2) {
		assert(V1.size() == V2.size()); 
		std::vector<float> result; 
		std::cout << "Size of branch (vector): " << V1.size() << std::endl; 
		for (unsigned int i=0; i<V1.size(); i++) 
		{
			auto V = V1.at(i) + V2.at(i); 
			result.push_back(V.M()); 
		}
		return result; 
	}; 

	//auto invMass_v2 = [invMass](std::vector<TLorentzVector>& V1, std::vector<TLorentzVector>& V2) {
	//	assert(V1.size() == V2.size()); 
	//	std::vector<float> result; 
	//	for (unsigned int i=0; i<V1.size(); i++) 
	//	{
	//		result.push_back(std::invoke(invMass, V1, V2)); 
	//	}
	//	return result; 
	//}; 

	auto mass = [](TLorentzVector& V) {
		return V.M(); 
	}; 

	// overloading lambdas: if constexpr (std::is_same_v<T, int>) 

	std::vector<TString> quantitiesToPlotFromTree = { "b_tau_rhomass1>>h(100, 0., 3.)", "b_tau_rhomass2>>h(100, 0., 3.)" }; // "BsDstarTauNu_B_mass>>h(100, 2., 6.)", "b_tau_mass>>h(100, 0., 2.)"


	TTree *tree = filemanager.GetItem<TTree*>("DstarDsMCfirst"); 

	std::vector<TString> formatsToPlot = {".pdf", ".png"}; 

	//TString campaignName = "PlotsDstarDsVsSignal/"; 

	bool normalise = true; 

	TString outfolder = "plots/"+campaignName; //plots/"+campaignName; //"plots/PlotsGenmactchedFinal/"

	bool webpublication = false; 

	//TString region = "SB"; 
	std::vector<TString> regions = {"SR", "CR", "SB"}; 

	if (gSystem->AccessPathName(outfolder)) gSystem->Exec("mkdir -p "+outfolder); 

	TString webfolder = "/eos/home-m/mhuwiler/www/Analysis/"+campaignName; // "/eos/home-m/mhuwiler/www/Analysis/DataMCplotsGenmatchedFinal/"
	if (webpublication) 
	{ 
		if (gSystem->AccessPathName(webfolder)) gSystem->Exec("mkdir -p "+webfolder); 
		TString webenginesource = "/eos/home-m/mhuwiler/software/php-plots/"; 
		gSystem->Exec("cp -r "+webenginesource+"res "+webfolder); 
		gSystem->Exec("cp "+webenginesource+"index.php "+webfolder); 
		// Get and edit the permission file
		ifstream accessFileSource(webenginesource+"example/htaccess"); 
		ofstream accessFileTarget(webfolder+".htaccess"); 
		std::string accessFileLine; 
		while (getline(accessFileSource, accessFileLine)) 
		{
			TString targetFileLine(accessFileLine); 
			targetFileLine.ReplaceAll("/<me>/<my-project>/", webfolder); 
			accessFileTarget << targetFileLine; 
		}
		accessFileSource.close(); 
		accessFileTarget.close(); 

		// Hack to have the file synchronised over the cernbox client(without hidden fyle sync), needs to copy it back from htaccess to .htaccess in the target directory. 
		gSystem->Exec("cp "+webfolder+".htaccess "+webfolder+"htaccess"); 

		// TODO: cleanup this hacky editing of variable! 
		webfolder = webfolder+"plots/"; 
		if (gSystem->AccessPathName(webfolder)) gSystem->Exec("mkdir -p "+webfolder); 
	}

	for (auto region : regions) 
	{
		for (auto quantity : quantitiesToPlotFromTree) 
		{
			auto strings = quantity.Tokenize(">>"); 

			//strings->Print(); 

			TString name = static_cast<TObjString*>(strings->At(0))->GetString(); 

			std::cout << "To be drawn: " << quantity << std::endl; 
			TCanvas *canvas = new TCanvas(name, name, 800, 600); 

			
			//tree->Draw(quantity); 

			if (not quantity.Contains(">>")) quantity+=">>h"; 

			TString quantityData = quantity; 
			quantityData.ReplaceAll("_genmatched", ""); 

			std::cout << "Quantity: " << quantity << std::endl; 
			std::cout << "Quantity data: " << quantityData << std::endl; 

			//if (quantity.Contains("track_genmatched_dR_Dstar")) quantityData.ReplaceAll("track_", "track_tracks_"); 
			//if (quantity.Contains("track_genmatched_pvAssociationQuality")) quantityData.ReplaceAll("track_", "track_track_"); 
			//if (quantity.Contains("track_dR_Dstar")) quantityData.ReplaceAll("track_", "track_tracks_"); 
			//if (quantity.Contains("track_pvAssociationQuality")) quantityData.ReplaceAll("track_", "track_track_"); 

			std::cout << "Quantity data: " << quantityData << std::endl; 

			//histo->GetXaxis()->SetRangeUser(0., 100.); 

			TString cut = "(BsDstarTauNu_D0_vprob>0.1) && (BsDstarTauNu_Ds_vprob>0.1)"; 
			Double_t normSig = 1.; 
			Double_t normDs = 1.; 
			Double_t norm3pi = 1.; 
			Double_t normDsstar = 1.; 
			Double_t factorTauWS = 1.84; 
			if (region == "SR") 
			{
				cut = "(BsDstarTauNu_D0_vprob>0.1) && (BsDstarTauNu_Ds_vprob>0.1) && (mvaScore >0.9)"; // && (mvaScore > 0.0) && (mvaScore <= 0.9)
				normSig = 638; 
				normDs = 615; 
				norm3pi = 1.1; 
				normDsstar = 1430; 
			}
			else if (region == "CR") 
			{
				cut = "(BsDstarTauNu_D0_vprob>0.1) && (BsDstarTauNu_Ds_vprob>0.1) && (mvaScore > 0.0) && (mvaScore <= 0.9)"; 
				normSig = 300.; 
				normDs = 607.; 
				norm3pi = 1.1; 
				normDsstar = 800.; 
			}
			else if (region == "SB") 
			{
				cut = "(BsDstarTauNu_D0_vprob>0.1) && (BsDstarTauNu_Ds_vprob>0.1) && (mvaScore > -0.5) && (mvaScore <= 0.0)"; // && (mvaScore > 0.0) && (mvaScore <= 0.9)
				normSig = 62.8; 
				normDs = 186.;
				norm3pi = 1.1; 
				normDsstar = 188.;  
			}
			else 
			{
				std::cout << "ERROR: no known region named " << region << std::endl; 
			}
			TString dummycut = "(BsDstarTauNu_D0_vprob>0.1) && (BsDstarTauNu_Ds_vprob>0.1)"; 

			std::cout << cut << std::endl; 

			filemanager.GetItem<TTree*>("Sig")->Draw(">>eventlist", "1", "goff"); 
			TEventList *eventlist = static_cast<TEventList*>(gDirectory->Get("eventlist")); 
			Int_t numberBeforeMC = eventlist->GetN(); 
			filemanager.GetItem<TTree*>("BkgDstarDs")->Draw(">>eventlist", "1", "goff"); 
			Int_t numberBeforeDs = eventlist->GetN();
			filemanager.GetItem<TTree*>("BkgDstar3pi")->Draw(">>eventlist", "1", "goff"); 
			Int_t numberBefore3pi = eventlist->GetN();
			filemanager.GetItem<TTree*>("BkgDstarDsstar")->Draw(">>eventlist", "1", "goff"); 
			Int_t numberBeforeDsstar = eventlist->GetN();
			filemanager.GetItem<TTree*>("dataD2WS")->Draw(">>eventlist", "1", "goff"); 
			Int_t numberBeforeWS = eventlist->GetN();
			std::cout << "Number of events: " << numberBeforeMC << std::endl; 

			filemanager.GetItem<TTree*>("dataD2")->Draw(quantityData.ReplaceAll(">>h", ">>h1"), cut); 
			TH1 *histoData = static_cast<TH1*>(canvas->GetPrimitive("h1")); 
			histoData->GetYaxis()->SetTitleOffset(0.9); 
			histoData->SetTitle(""); 
			filemanager.GetItem<TTree*>("Sig")->Draw(TString(quantity).ReplaceAll(">>h", ">>h2"), cut); 
			TH1 *histoMC = static_cast<TH1*>(canvas->GetPrimitive("h2")); 
			histoMC->SetTitle(""); 
			filemanager.GetItem<TTree*>("BkgDstarDs")->Draw(TString(quantity).ReplaceAll(">>h", ">>h3"), cut); 
			TH1 *histoBkgDs = static_cast<TH1*>(canvas->GetPrimitive("h3")); 
			histoBkgDs->SetTitle(""); 
			filemanager.GetItem<TTree*>("BkgDstarDsstar")->Draw(TString(quantity).ReplaceAll(">>h", ">>h4"), cut); 
			TH1 *histoBkgDsstar = static_cast<TH1*>(canvas->GetPrimitive("h4")); 
			histoBkgDsstar->SetTitle(""); 
			filemanager.GetItem<TTree*>("BkgDstar3pi")->Draw(TString(quantity).ReplaceAll(">>h", ">>h5"), cut); 
			TH1 *histoBkg3Pi = static_cast<TH1*>(canvas->GetPrimitive("h5")); 
			histoBkg3Pi->SetTitle(""); 
			filemanager.GetItem<TTree*>("dataD2WS")->Draw(TString(quantity).ReplaceAll(">>h", ">>h6"), cut); 
			TH1 *histoDataWS = static_cast<TH1*>(canvas->GetPrimitive("h6")); 
			histoDataWS->SetTitle(""); 
			filemanager.GetItem<TTree*>("dataD2TauWS")->Draw(TString(quantity).ReplaceAll(">>h", ">>h7"), cut); 
			TH1 *histoDataTauWS = static_cast<TH1*>(canvas->GetPrimitive("h7")); 
			histoDataTauWS->SetTitle(""); 

			// Plot roc curve here 
			//if (name == "track_genmatched_doca") 
			//{
			//	TString outname = outfolder+name+"ROC.root"; 
			//	std::cout << "making ROC curve for: " << name << " in file: " << outname << std::endl; 
			//	TCanvas *efficiencyCanvas = new TCanvas("efficiencyCanvas", "efficiencyCanvas", 800, 600); 
			//	filemanager.GetItem<TTree*>("MCgenmatched")->Draw(name+">>h(1000, -1.5, 0.5)"); 

	//		//	TH1 * effHisto = static_cast<TH1*>(efficiencyCanvas->GetPrimitive("h")); 
			//	Double_t denominator = effHisto->Integral(); 

	//		//	for (unsigned int i=0; i<1000; i++) {
			//		
			//	}
			//}

			Int_t numberAfterMC = histoMC->GetEntries(); 
			Int_t numberAfterDs = histoBkgDs->GetEntries(); 
			Int_t numberAfter3pi = histoBkg3Pi->GetEntries(); 
			Int_t numberAfterDsstar = histoBkgDsstar->GetEntries(); 
			Int_t numberAfterWS = histoDataWS->GetEntries(); 
			std::cout << "Efficiency of SR for MC: " << static_cast<float>(numberAfterMC)/static_cast<float>(numberBeforeMC) << std::endl; 
			std::cout << "Efficiency of SR for Ds bkg: " << static_cast<float>(numberAfterDs)/static_cast<float>(numberBeforeDs) << std::endl; 
			std::cout << "Efficiency of SR for 3pi bkg: " << static_cast<float>(numberAfter3pi)/static_cast<float>(numberBefore3pi) << std::endl; 
			std::cout << "Efficiency of SR for Dsstar bkg: " << static_cast<float>(numberAfterDsstar)/static_cast<float>(numberBeforeDsstar) << std::endl; 
			std::cout << "Efficiency of SR for WS bkg: " << static_cast<float>(numberAfterWS)/static_cast<float>(numberBeforeWS) << std::endl; 


			std::cout << "Number of MC events: " << histoMC->Integral() << std::endl; 

			//histoMC->Scale(histoData->Integral()/histoMC->Integral()); 
			if (normalise) 
			{
				//histoData->Scale(1./histoData->Integral()); 
				histoMC->Scale(normSig/histoMC->Integral()); 
				histoBkgDs->Scale(normDs/histoBkgDs->Integral());
				histoBkgDsstar->Scale(normDsstar/histoBkgDsstar->Integral()); 
				histoBkg3Pi->Scale(norm3pi/histoBkg3Pi->Integral()); 
				histoDataTauWS->Scale(factorTauWS); 
			}

			histoData->SetLineColor(kBlue); 
			histoMC->SetLineColor(kRed); 
			histoBkgDs->SetLineColor(kGreen); 
			histoBkgDsstar->SetLineColor(kGreen+3); 
			histoBkg3Pi->SetLineColor(kOrange+2); 
			histoDataWS->SetLineColor(kMagenta+3); 
			histoDataTauWS->SetLineColor(kCyan-7); 


			histoData->SetLineWidth(2); 
			histoMC->SetLineWidth(2); 
			histoBkgDs->SetLineWidth(2); 
			histoBkgDsstar->SetLineWidth(2); 
			histoBkg3Pi->SetLineWidth(2); 
			histoDataWS->SetLineWidth(2); 
			histoDataTauWS->SetLineWidth(2); 

			std::vector<double> maxes = { histoData->GetMaximum(), histoMC->GetMaximum(), histoBkgDs->GetMaximum() }; 



			TLegend *legend= new TLegend( canvas->GetLeftMargin()+0.35, 
	                                    1-canvas->GetTopMargin()-.2, 
	                                    //canvas->GetLeftMargin()+.4, 
	                                    canvas->GetLeftMargin()+(1.-(canvas->GetLeftMargin()+canvas->GetRightMargin())),
	                                    1-canvas->GetTopMargin() );
	      	legend->SetFillStyle(1);
	      	legend->AddEntry(histoData,"data","F");
	      	legend->AddEntry(histoMC,"signal MC (genmatched)","F");
	      	legend->AddEntry(histoBkgDs, "B^{0}#rightarrow D*D_{s} Inclusive"); 
	      	legend->AddEntry(histoBkgDsstar, "B^{0}#rightarrow D*D_{s}* Inclusive"); 
	      	legend->AddEntry(histoBkg3Pi, "B^{0}#rightarrow D*3#pi Nonresonant"); 
	      	legend->AddEntry(histoDataWS, "WS (wrong sign) data"); 
	      	legend->AddEntry(histoDataTauWS, "tau WS (|q| = 3) data"); 
	      	legend->SetBorderSize(1);
	      	legend->SetMargin( 0.3 );
	      	legend->SetTextSize(0.04);

	      	histoData->SetLineColor(kBlue); 
	      	histoMC->SetLineColor(kRed); 

	      	//histoData->SetTitle(name); 

	      	Float_t sc = 1.3;
	      	//if not (histo->GetMaximum())
	      	auto it = max_element(std::begin(maxes), std::end(maxes));

	      	histoData->SetMaximum( *it*static_cast<float>(sc) ); // TMath::Max(histoData->GetMaximum(), histoMC->GetMaximum())*sc
	      	std::cout << histoData->GetMaximum() << " " << histoMC->GetMaximum() << std::endl; 

	      	histoData->Draw("HIST"); 
	      	histoMC->Draw("HISTSAME"); 
	      	histoBkgDs->Draw("HISTSAME"); 
	      	histoBkgDsstar->Draw("HISTSAME"); 
	      	histoBkg3Pi->Draw("HISTSAME"); 
	      	histoDataWS->Draw("HISTSAME"); 
	      	histoDataTauWS->Draw("HISTSAME"); 

	      	legend->Draw(); 


			//histo->GetXaxis()->SetRangeUser(0., 100.); 

			canvas->Draw(); 

			for (auto format : formatsToPlot) 
			{
				canvas->Print(outfolder+name+region+format); 
				if (webpublication) canvas->Print(webfolder+name+region+format); 
			}

			delete canvas; 
		}
	}

	//histo1->Draw(); 

	//histo2->Draw(); 

	//canvas->Draw(); 

	//canvas->Print("Example.pdf"); // Works 

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


}

