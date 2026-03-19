#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>





using namespace ROOT; 


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


void PlotGenLevel(TString campaignName = "BkgGenParticleDefinitions/") 
{
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C+");
	//gSystem->Load("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.so"); 

	FileManager filemanager; 


	//filemanager.AddItem("reference", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/dummy.root", "GenEvents"); 
	//filemanager.AddItem("filters", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/RunDstarfilterConverted.root", "GenEvents"); 
	//filemanager.AddItem("filters", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/TestDstarFilterLargerConverted.root", "GenEvents"); 
	//filemanager.AddItem("baseline", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/CONVERTEDnoEtaCut.root", "GenEvents"); 
	//filemanager.AddItem("nocut", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/CONVERTEDsecond.root", "GenEvents"); 
	//filemanager.AddItem("restricted", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/CONVERTEDrestrictedPhaseSpace.root", "GenEvents"); 
	filemanager.AddItem("particles", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v6.8/convertedParticles.root", "GenEvents"); 
	filemanager.AddItem("definitions", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/v6.8/convertedDefinition.root", "GenEvents"); 


	filemanager.OpenAllItems(); 

	gStyle->SetOptStat(0); 

	std::vector<TString> samples = {"particles", "definitions"}; 

	vector<Int_t> colors = {4, 2, 3, 6, 7, 5, 9, 8, 15};



	std::vector<TString> quantitiesToPlot = {
		"B0_pt>>h(50, 0., 30.)", 
		"B0_eta>>h(50, -5., 5.)", 
		"B0_phi>>h(50, -5., 5.)", 
		"B0_m>>h(50, 4.5, 6.)", 
		"B0_q>>h(5, -2.5, 2.5)", 
		"Dstar_pt>>h(50, 0., 20.)", 
		"Dstar_eta>>h(50, -5., 5.)", 
		"Dstar_phi>>h(50, -5., 5.)", 
		"Dstar_m>>h(50, 1.9, 2.1)", 
		"Dstar_q>>h(5, -2.5, 2.5)", 
		"D0_pt>>h(50, 0., 15.)", 
		"D0_eta>>h(50, -5., 5.)", 
		"D0_phi>>h(50, -5., 5.)", 
		"D0_m>>h(50, 1.8, 1.95)", 
		"D0_q>>h(5, -2.5, 2.5)", 
		"K_pt>>h(50, 0., 10.)", 
		"K_eta>>h(50, -5., 5.)", 
		"K_phi>>h(50, -5., 5.)", 
		"K_m>>h(50, 0.45, 0.55)", 
		"K_q>>h(5, -2.5, 2.5)", 
		"pi_pt>>h(50, 0., 10.)", 
		"pi_eta>>h(50, -5., 5.)", 
		"pi_phi>>h(50, -5., 5.)", 
		"pi_m>>h(50, 0.1, 0.2)", 
		"pi_q>>h(5, -2.5, 2.5)", 
		"pis_pt>>h(50, 0., 5.)", 
		"pis_eta>>h(50, -5., 5.)", 
		"pis_phi>>h(50, -5., 5.)", 
		"pis_m>>h(50, 0.1, 0.2)", 
		"pis_q>>h(5, -2.5, 2.5)", 
		"tau_pt>>h(50, 0., 20.)", 
		"tau_eta>>h(50, -5., 5.)", 
		"tau_phi>>h(50, -5., 5.)", 
		"tau_m>>h(50, 1.5, 2.0)", 
		"tau_q>>h(5, -2.5, 2.5)", 
		"tauPi_pt>>h(50, 0., 10.)", 
		"tauPi_eta>>h(50, -5., 5.)", 
		"tauPi_phi>>h(50, -5., 5.)", 
		"tauPi_m>>h(50, 0.1, 0.2)", 
		"tauPi_q>>h(5, -2.5, 2.5)", 
		"tauPi0_pt>>h(50, 0., 10.)", 
		"tauPi0_eta>>h(50, -5., 5.)", 
		"tauPi0_phi>>h(50, -5., 5.)", 
		"tauPi0_m>>h(50, 0.1, 0.2)", 
		"tauPi0_q>>h(5, -2.5, 2.5)", 
		"nutau_pt>>h(50, 0., 10.)", 
		"nutau_eta>>h(50, -5., 5.)", 
		"nutau_phi>>h(50, -5., 5.)", 
		"nutau_m>>h(50, 0., 0.1)", 
		"nutau_q>>h(5, -2.5, 2.5)", 
		"nutaubar_pt>>h(50, 0., 10.)", 
		"nutaubar_eta>>h(50, -5., 5.)", 
		"nutaubar_phi>>h(50, -5., 5.)", 
		"nutaubar_m>>h(50, 0., 0.1)", 
		"nutaubar_q>>h(5, -2.5, 2.5)"
	}; 

	int totalnum = 10000000; 

	for (auto sample : samples) 
	{
		std::cout << "Summary: " << sample << " number: " << filemanager.GetItem<TTree*>(sample)->GetEntries() << " , eff: " << static_cast<double>(filemanager.GetItem<TTree*>(sample)->GetEntries())/static_cast<double>(totalnum) << std::endl; 
	}


	std::vector<TString> formatsToPlot = {".pdf", ".png"}; 

	TString outfolder = "plots/BkgGeneration/"+campaignName; 

	bool webpublication = true; 

	if (gSystem->AccessPathName(outfolder)) gSystem->Exec("mkdir -p "+outfolder); 

	TString webfolder = "/eos/home-m/mhuwiler/www/Analysis/"+campaignName; 
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

	TCanvas *dummycanvas = new TCanvas("dummycanvas", "dummycanvas", 800, 600); 


	for (auto quantity : quantitiesToPlot) 
	{
		auto strings = quantity.Tokenize(">>"); 

		//strings->Print(); 

		TString name = static_cast<TObjString*>(strings->At(0))->GetString(); 

		std::cout << "To be drawn: " << quantity << std::endl; 
		TCanvas *canvas = new TCanvas(name, name, 800, 600); 

		if (not quantity.Contains(">>")) quantity+=">>h";

		std::map<TString, TH1*> histocollection; 

		
		TLegend *legend= new TLegend( canvas->GetLeftMargin()+0.5, 
                                    1-canvas->GetTopMargin()-std::min(0.7*static_cast<double>(samples.size()), 0.31), 
                                    //canvas->GetLeftMargin()+.4, 
                                    canvas->GetLeftMargin()+(1.-(canvas->GetLeftMargin()+canvas->GetRightMargin())),
                                    1-canvas->GetTopMargin() );
      	legend->SetFillStyle(1);
      	legend->SetBorderSize(1);
      	legend->SetMargin( 0.3 );
      	legend->SetTextSize(0.04);

      	dummycanvas->cd(); 
		filemanager.GetItem<TTree*>(samples.at(0))->Draw(TString(quantity).ReplaceAll(">>h", ">>ref")); 
		TH1 *reference = static_cast<TH1*>(dummycanvas->GetPrimitive("ref")); 
		histocollection.insert(std::make_pair(samples.at(0), reference)); 
		reference->GetYaxis()->SetTitleOffset(0.9); 
		
		reference->SetLineColor(1); 
		reference->SetFillColor(1); 
		reference->SetFillStyle(3356);
		reference->SetLineWidth(2); 
		legend->AddEntry(reference,"baseline","F");
		reference->SetTitle(name); 
		canvas->cd(); 
		reference->Draw("HIST"); 

		//canvas->Draw(); 

		Double_t normalisation = reference->Integral(); 
		std::cout << "Number of events in reference: " << normalisation << std::endl; 

		Double_t previousMax = reference->GetMaximum(); //-999.; 
		std::cout << "Max reference: " << previousMax << std::endl; 

		Int_t index = 0; 

		for (int i=1; i<samples.size(); i++)
		{
			auto sample = samples.at(i); 
			TString variable = quantity; // Making a local copy for modification 
			variable.ReplaceAll(">>h", ">>"+sample); 
			std::cout << "Processing sample " << sample << std::endl; 
			TTree *tree = filemanager.GetItem<TTree*>(sample); 
			if (!tree) 
			{
				std::cout << "WARNING: Tree missing with reference: " << sample << std::endl; 
				continue; 
			}
			dummycanvas->cd(); 
			filemanager.GetItem<TTree*>(sample)->Draw(variable); 
			TH1 *histo = static_cast<TH1*>(dummycanvas->GetPrimitive(sample)); 
			histocollection.insert(std::make_pair(sample, histo)); 
		
			histo->Scale(normalisation/histo->Integral()); 
			previousMax = TMath::Max(previousMax, histo->GetMaximum()); 

			histo->SetLineColor(colors.at(index)); 
			histo->SetFillColor(colors.at(index)); 
			histo->SetFillStyle(3003); 

			histo->SetLineWidth(2); 

			legend->AddEntry(histo, sample, "F");

			canvas->cd(); 
			histo->Draw("HISTSAME"); 

			//canvas->Modified(); 

			index++; 
		

		}


		

      	

      	

      	Float_t sc = 1.3;
      	//if not (histo->GetMaximum())
      	for (auto it = histocollection.begin(); it != histocollection.end(); it++) 
      	{
      		//std::cout << "Scaling " << it->first << std::endl; 
      		it->second->SetMaximum( previousMax*sc );
      	}

      	
      	

      	legend->Draw(); 


		//histo->GetXaxis()->SetRangeUser(0., 100.); 

		canvas->Draw(); 
		//canvas->Update(); 

		//PauseUntilEnter(); 

		for (auto format : formatsToPlot) 
		{
			canvas->Print(outfolder+name+format); 
			if (webpublication) canvas->Print(webfolder+name+format); 
		}

		delete canvas; 
	}

	delete dummycanvas; 

	//histo1->Draw(); 

	//histo2->Draw(); 

	//canvas->Draw(); 

	//canvas->Print("Example.pdf"); // Works 

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


}

