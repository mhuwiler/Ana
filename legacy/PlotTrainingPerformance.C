#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"
ClassImp(FileManager)
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include "DrawTMVAHistogram.C"
#include "GetSeparation.C"


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

TCanvas* PlotROCStuff(TDirectory *directory, const TString& graphTitle = "Overtraining assessement") 
{
	TCanvas *canvas = new TCanvas("canvas", graphTitle, 800, 600); 
	// Setting the histograms 
	TH1 *sigTrain = static_cast<TH1*>(directory->Get("signalTrain")); 
	TH1 *bkgTrain = static_cast<TH1*>(directory->Get("backgroundTrain")); 
	TH1 *sigTest = static_cast<TH1*>(directory->Get("signalTest")); 
	TH1 *bkgTest = static_cast<TH1*>(directory->Get("backgroundTest")); 


	Int_t binning = 5; //10.; 

	gStyle->SetHatchesSpacing(0.7); 

	Float_t textsize = 0.04; 

	Bool_t displaylegend = true; 

	sigTrain->Draw("HIST"); 
	sigTrain->SetLineWidth(2); 
	sigTrain->SetLineColor(kBlue+2); 
	sigTrain->SetFillColorAlpha(kBlue+2, 0.17); //0.35
	//sigTrain->SetFillStyle(3003); //3356
	sigTrain->SetFillStyle(1001);
	sigTrain->Rebin(binning); 
	sigTrain->SetTitle(graphTitle); 
	sigTrain->GetXaxis()->SetTitle("MVA score"); 
	sigTrain->GetYaxis()->SetLabelSize(textsize);
    sigTrain->GetYaxis()->SetTitleSize(textsize); 
    sigTrain->GetXaxis()->SetTitleSize(textsize);
    sigTrain->GetXaxis()->SetLabelSize(textsize);

	bkgTrain->Draw("HIST SAME"); 
	bkgTrain->SetLineWidth(2); 
	bkgTrain->SetLineColor(kRed+2); 
	bkgTrain->SetFillColorAlpha(kRed+2, 0.17); 
	//bkgTrain->SetFillStyle(3003);
	bkgTrain->SetFillStyle(1001);
	bkgTrain->Rebin(binning); 

	sigTest->Draw("HIST SAME"); 
	sigTest->SetLineWidth(2); 
	sigTest->SetLineColor(kBlue+2); 
	sigTest->SetFillColor(kBlue+2); 
	sigTest->SetFillStyle(3356);
	sigTest->Rebin(binning); 

	bkgTest->Draw("HIST SAME"); 
	bkgTest->SetLineWidth(2); 
	bkgTest->SetLineColor(kRed+2); 
	bkgTest->SetFillColor(kRed+2); 
	bkgTest->SetFillStyle(3356);
	bkgTest->Rebin(binning); 

	TLegend *legend = nullptr; 
	if (displaylegend) 
	{
    	legend = new TLegend(canvas->GetLeftMargin() + .25, 
                             1-canvas->GetTopMargin() - 0.20, 
                             1-canvas->GetRightMargin() - .15, 
                             1-canvas->GetTopMargin());

      	legend->SetFillStyle(1);
      	legend->AddEntry(sigTrain,"signal train","F");
      	legend->AddEntry(bkgTrain,"background train","F");
      	legend->AddEntry(sigTest,"signal test","F");
      	legend->AddEntry(bkgTest,"background test","F");
      	legend->SetBorderSize(1);
      	legend->SetMargin( 0.3 );
      	legend->Draw("same");
    }

	Double_t maxY = 1.2*TMath::Max(TMath::Max(sigTrain->GetMaximum(), bkgTrain->GetMaximum()), TMath::Max(sigTest->GetMaximum(), bkgTest->GetMaximum())); 

	sigTrain->SetMaximum(maxY); 

	canvas->Draw(); 


	canvas->SetLogy(); 


	return canvas; 
}

// Method to add MVA variable


void PlotTrainingPerformance(TString campaignName = "MVAinputVariablesOfficial50M/") 
{
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C+");
	//gSystem->Load("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.so"); 

	FileManager filemanager; 


	filemanager.AddItem("prodlatest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/flatTupleDataLatest.root", "ntuplizer/tree"); 
	filemanager.AddItem("MCgenmatched", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/flatTupleGenmatchedAllSingleTauLatest.root", "ntuplizer/tree"); 

	filemanager.AddItem("MCOfficialSample", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/all.root", "ntuplizer/tree"); 
	filemanager.AddItem("DstarDsMCfirst", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/bkgDstarDsFirst.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataWS", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataWS.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataLatest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataLast.root", "ntuplizer/tree"); 
	filemanager.AddItem("MCofficial", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/MCFirstSubmission.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataLarge", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge.root", "ntuplizer/tree"); 

	// For ABCD estimation
	filemanager.AddItem("DataLargeMVAfirst", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_mva.root", "tree"); 
	filemanager.AddItem("DataLargeMVA", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_mvaxgb.root", "tree"); 
	filemanager.AddItem("DataLargeMVASimple", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_converted_mvaxgbsimple.root", "tree"); 
	filemanager.AddItem("MCvalidation", "/eos/home-m/mhuwiler/data/Analysis/v9/TauCutflowSample_mva.root", "tree"); 
	filemanager.AddItem("signalTrain", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "signalTrain"); 
	filemanager.AddItem("signalTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "signalTest"); 
	filemanager.AddItem("backgroundTrain", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "backgroundTrain"); 
	filemanager.AddItem("backgroundTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root", "backgroundTest"); 
	filemanager.AddItem("DataBackground", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge_converted.root", "tree"); 
	filemanager.AddItem("MCSignal", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/PrivateProductionGenDstar_converted.root", "tree"); 
	filemanager.AddItem("SignalOfficialMC50M", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M_converted.root", "tree"); 
	filemanager.AddItem("ParkingBPHAllRun2018B", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPHRun2018B_converted.root", "tree"); 


	filemanager.OpenAllItems(); 

	gStyle->SetOptStat(0); 


	auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>("DataLargeMVASimple")); // tree100k
	

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

	TCanvas *trainingperfcanvas = PlotROCStuff(TFile::Open("/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/newtest/plots.root")); 

	//auto frame2 = dataframe.Define("P_D0", P_v, {"BsDstarTauNu_D0_pt", "BsDstarTauNu_D0_eta", "BsDstarTauNu_D0_phi", "BsDstarTauNu_D0_mass"}) // auto frame2 = dataframe.Define("LV_D0", "TLorentzVector LV_D0; LV_D0.SetPtEtaPhiM(BsDstarTauNu_D0_pt, BsDstarTauNu_D0_eta, BsDstarTauNu_D0_phi, BsDstarTauNu_D0_mass); return LV_D0"); 
	//				.Define("P_Ds", P_v, {"BsDstarTauNu_Ds_pt", "BsDstarTauNu_Ds_eta", "BsDstarTauNu_Ds_phi", "BsDstarTauNu_Ds_mass"})
	//				.Define("P_tau", P_v, {"BsDstarTauNu_tau_pt", "BsDstarTauNu_tau_eta", "BsDstarTauNu_tau_phi", "BsDstarTauNu_tau_mass"})
	//				.Define("B_mass", invMass_v, {"P_Ds", "P_tau"}); 

	

	//auto histo1 = frame2.Histo1D("B_mass"); 

	//auto histo2 = frame2.Histo2D({"Bmass_vs_Dmass", "Correlation plot between B and D masses", 100, 0., 7000., 100, 0., 5000.}, "BsDstarTauNu_B_mass", "BsDstarTauNu_D0_unfit_mass"); 

	std::vector<std::pair<TString, TString> > quantitiesToPlot2D = {{"BsDstarTauNu_B_mass", "BsDstarTauNu_D0_unfit_mass"}, {"BsDstarTauNu_B_mass", "BsDstarTauNu_Ds_unfit_mass"}, {"BsDstarTauNu_B_mass", "BsDstarTauNu_tau_mass"}, {"BsDstarTauNu_tau_mass", "BsDstarTauNu_Ds_unfit_mass"}}; 


	std::vector<TString> quantitiesToPlotMVA = {
		"BsDstarTauNu_D0_pt>>h(50, 0., 30.)", "BsDstarTauNu_D0_eta>>h(50, -3., 3.)", "BsDstarTauNu_D0_phi>>h(50, -3.5, 3.5)", 
		"BsDstarTauNu_Ds_pt>>h(50, 0., 50.)", "BsDstarTauNu_Ds_eta>>h(50, -3., 3.)", "BsDstarTauNu_Ds_phi>>h(50, -3.5, 3.5)",  
		"BsDstarTauNu_tau_pt>>h(100, 0., 50.)", "BsDstarTauNu_tau_eta>>h(100, -3., 3.)", "BsDstarTauNu_tau_phi>>h(100, -3.5, 3.5)", 
		//"tau_doca>>h(100, -7., 2.)", "tau_docaerror>>h(100, 0., 1.)", "tau_docasigma>>h(100, -40., 20.)", "tau_tracks_dR_Dstar", 
		"BsDstarTauNu_D0_vprob>>h(50, -0.1, 1.1)", "BsDstarTauNu_Ds_vprob>>h(50, -0.1, 1.1)", "BsDstarTauNu_tau_vprob>>h(50, -0.1, 1.1)", "BsDstarTauNu_tau_alpha>>h(50, -0.1, 1.1)", 
		"BsDstarTauNu_tau_q>>h(5, -2., 2.)", 
		"BsDstarTauNu_D0_fl3d>>h(200, -2., 5.)", "BsDstarTauNu_D0_fls3d>>h(200, 0., 10.)", "BsDstarTauNu_D0_pvip>>h(200, -1., 1.)", //"BsDstarTauNu_D0_pvips>>h(200, -1., 10.)", 
		"BsDstarTauNu_D0_lip>>h(200, -1., 1.)", "BsDstarTauNu_D0_lips>>h(200, -1., 10.)", 
		"BsDstarTauNu_Ds_fl3d>>h(200, -2., 5.)", "BsDstarTauNu_Ds_fls3d>>h(100, -0.01, 10.)", "BsDstarTauNu_Ds_pvip>>h(200, -1., 1.)", //"BsDstarTauNu_Ds_pvips>>h(200, -1., 10.)", 
		"BsDstarTauNu_Ds_lip>>h(200, -1., 1.)", "BsDstarTauNu_Ds_lips>>h(200, -0.5, 0.5)", 
		"BsDstarTauNu_tau_fl3d>>h(200, -0.5, 0.5)", "BsDstarTauNu_tau_fls3d>>h(200, -0.5, 10.)", "BsDstarTauNu_tau_pvip>>h(200, -0.5, 0.5)", "BsDstarTauNu_tau_pvips>>h(200, -1., 20.)", "BsDstarTauNu_tau_lip>>h(200, -1., 1.)", "BsDstarTauNu_tau_lips>>h(200, -1., 10.)", 
		//"BsDstarTauNu_tau_sumofdnn>>h(200, -0.5, 3.5)", "max(max(BsDstarTauNu_tau_pfidx1, BsDstarTauNu_tau_pfidx2), BsDstarTauNu_tau_pfidx3)>>h(200, -0.5, 15.)", "min(min(BsDstarTauNu_tau_pfidx1, BsDstarTauNu_tau_pfidx2), BsDstarTauNu_tau_pfidx3)>>h(200, -0.5, 15.)",
		"BsDstarTauNu_tau_max_dr_3prong>>h(100, 0., 1.)", 
		"BsDstarTauNu_mu1_q>>h(5, -2., 2.)", "BsDstarTauNu_mu1_vx>>h(50, -10., 10.)", "BsDstarTauNu_mu1_vy>>h(50, -10., 10.)", "BsDstarTauNu_mu1_vz>>h(50, -10., 10.)", 
		"BsDstarTauNu_tau_pi1_pt>>h(50, 0., 15.)", "BsDstarTauNu_tau_pi1_eta>>h(50, -3., 3.)", "BsDstarTauNu_tau_pi1_phi>>h(50, -3.5, 3.5)", "BsDstarTauNu_tau_pi1_charge>>h(5, -2., 2.)", 
		"BsDstarTauNu_tau_pi2_pt>>h(50, 0., 15.)", "BsDstarTauNu_tau_pi2_eta>>h(50, -3., 3.)", "BsDstarTauNu_tau_pi2_phi>>h(50, -3.5, 3.5)", "BsDstarTauNu_tau_pi2_charge>>h(5, -2., 2.)", 
		"BsDstarTauNu_tau_pi3_pt>>h(50, 0., 15.)", "BsDstarTauNu_tau_pi3_eta>>h(50, -3., 3.)", "BsDstarTauNu_tau_pi3_phi>>h(50, -3.5, 3.5)", "BsDstarTauNu_tau_pi3_charge>>h(5, -2., 2.)", 
		"BsDstarTauNu_k_charge>>h(5, -2., 2.)", "BsDstarTauNu_pi_charge>>h(5, -2., 2.)", "BsDstarTauNu_spi_charge>>h(5, -2., 2.)"
		}; //"BsDstarTauNu_B_unfit_mass",  "BsDstarTauNu_Ds_unfit_mass-BsDstarTauNu_Ds_unfit_mass" 

	// generic options class as dictionary, scripts default overriden by provided options
	// Site for referencing plots made with jsroot 



	std::vector<TString> formatsToPlot = {".pdf", ".png"}; 

	//TString campaignName = "PlotsDstarDsVsSignal/"; 

	bool normalise = true; 

	TString outfolder = "plots/"+campaignName; //plots/"+campaignName; //"plots/PlotsGenmactchedFinal/"

	bool webpublication = true; 

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

	trainingperfcanvas->Print(outfolder+"OvertrainingAssessement.pdf"); 
	trainingperfcanvas->Print(webfolder+"OvertrainingAssessement.pdf"); 


	std::ofstream outputfile; 
    outputfile.open("VariableSeparationMVA.txt");

	for (auto quantity : quantitiesToPlotMVA) //auto quantity : quantitiesToPlotFromTree
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

		filemanager.GetItem<TTree*>("ParkingBPHAllRun2018B")->Draw(quantityData.ReplaceAll(">>h", ">>h1"), cut); 
		TH1 *histoData = static_cast<TH1*>(canvas->GetPrimitive("h1")); 
		histoData->GetYaxis()->SetTitleOffset(0.9); 
		histoData->SetTitle(""); 
		filemanager.GetItem<TTree*>("SignalOfficialMC50M")->Draw(quantity.ReplaceAll(">>h", ">>h2"), cut); 
		TH1 *histoMC = static_cast<TH1*>(canvas->GetPrimitive("h2")); 
		histoMC->SetTitle(""); 

		Double_t sep = GetSeparation(*static_cast<TH1D*>(histoMC), *static_cast<TH1D*>(histoData)); 
		std::cout << "Separation: " << name << ": " << sep << std::endl; 
		outputfile << name << ": " << sep << std::endl; 

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


		std::cout << "Number of MC events: " << histoMC->Integral() << std::endl; 

		//histoMC->Scale(histoData->Integral()/histoMC->Integral()); 
		if (normalise) histoData->Scale(histoMC->Integral()/histoData->Integral()); 

		histoData->SetLineColor(kBlue); 
		histoMC->SetLineColor(kRed); 

		histoData->SetLineWidth(2); 
		histoMC->SetLineWidth(2); 


		TLegend *legend= new TLegend( canvas->GetLeftMargin()+0.35, 
                                    1-canvas->GetTopMargin()-.15, 
                                    //canvas->GetLeftMargin()+.4, 
                                    canvas->GetLeftMargin()+(1.-(canvas->GetLeftMargin()+canvas->GetRightMargin())),
                                    1-canvas->GetTopMargin() );
      	legend->SetFillStyle(1);
      	legend->AddEntry(histoData,"data","F");
      	legend->AddEntry(histoMC,"signal MC (genmatched)","F");
      	legend->SetBorderSize(1);
      	legend->SetMargin( 0.3 );
      	legend->SetTextSize(0.04);

      	histoData->SetLineColor(kBlue); 
      	histoMC->SetLineColor(kRed); 

      	//histoData->SetTitle(name); 

      	Float_t sc = 1.3;
      	//if not (histo->GetMaximum())
      	histoData->SetMaximum( TMath::Max(histoData->GetMaximum(), histoMC->GetMaximum())*sc );
      	std::cout << histoData->GetMaximum() << " " << histoMC->GetMaximum() << std::endl; 

      	histoData->Draw("HIST"); 
      	histoMC->Draw("HISTSAME"); 

      	legend->Draw(); 


		//histo->GetXaxis()->SetRangeUser(0., 100.); 

		canvas->Draw(); 

		for (auto format : formatsToPlot) 
		{
			canvas->Print(outfolder+name+format); 
			if (webpublication) canvas->Print(webfolder+name+format); 
		}

		delete canvas; 
	}

	outputfile.close(); 

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


}

