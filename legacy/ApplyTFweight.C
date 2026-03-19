#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include "DrawTMVAHistogram.C"
#include "GetSeparation.C"
#include "FileFlow.h"
#include "Python.h"
#include "TPython.h"
#include <numpy/arrayobject.h>
#include "PythonInterface.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


using namespace ROOT; 
using namespace Ana; 


// To run this macro, the python class TFEvaluation.py needs to be loaded into root prior to execution. e.g.:
// root -e 'TPython::LoadMacro("TFEvaluation.py");' ApplyTFweight.C



TLorentzVector LV(double pt, double eta, double phi, double m) 
{
	TLorentzVector V; 
	V.SetPtEtaPhiM(pt, eta, phi, m); 
	return V; 
}

std::vector<float> EvaluateTFresponse(std::vector<float> pt, std::vector<float> eta, std::vector<float> phi, std::vector<float> q, std::vector<float> DOCA2D, std::vector<float> DOCA2DErr, std::vector<float> DOCA3D, std::vector<float> DOCA3DErr, std::vector<float> dzToPV, std::vector<float> dzToClosest, std::vector<float> isAssociate, std::vector<float> assocQualityToPV, std::vector<int> genmatch) 
{
	std::vector<float> response; 


	response.push_back(-999.); 
	return response; 
}

std::vector<float>  extendArray(const std::vector<float>& array, const int dim) 
{
	std::vector<float> result; 
	result.reserve(dim); 
	if (array.size() > dim) 
	{
		result = std::vector<float>(array.begin(), array.begin()+dim); 
	}
	else 
	{
		result = array; 
		while(result.size() < dim)
		{
			result.push_back(0.); 
		}
	}
	return result; 
}

// TODO: function that does the same inplace 

template<typename T>
void PrintArray(const std::vector<T>& vec) 
{
	std::cout << "Vector content: "; 
	for (auto i : vec) 
	{
		std::cout << i << ", "; 
	}
	std::cout << std::endl; 
}

std::vector<float> concatenateVectors(const std::vector<std::vector<float>* > vectors) 
{
	std::vector<float> result; // TODO: reserve the size 
	for (auto vec : vectors) 
	{
		result.insert(result.end(), vec->begin(), vec->end()); 
	}
	return result; 
}

void extendArray(std::vector<float>*& array, const int dim, std::vector<std::vector<float>* >& garbageCollector) 
{
	if (array->size() > dim) 
	{
		array = new std::vector<float>(array->begin(), array->begin()+dim); 
		garbageCollector.push_back(array); 
	}
	else 
	{
		array->insert(array->end(), dim - array->size(), 0.); 
	}
}

std::vector<double> castVector(std::vector<float> vec) 
{
	std::vector<double> result; 
	result.reserve(vec.size()); 
	for (auto element : vec)
	{
		result.push_back(static_cast<double>(element)); 
	}
	return result; 
}

 

void ApplyTFweight(const TString& identifier, const TString& version = "", const Int_t Nmax = 0, const Int_t Nmin = 0, const TString& destination = "") 
{
	//ROOT::EnableImplicitMT(); //ROOT::DisableImplicitMT(); 
	Init(version); 


	TString inIdentifier = identifier+"_ntuple"; 
	TString outIndentifier = identifier+"_tf"; 
	filemanager.AddItem("effInfo", filemanager.GetFile(inIdentifier), "ntuplizer/EffCalc"); 

	filemanager.OpenItem(inIdentifier); 

	gStyle->SetOptStat(0); 


	auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>(inIdentifier)); 


    //std::cout << "Before making class" << std::endl; 

	PythonInterface pyEvaluation("PyTFEval"); 

    //std::cout << "After making class" << std::endl; 

    pyEvaluation.Initialise(std::string(Ana::model[version.Data()]), 10, 20);

    int counter = 0; 

	//auto histo1 = frame2.Histo1D("B_mass"); 

	//auto histo2 = frame2.Histo2D({"Bmass_vs_Dmass", "Correlation plot between B and D masses", 100, 0., 7000., 100, 0., 5000.}, "BsDstarTauNu_B_mass", "BsDstarTauNu_D0_unfit_mass"); 

	auto TFresponse = [&pyEvaluation, &counter](float Dstarpt, float Dstareta, float Dstarphi, int Dstarcharge, std::vector<float> pt, std::vector<float> eta, std::vector<float> phi, std::vector<float> q, std::vector<float> DOCA2D, std::vector<float> DOCA2DErr, std::vector<float> DOCA3D, std::vector<float> DOCA3DErr, std::vector<float> dzToPV, std::vector<float> dzToClosest, std::vector<float> isAssociate, std::vector<float> assocQualityToPV, std::vector<float> doca2DToPV, std::vector<int> genmatch) 
	{
		//assert(Dstarpt.size() == 1); 
		//assert(Dstareta.size() == 1); 
		//assert(Dstarphi.size() == 1); 
		//assert(Dstarcharge.size() == 1); 

		const int initialSize = pt.size(); 

		pt.insert(pt.begin(), Dstarpt); 
		eta.insert(eta.begin(), Dstareta); 
		phi.insert(phi.begin(), Dstarphi); 
		//q.insert(q.begin(), Dstarcharge[0]); 

		if ((counter % 1000 == 0) || ((counter % 100 == 0) && (counter < 1000)) || (counter < 11)) std::cout << "Event no: " << counter << std::endl; 
		counter++; 

		// create vector saying whether it is a Dstar 
		std::vector<float> flag = {1.}; 

		std::vector<std::vector<float>* > vectors = {&eta, &phi, &pt}; // , &q

		// Hack to fit the trained model 
		//auto fakePVassoc = new std::vector<float>(dzToClosest.begin(), dzToClosest.begin()+dzToClosest.size()); 
		//auto fakeAssoc = new std::vector<float>(dzToClosest.begin(), dzToClosest.begin()+dzToClosest.size()); 

		std::vector<std::vector<float>* > additionalvectors = {&doca2DToPV, &assocQualityToPV, &DOCA3D, &DOCA2D, &DOCA3DErr, &DOCA2DErr, &dzToPV, &isAssociate, &dzToClosest}; 
		// End hack 

		for (auto vec : additionalvectors) 
		{
			vec->insert(vec->begin(), 0.); 
		}

		vectors.insert(vectors.end(), additionalvectors.begin(), additionalvectors.end()); 

		vectors.insert(vectors.end(), &flag); 

		std::vector<std::vector<float>* > garbageCollector; 
		for (auto& vec : vectors) 
		{
			extendArray(vec, 20, garbageCollector); 
		}

		// Hack to fit the trained model
		//garbageCollector.push_back(fakePVassoc); 
		//garbageCollector.push_back(fakeAssoc); 
		// End hack 

		auto concatenated = concatenateVectors(vectors); 

		//PrintArray(concatenated); 

		assert(concatenated.size() = 20*12); 

		//auto response = pyEvaluation.EvaluateArray(datavec);
		auto response = pyEvaluation.Evaluate(concatenated);

		//std::cout << "Response size: " << response.size() << std::endl; 

		response.erase(response.begin()); 

		response.resize(initialSize); 

		//std::cout << "Response size: " << response.size() << std::endl; 

    	//std::cout << "After evaluation" << std::endl; 

    	/*for (auto element : response) 
    	{
        	std::cout << element << ", "; 
    	}
    	std::cout << std::endl; */

    	for (auto element: garbageCollector) 
    	{
    		delete element; 
    	}

		return response; 
	};

	auto withWeight = dataframe.Range(0, Nmax).Define("TFscore", TFresponse, {"Dstar_pt", "Dstar_eta", "Dstar_phi", "Dstar_q", "track_pt", "track_eta", "track_phi", "track_charge", "track_doca2D", "track_doca2Derror", "track_doca", "track_docaerror", "track_dzToPV", "track_dzToClosestVertex", "track_isAssociatedToPV", "track_pvAssociationQuality", "track_PV_doca2D", "track_isgenmatched"}); 

	TString outfile = filemanager.GetFile(outIndentifier); 
	
	if (destination != "") 
	{
		auto tokens = outfile.Tokenize("/"); 
		TString outfilename = static_cast<TObjString*>(tokens->At(tokens->GetEntries()-1))->GetString(); 
		std::cout << "File name written out: " << outfilename << std::endl; 

		outfile = destination + outfilename; 
	}

	#include "stringbranches.gcf"

	for (auto branch : stringbranches) // Hack to fix string branche 
	{
		withWeight = withWeight.Redefine(branch, [](const ROOT::RVec<std::string> &v) {return std::vector<std::string>(v.begin(), v.end());}, {branch}); 
	}


	withWeight.Snapshot(filemanager.GetObject(outIndentifier), outfile.Data()); 

	std::cout << "Processed " << counter << " events. " << std::endl;

	// Update the created file with eff info
	TFile *output = TFile::Open(outfile.Data(), "UPDATE"); 
	output->cd("ntuplizer");
	filemanager.GetItem<TTree*>("effInfo", true)->CloneTree(); 
	output->Write(); 
	output->Close();

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


}

