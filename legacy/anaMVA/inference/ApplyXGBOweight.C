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
#include "../../FileFlow.h"
#include "Python.h"
#include "TPython.h"
#include <numpy/arrayobject.h>
#include "PythonInterface.h"
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


using namespace ROOT; 
using namespace Ana; 


// To run this macro, the python class TFEvaluation.py needs to be loaded into root prior to execution. e.g.:
// root -e 'TPython::LoadMacro("TFEvaluation.py");' ApplyTFweight.C


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

 

void ApplyXGBOweight(const TString& inIdentifier, const TString& outIdentifier, const TString& version = "", const Int_t Nmax = 0, const Int_t Nmin = 0, const TString& destination = "") 
{
	ROOT::DisableImplicitMT(); 
	Init(version); 

	// Hack to create an intermediate file for applications
	TString newInIdentifier = inIdentifier+"_tmp"; 
	TString tempFileName = TString(Ana::filemanager.GetFile(inIdentifier)).ReplaceAll("mva.root", "tmp.root"); 
	filemanager.AddItem(newInIdentifier, tempFileName, "tree"); 
	gSystem->Exec(TString::Format("cp %s %s", filemanager.GetFile(inIdentifier).c_str(), filemanager.GetFile(newInIdentifier).c_str())); 


	filemanager.OpenItem(newInIdentifier); 

	gStyle->SetOptStat(0); 


	auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>(newInIdentifier)); 


    //std::cout << "Before making class" << std::endl; 

	PythonInterface pyEvaluation("PyXGBEval"); 

    //std::cout << "After making class" << std::endl; 

    std::string TMVAweightFile = std::string(TString::Format("../../%s", Ana::MVA[version.Data()].c_str()).Data()); 

    // Hack to get the pkl file from the xml path 
    std::cout << "TMVA weight file: " << TMVAweightFile << std::endl; 

    std::string pklFile = TString(TMVAweightFile).ReplaceAll("/weights.xml", ".pkcl").ReplaceAll("model_", "models/model_").Data(); 

    std::cout << "pkl weight file: " << pklFile << std::endl; 

    pyEvaluation.Initialise(pklFile);

    int counter = 0; 

	//auto histo1 = frame2.Histo1D("B_mass"); 

	//auto histo2 = frame2.Histo2D({"Bmass_vs_Dmass", "Correlation plot between B and D masses", 100, 0., 7000., 100, 0., 5000.}, "BsDstarTauNu_B_mass", "BsDstarTauNu_D0_unfit_mass"); 

 //    struct expand_type 
 //    {
 //  		template<typename... T>
 //  		expand_type(T&&...) {}
	// };

	//template<typename... ArgTypes>
	//auto TFresponse = [&pyEvaluation, &counter](float D0pt, float D0eta, float D0phi, float D0vprob, float D0fl, float D0fsig, float Dstarpt, float Dstareta, "b_Ds_phi", "b_Ds_vprob", "b_Ds_fl", "b_Ds_fsig", "b_D0_lip", "b_D0_lips", "b_D0_pvip", "b_Ds_lip", "b_Ds_lips", "b_Ds_pvip", "b_tau_pt", "b_tau_eta", "b_tau_phi", "b_tau_fl", "b_tau_fsig", "b_tau_vprob", "b_tau_lip", "b_tau_pvip", "b_tau_pvipsig", "b_tau_alpha", "b_tau_legacyMaxdr", "b_tau_pi1pt", "b_tau_pi1eta", "b_tau_pi1phi", "b_tau_pi2pt", "b_tau_pi2eta", "b_tau_pi2phi", "b_tau_pi3pt", "b_tau_pi3eta", "b_tau_pi3phi", "b_tau_sumdnn") 
	auto MVAResponse = [&pyEvaluation, &counter](float D0_pt, float D0_eta, float D0_phi, float D0_vprob, float D0_fl, float D0_fsig, float Ds_pt, float Ds_eta, float Ds_phi, float Ds_vprob, float Ds_fl, float Ds_fsig, float D0_lip, float D0_lips, float D0_pvip, float Ds_lip, float Ds_lips, float Ds_pvip, float tau_pt, float tau_eta, float tau_phi, float tau_fl, float tau_fsig, float tau_vprob, float tau_lip, float tau_pvip, float tau_pvipsig, float tau_alpha, float tau_legacyMaxdr, float tau_pi1pt, float tau_pi1eta, float tau_pi1phi, float tau_pi2pt, float tau_pi2eta, float tau_pi2phi, float tau_pi3pt, float tau_pi3eta, float tau_pi3phi, float tau_sumdnn) 
	{
		std::cout << "Event no: " << counter << std::endl; 
		counter++; 

		std::vector<float> variables = {D0_pt, D0_eta, D0_phi, D0_vprob, D0_fl, D0_fsig, Ds_pt, Ds_eta, Ds_phi, Ds_vprob, Ds_fl, Ds_fsig, D0_lip, D0_lips, D0_pvip, Ds_lip, Ds_lips, Ds_pvip, tau_pt, tau_eta, tau_phi, tau_fl, tau_fsig, tau_vprob, tau_lip, tau_pvip, tau_pvipsig, tau_alpha, tau_legacyMaxdr, tau_pi1pt, tau_pi1eta, tau_pi1phi, tau_pi2pt, tau_pi2eta, tau_pi2phi, tau_pi3pt, tau_pi3eta, tau_pi3phi, tau_sumdnn}; 

		//for (auto element : variables) std::cout << element << std::endl; 

		// create vector saying whether it is a Dstar 
		/*std::vector<float> flag = {1.}; 

		std::vector<std::vector<float>* > vectors = {&eta, &phi, &pt, &q}; 

		// Hack to fit the trained model 
		auto fakePVassoc = new std::vector<float>(dzToClosest.begin(), dzToClosest.begin()+dzToClosest.size()); 
		auto fakeAssoc = new std::vector<float>(dzToClosest.begin(), dzToClosest.begin()+dzToClosest.size()); 

		std::vector<std::vector<float>* > additionalvectors = {fakePVassoc, &DOCA3D, &DOCA2D, &DOCA3DErr, &DOCA2DErr, &dzToPV, fakeAssoc, &dzToClosest}; 
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
		garbageCollector.push_back(fakePVassoc); 
		garbageCollector.push_back(fakeAssoc); 
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

    	/*for (auto element: garbageCollector) 
    	{
    		delete element; 
    	}

		return response; */

		auto response = pyEvaluation.Evaluate(variables); 

		assert(response.size() == 1); 

		return response; 
	};

	auto withWeight = dataframe.Range(0, Nmax).Define("mvaScoreNew", MVAResponse, {"D0_pt", "D0_eta", "D0_phi", "D0_vprob", "D0_fl", "D0_fsig", "Dstar_pt", "Dstar_eta", "Dstar_phi", "Dstar_vprob", "Dstar_fl", "Dstar_fsig", "D0_lip", "D0_lipsig", "D0_pvip", "Dstar_lip", "Dstar_lipsig", "Dstar_pvip", "b_tau_pt", "b_tau_eta", "b_tau_phi", "b_tau_fl", "b_tau_fsig", "b_tau_vprob", "b_tau_lip", "b_tau_pvip", "b_tau_pvipsig", "b_tau_alpha", "b_tau_legacyMaxdr", "b_tau_pi1pt", "b_tau_pi1eta", "b_tau_pi1phi", "b_tau_pi2pt", "b_tau_pi2eta", "b_tau_pi2phi", "b_tau_pi3pt", "b_tau_pi3eta", "b_tau_pi3phi", "b_tau_sumdnn"}); 

	TString outfile = filemanager.GetFile(inIdentifier); 
	
	if (destination != "") 
	{
		auto tokens = outfile.Tokenize("/"); 
		TString outfilename = static_cast<TObjString*>(tokens->At(tokens->GetEntries()-1))->GetString(); 
		std::cout << "File name written out: " << outfilename << std::endl; 

		outfile = destination + outfilename; 
	}

	#include "../../stringbranches.gcf"

	for (auto branch : stringbranches) // Hack to fix string branches 
	{
		withWeight = withWeight.Redefine(branch, [](const ROOT::RVec<std::string> &v) {return std::vector<std::string>(v.begin(), v.end());}, {branch}); 
	}

	withWeight = withWeight.Redefine("v_taucandidates", [](const ROOT::RVec<Tau> &v) {return std::vector<Tau>(v.begin(), v.end());}, {"v_taucandidates"});

	withWeight.Snapshot(filemanager.GetObject(inIdentifier), outfile.Data()); 

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	//gSystem->Exec(TString::Format("rm %s", filemanager.GetFile(newInIdentifier).c_str())); 

	filemanager.CloseAll(); 


}

