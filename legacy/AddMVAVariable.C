#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TRandom3.h>
#include <TSystem.h>
#include <TPRegexp.h>
#include <TMVA/Reader.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>



struct TMVAWeightfileVariables 
{
    std::vector<std::pair<TString, TString> > variables; 

    void load(const std::string& filename) 
    {
        boost::property_tree::ptree propertyTree; 

        boost::property_tree::read_xml(filename, propertyTree); 

        std::cout << "File name: " << filename << std::endl; 


        for (auto variable : propertyTree.get_child("MethodSetup.Variables")) 
        {
            if (variable.first == "<xmlattr>") continue; 
            // Loading the xml attributes Expression (variable name) and Type into the vector of pairs 
            variables.push_back(std::make_pair<TString, TString>(variable.second.get("<xmlattr>.Label", "").data(), variable.second.get("<xmlattr>.Type", "kNoType").data())); 
        }


    }
};


void AddMVAVariable(const TString infile, const TString trees, const TString weightfile = "newtest/model_optimized/weights.xml", TString branchName = "mvaScore", TString suffix = "_mva") 
{
	TFile *inFile = TFile::Open(infile.Data(), "READ"); 

	TTree *tree = static_cast<TTree*>(inFile->Get(trees)); 


	std::vector<TString> variables = {"BsDstarTauNu_mu1_q/I", "BsDstarTauNu_D0_pt/F", "BsDstarTauNu_D0_eta/F", "BsDstarTauNu_D0_phi/F", "BsDstarTauNu_D0_vprob/F", "BsDstarTauNu_D0_fl3d/F", "BsDstarTauNu_D0_fls3d/F",
	"BsDstarTauNu_Ds_pt/F", "BsDstarTauNu_Ds_eta/F", "BsDstarTauNu_Ds_phi/F", "BsDstarTauNu_Ds_vprob/F", "BsDstarTauNu_Ds_fl3d/F", "BsDstarTauNu_Ds_fls3d/F",
	"BsDstarTauNu_D0_lip/F", "BsDstarTauNu_D0_lips/F", "BsDstarTauNu_D0_pvip/F", "BsDstarTauNu_Ds_lip/F", "BsDstarTauNu_Ds_lips/F", "BsDstarTauNu_Ds_pvip/F", 
	"BsDstarTauNu_tau_pt/F", "BsDstarTauNu_tau_eta/F", "BsDstarTauNu_tau_phi/F", "BsDstarTauNu_tau_q/I", "BsDstarTauNu_tau_fl3d/F", "BsDstarTauNu_tau_fls3d/F", "BsDstarTauNu_tau_vprob/F", 
	"BsDstarTauNu_tau_lip/F", "BsDstarTauNu_tau_lips/F", "BsDstarTauNu_tau_pvip/F", "BsDstarTauNu_tau_pvips/F", "BsDstarTauNu_tau_alpha/F", "BsDstarTauNu_tau_max_dr_3prong/F", 
	"BsDstarTauNu_tau_pi1_pt/F", "BsDstarTauNu_tau_pi1_eta/F", "BsDstarTauNu_tau_pi1_phi/F", "BsDstarTauNu_tau_pi1_charge/F", 
	"BsDstarTauNu_tau_pi2_pt/F", "BsDstarTauNu_tau_pi2_eta/F", "BsDstarTauNu_tau_pi2_phi/F", "BsDstarTauNu_tau_pi2_charge/F", 
	"BsDstarTauNu_tau_pi3_pt/F", "BsDstarTauNu_tau_pi3_eta/F", "BsDstarTauNu_tau_pi3_phi/F", "BsDstarTauNu_tau_pi3_charge/F", 
	"BsDstarTauNu_k_charge/F", "BsDstarTauNu_pi_charge/F", "BsDstarTauNu_spi_charge/F", 
	"BsDstarTauNu_mu1_vx/F", "BsDstarTauNu_mu1_vy/F", "BsDstarTauNu_mu1_vz/F"}; 

	std::unique_ptr<TMVA::Reader> reader = std::make_unique<TMVA::Reader>("!Color:Silent"); 

	Int_t numVars = variables.size(); 

	Float_t floatVars[numVars]; 
    Int_t intVars[numVars]; 

    Float_t readerVars[numVars]; 

    std::vector<std::vector<float>* > floatVecs(numVars); 
    std::vector<std::vector<Int_t>* > intVecs(numVars); 

    std::vector<float> *var1 = nullptr; 

    std::vector<int> varsToCast; 

    tree->SetBranchAddress("BsDstarTauNu_D0_pt", &var1); 

	for (Long64_t i=0; i<numVars; i++) //auto variable : variables 
	{
		auto variable = variables.at(i); 

		auto components = variable.Tokenize("/"); 

		assert(components.GetEntries() == 2); 

		TString name = static_cast<TObjString*>(components->At(0))->GetString(); 
		TString type = static_cast<TObjString*>(components->At(1))->GetString(); 

		std::cout << "Variable: " << name << " of type: " << type << std::endl; 

		if (type == "F") 
		{ 
			floatVecs[i] = new std::vector<Float_t>; 
			floatVecs[i]->push_back(-999.); // Hack: To have at least a variable
			tree->SetBranchAddress(name, &floatVecs[i]); //&floatVars[i] &readerVars[i]
		}
		else if (type == "I") 
		{
			//auto vec = new std::vector<Int_t>(); 
			//intVecs.push_back(vec); 
			intVecs[i] = new std::vector<Int_t>; 
			tree->SetBranchAddress(name, &intVecs[i]); //&intVars[i]
			varsToCast.push_back(i); 
		}
		else 
		{
			std::cerr << "ERROR: Unsupported type: " << type << " for variable " << name << std::endl; 
			return; 
		}



		reader->AddVariable(name, &readerVars[i]); 

	}

	reader->BookMVA("BDT", weightfile); 

	TFile *outFile = TFile::Open(TString(infile).ReplaceAll(".root", suffix+".root"), "RECREATE"); 

	tree->SetBranchStatus("*", 1); 

	TTree *newTree = tree->CloneTree(0); 

	//newTree->SetDirectory(0); // Leaving it in memory 

    Float_t weightVar = 0.; 
    TBranch *branch = newTree->Branch(branchName, &weightVar, branchName+"/F"); 

	for (Long64_t i=0; i<tree->GetEntries(); i++) 
	{
		tree->GetEntry(i); 
		//for (auto index : varsToCast) 
		//{
		//	readerVars[i] = static_cast<Float_t>(intVars[i]); 
		//}
		for (auto index : varsToCast) 
		{
			auto vec = intVecs.at(index); 
			floatVecs[index] = new std::vector<Float_t>(vec->begin(), vec->end()); 
		}
		std::cout << "Size of vector: " << var1->size() << std::endl; 
		//std::cout << "Var 1: " << var1->at(0) << std::endl; 

		std::cout << "Size of branch vector: " << floatVecs.at(0)->size() << std::endl; 
		std::cout << "Value: " << floatVecs.at(0)->at(0) << std::endl; 

		/*for (int idx=0; idx<varsToCast.size(); idx++) 
		{
			auto index = varsToCast.at(idx); 
			auto vec = intVecs.at(idx); 
			//delete floatVecs[index]; 
			//floatVecs[index] = std::vector<Float_t>(vec->begin(), vec->end()); // Cast the vector of int to vector of float 
			std::cout << "Size of int vector: " << vec->size() << std::endl; 
		}*/
		//readerVars[i] = floatVars[i]; 

		// CRAZYYYYY hack fix this bullcrap! 
		assert(floatVecs.size() == numVars); 
		std::cout << "Vector size: " << floatVecs.size() << " num float: " << numVars << std::endl; 
		for (int i=0; i < floatVecs.size(); i++) 
		{
			assert(floatVecs.at(i)->size()>0); 
			std::cout << "Size of variable vector: " << floatVecs.at(i)->size() << std::endl; 
			readerVars[i] = floatVecs.at(i)->at(0); 
		}

		weightVar = reader->EvaluateMVA("BDT"); 

		newTree->Fill(); 

	}

	outFile->Write(); 
	outFile->Close(); 

	inFile->Close(); 

	// Clean up memory 
	for (int i=0; i<floatVecs.size(); i++) 
	{
		delete floatVecs[i]; 
	}
	for (int i=0; i<intVecs.size(); i++) 
	{
		delete intVecs[i]; 
	}
}

