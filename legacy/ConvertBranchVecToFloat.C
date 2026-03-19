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

std::vector<std::pair<TString, TString> > GetBranchTypes(const TTree* tree) 
{
	std::vector<std::pair<TString, TString> > result; 

	TObjArray* listOfBranches = const_cast<TTree*>(tree)->GetListOfBranches();   // Returns a TObjArray of branches 

	std::vector<TString> types; 

	for (unsigned int i=0; i< listOfBranches->GetSize() ; i++) 
	{
		std::cout << "Calss name: " << listOfBranches->At(i)->ClassName() << std::endl; 
		TString branchType = listOfBranches->At(i)->ClassName(); 
		if (branchType == "TBranchElement") 
		{
			std::cout << "Branch of type 'TBranchElement', most probably a vector type... " << std::endl; 
			types.push_back(static_cast<TBranchElement*>(listOfBranches->At(i))->GetTypeName()); 
		}
		else if (branchType == "TBranch") 
		{
			std::cout << "Branch of type 'TBranch', most probably a normal type... " << std::endl; 
			types.push_back(static_cast<TLeafObject*>(listOfBranches->At(i))->GetTypeName()); 
		}
		else 
		{
			std::cout << "ERROR (GetBranchTypes): Unkown branch type: " << branchType << std::endl; 
		}

		std::cout << "Branch types: " << types.at(i) << std::endl; 

		TString varType = "", vecType = ""; 
		varType = static_cast<TLeafObject*>(listOfBranches->At(i))->GetTypeName(); 
		vecType = static_cast<TBranchElement*>(listOfBranches->At(i))->GetTypeName(); 
		assert(varType || vecType); 
		std::cout << "Is of a known type: " << (varType || vecType) << std::endl; 
	  	//std::cout << "Branch name: " << listOfBranches->At(i)->GetName() << " branch type: " << static_cast<TBranchElement*>(listOfBranches->At(i))->GetTypeName() << std::endl; 
	  	result.push_back(std::make_pair<TString, TString>(listOfBranches->At(i)->GetName(), static_cast<TLeafObject*>(listOfBranches->At(i))->GetTypeName())); 
	}

	//delete listOfBranches; 
	return result; 
}


void ConvertBranchVecToFloat(const TString& infile, const TString& trees, const TString& outfile) 
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
	"BsDstarTauNu_mu1_vx/F", "BsDstarTauNu_mu1_vy/F", "BsDstarTauNu_mu1_vz/F", "BsDstarTauNu_B_q2/F", 
	// Additional spectator branches 
	"BsDstarTauNu_tau_mass/F", 
	"BsDstarTauNu_B_mass/F", 
	"BsDstarTauNu_Ds_mass/F", "BsDstarTauNu_D0_mass/F", "BsDstarTauNu_Ds_unfit_mass/F", "BsDstarTauNu_D0_unfit_mass/F", 
	"BsDstarTauNu_tau_rhomass1/F", "BsDstarTauNu_tau_rhomass2/F", 
	//"BsDstarTauNu_K_pt/F", "BsDstarTauNu_K_eta/F", "BsDstarTauNu_K_phi/F", "BsDstarTauNu_pi_pt/F", "BsDstarTauNu_pi_eta/F", "BsDstarTauNu_pi_phi/F", 
	//"BsDstarTauNu_pis_pt/F", "BsDstarTauNu_pis_eta/F", "BsDstarTauNu_pis_phi/F", 
	"BsDstarTauNu_tau_sumofdnn/F"}; 

	//auto branchTypes = GetBranchTypes(tree); 

	Int_t numVars = variables.size(); 

	Float_t floatVars[numVars]; 
    Int_t intVars[numVars]; 

    std::vector<std::vector<Float_t>* > floatVecs(numVars); 
    std::vector<std::vector<Int_t>* > intVecs(numVars); 

    Int_t numFloat = 0, numInt = 0; 


    TFile *outFile = TFile::Open(TString(infile).ReplaceAll(".root", "_converted.root"), "RECREATE"); 

    //outFile->mkdir("ntuplizer"); 
    //outFile->cd("ntuplizer"); 

    TTree *newTree = new TTree(tree->GetName(), tree->GetTitle()); 

//    for (auto element : branchTypes) {
//    	std::cout << "Branch: " << element.first << " of type: " << element.second << std::endl; 
//    }


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
			floatVecs[numFloat] = new std::vector<Float_t>; 
			//floatVecs[numFloat]->push_back(-999.); // Hack: To have at least a variable
			tree->SetBranchAddress(name, &floatVecs[numFloat]); //&floatVars[i] &readerVars[i]
			newTree->Branch(name, &floatVars[numFloat], variable); 
			numFloat++; 
		}
		else if (type == "I") 
		{
			//auto vec = new std::vector<Int_t>(); 
			//intVecs.push_back(vec); 
			intVecs[numInt] = new std::vector<Int_t>; 
			tree->SetBranchAddress(name, &intVecs[numInt]); //&intVars[i]
			newTree->Branch(name, &intVars[numInt], variable); 
			numInt++; 
		}
		else 
		{
			std::cerr << "ERROR: Unsupported type: " << type << " for variable " << name << std::endl; 
			return; 
		}

	}

	tree->SetBranchStatus("*", 1); 


	//newTree->SetDirectory(0); // Leaving it in memory 


	assert(floatVecs.size() == numFloat); 
	assert(intVecs.size() == numInt); 

	std::cout << "Variable sizes: (" << floatVecs.size() << ", " << numFloat << "), (" << intVecs.size() << ", " << numInt << ") " << std::endl; 

	for (Long64_t i=0; i<tree->GetEntries(); i++) 
	{
		tree->GetEntry(i); 
		
		// TODO: decide what to do, use all indices or only the first 
		for (int idx = 0; idx<numFloat; idx++) 
		{
			assert(floatVars[idx]); 
			floatVars[idx] = floatVecs.at(idx)->at(0); 
		}
		for (int idx = 0; idx < numInt; idx++) 
		{
			assert(intVars[idx]); 
			intVars[idx] = intVecs.at(idx)->at(0); 
		}

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

