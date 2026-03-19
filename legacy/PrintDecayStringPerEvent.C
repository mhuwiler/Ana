#include "TChain.h"
#include "TFile.h"



void PrintDecayString(TTree *tree, bool fulldecay = false, int maxnum = -1, bool skipincomplete = false) 
{

	std::string *decaystring = new std::string(); 
	std::string *decaystring1 = new std::string(); 
	std::string *decaystring2 = new std::string(); 
	std::string *decaystring3 = new std::string(); 
	std::string *decaystringK = new std::string(); 
	std::string *decaystringpi = new std::string(); 
	std::string *decaystringspi = new std::string(); 

	tree->SetBranchAddress("pttau_tau_gen1str", &decaystring1); 
	tree->SetBranchAddress("pttau_tau_gen2str", &decaystring2); 
	tree->SetBranchAddress("pttau_tau_gen3str", &decaystring3); 
	tree->SetBranchAddress("D0_genkstr", &decaystringK); 
	tree->SetBranchAddress("D0_genpistr", &decaystringpi); 
	tree->SetBranchAddress("Dstar_genpistr", &decaystringspi); 

	tree->SetBranchAddress("genstring", &decaystring); 

	int max = static_cast<int>(tree->GetEntries()); 
	if (maxnum > 0) 
	{
		max = min(max, maxnum); 
	}
	for (int i=0; i<max; i++) 
	{
		tree->GetEntry(i);

		if (skipincomplete && (*decaystringK == "" || *decaystringpi == "" || *decaystringspi == "" || *decaystring1 == "" || *decaystring2 == "" || *decaystring3 == "" || *decaystring1 == "|" || *decaystring2 == "|" || *decaystring3 == "|" )) continue; 

		for (unsigned int i=0; i<45; i++) 
		{
			std::cout << "-"; 
		}
		std::cout << endl; 

		if (fulldecay) std::cout << *decaystring << std::endl; 

		std::cout << "K: " << *decaystringK << std::endl; 
		std::cout << "pi: " << *decaystringpi << std::endl; 
		std::cout << "spi: " << *decaystringspi << std::endl; 
		std::cout << "pi1: " << *decaystring1 << std::endl; 
		std::cout << "pi2: " << *decaystring2 << std::endl; 
		std::cout << "pi3: " << *decaystring3 << std::endl; 
	}
	std::cout << std::endl << "Analyzed " << max << " events. " << std::endl; 
}


void PrintDecayStringPerEvent() 
{
	TFile *file = TFile::Open("../../data/v5/BkgDstarDsPart.root", "READ"); 
	TTree *tree = static_cast<TTree*>(file->Get("ntuplizer/tree")); 

	PrintDecayString(tree); 


	file->Close(); 
}

