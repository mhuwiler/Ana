#ifndef FileFlow_h
#define FileFlow_h
#include "plugins/FileManager/CFileManager.C"
#include <string>
#include "TString.h"
#include "TCut.h"
#include <unordered_map>
#include "SampleData.C"
#include "ROOT/RDataFrame.hxx"
#define BOOST_JSON_STACK_BUFFER_SIZE 1024
#include <boost/json/src.hpp> // Currently not used, remove?
#include <json/value.h>
#include <json/json.h>
//R__ADD_LIBRARY_PATH($FOODIR) // if needed
R__LOAD_LIBRARY(/opt/local/lib/libjsoncpp.dylib) // Load the library
#include <fstream>


//gROOT->Load("/opt/local/lib/libjsoncpp.dylib")


namespace Ana 
{

	FileManager filemanager; 

	std::unordered_map<std::string, Int_t> colorold; 

	std::unordered_map<std::string, ROOT::RDF::TH1DModel> binning; 

	std::unordered_map<std::string, SampleData> samples; 

	std::unordered_map<std::string, std::string> FinalBDT; 


	std::vector<TString> LoadRegions(const TString& path) 
	{
		std::ifstream file(path+"/Info.json", std::ifstream::binary);
		Json::Value regions;
		file >> regions;

		std::vector<TString> loaded; 

		for (auto item : regions.getMemberNames()) 
		{
			//std::cout << item << std::endl; 

			for (auto region : regions[item].getMemberNames()) 
			{
				auto& regiondata = regions[item][region]; 
				assert(regiondata.size() == 2); 

				/*for (auto element : regiondata) 
				{
					std::cout << element << std::endl; 
				}*/

				TString name = TString::Format("%s_%s", item.c_str(), region.c_str());
				filemanager.AddItem(name, regiondata[0].asString(), regiondata[1].asString()); 
				loaded.push_back(name);
			}
		}

		/*for (auto item : regions)
		{
			std::cout << item << " " << std::endl; 
		}*/


		return std::move(loaded); 
	}


	void CommonInitialisation()
	{
	
		colorold = {{"Sig", 2}, {"BkgDstarDs", 3}, {"BkgDstarDsstar", 8}, {"BkgDstara1", 4}, {"dataD2WS", 6}, {"dataD2TauWS", 7}, {"other", 9}, {"yetanother", 1}}; // Legacy color scheme 


		Int_t nBins = 20; 

		binning = {{"b_tau_rhomass1", {"", "#rho_{12} mass;Invariant m_{#rho} [GeV];Counts", nBins, 0., 1.5}},
				{"b_tau_rhomass2", {"", "rho_{23} mass;Invariant m_{#rho} [GeV];Counts", nBins, 0., 1.5}},
				{"b_B_m", {"", "B mass;Reconstructed m_{B} [GeV];Counts", nBins, 2., 6.}},
				{"b_B_q2", {"", "q2;q^{2} [GeV];Counts", nBins, 0., 12.}},
				//{"B_m", {"", ";B mass [GeV];Counts", nBins, 0., 6.}},
				//{"B_q2", {"", ";B mass [GeV];Counts", nBins, 0., 12.}},
				//{"tau_rhomass1", {"", ";#rho_{12} mass [GeV];Counts", nBins, 0., 1.5}},
				//{"tau_rhomass2", {"", ";#rho_{12} mass [GeV];Counts", nBins, 0., 1.5}},
				{"b_B_proper_xi_rho1", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
				{"b_B_proper_xi_rho2", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
				{"b_tau_proper_alpha_rho1_pi", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
				{"b_tau_proper_alpha_rho2_pi", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
				{"b_tau_proper_theta_rho1", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
				{"b_tau_proper_theta_rho2", {"", "#tau mass;Reconstructed m_{#tau} [GeV];Counts", nBins, -1.1, 1.1}},
		}; 

		// Complete binning with default
		std::unordered_map<std::string, ROOT::RDF::TH1DModel> defaultbinning = { 
			#include "VariableDefinitions.ccf" 
		}; 
		for (auto item : defaultbinning) 
		{
			if (binning.find(item.first) == binning.end()) 
			{
				binning[item.first] = item.second; 
			}
		}

		samples = InitSamples(); 

		// Filling the versions of final MVA discriminator
		FinalBDT = {	
			{"v6.9", "./anaMVA/MVAinSBdata/model_optimized/weights.xml"}, //LatestVsData
			{"v6.95", "./anaMVA/MVAinSBdata/model_optimized/weights.xml"} //LatestVsData
		}; 

		//if (regions != "") LoadRegions(regions); 

	}


	std::vector<TString> Load(const TString& regions) 
	{
		CommonInitialisation(); 

		auto loaded = LoadRegions(regions); 

		std::cout << "Size: " << loaded.size() << std::endl; 
		for (auto item : loaded) 
		{
			std::cout << "item: "; 
			std::cout << item << std::endl; 
		}

		return loaded; 
	}

}




#endif

