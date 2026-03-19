#ifndef FileFlow_h
#define FileFlow_h
#include "plugins/FileManager/CFileManager.C"
#include <string>
#include "TString.h"
#include "TObjString.h"
#include "TCut.h"
#include <unordered_map>
#include "SampleData.C"
#include "ROOT/RDataFrame.hxx"
//#include <json/value.h>
//#include <json/json.h>
//R__ADD_LIBRARY_PATH($FOODIR) // if needed
//R__LOAD_LIBRARY(/opt/local/lib/libjsoncpp.dylib) // Load the library
//#include "external/jsoncpp/dist/jsoncpp.cpp"
#include <fstream>
#include "AnaBuildingBlocks.C"
#include "CutFlowMonitor.h"
//#include "HHbbtautauAnaElements.C"


namespace Ana 
{

	FileManager filemanager; 

	TString folder = "";

	TString folderbase = "";

	std::map<std::string, TCut> cut; //std::map<std::string, std::map<std::sting, TCut> > cuts; 

	//std::vector<TColor*> colors; 

	//std::vector<Int_t> mycolors; 

	//std::unordered_map<std::string, Int_t> color; 

	std::unordered_map<std::string, Int_t> colorold; 

	std::unordered_map<std::string, std::string> model; 

	std::unordered_map<std::string, std::string> MVA; 

	std::unordered_map<std::string, std::string> FinalBDT; 

	//std::unordered_map<std::string, std::string> legends; 

	std::unordered_map<std::string, ROOT::RDF::TH1DModel> binning; 

	//std::unordered_map<std::string, std::string> labels; 

	RedirectingMap<std::string, SampleData> samples; 


	std::map<std::string, TCut> cutstandalone;

	//ABCDcuts ABCDconfig;  



	void CommonInitialisation(const bool flag)
	{
	
		colorold = {{"Sig", 2}, {"BkgDstarDs", 3}, {"BkgDstarDsstar", 8}, {"BkgDstara1", 4}, {"dataD2WS", 6}, {"dataD2TauWS", 7}, {"other", 9}, {"yetanother", 1}}; // Legacy color scheme 


		Int_t nBins = 50; 

		// binning = {}; // Can add here custom binning that override VariableDefinitions.ccf 

		// Complete binning with default
		std::unordered_map<std::string, ROOT::RDF::TH1DModel> defaultbinning = { 
			#include "VariableDefinitions.ccf" 
		}; 

		std::vector<std::string> excluded; // Can exlude variables from VariableDefinitions.ccf to be added to binning

		
		for (auto item : defaultbinning) 
		{
			if ((binning.find(item.first) == binning.end()) && (std::find(excluded.begin(), excluded.end(), item.first) == excluded.end())) 
			{
				binning[item.first] = item.second; 
			}
		}

		samples = InitSamples(); 
	}


	void Init(const TString& cycle = "", const bool flag = false) 
	{
		if (cycle != "") // For legacy purpose 
		{

			folder = "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/Run3_2023_BoostedPrivate/"; // TODO: set this from environment variable 
			folderbase = "/eos/home-m/mhuwiler/software/rh9/AnaBoosted/data/";

			filemanager.AddItem("ggFofficial", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv12/Run3Summer22NanoAODv12_1-1.root", "Events"); 
			filemanager.AddItem("ggFv15", "/eos/home-m/mhuwiler/data/HHtobbtautau/NanoAODv15/signalggF.root", "Events"); 


			filemanager.AddItem("ggf", folder+"ggf.root", "Events"); 
			filemanager.AddItem("ggfall", folder+"ggFHHtobbtautau.root", "Events"); 
			filemanager.AddItem("ggfhh", folder+"ggFHHtobbtautau_hh", "Events"); 
			filemanager.AddItem("ggfhm", folder+"ggFHHtobbtautau_hm", "Events"); 
			filemanager.AddItem("ggfhe", folder+"ggFHHtobbtautau_he", "Events"); 

			filemanager.AddItem("VBF", folder+"VBF_SM.root", "Events"); 
			filemanager.AddItem("VBFall", folder+"VBFHHtobbtautau.root", "Events"); 
			filemanager.AddItem("VBFhh", folder+"VBFHHtobbtautau_hh.root", "Events"); 
			filemanager.AddItem("VBFhm", folder+"VBFHHtobbtautau_hm.root", "Events"); 
			filemanager.AddItem("VBFhe", folder+"VBFHHtobbtautau_he.root", "Events"); 

			filemanager.AddItem("QCD", folder+"qcd_HT_100-1200.root", "Events"); 
			filemanager.AddItem("QCDall", folder+"QCD.root", "Events"); 
			filemanager.AddItem("QCDhh", folder+"QCD_hh.root", "Events"); 
			filemanager.AddItem("QCDhm", folder+"QCD_hm.root", "Events"); 
			filemanager.AddItem("QCDhe", folder+"QCD_he.root", "Events"); 

			filemanager.AddItem("data", folder+"jetmet.root", "Events"); 
			filemanager.AddItem("dataall", folder+"data.root", "Events"); 
			filemanager.AddItem("datahh", folder+"data_hh.root", "Events"); 
			filemanager.AddItem("datahm", folder+"data_hm.root", "Events"); 
			filemanager.AddItem("datahe", folder+"data_he.root", "Events"); 


		}

		else 
		{
			
		}

		//ABCDconfig = {"(b_tau_sumdnn>1.9)", "(b_tau_sumdnn<1.5)", "b_B_nmu<1&&b_B_ne<1&&b_B_nh<1", "(b_B_nmu<1&&b_B_ne<1)&&b_B_nh>1"}; 

		
		CommonInitialisation(flag); 


	}


	/*std::vector<TString> LoadRegions(const TString& path) 
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
				}

				TString name = TString::Format("%s_%s", item.c_str(), region.c_str());
				filemanager.AddItem(name, regiondata[0].asString(), regiondata[1].asString()); 
				loaded.push_back(name);
			}
		}

		/*for (auto item : regions)
		{
			std::cout << item << " " << std::endl; 
		}


		return std::move(loaded); 
	}*/



	/*std::vector<TString> Load(const TString& regions, const bool denominator = false) 
	{
		CommonInitialisation(denominator); 

		auto loaded = LoadRegions(regions); 

		std::cout << "Size: " << loaded.size() << std::endl; 
		for (auto item : loaded) 
		{
			std::cout << "item: "; 
			std::cout << item << std::endl; 
		}

		return loaded; 
	}*/


	TChain* GetSample(const TString& name) 
	{
		TString samplename = name; 
		TString suffix = "";
		
		TChain *chain = filemanager.GetItem<TChain*>(name, true); 
		if (!chain) 
		{
			if (name.Contains("_")) 
			{
				// Split and get back after
				auto tokens = name.Tokenize("_"); 
				samplename = static_cast<TObjString*>(tokens->At(0))->GetString(); 
				suffix = TString(name).ReplaceAll(samplename, ""); 
			}
			std::cout << samplename << " " << suffix << std::endl;
		
			chain = new TChain(name, name); 
			for (auto item : samples.at(samplename.Data()).fileRefs) 
			{
				TString itemname = TString(item)+suffix; 
				std::cout << itemname << std::endl;
				chain->AddFile(TString::Format("%s/%s", filemanager.GetFile(itemname).c_str(), filemanager.GetObject(itemname).c_str())); 
			}
		}
		std::cout << "Number of events " << chain->GetEntries() << std::endl;

		return chain; 
	}

}




#endif

