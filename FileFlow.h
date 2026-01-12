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


	// Constants being defined centrally 
	Double_t mvaCutSR = 0.7; // TODO: make cuts per version 
	Double_t mvaCutSB = -0.1; 
	Double_t mvaCutCR = -0.5; 


	void CommonInitialisation(const bool denominator)
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


	void Init(const TString& cycle = "", const bool denominator = false) 
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

		if (!denominator) 
		{
			cut.emplace(std::make_pair("training", "(b_B_fsig>2.)"));
			cut.emplace(std::make_pair("bare", "(b_B_fsig>2.)&&(b_tau_fsig>3.)&&b_tau_vprob>0.1&&(b_tau_alpha>0.)&&(b_B_r<1.)&&(b_B_mu_alpha>1.)&&(b_tau_rhomass_max>0.5)&&(b_tau_m<1.7)&&(b_B_m>3.)&&(b_B_m<5.28)")); // (b_tau_m<1.7)&&(b_B_m>3.)&&(b_B_m<5.5)&&&&(b_tau_min_dr_mu>0.5)&&mvaScore>0.7
			//cut.emplace(std::make_pair("base", cut["bare"]+TCut(TString::Format("(%s)&&%s", ABCDconfig.cutYh.c_str(), ABCDconfig.cutXh.c_str())))); // (b_tau_sumdnn>2.)b_B_nmu<1&&b_B_ne<1&&b_B_nh<1&&b_tau_vprob>0.1(b_tau_m<1.7)&&(b_B_m>3.)&&(b_B_m<5.5)&&&&(b_tau_min_dr_mu>0.5)&&mvaScore>0.7
			// &&b_B_nmu<1&&b_B_ne<1&&b_B_nh<3 &&(b_tau_min_dr_mu>0.5)&&(b_tau_alpha>0)Dstar_vprob>0.1 (b_tau_m<1.7)&&(mvaScore>0.3)&&(b_tau_min_dr_mu>0.5)&&(b_tau_alpha>0)&&(b_B_fsig>2.)&&b_B_m>3.&&b_B_nmu<1&&b_B_nh<1&&b_B_ne<1&&(b_B_npi0+b_B_ngamma)<2&&Dstar_vprob>0.1       (b_tau_m<1.7)&&(mvaScore>0.7)&&(b_tau_min_dr_mu>0.5)&&(b_tau_alpha>0)&&(D0_fsig>3.)&&(b_B_m>3.)     (b_tau_m<1.7)&&(b_tau_min_dr_mu>0.5)&&(b_tau_alpha>0)&&(D0_fsig>3.)&&(b_B_m>3.) //&&(b_tau_min_dr_mu>0.5)&&(b_tau_min_dr_e>0.5)//&& (b_B_q2 > 6.) // (b_Ds_vprob>0.1) && (b_D0_vprob>0.1) && (mvaScore>-2.) && (b_B_mu_alpha > 1.) &&((b_tau_m_KKpi1>2.)||(b_tau_m_KKpi1<1.9))&&((b_tau_m_KKpi2>2.)||(b_tau_m_KKpi2<1.9)) &&(b_tau_min_dr_mu>0.5)&&(b_tau_min_dr_e>0.5)&&(b_tau_alpha>0) (b_tau_m<1.7)&&(b_tau_alpha>0)&&(D0_fsig>3.)&&(b_B_m>3.)&&(b_tau_vtx4trkProb<0.2)
			cut.emplace(std::make_pair("SR", cut["base"]+TCut(TString::Format("mvaScore>=%f", mvaCutSR)))); 
			cut.emplace(std::make_pair("SB", cut["base"]+TCut(TString::Format("(mvaScore >= %f) && (mvaScore < %f)", mvaCutSB, mvaCutSR)))); 
			cut.emplace(std::make_pair("CR", cut["base"]+TCut(TString::Format("(mvaScore >= %f) && (mvaScore < %f)", mvaCutCR, mvaCutSB)))); 
			// (b_B_fsig>2.)&&(b_tau_fsig>3.)&&b_tau_vprob>0.1&&(b_tau_alpha>0.)&&(b_B_r<1.)&&(b_B_mu_alpha>1.)&&(b_tau_rhomass_max>0.5)&&(b_tau_m<1.7)&&(b_B_m>3.)&&(b_B_m<5.28)
		}
		else 
		{
			cut.emplace(std::make_pair("base", "b_B_m>0.&&b_B_nmu<1&&b_B_nh<3&&b_B_ne<1&&b_B_fsig>2")); //&&b_B_npi0<15 &&b_B_mm2<1.&&pttau_B_fl>0.2 &&b_B_mm2<1.&&pttau_B_fsig>3
			cut.emplace(std::make_pair("mw", "b_B_m>5.15&&b_B_m<5.4&&b_B_fsig>2.&&b_B_nmu<1&&b_B_ne<1&&b_B_nh<1")); //&&b_B_npi0<3&&pttau_B_ngamma<5
			// superclean: pttau_B_m>0.&&pttau_B_nmu<1&&pttau_B_nh<1&&pttau_B_ne<1&&(pttau_B_npi0+pttau_B_ngamma)<2&&pttau_B_mm2<1.&&pttau_B_fl>0.2 (or with fsig > 7.)
		}

		cutstandalone.emplace(std::make_pair("SR", TCut("base", "(b_Ds_vprob>0.1) && (b_D0_vprob>0.1) && (mvaScore>-2.)")+TCut(TString::Format("mvaScore>=%f", mvaCutSR)))); 
		cutstandalone.emplace(std::make_pair("SB", TCut("base", "(b_Ds_vprob>0.1) && (b_D0_vprob>0.1) && (mvaScore>-2.)")+TCut(TString::Format("(mvaScore >= %f) && (mvaScore < %f)", mvaCutSB, mvaCutSR)))); 
		cutstandalone.emplace(std::make_pair("CR", TCut("base", "(b_Ds_vprob>0.1) && (b_D0_vprob>0.1) && (mvaScore>-2.)")+TCut(TString::Format("(mvaScore >= %f) && (mvaScore < %f)", mvaCutCR, mvaCutSB)))); 

		// https://colorbrewer2.org/?type=diverging&scheme=RdYlBu&n=7 
		//colors = {new TColor(TColor::GetFreeColorIndex(), 165,0,38), new TColor(TColor::GetFreeColorIndex(), 215,48,39), new TColor(TColor::GetFreeColorIndex(), 244,109,67), new TColor(TColor::GetFreeColorIndex(), 253,174,97), new TColor(TColor::GetFreeColorIndex(), 254,224,144), new TColor(TColor::GetFreeColorIndex(), 255,255,191), new TColor(TColor::GetFreeColorIndex(), 224,243,248), new TColor(TColor::GetFreeColorIndex(), 171,217,233), new TColor(TColor::GetFreeColorIndex(), 116,173,209), new TColor(TColor::GetFreeColorIndex(), 69,117,180), new TColor(TColor::GetFreeColorIndex(), 49,54,149)}; 
			// {new TColor(TColor::GetFreeColorIndex(), 215,48,39), new TColor(TColor::GetFreeColorIndex(), 252,141,89), new TColor(TColor::GetFreeColorIndex(), 254,224,144), new TColor(TColor::GetFreeColorIndex(), 255,255,191), new TColor(TColor::GetFreeColorIndex(), 224,243,248), new TColor(TColor::GetFreeColorIndex(), 145,191,219), new TColor(TColor::GetFreeColorIndex(), 69,117,180)}; 

		//mycolors.reserve(colors.size()); 
		//for (auto color : colors) 
		//{
		//	mycolors.push_back(color->GetNumber()); 
		//}

		//color = {{"Sig", mycolors[0]}, {"B0toDstarDs", mycolors[1]}, {"BkgDstarDs", mycolors[1]}, {"B0toDstarDsstar", mycolors[2]}, {"BkgDstarDsstar", mycolors[2]}, {"B0toDstarD", mycolors[3]}, {"B0toDstarD0K", mycolors[5]}, {"ButoDstarDK", mycolors[4]}, {"B0toDstar3pi", mycolors[7]}, {"BkgDstara1", mycolors[6]},{"WS", mycolors[10]}, {"WSTau", mycolors[10]}, {"dataD2WS", mycolors[10]}, {"dataD2TauWS", mycolors[10]}, }; 
			// {{"Sig", mycolors[0]}, {"SigPart", mycolors[1]}, {"B0toDstarDs", mycolors[2]}, {"BkgDstarDs", mycolors[2]}, {"B0toDstarDsstar", mycolors[3]}, {"BkgDstarDsstar", mycolors[3]}, {"BkgDstara1", mycolors[4]},{"WS", mycolors[6]}, {"WSTau", mycolors[5]}, {"dataD2WS", mycolors[6]}, {"dataD2TauWS", mycolors[5]}, {"B0toDstarD0K", mycolors[5]}}; 

		CommonInitialisation(denominator); 

		model = {	
			{"v1", "./data/tautagger/batchsize_10/serialized"}, 
			{"v2", "./data/tautagger/FirstTopUp/serialized"}, 
			{"v3", "./data/tautagger/FirstTopUp/serialized"},  
			{"v3.5", "./data/tautagger/FirstTopUp/serialized"},  
			{"v4", "./data/tautagger/NewSelTopUpNoOverlap/serialized"}, 
			{"v6.7", "./data/tautagger/FlightSigCorrNoCharge/serialized"}, 
			{"v6.8", "./data/tautagger/FlightSigCorrNoCharge/serialized"},
			{"v6.9", "./data/tautagger/FlightSigCorrNoCharge/serialized"},
			{"v6.945", "./data/tautagger/FlightSigCorrNoCharge/serialized"},
			{"v6.95", "./data/tautagger/FlightSigCorrNoCharge/serialized"},
			{"v7", "./data/tautagger/trainingv18/serialized"}, 
			{"v7.1", "./data/tautagger/trainingv18/serialized"},
			{"v7.2", "./data/tautagger/trainingv18/serialized"}
		}; 

		MVA = {	
			{"v1", "./anaMVA/secondtraining/model_optimized/weights.xml"}, 
			{"v2", "./anaMVA/secondtraining/model_optimized/weights.xml"}, 
			{"v3", "./anaMVA/secondtraining/model_optimized/weights.xml"},  
			{"v3.5", "./anaMVA/NewSelection/model_optimized/weights.xml"},  
			{"v4", "./anaMVA/NewSelection/model_optimized/weights.xml"}, 
			{"v6.7", "./anaMVA/NewFixTauFL/model_optimized/weights.xml"}, 
			{"v6.8", "./anaMVA/NewIsoWithCutsvsdata/model_optimized/weights.xml"}, //./anaMVA/NewIsoVariablesAgainstWS/model_optimized/weights.xml" ./anaMVA/NewIsoVariables/model_optimized/weights.xml
			{"v6.9", "./anaMVA/NewFixTauFL/model_optimized/weights.xml"}, //LatestVsData
			{"v6.95", "./anaMVA/TrimmedVariablesWS/model_optimized/weights.xml"}, //LatestVsData
			{"v7", "./anaMVA/NewFixTauFL/model_optimized/weights.xml"}, //trainingv7onD1
			{"v7.1", "./anaMVA/NewFixTauFL/model_optimized/weights.xml"}, //trainingv7onD1
			{"v7.2", "./anaMVA/NewFixTauFL/model_optimized/weights.xml"} //trainingv7onD1
		}; 

		//legends = {{"Sig", "signal"}, {"SigPart", "part. signal"}, {"B0toDstarDs", "B^{0}#rightarrowD*D_{s}"}, {"BkgDstarDs", "B^{0}#rightarrowD*D_{s}"}, {"B0toDstarDsstar", "B^{0}#rightarrowD*D*_{s}"}, {"BkgDstarDsstar", "B^{0}#rightarrowD*D*_{s}"}, {"BkgDstara1", "B^{0}#rightarrowD*a_{1}"}, {"B0toDstar3pi", "B^{0}#rightarrowD^{*}3pi"}, {"B0toDstarD", "B^{0}#rightarrowD^{*}D"}, {"ButoDstarDK", "B^{+}#rightarrowD^{*}DK"}, {"B0toDstarD0K", "B^{0}#rightarrowD^{*}D^{0}K"},{"WS", "|q_{B}| = 2  WS"}, {"WSTau", "|q_{#tau}| = 3  WS"}, {"dataD2WS", "|q_{B}|=2 WS"}, {"dataD2TauWS", "|q_{#tau}|=3 WS"}}; 

		//labels = {{"b_tau_rhomass1", "Invariant m_{#rho}"}, {"b_tau_rhomass2", "Invariant m_{#rho}"}, {"b_B_q2", "q2"}, {"b_B_m", "Reconstructed m_{B}"}}; 

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

