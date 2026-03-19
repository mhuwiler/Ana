#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
//#include "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/Ana/plugins/FileManager/CFileManager.C"
ClassImp(FileManager)
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include "FileFlow.h"
#include "DrawTMVAHistogram.C"
#include "GetSeparation.C"
//#include "Python.h"
//#include "TPython.h"
#include "Tau.h"


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


ROOT::VecOps::RVec<float> FillTauDNNscore(ROOT::VecOps::RVec<int> trackIndices, std::vector<float> trackScores) 
{
	std::vector<float> responses; 
	for (auto index : trackIndices) 
	{
		float dnnScore = -1; 
		if (index < trackScores.size()) 
		{
			dnnScore = trackScores.at(index); 
		}
		else 
		{
			std::cout << "WARNING: index out of range! " << std::endl; 
		}
		responses.push_back(dnnScore); 
	}
	assert(responses.size() == trackIndices.size()); 
	return responses; 
}


std::vector<float> FillSumDNN(ROOT::VecOps::RVec<float> pion1_dnn, ROOT::VecOps::RVec<float> pion2_dnn, ROOT::VecOps::RVec<float> pion3_dnn) 
{
	std::vector<float> response; 
	const int dim = pion1_dnn.size(); 
	assert(pion2_dnn.size() == dim); 
	assert(pion3_dnn.size() == dim); 
	for (unsigned int i=0; i<dim; i++) 
	{
		response.push_back(pion1_dnn.at(i)+pion2_dnn.at(i)+pion3_dnn.at(i)); 
	}
	return response; 
}


struct Basictau 
{
	float pt; 
	float eta; 
	float phi; 
	int q; 
	float m; 
	int idx1;
	int idx2; 
	int idx3; 
	float sumdnn; 
	float dnn1; 
	float dnn2; 
	float dnn3; 

}; 


std::vector<Tau> BuildTauCandidates(//ROOT::VecOps::RVec<float> taupt, ROOT::VecOps::RVec<float> taueta, ROOT::VecOps::RVec<float> tauphi, ROOT::VecOps::RVec<int> taucharge, ROOT::VecOps::RVec<float> taumass, ROOT::VecOps::RVec<float> tauVprob, ROOT::VecOps::RVec<float> taufsig, ROOT::VecOps::RVec<float> taulip, ROOT::VecOps::RVec<int> idx1, ROOT::VecOps::RVec<int> idx2, ROOT::VecOps::RVec<int> idx3, ROOT::VecOps::RVec<float> dnn1, ROOT::VecOps::RVec<float> dnn2, ROOT::VecOps::RVec<float> dnn3, std::vector<float> sumdnn, 
										//ROOT::VecOps::RVec<float> alpha, ROOT::VecOps::RVec<float> maxDr, ROOT::VecOps::RVec<float> taufl, ROOT::VecOps::RVec<float> pvip, ROOT::VecOps::RVec<float> pvips, ROOT::VecOps::RVec<float> dau1pt, ROOT::VecOps::RVec<float> dau1eta, ROOT::VecOps::RVec<float> dau1phi, ROOT::VecOps::RVec<float> dau2pt, ROOT::VecOps::RVec<float> dau2eta, ROOT::VecOps::RVec<float> dau2phi, ROOT::VecOps::RVec<float> dau3pt, ROOT::VecOps::RVec<float> dau3eta, ROOT::VecOps::RVec<float> dau3phi, ROOT::VecOps::RVec<float> rhomass1, ROOT::VecOps::RVec<float> rhomass2
		#include "tauvariableargs.gcf"
) 
{
	std::vector<Tau> mytaus; 
	for (unsigned int i=0; i<init_pt.size(); i++) 
	{
		Tau tau; 
		//Tau tau(taupt.at(i), taueta.at(i), tauphi.at(i), taucharge.at(i), taumass.at(i)); 
		////tau.pt = taupt.at(i); 

		//tau.SetKinematics(tauVprob.at(i), taufsig.at(i), taulip.at(i)); 
		//tau.SetIndices(idx1.at(i), idx2.at(i), idx3.at(i)); 
		//tau.SetDNN(dnn1.at(i), dnn2.at(i), dnn3.at(i), sumdnn.at(i)); 

		//tau.SetEventKinematics(alpha.at(i), maxDr.at(i), taufl.at(i), pvip.at(i), pvips.at(i)); 
		//tau.SetDau1Kin(dau1pt.at(i), dau1eta.at(i), dau1phi.at(i)); 
		//tau.SetDau2Kin(dau2pt.at(i), dau2eta.at(i), dau2phi.at(i)); 
		//tau.SetDau3Kin(dau3pt.at(i), dau3eta.at(i), dau3phi.at(i)); 
		//tau.SetRhoMasses(rhomass1.at(i), rhomass2.at(i)); 

		#include "tauvariableaffectation.gcf"

		mytaus.push_back(tau); 
	}
	std::cout << "Built tau candidates" << std::endl; 
	assert(mytaus.size() == init_pt.size()); 
	return mytaus; 
}

struct 
{
    bool operator()(const Tau& B1, const Tau& B2) const 
    { 
    	return  measure(B1) > measure(B2); 
    }

    float measure(const Tau& B) const 
    {
    	return B.pt; 
    }

} SortBCandidates; 



void WriteEffInfo(std::string directory, const Double_t Initial, Double_t final = 0., const TString structure = "EffUpdateTau") 
{
	TFile *outfile = TFile::Open(filemanager.GetFile(directory).c_str(), "UPDATE"); 

	TDirectory *dir = static_cast<TDirectory*>(outfile->Get("ntuplizer")); // TODO: tokenize this from the object

	if (final == 0) 
	{
		RDataFrame frame(*static_cast<TTree*>(outfile->Get(filemanager.GetObject(directory).c_str()))); 
		final = frame.Count().GetValue(); 
		std::cout << final << std::endl; 
	}

	dir->cd(); 
	//TDirectory *newdir = directory->mkdir("SelectTau"); 
	dir->cd(); 
	TVectorD vec(2);
	vec[0] = Initial; 
	vec[1] = final;  
	vec.Write(structure); 


	outfile->Write(); 

	outfile->Close(); 
}

void UpdateEffInfo(const TString& outIdentifier, const TString& effTreeName = "ntuplizer/EffCalc", const TString& effInfoObject = "ntuplizer/EffInfo") 
{
	TString filename = filemanager.GetFile(outIdentifier);
 
	//if (!selEff) std::cerr << "ERROR: No selection efficiency information found in file. Please make sure it was added or copied to the file under: " << effTreeName << std::endl;

	//TDirectory *dir = static_cast<TDirectory*>(file->Get("ntuplizer"));

	auto selEffFrame = RDataFrame(*filemanager.GetItem<TTree*>("effInfo", true));

	Double_t n = selEffFrame.Sum("numSelected").GetValue();

	Double_t N = selEffFrame.Sum("numTotal").GetValue();

	RDataFrame frame(*filemanager.GetItem<TTree*>(outIdentifier, true)); 

	Double_t nsel = frame.Count().GetValue();

	TFile *file = TFile::Open(filename, "UPDATE");

	auto tokens = effInfoObject.Tokenize("/"); 
	TString effInfoName = static_cast<TObjString*>(tokens->At(tokens->GetEntries()-1))->GetString(); 
	TString effInfoFolder = effInfoObject; 
	effInfoFolder.ReplaceAll("/"+effInfoName, "");

	TDirectory *dir = static_cast<TDirectory*>(file->Get(effInfoFolder));

	//std::cout << "Yields: " << N << " " << n << " " << nsel << std::endl;

	dir->cd(); 
	TVectorD vec(3);
	vec[0] = N; 
	vec[2] = nsel; // Always the first two elements give the efficiency
	vec[1] = n;
	vec.Write(effInfoName); 

	file->Close(); 
}

Tau SelectTauCandidate(std::vector<Tau> collection, const int q) 
{
	std::sort(collection.begin(), collection.end(), SortBCandidates); //std::greater<>()

	for (auto candidate : collection) 
	{
		if (true) return candidate; // candidate.m < 1.7 candidate.q == q
	}

	return Tau(); // This wil return an invalid tau
}

Tau SelectGenmatchedTauCandidate(std::vector<Tau> collection) // For taking explicitly candidates in the list that are genmatched with a certain number of pions 
{
	// Sort the taus by pT 
	std::sort(collection.begin(), collection.end(), SortBCandidates); //std::greater<>()
	int matchThreshold = 3; 

	for (auto cand : collection) 
	{
		int numMatched = static_cast<int>(cand.match1) + static_cast<int>(cand.match2) + static_cast<int>(cand.match3); 
		if (numMatched == matchThreshold) 
		{
			std::cout << "Num matched: " << numMatched << std::endl; 
			return cand; // Returning the first in the list which matches the required number of matched pions
		}
	}

	return Tau(); 
}

double deltaR(double eta1, double phi1, double eta2, double phi2) 
{
	double deta = eta1 - eta2;
    double dphi = std::abs(phi1 - phi2);
    if (dphi > double(M_PI))
      dphi -= double(2 * M_PI);
    double dr2 = deta * deta + dphi * dphi;
	return std::sqrt(dr2); 
}

std::vector<double> ComputeAlpha(ROOT::VecOps::RVec<float> muEta, ROOT::VecOps::RVec<float> muPhi, ROOT::VecOps::RVec<float> candidateEta, ROOT::VecOps::RVec<float> candidatePhi) // TODO: change to double 
{
	std::vector<double> result; 
	result.reserve(candidateEta.size()); 
	assert(muEta.size() > 0); 
	assert(muEta.size() == muPhi.size()); 
	assert(candidateEta.size() == candidatePhi.size()); 
	double refEta = muEta.at(0); 
	double refPhi = muPhi.at(0); 

	for (unsigned int i=0; i<candidateEta.size(); i++) 
	{
		result.push_back(deltaR(refEta, refPhi, candidateEta.at(i), candidateEta.at(i))); 
	}

	return result; 
}


float findMin(const float& pt1, const float& pt2, const float& pt3)
{
	return min(min(pt1, pt2), pt3); 
}


float findMax(const float& pt1, const float& pt2, const float& pt3)
{
	return max(max(pt1, pt2), pt3); 
}

float getMin(const float& pt1, const float& pt2)
{
	return min(pt1, pt2); 
}


float getMax(const float& pt1, const float& pt2)
{
	return max(pt1, pt2); 
}



float extractFirstElement(const ROOT::VecOps::RVec<float>& vec) 
{
	assert(vec.size() > 0); 
	return vec.at(0); 
}


float ratioOfRadius(const float r1, const float r2) 
{
	return r1/r2; 
}

std::map<unsigned int, std::map<unsigned int, std::vector<unsigned int> > > LoadEventList(const TString& path) 
{
		std::ifstream file(path, std::ifstream::binary);
		Json::Value events;
		file >> events;

		std::map<unsigned int, std::map<unsigned int, std::vector<unsigned int> > > loaded; 

		for (auto item : events.getMemberNames()) 
		{
			//std::cout << item << std::endl; 
			unsigned int run = std::stol(item);

			//std::cout << run << std::endl;

			std::map<unsigned int, std::vector<unsigned int> > lumis; 

			for (auto block : events[item].getMemberNames()) 
			{
				unsigned int lumiblock = std::stol(block); 
				auto& eventlist = events[item][block]; 
				//assert(eventlist.size() == 2); 

				/*for (auto element : eventlist) 
				{
					std::cout << element << std::endl; 
				}*/

				//std::cout << eventlist.size() << std::endl; 

				std::vector<unsigned int> vec; 
				vec.reserve(eventlist.size()); 
				std::transform(eventlist.begin(), eventlist.end(), std::back_inserter(vec), [](const auto& e) { return e.asDouble(); });
				lumis.emplace(std::piecewise_construct, std::make_tuple(lumiblock), std::make_tuple(vec));
			}

			loaded.insert(std::make_pair(run, std::move(lumis))); 
		}

		/*for (auto item : events)
		{
			std::cout << item << " " << std::endl; 
		}*/


		file.close();
		return std::move(loaded); 
}

 

void ProcessingDenominatorOffline(const TString& identifier, const TString& cycle, const bool ws = false) 
{
	ROOT::DisableImplicitMT(); 
	//ROOT::EnableImplicitMT(); //ROOT::DisableImplicitMT(); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C+");
	//gSystem->Load("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.so"); 
	Init(cycle); 

	std::cout << "Starting processing" << std::endl; 


	TString inIdentifier = identifier+"_ntuple";
	TString outIndentifier = identifier+"_DNN";


	filemanager.OpenItem(inIdentifier); 

	gStyle->SetOptStat(0); 


	auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>(inIdentifier)); // tree100k

	filemanager.AddItem("effInfo", filemanager.GetFile(inIdentifier), "ntuplizer/EffCalc"); 

	Double_t Ninitial = (*filemanager.GetItem<TTree*>(inIdentifier)).GetEntries(); 
	

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

	//auto frame2 = dataframe.Define("P_D0", P_v, {"BsDstarTauNu_D0_pt", "BsDstarTauNu_D0_eta", "BsDstarTauNu_D0_phi", "BsDstarTauNu_D0_mass"}) // auto frame2 = dataframe.Define("LV_D0", "TLorentzVector LV_D0; LV_D0.SetPtEtaPhiM(BsDstarTauNu_D0_pt, BsDstarTauNu_D0_eta, BsDstarTauNu_D0_phi, BsDstarTauNu_D0_mass); return LV_D0"); 
	//				.Define("P_Ds", P_v, {"BsDstarTauNu_Ds_pt", "BsDstarTauNu_Ds_eta", "BsDstarTauNu_Ds_phi", "BsDstarTauNu_Ds_mass"})
	//				.Define("P_tau", P_v, {"BsDstarTauNu_tau_pt", "BsDstarTauNu_tau_eta", "BsDstarTauNu_tau_phi", "BsDstarTauNu_tau_mass"})
	//				.Define("B_mass", invMass_v, {"P_Ds", "P_tau"}); 


	//auto histo1 = frame2.Histo1D("B_mass"); 

	//auto histo2 = frame2.Histo2D({"Bmass_vs_Dmass", "Correlation plot between B and D masses", 100, 0., 7000., 100, 0., 5000.}, "BsDstarTauNu_B_mass", "BsDstarTauNu_D0_unfit_mass"); 

	auto TFresponse = [](std::vector<float> Dstarpt, std::vector<float> Dstareta, std::vector<float> Dstarphi, std::vector<float> Dstarcharge, std::vector<float> pt, std::vector<float> eta, std::vector<float> phi, std::vector<float> q, std::vector<float> DOCA2D, std::vector<float> DOCA2DErr, std::vector<float> DOCA3D, std::vector<float> DOCA3DErr, std::vector<float> dzToPV, std::vector<float> dzToClosest, std::vector<float> isAssociate, std::vector<float> assocQualityToPV, std::vector<int> genmatch) 
	{
		std::vector<float> response; 

		return response; 
	};

	auto veto = LoadEventList(Ana::folder+"/numveto.json"); 

	auto numeratorVeto = [&veto](unsigned int run, unsigned int block, unsigned long long event) 
	{
		bool isFound = false; 
		const auto& evt = veto.find(run);
		if (!(evt == veto.end())) 
		{
			auto lumiblock = (*evt).second.find(block); 
			if (!(lumiblock == (*evt).second.end())) 
			{
				auto& vec = (*lumiblock).second; 
				auto foundevt = std::find(vec.begin(), vec.end(), event); 
				if (foundevt != vec.end()) isFound = true; 
			}
		} 
		if (isFound) std::cout << "Event in numerator!" << std::endl; 
		return isFound; 
	};

	int count = 0; 

	auto SelectBCandidate = [&count](//ROOT::VecOps::RVec<float> taupt, ROOT::VecOps::RVec<float> taueta, ROOT::VecOps::RVec<float> tauphi, ROOT::VecOps::RVec<int> taucharge, ROOT::VecOps::RVec<float> taumass, ROOT::VecOps::RVec<float> tauVprob, ROOT::VecOps::RVec<float> taufsig, ROOT::VecOps::RVec<float> taulip, ROOT::VecOps::RVec<int> idx1, ROOT::VecOps::RVec<int> idx2, ROOT::VecOps::RVec<int> idx3, ROOT::VecOps::RVec<float> dnn1, ROOT::VecOps::RVec<float> dnn2, ROOT::VecOps::RVec<float> dnn3, std::vector<float> sumdnn, 
										//ROOT::VecOps::RVec<float> alpha, ROOT::VecOps::RVec<float> maxDr, ROOT::VecOps::RVec<float> taufl, ROOT::VecOps::RVec<float> pvip, ROOT::VecOps::RVec<float> pvips, ROOT::VecOps::RVec<float> dau1pt, ROOT::VecOps::RVec<float> dau1eta, ROOT::VecOps::RVec<float> dau1phi, ROOT::VecOps::RVec<float> dau2pt, ROOT::VecOps::RVec<float> dau2eta, ROOT::VecOps::RVec<float> dau2phi, ROOT::VecOps::RVec<float> dau3pt, ROOT::VecOps::RVec<float> dau3eta, ROOT::VecOps::RVec<float> dau3phi, ROOT::VecOps::RVec<float> rhomass1, ROOT::VecOps::RVec<float> rhomass2, 
										//ROOT::VecOps::RVec<int> match1, ROOT::VecOps::RVec<int> match2, ROOT::VecOps::RVec<int> match3, ROOT::VecOps::RVec<float> Bm, ROOT::VecOps::RVec<float> Bq2
			#include "tauvariableargs.gcf"

	) // TODO: set to int  
	{
		std::vector<Tau> mytaus; 
		for (unsigned int i=0; i<init_pt.size(); i++) 
		{
			Tau tau;

			#include "tauvariableaffectation.gcf"

			mytaus.push_back(tau); 
		}
		count++; 
		if ((count <= 10)) std::cout << "Built tau candidates for event " << count << std::endl;
		if (count == 10) std::cout << "..." << std::endl;
		if ((count % 100 == 0) && (count < 1000)) std::cout << "Built tau candidates for event " << count << std::endl;
		if (count == 1001) std::cout << "..." << std::endl;
		if (count % 1000 == 0) std::cout << "Built tau candidates for event " << count << std::endl; 
		assert(mytaus.size() == taupt.size()); 
		return mytaus; 
	};

	auto SelectBCandidateRegion = [&ws](std::vector<Tau> collection, const int q) 
	{
		std::sort(collection.begin(), collection.end(), SortBCandidates); //std::greater<>()

		int charge = ws ? q : -q; 

		for (auto candidate : collection) 
		{
			if (candidate.q == charge) return candidate; // candidate.m < 1.7 candidate.q == q
		}

		return Tau(); // This wil return an invalid tau
	};
	// Could use the following type: std::function<Tau (std::vector<Tau>, const int)>



	auto withB = dataframe.Define("v_Bcandidates", SelectBCandidate, {//"v_tau_pt", "v_tau_eta", "v_tau_phi", "v_tau_q", "v_tau_m", "v_tau_vprob", "v_tau_fsig", "v_tau_lip", "v_tau_idx1", "v_tau_idx2", "v_tau_idx3", "v_tau_dnn_1", "v_tau_dnn_2", "v_tau_dnn_3", "v_tau_sumdnn", 
																		//"v_tau_alpha", "v_tau_fl", "v_tau_pvip", "v_tau_pvipsig", "v_tau_legacyMaxdr", "v_tau_pi1pt", "v_tau_pi1eta", "v_tau_pi1phi", "v_tau_pi2pt", "v_tau_pi2eta", "v_tau_pi2phi", "v_tau_pi3pt", "v_tau_pi3eta", "v_tau_pi3phi", "v_tau_rhomass1", "v_tau_rhomass2", "v_tau_match1", "v_tau_match2", "v_tau_match3", "v_B_m", "v_B_q2"});
																		////"v_tau_alpha", "v_tau_fl3d", "v_tau_pvip", "v_tau_pvips", "v_tau_max_dr_3prong", "v_tau_pi1_pt", "v_tau_pi1_eta", "v_tau_pi1_phi", "v_tau_pi2_pt", "v_tau_pi2_eta", "v_tau_pi2_phi", "v_tau_pi3_pt", "v_tau_pi3_eta", "v_tau_pi3_phi", "v_tau_rhomass1", "v_tau_rhomass2"
		#include "taubranchreading.gcf"

	});
	
	withB = withB.Define("b_tau", SelectBCandidateRegion, {"v_Bcandidates", "Dstar_q"})
				// .Define("b_tau_pt", Tau::WritePt, {"b_tau"})
				// .Define("b_tau_eta", Tau::WriteEta, {"b_tau"})
				// .Define("b_tau_phi", Tau::WritePhi, {"b_tau"})
				// .Define("b_tau_q", Tau::WriteCharge, {"b_tau"})
				// .Define("b_tau_m", Tau::WriteMass, {"b_tau"})
				// .Define("b_tau_vprob", Tau::WriteVprob, {"b_tau"})
				// .Define("b_tau_fsig", Tau::WriteFsig, {"b_tau"})
				// .Define("b_tau_lip", Tau::WriteLip, {"b_tau"}) // Add lips
				// .Define("b_tau_dnn1", Tau::WriteDNN1, {"b_tau"})
				// .Define("b_tau_dnn2", Tau::WriteDNN2, {"b_tau"})
				// .Define("b_tau_dnn3", Tau::WriteDNN3, {"b_tau"})
				// .Define("b_tau_sumdnn", Tau::WriteSumDNN, {"b_tau"})
				// .Define("b_tau_idx1", Tau::WriteIdx1, {"b_tau"})
				// .Define("b_tau_idx2", Tau::WriteIdx2, {"b_tau"})
				// .Define("b_tau_idx3", Tau::WriteIdx3, {"b_tau"})
				// .Define("b_tau_fl", Tau::WriteFl, {"b_tau"})
				// .Define("b_tau_alpha", Tau::WriteAlpha, {"b_tau"})
				// .Define("b_tau_pvip", Tau::WritePVIP, {"b_tau"})
				// .Define("b_tau_pvips", Tau::WritePVIPsig, {"b_tau"})
				// .Define("b_tau_maxdr", Tau::WriteDr, {"b_tau"})
				// .Define("b_tau_pi1pt", Tau::WriteDau1Pt, {"b_tau"})
				// .Define("b_tau_pi1eta", Tau::WriteDau1Eta, {"b_tau"})
				// .Define("b_tau_pi1phi", Tau::WriteDau1Phi, {"b_tau"})
				// .Define("b_tau_pi2pt", Tau::WriteDau2Pt, {"b_tau"})
				// .Define("b_tau_pi2eta", Tau::WriteDau2Eta, {"b_tau"})
				// .Define("b_tau_pi2phi", Tau::WriteDau2Phi, {"b_tau"})
				// .Define("b_tau_pi3pt", Tau::WriteDau3Pt, {"b_tau"})
				// .Define("b_tau_pi3eta", Tau::WriteDau3Eta, {"b_tau"})
				// .Define("b_tau_pi3phi", Tau::WriteDau3Phi, {"b_tau"})
				// .Define("b_tau_rhomass1", Tau::WriteRhomass1, {"b_tau"})
				// .Define("b_tau_rhomass2", Tau::WriteRhomass2, {"b_tau"})
				// .Define("b_tau_match", Tau::WriteMatch, {"b_tau"})
				// .Define("b_tau_sumMatch", Tau::WriteSumMatch, {"b_tau"})
				// .Define("b_B_m", Tau::WriteBmass, {"b_tau"})
				// .Define("b_B_q2", Tau::WriteBq2, {"b_tau"});
				#include "taubranchcreation.gcf" 
				.Define("b_B_rrh", ratioOfRadius, {"b_B_min_dr_h", "b_B_r"})
				.Define("b_B_rre", ratioOfRadius, {"b_B_min_dr_e", "b_B_r"})
				.Define("b_B_rrmu", ratioOfRadius, {"b_B_min_dr_mu", "b_B_r"})
				.Define("b_B_rrh_first", ratioOfRadius, {"b_B_min_dr_h", "b_B_r_first"})
				.Define("b_B_rre_first", ratioOfRadius, {"b_B_min_dr_e", "b_B_r_first"})
				.Define("b_B_rrmu_first", ratioOfRadius, {"b_B_min_dr_mu", "b_B_r_first"})
				.Define("b_B_rrh_second", ratioOfRadius, {"b_B_min_dr_h", "b_B_r_second"})
				.Define("b_B_rre_second", ratioOfRadius, {"b_B_min_dr_e", "b_B_r_second"})
				.Define("b_B_rrmu_second", ratioOfRadius, {"b_B_min_dr_mu", "b_B_r_second"});
	;

	/*withB = withB.Define("b_D0_pt", extractFirstElement, {"BsDstarTauNu_D0_pt"})
				.Define("b_D0_eta", extractFirstElement, {"BsDstarTauNu_D0_eta"})
				.Define("b_D0_phi", extractFirstElement, {"BsDstarTauNu_D0_phi"})
				.Define("b_D0_vprob", extractFirstElement, {"BsDstarTauNu_D0_vprob"})
				.Define("b_D0_fl", extractFirstElement, {"BsDstarTauNu_D0_fl3d"})
				.Define("b_D0_fsig", extractFirstElement, {"BsDstarTauNu_D0_fls3d"})
				.Define("b_D0_lip", extractFirstElement, {"BsDstarTauNu_D0_lip"})
				.Define("b_D0_lips", extractFirstElement, {"BsDstarTauNu_D0_lips"})
				.Define("b_D0_pvip", extractFirstElement, {"BsDstarTauNu_D0_pvip"})
				.Define("b_Ds_pt", extractFirstElement, {"BsDstarTauNu_Ds_pt"})
				.Define("b_Ds_eta", extractFirstElement, {"BsDstarTauNu_Ds_eta"})
				.Define("b_Ds_phi", extractFirstElement, {"BsDstarTauNu_Ds_phi"})
				.Define("b_Ds_vprob", extractFirstElement, {"BsDstarTauNu_Ds_vprob"})
				.Define("b_Ds_fl", extractFirstElement, {"BsDstarTauNu_Ds_fl3d"})
				.Define("b_Ds_fsig", extractFirstElement, {"BsDstarTauNu_Ds_fls3d"})
				.Define("b_Ds_lip", extractFirstElement, {"BsDstarTauNu_Ds_lip"})
				.Define("b_Ds_lips", extractFirstElement, {"BsDstarTauNu_Ds_lips"})
				.Define("b_Ds_pvip", extractFirstElement, {"BsDstarTauNu_Ds_pvip"}); 
	*/

	auto filtered = withB.Filter("(b_B_m > 0.)"); //&&(b_tau_min_dr_mu>0.5)&&(b_tau_min_dr_e>0.5)

	filtered = filtered.Define("isInNumerator", numeratorVeto, {"EVENT_run", "EVENT_lumiBlock", "EVENT_event"}).Filter("!isInNumerator"); 

	auto pimped = filtered.Define("b_tau_minpipt", findMin, {"b_tau_pi1pt", "b_tau_pi2pt", "b_tau_pi3pt"}).Define("b_tau_maxpipt", findMax, {"b_tau_pi1pt", "b_tau_pi2pt", "b_tau_pi3pt"})
						.Define("b_tau_minpieta", findMin, {"b_tau_pi1eta", "b_tau_pi2eta", "b_tau_pi3eta"}).Define("b_tau_maxpieta", findMax, {"b_tau_pi1eta", "b_tau_pi2eta", "b_tau_pi3eta"})
						.Define("b_tau_minpiphi", findMin, {"b_tau_pi1phi", "b_tau_pi2phi", "b_tau_pi3phi"}).Define("b_tau_maxpiphi", findMax, {"b_tau_pi1phi", "b_tau_pi2phi", "b_tau_pi3phi"})
						.Define("b_tau_rhomass_min", getMin, {"b_tau_rhomass1", "b_tau_rhomass2"}).Define("b_tau_rhomass_max", getMax, {"b_tau_rhomass1", "b_tau_rhomass2"}); 

	//pimped = pimped.Define("b_tau_rho_m_unrolled", computeUnrolledMass, {"b_tau_rhomass1", "b_tau_rhomass2"}); 
	pimped = pimped.Define("b_tau_rho_m_unrolled", "int((min(b_tau_rhomass2, float(1.3)) - 0.2)/0.22) + 6*int((min(b_tau_rhomass1, float(1.3)) - 0.2)/0.22)");

	#include "stringbranches.gcf"

	for (auto branch : stringbranches) // Hack to fix string branche 
	{
		pimped = pimped.Redefine(branch, [](const ROOT::RVec<std::string> &v) {return std::vector<std::string>(v.begin(), v.end());}, {branch}); 
	}

	const std::vector<std::string> blacklist = {"v_Bcandidates", "b_tau"};

	pimped.Snapshot(filemanager.GetObject(outIndentifier), filemanager.GetFile(outIndentifier), purgeColumns(pimped.GetColumnNames(), blacklist)); 

	std::cout << "Processed " << count << " events. " << std::endl;

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


	// Writing out eff info 
	//filemanager.OpenItem(outIndentifier, "READ"); 

	//auto effFrame = RDataFrame(*filemanager.GetItem<TTree*>(outIndentifier)); 

	//Double_t Nfinal = effFrame.Count().GetValue(); 


	//std::cout << "N initial: " << Ninitial << ", N final: " << Nfinal << std::endl; 

	WriteEffInfo(outIndentifier.Data(), Ninitial); // Todo: remove

	UpdateEffInfo(outIndentifier); 
	

	//filemanager.CloseAll(); 


}

