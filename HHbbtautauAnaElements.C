#ifndef HHbbtautauAnaElements_hxx
#define HHbbtautauAnaElements_hxx
#include <string>
#include <vector>
#include "TChain.h"
#include "ROOT/RDataFrame.hxx"
#include "TLorentzVector.h"
#include "Math/Vector4D.h"
#include <thread>
#include "Particle.h"
#include "FileFlow.h"


constexpr double Pion_Mass = 0.13957; // The pion mass from the PDG (used as default mass hypothesis)
constexpr double Muon_Mass = -1.; // TODO: set value from PDG, perhaps also for electron and Kaon


constexpr double drThres = 0.05; // The maximal dR separation when genmatching 

// TODO: 
// - implement autoblacklist
// - do not hardcode names, write all initial names and final names for collections
// - test templated P4 TLorentzVector functions
// - separate into NanoBuildingBlocks.C and HHbbtautauAnaElements.C


//typedef ROOT::Math::PtEtaPhiM4D<double> R4Vec; 
template <typename T>
using R4Vec = ROOT::Math::PtEtaPhiM4D<T>; 


namespace Ana 
{
	extern std::unordered_map<std::string, int> autoblacklist; 

	template <typename T>
	T defaultValue();	

	// ---- Specializations ----	

	// int
	template <>
	inline int defaultValue<int>() {
	    return -999;
	}	

	// float
	template <>
	inline float defaultValue<float>() {
	    return -999.f;
	}	

	// double
	template <>
	inline double defaultValue<double>() {
	    return -999.0;
	}	

	// std::string
	template <>
	inline std::string defaultValue<std::string>() {
	    return "";
	}


	std::unordered_map<std::string, int> PDGid = { 
		{ "Muon", 13 }, 
		{ "Electron", 11 }, 
		{ "Tau", 15 }, 
		{ "Higgs", 25 }, 
		{ "b", 5 }, 
		{"Pi", 211}, 
	}; 

	std::string RevertPDGid(int id) 
	{
		std::string response = ""; 
		for (auto it = PDGid.begin(); it != PDGid.end(); it++) 
		{
			if (it->second == abs(id)) 
			{
				response = it->first; 
				break; 
			}
		}
		return response; 
	}

	enum GenDecay 
	{
		TauhTauh = 1, 
		TauhTaumu = 2, 
		TauhTaue = 3, 
		None = 0
	}; 

	class GenMatchingResult 
	{
		public: 
			int Htotau = -999; 
			int Htob = -999; 
			int b1 = -999; 
			int b2 = -999; 
			int tau1 = -999; 
			int tau2 = -999; 
			int mu = -999; 
			int e = -999; 
			int decayType = None; 
			int VBFjet1 = -999.; 
			int VBFjet2 = -999.; 
			int VBFgenJet1 = -999.; 
			int VBFgenJet2 = -999.; 
	};

	template<typename T>
	int VecSize(ROOT::VecOps::RVec<T> vec)
	{
		return vec.size(); 
	}


	Particle computeP4(const double pt, const double eta, const double phi, const double m = Pion_Mass)   
	{
		Particle P(pt, eta, phi, m, 15); 
		return P; 
	}


	ROOT::VecOps::RVec<Particle> computeP4Vec(const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& m, const ROOT::VecOps::RVec<int>& pdgid = {-999})   
	{
		// Return a vector of P4 for collections with vector branches
		ROOT::VecOps::RVec<Particle> P; 
		int n = pt.size(); 
		P.reserve(n); 

		assert(eta.size() == n); 
		assert(phi.size() == n); 
		assert(m.size() == n); 

		// Loop over the elements
		for (unsigned int i=0; i<n; i++) 
		{
			P.emplace_back(pt[i], eta[i], phi[i], m[i], pdgid[i]); 
		}

		return P; 
	}


	//template<typename T>
	ROOT::RDF::RNode GetP4(ROOT::RDF::RNode& frame, TString prefix, TString name = "") 
	{
		//TLorentzVector P4; 
		std::string pfx = prefix.Data(); 
		//ROOT::RDF::RNode *extended 
		frame = frame.Define(pfx+"_P4", computeP4Vec, {pfx+"_pt", pfx+"_eta", pfx+"_phi", pfx+"_mass", pfx+"_pdgId"});
		//*frame = frame->Define(pfx+"_P4", computeP4Vec<T>, {pfx+"_pt", pfx+"_eta", pfx+"_phi", pfx+"_mass"});
		return frame; 
	}


	/*std::vector<std::string>& purgeColumns(std::vector<std::string> &&columns, const std::vector<std::string>& blacklist)
	{
   			// a lambda that checks if `s` is in the blacklist
   			auto is_blacklisted = [&blacklist](const std::string &s)  { return std::find(blacklist.begin(), blacklist.end(), s) != blacklist.end(); };

   			// removing elements from std::vectors is not pretty, see https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
   			columns.erase(std::remove_if(columns.begin(), columns.end(), is_blacklisted), columns.end());

   			return columns; 
	}*/

	
	ROOT::RDF::RNode GetGenParticleCollection(ROOT::RDF::RNode& frame, const std::string genprefix = "GenPart") 
	{
		auto columnNames = frame.GetColumnNames(); 
		if (std::find(columnNames.begin(), columnNames.end(), genprefix+"_Particle") == columnNames.end()) 
		frame = frame.Define(genprefix+"_Particle", computeP4Vec, {genprefix+"_pt", genprefix+"_eta", genprefix+"_phi", genprefix+"_mass", genprefix+"_pdgId"}); 
		autoblacklist[genprefix+"_Particle"]++; 
		return frame; 
	}


	ROOT::RDF::RNode GetGenParticles(ROOT::RDF::RNode& frame, const int id, const std::string name, const std::string genprefix = "GenPart") 
	{
		frame = GetGenParticleCollection(frame, genprefix); 

		auto FilterParticles = [id](ROOT::VecOps::RVec<Particle> particles) 
		{
			ROOT::VecOps::RVec<Particle> result; 
			result.reserve(particles.size()); 
			std::copy_if(particles.begin(), particles.end(), std::back_inserter(result), [&](const Particle& particle){ return abs(particle.pdgid) == id; }); // Copy particles where the pdgid matches the requirement
			return result; 
		}; 

		frame = frame.Define(name, FilterParticles, {genprefix+"_Particle"}); 
		return frame; 
	}


	ROOT::RDF::RNode GetGenParticles(ROOT::RDF::RNode& frame, const std::string prefix, const std::string genprefix = "GenPart") 
	{
		const int id = PDGid[prefix]; 
		const std::string name = "Gen"+prefix; 
		autoblacklist[name]++; 
		return GetGenParticles(frame, id, name, genprefix); 
	}

	/*Particle IdentifyGenMuon(ROOT::RDF::RNode *frame, const std::string value = "") 
	{
		/*for (auto particle : genParticles) 
		{
			if (particle.pdgid == PDGid["Muon"]) 
			{
				// Could be the gen muon from Taumu

			}
		}*/
	//}


	int getMother(const ROOT::VecOps::RVec<float>& mother, int particle) 
	{
		return mother[particle]; 
	}


	bool isDescendantOf(int descendant, int ancestor, const ROOT::VecOps::RVec<float>& mothers) 
	{

		bool particleFound = false; 

		int particle = descendant; 

		while (particle != ancestor)
		{
			particle = mothers[particle]; 

			if (particle < 0) return false; 
		}

		return true; 
	}


	bool isAncestor(const int ancestor, const int descendant, const ROOT::VecOps::RVec<float>& mothers)
	{
	  	return isDescendantOf(descendant, ancestor, mothers); 
	}


	bool isLastCopy(int flags) 
	{
    	return flags & (1u << 13);
	}

	bool isHardProcess(int flag) 
	{
		return flag & (1u << 7); 
	}

	bool fromHardProcess(int flag) 
	{
		return flag & (1u << 8); 
	}


	std::vector<int> findMothers(int particle, int motherId, const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mothers, const ROOT::VecOps::RVec<int>& statusFlags, const long int flag = -999.) 
	{
		int currentId = 0; 

		std::vector<int> results; 

		while (particle >= 0 ) // put here >= 0 ?
		{
			particle = mothers[particle]; 

			bool flagOK = isLastCopy(statusFlags[particle]); 
			if (flag > 0 ) flagOK = flagOK && (statusFlags[particle] & flag); 
			//std::cout << (statusFlags[particle] & (1u << 8)) << " " << (statusFlags[particle] & flag) << " " << isLastCopy(statusFlags[particle]) << " " << flag << " " << flagOK << std::endl; 
			if ((abs(id[particle]) == motherId) && flagOK) results.push_back(particle); 
		}

		return results; 
	}

	std::vector<int> findMothers(int particle, std::string motherType, const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mothers, const ROOT::VecOps::RVec<int>& statusFlags, const long int flag = -999.) 
	{
		int motherId = PDGid[motherType]; 

		return findMothers(particle, motherId, id, mothers, statusFlags, flag); 
	}


	std::vector<int> findDescendants(int particle, int descendantType, const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mothers, const ROOT::VecOps::RVec<int>& statusFlags, const long int flag = -999.) 
	{
		std::vector<int> descendants; 

		//std::cout << "Sizes: " << id.size() << " " << mothers.size() << " " << statusFlags.size() << std::endl; 
		assert(id.size() == mothers.size()); 
		assert(id.size() == statusFlags.size()); 
		for (unsigned int i=0; i<id.size(); i++) 
		{
			//std::cout << "mother: " << mothers.at(i) << std::endl; 
			if (abs(id[i]) != descendantType) continue; 
			bool flagOK = isLastCopy(statusFlags[i]); 
			if (flag > 0 ) flagOK = flagOK && (statusFlags[i] & flag); 
			if (!flagOK) continue; 
			auto possibleMothers = findMothers(i, id[particle], id, mothers, statusFlags); 
			if (std::find(possibleMothers.begin(), possibleMothers.end(), particle) != possibleMothers.end()) 
			{
				descendants.push_back(i); 
			}
		}

		return descendants; 
	}


	std::vector<int> findDescendants(int particle, std::string descendantType, const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mothers, const ROOT::VecOps::RVec<int>& statusFlags, const long int flag = -999.) 
	{
		int descendantId = PDGid[descendantType]; 

		return findDescendants(particle, descendantId, id, mothers, statusFlags, flag); 
	}


	/* ChatGPT on status flags
	0  : isPrompt
	1  : isDecayedLeptonHadron
	2  : isTauDecayProduct
	3  : isPromptTauDecayProduct
	4  : isDirectTauDecayProduct
	5  : isDirectPromptTauDecayProduct
	6  : isDirectHadronDecayProduct
	7  : isHardProcess
	8  : fromHardProcess
	9  : isHardProcessTauDecayProduct
	10 : isDirectHardProcessTauDecayProduct
	11 : fromHardProcessBeforeFSR
	12 : isFirstCopy
	13 : isLastCopy     ← THIS ONE
	14 : isFirstCopyBeforeFSR
	*/


	GenMatchingResult DecayGenMatching(const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mother, const ROOT::VecOps::RVec<int>& statusFlag) 
	{

		// statusFlags bit helpers (bit numbers, zero-indexed)
    	//const unsigned int BIT_isLastCopy = (1u << 13);         // 13 => isLastCopy
    	const unsigned int BIT_isDirectTauDecayProduct = (1u << 5); // 5 => isDirectTauDecayProduct (useful)


    	std::vector<int> muons; 
    	muons.reserve(1); 
    	std::vector<int> taus; 
    	taus.reserve(2); 
    	std::vector<int> Higgses; 
    	Higgses.reserve(2); 
    	std::vector<int> bs; 
    	bs.reserve(2); 
    	std::vector<int> electrons; 
    	electrons.reserve(1); 
    	std::vector<int> Higgstob; 
    	Higgstob.reserve(1); 


    	GenMatchingResult result; 


    	const unsigned int hardProcess = (1u << 8); 

    	for (unsigned int i=0; i<id.size(); i++) 
    	{
    		if ((abs(id[i]) == PDGid["Higgs"]) &&  (isLastCopy(statusFlag[i]))) 
    		{
    			auto taudaughters = findDescendants(i, "Tau", id, mother, statusFlag, hardProcess); 
    			auto bdaughters = findDescendants(i, "b", id, mother, statusFlag, hardProcess); 

    			//if (taudaughters.size() && bdaughters.size()) continue;

    			if (bdaughters.size() < 2) continue; // making sure we have 2 b

    			bs.insert(bs.end(), bdaughters.begin(), bdaughters.end()); 
    			Higgstob.push_back(i); 

    		}
    	}

    	std::cout << "N bs: " << bs.size() << std::endl; 
    	if (bs.size() > 2) return result; // if more than 2 b, probably 4b or other weird stuff

    	//std::cout << "Sizes: " << id.size() << " " << mother.size() << " " << statusFlag.size() << std::endl; 


		for (unsigned int i=0; i<id.size(); i++) 
		{
			if ((abs(id[i]) == PDGid["Muon"]) && (isLastCopy(statusFlag[i]))) 
			{
				// Might be the muon
				auto localtaus = findMothers(i, "Tau", id, mother, statusFlag); 

				if (localtaus.size() < 1 ) continue; 

				//std::cout << "Flag: " << (1u << 7) << " " << hardProcess << std::endl; 
				auto localHiggses = findMothers(i, "Higgs", id, mother, statusFlag, hardProcess); 

				std::cout << "N taus: " << localtaus.size() << ", N Higgses: " << localHiggses.size() << std::endl; 

				if (localHiggses.size() < 1) continue; 


				std::vector<int> otherTaus = findDescendants(localHiggses[0], "Tau", id, mother, statusFlag, hardProcess); 
				std::cout << "N taus: " << otherTaus.size() << std::endl; 
				for (auto element : otherTaus) 
				{
					std::string text = RevertPDGid(id[element]); 
					/*for (auto it = PDGid.begin(); it != PDGid.end(); it++) 
					{
						if (it->second == element) text = it->first; 
					}*/
					std::cout << text << ": " << id[element] << " (id), " << RevertPDGid(id[mother[element]]) << " (mother = " << id[mother[element]] << "), " << statusFlag[element] << " (status)" << std::endl; 
				}
				otherTaus.erase(std::remove(otherTaus.begin(), otherTaus.end(), localtaus.at(0)), otherTaus.end()); // Remove the muonic tau
				bool notTauh = false; 
				for (unsigned int j=0; j<otherTaus.size(); j++) // Make sure the other tau decay is not electronic
				{
					if (findDescendants(otherTaus[j], "Electron", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
					{
						notTauh = true; 
					}
				}

				if (notTauh) continue; 




				muons.push_back(i); 
				Higgses.insert(Higgses.end(), localHiggses.begin(), localHiggses.end()); 
				taus.insert(taus.end(), localtaus.begin(), localtaus.end()); 
				taus.insert(taus.end(), otherTaus.begin(), otherTaus.end()); 

				result.decayType = TauhTaumu; 


			}

			if ((abs(id[i]) == PDGid["Electron"]) && (isLastCopy(statusFlag[i]))) 
			{
				// Might be the muon
				auto localtaus = findMothers(i, "Tau", id, mother, statusFlag); 

				if (localtaus.size() < 1 ) continue; 

				//std::cout << "Flag: " << (1u << 7) << " " << hardProcess << std::endl; 
				auto localHiggses = findMothers(i, "Higgs", id, mother, statusFlag, hardProcess); 

				std::cout << "N taus: " << localtaus.size() << ", N Higgses: " << localHiggses.size() << std::endl; 

				if (localHiggses.size() < 1) continue; 


				std::vector<int> otherTaus = findDescendants(localHiggses[0], "Tau", id, mother, statusFlag, hardProcess); 
				std::cout << "N taus: " << otherTaus.size() << std::endl; 
				for (auto element : otherTaus) 
				{
					std::string text = RevertPDGid(id[element]); 
					/*for (auto it = PDGid.begin(); it != PDGid.end(); it++) 
					{
						if (it->second == element) text = it->first; 
					}*/
					std::cout << text << ": " << id[element] << " (id), " << RevertPDGid(id[mother[element]]) << " (mother = " << id[mother[element]] << "), " << statusFlag[element] << " (status)" << std::endl; 
				}
				otherTaus.erase(std::remove(otherTaus.begin(), otherTaus.end(), localtaus.at(0)), otherTaus.end()); // Remove the muonic tau
				bool notTauh = false; 
				for (unsigned int j=0; j<otherTaus.size(); j++) // Make sure the other tau decay is not electronic
				{
					if (findDescendants(otherTaus[j], "Muon", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
					{
						notTauh = true; 
					}
				}

				if (notTauh) continue; 




				electrons.push_back(i); 
				Higgses.insert(Higgses.end(), localHiggses.begin(), localHiggses.end()); 
				taus.insert(taus.end(), localtaus.begin(), localtaus.end()); 
				taus.insert(taus.end(), otherTaus.begin(), otherTaus.end()); 

				result.decayType = TauhTaue; 


			}
		}


		assert(Higgstob.size() == 1); 
		assert(bs.size() == 2); 


		if (!((result.decayType == TauhTaumu) || (result.decayType == TauhTaue))) 
		{
			for (unsigned int i=0; i<id.size(); i++) 
	    	{
	    		if ((abs(id[i]) == PDGid["Higgs"]) &&  (isLastCopy(statusFlag[i])) && (i != Higgstob[0])) 
	    		{
	    			auto taudaughters = findDescendants(i, "Tau", id, mother, statusFlag, hardProcess); 

	    			//if (taudaughters.size() && bdaughters.size()) continue;

	    			if (taudaughters.size() < 2) continue; // making sure we have 2 taus

	    			bool notTauh = false; 
	    			for (unsigned int j=0; j<taudaughters.size(); j++) // Make sure the tau decays are not muonic or electronic
					{
						if (findDescendants(taudaughters[j], "Muon", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
						{
							notTauh = true; 
						}
						if (findDescendants(taudaughters[j], "Electron", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
						{
							notTauh = true; 
						}
					}

					if (notTauh) continue; 

	    			taus.insert(taus.end(), taudaughters.begin(), taudaughters.end()); 
	    			Higgses.push_back(i); 

	    		}
	    	}
	    	if (taus.size() == 2) result.decayType = TauhTauh; 
		}


		if (result.decayType == TauhTauh) 
		{
			assert((taus.size() == 2) && (muons.size() == 0) && (electrons.size() == 0)); 
			result.tau1 = taus.at(0); 
			result.tau2 = taus.at(1); 
		}
		if (result.decayType == TauhTaumu) 
		{
			assert((taus.size() == 2) && (muons.size() == 1) && (electrons.size() == 0)); 
			result.tau1 = taus[0]; 
			result.tau2 = taus[1]; 
			result.mu = muons[0]; 
		}
		if (result.decayType == TauhTaue) 
		{
			assert((taus.size() == 2) && (muons.size() == 0) && (electrons.size() == 1)); 
			result.tau1 = taus[0]; 
			result.tau2 = taus[1]; 
			result.e = electrons[0]; 
		}
		assert(Higgses.size() == 1); 


		// Filling the gen particle indices
		result.Htotau = Higgses[0]; 
		result.Htob = Higgstob[0]; 
		result.b1 = bs[0]; 
		result.b2 = bs[1]; 

		

		return result; 
	}


	GenMatchingResult DecayGenMatchingVBF(const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& mother, const ROOT::VecOps::RVec<int>& statusFlag) 
	{

		// statusFlags bit helpers (bit numbers, zero-indexed)
    	//const unsigned int BIT_isLastCopy = (1u << 13);         // 13 => isLastCopy
    	const unsigned int BIT_isDirectTauDecayProduct = (1u << 5); // 5 => isDirectTauDecayProduct (useful)


    	std::vector<int> muons; 
    	muons.reserve(1); 
    	std::vector<int> taus; 
    	taus.reserve(2); 
    	std::vector<int> Higgses; 
    	Higgses.reserve(2); 
    	std::vector<int> bs; 
    	bs.reserve(2); 
    	std::vector<int> electrons; 
    	electrons.reserve(1); 
    	std::vector<int> Higgstob; 
    	Higgstob.reserve(1); 


    	GenMatchingResult result; 


    	const unsigned int hardProcess = (1u << 8); 

    	for (unsigned int i=0; i<id.size(); i++) 
    	{
    		if ((abs(id[i]) == PDGid["Higgs"]) &&  (isLastCopy(statusFlag[i]))) 
    		{
    			auto taudaughters = findDescendants(i, "Tau", id, mother, statusFlag, hardProcess); 
    			auto bdaughters = findDescendants(i, "b", id, mother, statusFlag, hardProcess); 

    			//if (taudaughters.size() && bdaughters.size()) continue;

    			if (bdaughters.size() < 2) continue; // making sure we have 2 b

    			bs.insert(bs.end(), bdaughters.begin(), bdaughters.end()); 
    			Higgstob.push_back(i); 

    		}
    	}

    	std::cout << "N bs: " << bs.size() << std::endl; 
    	if (bs.size() > 2) return result; // if more than 2 b, probably 4b or other weird stuff

    	//std::cout << "Sizes: " << id.size() << " " << mother.size() << " " << statusFlag.size() << std::endl; 


		for (unsigned int i=0; i<id.size(); i++) 
		{
			if ((abs(id[i]) == PDGid["Muon"]) && (isLastCopy(statusFlag[i]))) 
			{
				// Might be the muon
				auto localtaus = findMothers(i, "Tau", id, mother, statusFlag); 

				if (localtaus.size() < 1 ) continue; 

				//std::cout << "Flag: " << (1u << 7) << " " << hardProcess << std::endl; 
				auto localHiggses = findMothers(i, "Higgs", id, mother, statusFlag, hardProcess); 

				std::cout << "N taus: " << localtaus.size() << ", N Higgses: " << localHiggses.size() << std::endl; 

				if (localHiggses.size() < 1) continue; 


				std::vector<int> otherTaus = findDescendants(localHiggses[0], "Tau", id, mother, statusFlag, hardProcess); 
				std::cout << "N taus: " << otherTaus.size() << std::endl; 
				for (auto element : otherTaus) 
				{
					std::string text = RevertPDGid(id[element]); 
					/*for (auto it = PDGid.begin(); it != PDGid.end(); it++) 
					{
						if (it->second == element) text = it->first; 
					}*/
					std::cout << text << ": " << id[element] << " (id), " << RevertPDGid(id[mother[element]]) << " (mother = " << id[mother[element]] << "), " << statusFlag[element] << " (status)" << std::endl; 
				}
				otherTaus.erase(std::remove(otherTaus.begin(), otherTaus.end(), localtaus.at(0)), otherTaus.end()); // Remove the muonic tau
				bool notTauh = false; 
				for (unsigned int j=0; j<otherTaus.size(); j++) // Make sure the other tau decay is not electronic
				{
					if (findDescendants(otherTaus[j], "Electron", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
					{
						notTauh = true; 
					}
				}

				if (notTauh) continue; 




				muons.push_back(i); 
				Higgses.insert(Higgses.end(), localHiggses.begin(), localHiggses.end()); 
				taus.insert(taus.end(), localtaus.begin(), localtaus.end()); 
				taus.insert(taus.end(), otherTaus.begin(), otherTaus.end()); 

				result.decayType = TauhTaumu; 


			}

			if ((abs(id[i]) == PDGid["Electron"]) && (isLastCopy(statusFlag[i]))) 
			{
				// Might be the muon
				auto localtaus = findMothers(i, "Tau", id, mother, statusFlag); 

				if (localtaus.size() < 1 ) continue; 

				//std::cout << "Flag: " << (1u << 7) << " " << hardProcess << std::endl; 
				auto localHiggses = findMothers(i, "Higgs", id, mother, statusFlag, hardProcess); 

				std::cout << "N taus: " << localtaus.size() << ", N Higgses: " << localHiggses.size() << std::endl; 

				if (localHiggses.size() < 1) continue; 


				std::vector<int> otherTaus = findDescendants(localHiggses[0], "Tau", id, mother, statusFlag, hardProcess); 
				std::cout << "N taus: " << otherTaus.size() << std::endl; 
				for (auto element : otherTaus) 
				{
					std::string text = RevertPDGid(id[element]); 
					/*for (auto it = PDGid.begin(); it != PDGid.end(); it++) 
					{
						if (it->second == element) text = it->first; 
					}*/
					std::cout << text << ": " << id[element] << " (id), " << RevertPDGid(id[mother[element]]) << " (mother = " << id[mother[element]] << "), " << statusFlag[element] << " (status)" << std::endl; 
				}
				otherTaus.erase(std::remove(otherTaus.begin(), otherTaus.end(), localtaus.at(0)), otherTaus.end()); // Remove the muonic tau
				bool notTauh = false; 
				for (unsigned int j=0; j<otherTaus.size(); j++) // Make sure the other tau decay is not electronic
				{
					if (findDescendants(otherTaus[j], "Muon", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
					{
						notTauh = true; 
					}
				}

				if (notTauh) continue; 




				electrons.push_back(i); 
				Higgses.insert(Higgses.end(), localHiggses.begin(), localHiggses.end()); 
				taus.insert(taus.end(), localtaus.begin(), localtaus.end()); 
				taus.insert(taus.end(), otherTaus.begin(), otherTaus.end()); 

				result.decayType = TauhTaue; 


			}
		}


		assert(Higgstob.size() == 1); 
		assert(bs.size() == 2); 


		if (!((result.decayType == TauhTaumu) || (result.decayType == TauhTaue))) 
		{
			for (unsigned int i=0; i<id.size(); i++) 
	    	{
	    		if ((abs(id[i]) == PDGid["Higgs"]) &&  (isLastCopy(statusFlag[i])) && (i != Higgstob[0])) 
	    		{
	    			auto taudaughters = findDescendants(i, "Tau", id, mother, statusFlag, hardProcess); 

	    			//if (taudaughters.size() && bdaughters.size()) continue;

	    			if (taudaughters.size() < 2) continue; // making sure we have 2 taus

	    			bool notTauh = false; 
	    			for (unsigned int j=0; j<taudaughters.size(); j++) // Make sure the tau decays are not muonic or electronic
					{
						if (findDescendants(taudaughters[j], "Muon", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
						{
							notTauh = true; 
						}
						if (findDescendants(taudaughters[j], "Electron", id, mother, statusFlag).size() > 0) // If we find an electron in the other tau decay
						{
							notTauh = true; 
						}
					}

					if (notTauh) continue; 

	    			taus.insert(taus.end(), taudaughters.begin(), taudaughters.end()); 
	    			Higgses.push_back(i); 

	    		}
	    	}
	    	if (taus.size() == 2) result.decayType = TauhTauh; 
		}


		if (result.decayType == TauhTauh) 
		{
			assert((taus.size() == 2) && (muons.size() == 0) && (electrons.size() == 0)); 
			result.tau1 = taus.at(0); 
			result.tau2 = taus.at(1); 
		}
		if (result.decayType == TauhTaumu) 
		{
			assert((taus.size() == 2) && (muons.size() == 1) && (electrons.size() == 0)); 
			result.tau1 = taus[0]; 
			result.tau2 = taus[1]; 
			result.mu = muons[0]; 
		}
		if (result.decayType == TauhTaue) 
		{
			assert((taus.size() == 2) && (muons.size() == 0) && (electrons.size() == 1)); 
			result.tau1 = taus[0]; 
			result.tau2 = taus[1]; 
			result.e = electrons[0]; 
		}
		assert(Higgses.size() == 1); 


		// Filling the gen particle indices
		result.Htotau = Higgses[0]; 
		result.Htob = Higgstob[0]; 
		result.b1 = bs[0]; 
		result.b2 = bs[1]; 

		// TODO: add VBF gen particle selection 

		

		return result; 
	}


	template<typename T>
	TLorentzVector getP4(int idx, const ROOT::VecOps::RVec<T>& pt, const ROOT::VecOps::RVec<T>& eta, const ROOT::VecOps::RVec<T>& phi, const ROOT::VecOps::RVec<T>& m) 
	{
		TLorentzVector P4; 

		unsigned int n = pt.size(); 
		assert(eta.size() == n); 
		assert(phi.size() == n); 
		assert(m.size() == n); 
		// TODO: turn off the above checks when running for speed

		if ((idx > pt.size() -1) || (idx < 0)) return P4; 

		P4.SetPtEtaPhiM(pt[idx], eta[idx], phi[idx], m[idx]); 
		return P4;  
	}


	GenMatchingResult MatchVBFJets(GenMatchingResult& genDecay, const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& mass, const ROOT::VecOps::RVec<float>& genPt, const ROOT::VecOps::RVec<float>& genEta, const ROOT::VecOps::RVec<float>& genPhi, const ROOT::VecOps::RVec<float>& genMass) 
	{
		GenMatchingResult result(genDecay); 

		int nJets = pt.size(); 
		assert(eta.size() == nJets); 
		assert(phi.size() == nJets); 
		assert(mass.size() == nJets); 

		// VBF selection thesholds
		double ptThres = 20.; 
		double etaThres = 5.; 

		double isoThres = 0.4; 

		double diJetMassThres = 300.; 

		double deltaEtaThres = 3.;

		TLorentzVector b1 = getP4(genDecay.b1, pt, eta, phi, mass); 
		TLorentzVector b2 = getP4(genDecay.b2, pt, eta, phi, mass); 
		TLorentzVector tau1 = getP4(genDecay.tau1, pt, eta, phi, mass); 
		TLorentzVector tau2 = getP4(genDecay.tau2, pt, eta, phi, mass); 


		TLorentzVector jet1, jet2; 


		double diJetMass = -999.; 

		int jetIdx1 = -999; 
		int jetIdx2 = -999; 
		

		for (unsigned int i=0; i<nJets; i++) 
		{
			// VBF jet selection requirements 
			if (pt[i] < ptThres) continue; 
			if (eta[i] > etaThres) continue; 

			jet1.SetPtEtaPhiM(pt[i], eta[i], phi[i], mass[i]); 

			if (jet1.DeltaR(b1) < isoThres) continue; 
			if (jet1.DeltaR(b2) < isoThres) continue; 
			if (jet1.DeltaR(tau1) < isoThres) continue; 
			if (jet1.DeltaR(tau2) < isoThres) continue; 


			for (unsigned int j=0; j<nJets; j++) 
			{
				if (j == i) continue; 

				// VBF jet selection requirements 
				if (pt[i] < ptThres) continue; 
				if (eta[i] > etaThres) continue; 

				jet2.SetPtEtaPhiM(pt[j], eta[j], phi[j], mass[j]); 

				if (jet2.DeltaR(b1) < isoThres) continue; 
				if (jet2.DeltaR(b2) < isoThres) continue; 
				if (jet2.DeltaR(tau1) < isoThres) continue; 
				if (jet2.DeltaR(tau2) < isoThres) continue; 

				// Both jets pass the VBF jet preselection requirements and do not overlap wth the Higgs decay products 


				double currentDiJetMass = (jet1 + jet2).M(); 

				// Selection on jet pair
				if (currentDiJetMass < diJetMassThres) continue; 
				if (abs(jet1.Eta() - jet2.Eta()) < deltaEtaThres) continue;
				if (jet1.Eta()*jet2.Eta() > 0.) continue;

				if (currentDiJetMass > diJetMass) 
				{
					diJetMass = currentDiJetMass; 

					jetIdx1 = i; 
					jetIdx2 = j; 
				}


			}

		}

		// Final selection on jet pair
		//bool finalSel = (diJetMass > diJetMassThres) && (abs(eta[jetIdx1] - eta[jetIdx2]) > deltaEtaThres) && (eta[jetIdx1]*eta[jetIdx2] < 0.); 
		//if (!finalSel ) return result; 

		result.VBFgenJet1 = jetIdx1; 
		result.VBFgenJet2 = jetIdx2; 

		return result;
	}


	template<typename T>
	T overflowProtected(const ROOT::VecOps::RVec<T>& collection, const long int index) 
	{
		if ((index > collection.size() -1) || (index < 0)) return defaultValue<T>(); 
		return collection[index]; 
	}


	float deltaR(int genParticle, const ROOT::VecOps::RVec<float>& genPt, const ROOT::VecOps::RVec<float>& genEta, const ROOT::VecOps::RVec<float>& genPhi, const ROOT::VecOps::RVec<float>& genMass, const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& m, const ROOT::VecOps::RVec<float>& pdgId = std::vector<float>(), const std::string particleType = "") 
	{
		double genpt = overflowProtected(genPt, genParticle); 
		if (genpt < 0) return genpt; // If gen particle is out of bounds, return default overflow value
		TLorentzVector gen; 
		gen.SetPtEtaPhiM(genpt, overflowProtected(genEta, genParticle), overflowProtected(genPhi, genParticle), overflowProtected(genMass, genParticle)); 

		int type = PDGid[particleType]; 


		float dR = 999.; 

		TLorentzVector P4; 
		for (unsigned int i = 0; i<pt.size(); i++) 
		{
			if ((pdgId.size()) && (abs(pdgId[i])) != type) continue; 
			P4.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]); 

			float currentdR = P4.DeltaR(gen); 

			if (currentdR < dR) dR = currentdR; 
		}

		return dR; 
	}


	int closestMatch(int genParticle, const ROOT::VecOps::RVec<float>& genPt, const ROOT::VecOps::RVec<float>& genEta, const ROOT::VecOps::RVec<float>& genPhi, const ROOT::VecOps::RVec<float>& genMass, const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& m, const ROOT::VecOps::RVec<float>& pdgId = std::vector<float>(), const std::string particleType = "") 
	{
		if (overflowProtected(genPt, genParticle) < 0.) return defaultValue<int>(); 
		TLorentzVector gen; 
		gen.SetPtEtaPhiM(overflowProtected(genPt, genParticle), overflowProtected(genEta, genParticle), overflowProtected(genPhi, genParticle), overflowProtected(genMass, genParticle)); 

		int type = PDGid[particleType]; 


		float dR = 999.; 
		int closest = -999; 

		TLorentzVector P4; 
		for (unsigned int i = 0; i<pt.size(); i++) 
		{
			if ((pdgId.size() != pt.size()) || (abs(pdgId[i]) != type)) continue; 
			P4.SetPtEtaPhiM(pt[i], eta[i], phi[i], m[i]); 

			float currentdR = P4.DeltaR(gen); 

			if ((currentdR < dR) && (currentdR < drThres)) 
			{
				dR = currentdR; 
				closest = i; 
			}
		}

		return closest; 
	}


	float deltaR(int genParticle, const ROOT::VecOps::RVec<float>& genPt, const ROOT::VecOps::RVec<float>& genEta, const ROOT::VecOps::RVec<float>& genPhi, const ROOT::VecOps::RVec<float>& genMass, const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& m, const int cand) 
	{
		double genpt = overflowProtected(genPt, genParticle); 
		double candpt = overflowProtected(pt, cand); 
		if (genpt < 0. || candpt < 0.) return genpt; // If gen particle is out of bounds, return default overflow value
		TLorentzVector gen, candidate; 
		gen.SetPtEtaPhiM(genpt, overflowProtected(genEta, genParticle), overflowProtected(genPhi, genParticle), overflowProtected(genMass, genParticle)); 
		candidate.SetPtEtaPhiM(candpt, eta[cand], phi[cand], m[cand]); 

		double dR = gen.DeltaR(candidate); 
		

		return dR; 
	}


	int RecoMuon(const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& id, const ROOT::VecOps::RVec<float>& dz, const ROOT::VecOps::RVec<float>& dxy) 
	{
		int n = pt.size(); 
		assert(eta.size() == n); 
		assert(phi.size() == n); 
		assert(mass.size() == n); 

		// Muon selection requirements 
		double ptThres = 20.; 
		double etaThres = 2.4; 
		double dzThres = 0.2; 
		double dxyThres = 0.045; 
		double idThres = 0.2; 


		int muon = -999.; 

		for (unsigned int i=0; i<n; i++) 
		{
			if (pt[i] < ptThres) continue; 
			if (eta[i] > etaThres) continue; 
			if (dz[i] > dzThres) continue; 
			if (dxy[i] > dxyThres) continue; 
			if (id[i] < idThres) continue; 

			muon = i;
			break; 
		}

		return muon; 
	}


	int RecoTauJet(const ROOT::VecOps::RVec<float>& pt, const ROOT::VecOps::RVec<float>& eta, const ROOT::VecOps::RVec<float>& phi, const ROOT::VecOps::RVec<float>& id) 
	{
		int nJets = pt.size(); 
		assert(eta.size() == nJets); 
		assert(phi.size() == nJets); 
		assert(mass.size() == nJets); 

		// Muon selection requirements 
		double ptThres = 20.; 
		double etaThres = 2.4; 
		double idThres = 0.2; 
		


		int tauJet = -999.; 

		for (unsigned int i=0; i<nJets; i++) 
		{
			if (pt[i] < ptThres) continue; 
			if (eta[i] > etaThres) continue; 
			if (id[i] < idThres) continue; 

			tauJet = i;
			break; 
		}

		return tauJet; 
	}



	void AddColumn(ROOT::RDF::RNode* df, const std::string &newColName) {
    	*df = df->Define(newColName, [](){ return 42; });
	}


	// Making non-templated wrappwrs
  	//ROOT::VecOps::RVec< R4Vec<float> > computeP4Vec_f(const ROOT::VecOps::RVec<float>& pt, ...)
  	//{ return computeP4Vec<float>(pt, ...); }

  	//ROOT::VecOps::RVec< R4Vec<double> > computeP4Vec_d(const ROOT::VecOps::RVec<double>& pt, ...)
  	//{ return computeP4Vec<double>(pt, ...); }


}


// // --- Explicit instantiations for the types you will use ---
// // This ensures the symbol for computeP4Vec<float> and computeP4Vec<double> exists
// // (only necessary if you compile this into a shared library; harmless otherwise)	

// template ROOT::VecOps::RVec< R4Vec<float> > Ana::computeP4Vec<float>(const ROOT::VecOps::RVec<float>&,
// const ROOT::VecOps::RVec<float>&,
// const ROOT::VecOps::RVec<float>&,
// const ROOT::VecOps::RVec<float>&);	

// template ROOT::VecOps::RVec< R4Vec<double> > Ana::computeP4Vec<double>(const ROOT::VecOps::RVec<double>&,
// const ROOT::VecOps::RVec<double>&,
// const ROOT::VecOps::RVec<double>&,
// const ROOT::VecOps::RVec<double>&);	

// // Also instantiate GetP4 for float/double
// template ROOT::RDF::RNode* Ana::GetP4<float>(ROOT::RDF::RNode*, TString, TString);
// template ROOT::RDF::RNode* Ana::GetP4<double>(ROOT::RDF::RNode*, TString, TString);



#endif


