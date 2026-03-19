#ifndef AnaBuildingBlocks_hxx
#define AnaBuildingBlocks_hxx
#include <string>
#include <vector>
#include "TChain.h"
#include "ROOT/RDataFrame.hxx"
#include <thread>


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
	//#include <chrono>
	//#include <thread>
	//std::this_thread::sleep_for(static_cast<std::chrono::seconds>(timeInSec));
	sleep(timeInSec); 
}


namespace Ana 
{
	static std::unordered_map<std::string, int> autoblacklist; 

	std::vector<std::string>& purgeColumns(std::vector<std::string> &&columns, const std::vector<std::string>& blacklist = {})
	{
			auto list(autoblacklist); 
			list.reserve(list.size()+blacklist.size()); 
			for (auto item : blacklist) 
			{
				list[item]++;
			}

			//char_separator<char> sep("*");
   			char sep = '*'; 

			std::vector<std::vector<std::string> > tokenizedBlacklistItems; 
			for (auto item : list) 
			{
				if (item.first.find(sep) != std::string::npos) // The blacklisted item contains at least one wildcard
				{
					const char *str = item.first.data(); 
			
					std::vector<std::string> chunks;
				    do
				    {
				        const char *begin = str;					

				        while(*str != sep && *str) str++;					

				        if (begin != str) chunks.push_back(std::string(begin, str)); 
				    } 
				    while (0 != *str++);
				    tokenizedBlacklistItems.push_back(chunks); 
				}
			}

   			// a lambda that checks if `s` is in the blacklist
   			auto is_blacklisted = [&list, &tokenizedBlacklistItems](const std::string &s)  { 
   				bool contained = (list.find(s) != list.end()); 

   				
   				bool matched = false; 
   				for (auto pattern : tokenizedBlacklistItems) 
   				{
   					int pos = 0; 
   					bool foundAll = true; 
   					for (auto chunk : pattern) 
   					{
   						int i = s.find(chunk); 
   						if ((i < pos) || (i == std::string::npos)) foundAll = false; 
   						//std::cout << s << " " << chunk << " " << i << " " << foundAll << std::endl; 
   					}

   					if (foundAll) 
   					{
   						matched = true; 
   						break; 
   					}
   				}
   				//if (matched) std::cout << s << std::endl; 
   				
   				return (contained || matched); 
   			};

   			// removing elements from std::vectors is not pretty, see https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
   			columns.erase(std::remove_if(columns.begin(), columns.end(), is_blacklisted), columns.end());

   			return columns; 
	}

	// Workaround to save RDataFrames containing string branches
	template<typename T>
	void fixStringVariables(T &dataframe)
	{
		#include "stringbranches.gcf"

		for (auto branch : stringbranches) // Hack to fix string branche 
		{
			dataframe = dataframe.Redefine(branch, [](const ROOT::RVec<std::string> &v) {return std::vector<std::string>(v.begin(), v.end());}, {branch}); 
		}
	}

	void normaliseBinContent(TH1* hist) 
	{
	  	for (Int_t i=0; i<hist->GetNbinsX(); i++) 
	  	{
	    	hist->SetBinContent(i, hist->GetBinContent(i)/hist->GetBinWidth(i)); 
	  	}
	}



	// Kept for legacy purposes
	namespace old 
	{

		std::vector<std::string>& purgeColumns(std::vector<std::string> &&columns)
		{
			const std::vector<std::string> blacklist = {"v_taucandidates", "b_tau"};
				// a lambda that checks if `s` is in the blacklist
				auto is_blacklisted = [&blacklist](const std::string &s)  { return std::find(blacklist.begin(), blacklist.end(), s) != blacklist.end(); };

				// removing elements from std::vectors is not pretty, see https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
				columns.erase(std::remove_if(columns.begin(), columns.end(), is_blacklisted), columns.end());

				return columns; 
		}

	}

}

#endif


