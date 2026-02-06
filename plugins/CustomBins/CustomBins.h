#ifndef customBins_class
#define customBins_class
#include <vector>

#include "TROOT.h"


class CustomBins 
{
	public:
	CustomBins(Int_t nBins, Double_t rangeMin, Double_t rangeMax); 

	CustomBins(std::vector<Double_t> binEdges); 

	CustomBins(Int_t size, Double_t* binEdges); 
	
	CustomBins() = delete; 	// Make sure there is no default constructor generated. 

	inline Double_t* GetBinEdges() { return bins.data(); }

	inline Int_t GetNBins() { return bins.size()-1; }

	inline std::vector<Double_t> GetBinVector() { return bins; }

	private:
	void ErrorTooSmallBinVector(); 
	std::vector<Double_t> bins; 
}; 


#endif // Define customBins_class 
