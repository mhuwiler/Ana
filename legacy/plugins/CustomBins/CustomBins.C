#include <cstdlib>
#include <vector>
#include <iostream>
#include "CustomBins.h"


CustomBins::CustomBins(Int_t nBins, Double_t rangeMin, Double_t rangeMax) 
{
	if ((rangeMin > rangeMax) || (nBins < 0)) 
	{
		std::cerr << "ERROR: The parameters provided to the function are odd. Please provide them in the following order :\n   CustomBins(  nBins (>0), rangeMin, rangeMax. )\nThe code returns an empty vector. " << std::endl; 
		return; 
	}
	Double_t step = (rangeMax - rangeMin)/static_cast<Double_t>(nBins + 1); 
	for (unsigned int i=0; i<=nBins+1; i++)
	{
		bins.push_back(rangeMin + static_cast<Double_t>(i) * step); 
	}
}

CustomBins::CustomBins(std::vector<Double_t> binEdges) 
{	
	if (binEdges.size() > 2) 
	{
		bins = binEdges; 
	}
	else 
	{
		ErrorTooSmallBinVector(); 
		return; 
	}
}

CustomBins::CustomBins(Int_t size, Double_t* binEdges) 
{	
	if (size <= 2) 
	{
		ErrorTooSmallBinVector(); 
		return; 
	}
	for (unsigned int i=0; i<size; i++) {
		bins.push_back(binEdges[i]); 
	}
}

void CustomBins::ErrorTooSmallBinVector() {
	std::cerr << "ERROR: The vector contains less than 2 bin edges (necessary for a single bin histogram), it is impossible to create the bin vector. " << std::endl; 
	return; 
}
