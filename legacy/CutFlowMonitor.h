#ifndef CutFlowMonitor_h
#define CutFlowMonitor_h


#include "TH1.h"
#include <map>
#include <vector>
#include "TDirectory.h"
#include <iostream>


class CutFlow 
{
public: 
    CutFlow(const TString& name, const TString& title) : histName(name), histTitle(title) {}; // We do not want to set a default directory. If it is not set, the histogram will be written to the gDirectory at this time. 
    CutFlow(const TString& name, const TString& title, TDirectory *const aDirectory) : CutFlow(name, title) 
    {
        SetDirectory(aDirectory); 
    }; 
    CutFlow(const TString& name, const TString& title, TDirectory *const aDirectory, bool wantRelEff) : CutFlow(name, title, aDirectory) 
    {
        doRelEff = wantRelEff; 
    }
    CutFlow(const TString& name, const TString& title, bool wantRelEff) : CutFlow(name, title) 
    {
        doRelEff = wantRelEff; 
    }
    CutFlow() 
    {
        TString name = TString::Format("cutflow_%d", no); 
        histName = name; 
        histTitle = name; 
        no++; 
    }
    void Add(const TString& cutname, Int_t value) 
    {
        cutyields.push_back(std::pair<TString, Int_t>(cutname, value)); 
    }; 

    void sort() 
    {
        std::sort(cutyields.begin(), cutyields.end(), sortbysecond); 
    }; 

    void SetDirectory(TDirectory *const aDirectory) 
    {
        directory = aDirectory; 
        wasSet = true; 
    }; 

    void Increment(TString cutname, Float_t step = 1.) 
    {
        auto it = std::find_if(cutyields.begin(), cutyields.end(), keymatch(cutname)); 
        if (it != cutyields.end()) 
        {
            it->second += step; 
        }
        else 
        {
            cutyields.push_back(std::make_pair(cutname, step)); 
        }
    }; 

    bool WriteHistogram() 
    {
        auto cachedir = gDirectory; 

        if (wasSet) 
        {
            directory->cd(); // TODO: cache the previous directory and restore it after writing 
        }
        TH1 *hist = GenerateHistogram(); 
        hist->Write(); 

        if (doRelEff) 
        {
            TH1 *effhist = MakeEffHistogram(); 
            effhist->Write(); 
        }

        cachedir->cd(); 

        return true; 
    }; 

    TH1* GenerateHistogram() 
    {
        TH1* histo = new TH1F(histName, histTitle, cutyields.size(), -0.5, cutyields.size()-0.5); // TODO: delete object
        for (unsigned int i=0; i<cutyields.size(); i++) 
        {
            histo->SetBinContent(i+1, cutyields.at(i).second); 
            histo->GetXaxis()->SetBinLabel(i+1, cutyields.at(i).first); 
        }
        return histo; 
    }; 

    TH1* MakeEffHistogram() 
    {
        // We make each step's relative efficiency
        TH1* histo = new TH1F(histName+"_eff", histTitle+" relative efficiency", cutyields.size()-1, -0.5, cutyields.size()-1.5); //TODO: delete 
        for (unsigned int i=1; i<cutyields.size(); i++) 
        {
            histo->SetBinContent(i, static_cast<Double_t>(cutyields.at(i).second)/static_cast<Double_t>(cutyields.at(i-1).second)); 
            histo->GetXaxis()->SetBinLabel(i, cutyields.at(i).first); 
        }
        return histo;
    }

    const std::vector<std::pair<TString, Int_t> >& Content() const
    {
        return cutyields; 
    }

    void Merge(const CutFlow& other) 
    {
        for (auto element : other.Content()) 
        {
            Increment(element.first, element.second); 
        }

    }

    int Evaluate(const TString& cutname) 
    {
        auto it = std::find_if(cutyields.begin(), cutyields.end(), keymatch(cutname)); 
        if (it != cutyields.end()) 
        {
            return it->second; 
        }
        else 
        {
            return 0; 
        }
    }

    void Print() 
    {
        std::cout << "CutFlow: " << histTitle << std::endl; 
        for (unsigned int i=0; i<cutyields.size(); i++) 
        {
            std::cout << "\t"+cutyields.at(i).first << ": " << cutyields.at(i).second << std::endl; 
        }
        if (cutyields.size()==0) std::cout << "[empty]" << std::endl; 
        std::cout << std::endl; 
    }

    static bool sortbysecond(const std::pair<TString, Int_t>& a, const std::pair<TString, Int_t>& b) 
    {
        return a.second > b.second; 
    }

    void reset() 
    {
        for (unsigned int i=0; i<cutyields.size(); i++) 
        {
            cutyields[i].second = 0; 
        }
    }; 

    void clear() 
    {
        cutyields.clear(); 
    }; 



  Float_t nCuts; 

private: 
    // Matching function for use in std::find_if 
    struct keymatch 
    {
        keymatch(TString key) : keyToMatch(key) {};  
        bool operator()(std::pair<TString, Int_t> element) 
        {
            return element.first == keyToMatch; 
        }; 
        TString keyToMatch; 
    }; 

    std::vector<std::pair<TString, Int_t> > cutyields; 

    TDirectory *directory = nullptr; 
    bool wasSet = false; 

    TString histName; 
    TString histTitle; 
    bool doRelEff = false; 

    inline static int no; 


}; 





#endif