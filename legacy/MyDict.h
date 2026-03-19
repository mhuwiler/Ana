// MyDict.h
#ifndef MYDICT_H
#define MYDICT_H

#include <vector>
#include "TLorentzVector.h"
#include "ROOT/RVec.hxx"

#pragma link C++ class TLorentzVector+;
#pragma link C++ class ROOT::VecOps::RVec<TLorentzVector>+;
#pragma link C++ class std::vector<TLorentzVector>+;
#pragma link C++ class ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiM4D<float>+; 
#pragma link C++ class ROOT::VecOps::RVec<ROOT::Math::PtEtaPhiM4D<double>+; 

#endif
