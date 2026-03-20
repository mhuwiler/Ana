#ifndef TAU_H
#define TAU_H
#ifdef __MAKECINT__ 
#pragma link C++ class vector<Particle>; 
#endif // __MAKECINT__

#include <TObject.h>
#include <vector>
#include <map>
#include "TLorentzVector.h"


class Particle : public TObject 
{
  public: 
	
	double pt; 
	double eta; 
	double phi; 
	double m; 
	int pdgid; 
	TLorentzVector P4; 

	//#import "tauvariables.gcf"


	ClassDef(Particle, 0);

	Particle() = default; 

	inline Particle(const float& initPt, const float& initEta, const float& initPhi, const double mass, const int initpdgid = 0) 
	  : pt(initPt), eta(initEta), phi(initPhi), m(mass), pdgid(initpdgid) {
	  		P4.SetPtEtaPhiM(pt, eta, phi, m); 
	  }

	


typedef std::vector<Particle> ParticleCollection; 
//typedef std::map<Tau, float> PtOrderedTauCollection; 

}; 


#endif 
