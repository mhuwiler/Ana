#ifndef TAU_H
#define TAU_H
#ifdef __MAKECINT__ 
#pragma link C++ class vector<Tau>; 
#endif // __MAKECINT__

#include <TObject.h>
#include <vector>
#include <map>


class Tau : public TObject 
{
  public: 
	
	double dnn1 = -999.; 
	double dnn2 = -000.; 
	double dnn3 = -999.; 
	double sumdnn = -999.; 

	#import "tauvariables.gcf"


	ClassDef(Tau, 0);

	Tau() = default; 

	Tau(const float& initPt, const float& initEta, const float& initPhi, const int charge, const float& mass) 
	  : pt(initPt), eta(initEta), phi(initPhi), q(charge), m(mass) {}

	void SetIndices(const int myIdx1, const int myIdx2, const int myIdx3) 
	{
		idx1 = myIdx1; 
		idx2 = myIdx2; 
		idx3 = myIdx3; 
	}

	/*void SetDNN(const float& myDNN1, const float& myDNN2, const float&myDNN3, const float& mySumOfDNN) 
	{
		dnn1 = myDNN1; 
		dnn2 = myDNN2; 
		dnn3 = myDNN3; 
		sumdnn = mySumOfDNN; 
	}*/

	void SetKinematics(const float& vertexProb, const float& flightSig, const float& longIP) 
	{
		vprob = vertexProb; 
		fsig = flightSig; 
		lip = longIP; 
	}


	bool operator>(const Tau& other) const 
	{ 
    	return pt > other.pt;
  	}

  	bool operator<(const Tau& other) const 
	{ 
    	return pt < other.pt;
  	}


  	/*void SetEventKinematics(const float& angle, const float& deltaR, const float& flightLenght, const float& IPtoPV, const float& sigIPtoPV) 
  	{
  		alpha = angle; 
  		dr = deltaR; 
  		fl = flightLenght; 
  		pvip = IPtoPV; 
  		pvips = sigIPtoPV; 
  	}*/


  	void SetDau1Kin(const float& dauPt, const float& dauEta, const float& dauPhi) 
  	{
  		pi1pt = dauPt; 
  		pi1eta = dauEta; 
  		pi1phi = dauPhi; 
  	}


  	void SetDau2Kin(const float& dauPt, const float& dauEta, const float& dauPhi) 
  	{
  		pi2pt = dauPt; 
  		pi2eta = dauEta; 
  		pi2phi = dauPhi; 
  	}


  	void SetDau3Kin(const float& dauPt, const float& dauEta, const float& dauPhi) 
  	{
  		pi3pt = dauPt; 
  		pi3eta = dauEta; 
  		pi3phi = dauPhi; 
  	}

  	void SetRhoMasses(const float& m12, const float& m23) 
  	{
  		rhomass1 = m12; 
  		rhomass2 = m23; 
  	}

  	void SetMatch(const bool genMatch1, const bool genMatch2, const bool genMatch3) 
  	{
  		match1 = genMatch1; 
  		match2 = genMatch2; 
  		match3 = genMatch3; 
  	}

  	void SetMatch(const int genMatch1, const int genMatch2, const int genMatch3) 
  	{
  		match1 = static_cast<bool>(genMatch1); 
  		match2 = static_cast<bool>(genMatch2); 
  		match3 = static_cast<bool>(genMatch3); 
  	}

  	void SetBQuantities(const float mass, const float q2) 
  	{
  		B_m = mass; 
  		B_q2 = q2; 
  	}





	#import "taubranchwriting.gcf"

  	// Non standard writing 
  	static float WriteDNN1(const Tau& tau) 
	{
		return tau.dnn1; 
	}; 

	static float WriteDNN2(const Tau& tau) 
	{
		return tau.dnn2; 
	}; 

	static float WriteDNN3(const Tau& tau) 
	{
		return tau.dnn3; 
	}; 

	static float WriteSumDNN(const Tau& tau) 
	{
		return tau.sumdnn; 
	}; 

	/*static int WriteMatch(const Tau& tau) 
	{
		return static_cast<int>(tau.match1 || tau.match2 || tau.match3); 
	}; */

	static int WriteSumMatch(const Tau& tau) 
	{
		return static_cast<int>(tau.match1) + static_cast<int>(tau.match2) + static_cast<int>(tau.match3); 
	}; 

}; 


typedef std::vector<Tau> TauCollection; 
//typedef std::map<Tau, float> PtOrderedTauCollection; 


#endif 
