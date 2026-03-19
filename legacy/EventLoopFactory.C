#include "TChain.h"
#include "TFile.h"
#include "AnaBuildingBlocks.C"



class EventLoopFactory 
{
	public: 

	EventLoopFactory(TTree *tree, bool fulldecay = false, int maxnum = -1) : fTree(tree), fFullDecay(fulldecay), fMaxNum(maxnum) { Initialise(); }; 

	EventLoopFactory() = delete; 

	void Initialise() 
	{
		canvas = new TCanvas("canvas", "canvas", 800, 600);
		etaphi = new TH2D("etaphi", "Phase space of the decay;#eta;#phi", 100, -2.5, 2.5, 100, -3.5, 3.5); 

		fTree->SetBranchAddress("D0_keta", &Keta); 
		fTree->SetBranchAddress("D0_kphi", &Kphi); 
		fTree->SetBranchAddress("D0_pieta", &pieta); 
		fTree->SetBranchAddress("D0_piphi", &piphi); 
		fTree->SetBranchAddress("Dstar_piseta", &piseta); 
		fTree->SetBranchAddress("Dstar_pisphi", &pisphi); 
		fTree->SetBranchAddress("b_tau_pi1eta", &pi1eta); 
		fTree->SetBranchAddress("b_tau_pi1phi", &pi1phi); 
		fTree->SetBranchAddress("b_tau_pi2eta", &pi2eta); 
		fTree->SetBranchAddress("b_tau_pi2phi", &pi2phi); 
		fTree->SetBranchAddress("b_tau_pi3eta", &pi3eta); 
		fTree->SetBranchAddress("b_tau_pi3phi", &pi3phi); 

		fTree->SetBranchAddress("b_B_eta", &Beta); 
		fTree->SetBranchAddress("b_B_phi", &Bphi); 
		fTree->SetBranchAddress("Dstar_eta", &Dstareta); 
		fTree->SetBranchAddress("Dstar_phi", &Dstarphi); 
		fTree->SetBranchAddress("b_tau_eta", &taueta); 
		fTree->SetBranchAddress("b_tau_phi", &tauphi); 

		etaphi->Draw();
		canvas->Draw();

		K = TMarker(0., 0., 20);
		K.SetMarkerColor(kOrange+7);
		pi = TMarker(0., 0., 20);
		pi.SetMarkerColor(kGreen+2);
		pis = TMarker(0., 0., 20);
		pis.SetMarkerColor(kGreen+4);
		pi1 = TMarker(0., 0., 20);
		pi1.SetMarkerColor(kBlue+1);
		pi2 = TMarker(0., 0., 20);
		pi2.SetMarkerColor(kBlue+1);
		pi3 = TMarker(0., 0., 20);
		pi3.SetMarkerColor(kBlue+1);
		B = TMarker(0., 0., 5);
		B.SetMarkerColor(kRed);
		Dstar = TMarker(0., 0., 4);
		Dstar.SetMarkerColor(kOrange+4);
		tau = TMarker(0., 0., 4);
		tau.SetMarkerColor(kRed+3);

		

		// Event string printing
		decaystring = new std::string(); 
		decaystring1 = new std::string(); 
		decaystring2 = new std::string(); 
		decaystring3 = new std::string(); 
		decaystringK = new std::string(); 
		decaystringpi = new std::string(); 
		decaystringspi = new std::string(); 	

		fTree->SetBranchAddress("b_tau_gen1str", &decaystring1); 
		fTree->SetBranchAddress("b_tau_gen2str", &decaystring2); 
		fTree->SetBranchAddress("b_tau_gen3str", &decaystring3); 
		fTree->SetBranchAddress("D0_genkstr", &decaystringK); 
		fTree->SetBranchAddress("D0_genpistr", &decaystringpi); 
		fTree->SetBranchAddress("Dstar_genpistr", &decaystringspi); 	

		fTree->SetBranchAddress("genstring", &decaystring); 	
	}

	~EventLoopFactory() {
		delete decaystring;
		delete decaystring1;
		delete decaystring2;
		delete decaystring3;
		delete decaystringK;
		delete decaystringpi;
		delete decaystringspi;	

		delete canvas; 
	}

	void PrintEvents() // Event loop 
	{	
		int max = static_cast<int>(fTree->GetEntries()); 
		if (fMaxNum > 0) 
		{
			max = min(max, fMaxNum); 
		}
		for (int i=0; i<max; i++) 
		{
			fTree->GetEntry(i);	

			for (unsigned int i=0; i<45; i++) 
			{
				std::cout << "-"; 
			}
			std::cout << endl << endl; 	

			// Start doing stuff per event
			PrintDecayString(); 
			PrintEtaPhiLive(); 
			//plotting = std::thread(&EventLoopFactory::PrintEtaPhiLive, this); 
			//plotting.join(); 

			pause = std::thread(&PauseUntilEnter); 
			pause.join();
		}
		std::cout << std::endl << "Analyzed " << max << " events. " << std::endl; 
	}

	void PrintDecayString() 
	{
		if (fFullDecay) std::cout << *decaystring << std::endl; 	

		std::cout << "K: " << *decaystringK << std::endl; 
		std::cout << "pi: " << *decaystringpi << std::endl; 
		std::cout << "spi: " << *decaystringspi << std::endl; 
		std::cout << "pi1: " << *decaystring1 << std::endl; 
		std::cout << "pi2: " << *decaystring2 << std::endl; 
		std::cout << "pi3: " << *decaystring3 << std::endl; 
	}

	void PrintEtaPhiLive() 
	{
		//std::cout << "Plotting eta phi " << Keta << ", " << Kphi << std::endl;
		//etaphi->Reset();
		//etaphi->Fill(Keta, Kphi); 
		//auto  = new TMarker(Keta, Kphi, 20); 
		//K->SetMarkerColor(kOrange+7);
		K.SetX(Keta);
		K.SetY(Kphi);
		K.Draw(); 
		pi.SetX(pieta);
		pi.SetY(piphi);
		pi.Draw(); 
		pis.SetX(piseta);
		pis.SetY(pisphi);
		pis.Draw(); 
		pi1.SetX(pi1eta);
		pi1.SetY(pi1phi);
		pi1.Draw(); 
		pi2.SetX(pi2eta);
		pi2.SetY(pi2phi);
		pi2.Draw(); 
		pi3.SetX(pi3eta);
		pi3.SetY(pi3phi);
		pi3.Draw(); 
		B.SetX(Beta);
		B.SetY(Bphi);
		B.Draw(); 
		Dstar.SetX(Dstareta);
		Dstar.SetY(Dstarphi);
		Dstar.Draw(); 
		tau.SetX(taueta);
		tau.SetY(tauphi);
		tau.Draw(); 
		canvas->Modified();
    	canvas->Update();
    	gSystem->ProcessEvents();
	}


	void Test() 
	{
		TFile *file = TFile::Open("../../data/v7/BkgDstarDsPart.root", "READ"); 
		fTree = static_cast<TTree*>(file->Get("ntuplizer/tree")); 	

		PrintDecayString(); 	
	

		file->Close(); 
	}


	TTree* fTree = nullptr; 
	bool fFullDecay = false; 
	int fMaxNum = -1; 

	// Variables for decay string printing
	std::string *decaystring = nullptr; 
	std::string *decaystring1 = nullptr; 
	std::string *decaystring2 = nullptr; 
	std::string *decaystring3 = nullptr; 
	std::string *decaystringK = nullptr; 
	std::string *decaystringpi = nullptr; 
	std::string *decaystringspi = nullptr; 	

	// Variables for eta phi printing
	TCanvas *canvas = nullptr; 
	TH2D *etaphi = nullptr; 

	TMarker K; 
	TMarker pi;
	TMarker pis;
	TMarker pi1;
	TMarker pi2;
	TMarker pi3;
	TMarker B;
	TMarker Dstar;
	TMarker tau;

	float Keta; 
	float Kphi; 
	float pieta; 
	float piphi; 
	float piseta; 
	float pisphi; 
	float pi1eta; 
	float pi1phi; 
	float pi2eta; 
	float pi2phi; 
	float pi3eta; 
	float pi3phi; 
	float Beta; 
	float Bphi;
	float Dstareta;
	float Dstarphi;
	float taueta;
	float tauphi;

	std::thread plotting; 
	std::thread pause; 
};

