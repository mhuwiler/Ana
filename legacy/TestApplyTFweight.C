#include "ROOT/RDataFrame.hxx"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TChain.h"
#include "/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"
ClassImp(FileManager)
#include "TLorentzVector.h"
#include "TGraph2D.h"
#include "TH2D.h"
#include "TLegend.h"
#include <iostream>
#include "DrawTMVAHistogram.C"
#include "GetSeparation.C"
#include "Python.h"
#include "TPython.h"
#include <numpy/arrayobject.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


using namespace ROOT; 


// To run this macro, the python class TFEvaluation.py needs to be loaded into root prior to execution. e.g.:
// root -e 'TPython::LoadMacro("TFEvaluation.py");' ApplyTFweight.C


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
	#include <chrono>
	#include <thread>
	//std::this_thread::sleep_for(static_cast<std::chrono::seconds>(timeInSec));
	sleep(timeInSec); 
}

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

class PythonInterface 
{
	public: 
	PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *python_class, *object, *result;

	PythonInterface(const std::string& moduleName) 
	{
		setenv("PYTHONPATH",".",1);
    	Py_Initialize ();
    	pName = PyUnicode_FromString ((char*)moduleName.data());

    	pModule = PyImport_Import(pName);

    	//pDict = PyModule_GetDict(pModule);

    	import_array ();                   

        if (pModule == nullptr) 
        {
            PyErr_Print();
            std::cerr << "Fails to import the module.\n";
            return;
        }

        // dict is a borrowed reference.
        pDict = PyModule_GetDict(pModule);
        if (pDict == nullptr) 
        {
            PyErr_Print();
            std::cerr << "Fails to get the dictionary.\n";
            return;
        }

        // Builds the name of a callable class
        python_class = PyDict_GetItemString(pDict, (char*)moduleName.data());
        if (python_class == nullptr) 
        {
            PyErr_Print();
            std::cerr << "Fails to get the Python class.\n";
            return;
        }

        // Creates an instance of the class
        if (PyCallable_Check(python_class)) 
        {
            object = PyObject_CallObject(python_class, nullptr);
            Py_DECREF(python_class);
        } 
        else 
        {
            std::cout << "Cannot instantiate the Python class" << std::endl;
            return;
        }
	}

	~PythonInterface() 
	{
		Py_DECREF(pName);                
    	Py_DECREF (pModule);
    	Py_DECREF (pDict);
        Py_DECREF(python_class);

    	Py_Finalize ();    
	}

	std::vector<double> EvaluateArray(const vector<double>& data)
	{
	    double *ptr = const_cast<double*>(data.data());
	    npy_intp dims[1] = { static_cast<npy_intp>(data.size()) };
	    PyObject *py_array;

	    

	    py_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, ptr);
	    

	    pArgs = PyTuple_New (1);
	    PyTuple_SetItem (pArgs, 0, py_array);

	    pFunc = PyObject_GetAttrString (object, (char*)"pyArray"); 

	    if (PyCallable_Check (pFunc))
	    {
	        result = PyObject_CallObject(pFunc, pArgs);
	    } 
        else
	    {
	        cout << "Function is not callable !" << endl;
	    }

        double *response = static_cast<double*>(PyArray_DATA((PyArrayObject*)result)); 

        std::vector<double> returnvec; 
        returnvec.reserve(data.size()); 

        for (int i=0; i<data.size(); i++) 
        {
            std::cout << *(response + i) << ", "; 
            returnvec.push_back(*(response + i)); 
            //response++; 
        }
        std::cout << std::endl; 

        //PyObject* myResult = PyObject_CallMethod(object, "Add2toNumber", "(d)", a); 

	    Py_DECREF (py_array);                             
	    Py_DECREF (pFunc);


	    return returnvec;
	}

};
 

void TestApplyTFweight(TString campaignName = "ApplyTFweight/") 
{
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C"); 
	//gROOT->LoadMacro("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.C+");
	//gSystem->Load("/Users/mhuwiler/coding/plugins/FileManager/CFileManager.so"); 

	FileManager filemanager; 


	filemanager.AddItem("prodlatest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/flatTupleDataLatest.root", "ntuplizer/tree"); 
	filemanager.AddItem("MCgenmatched", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/flatTupleGenmatchedAllSingleTauLatest.root", "ntuplizer/tree"); 

	filemanager.AddItem("MCOfficialSample", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/all.root", "ntuplizer/tree"); 
	filemanager.AddItem("DstarDsMCfirst", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/bkgDstarDsFirst.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataWS", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataWS.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataLatest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataLast.root", "ntuplizer/tree"); 
	filemanager.AddItem("MCofficial", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/MCFirstSubmission.root", "ntuplizer/tree"); 
	filemanager.AddItem("DataLarge", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge.root", "ntuplizer/tree"); 

	// For ABCD estimation
	filemanager.AddItem("DataLargeMVAfirst", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_mva.root", "tree"); 
	filemanager.AddItem("DataLargeMVA", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_mvaxgb.root", "tree"); 
	filemanager.AddItem("DataLargeMVASimple", "/eos/home-m/mhuwiler/data/Analysis/v9/DataVeryLarge_converted_mvaxgbsimple.root", "tree"); 
	filemanager.AddItem("MCvalidation", "/eos/home-m/mhuwiler/data/Analysis/v9/TauCutflowSample_mva.root", "tree"); 
	filemanager.AddItem("signalTrain", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "signalTrain"); 
	filemanager.AddItem("signalTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "signalTest"); 
	filemanager.AddItem("backgroundTrain", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "backgroundTrain"); 
	filemanager.AddItem("backgroundTest", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/trainingLarge/plots.root", "backgroundTest"); 
	filemanager.AddItem("DataBackground", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/DataVeryLarge_converted.root", "tree"); 
	filemanager.AddItem("MCSignal", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/PrivateProductionGenDstar_converted.root", "tree"); 
	filemanager.AddItem("SignalOfficialMC50M", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/SignalOfficialMC50M.root", "ntuplizer/tree"); 
	filemanager.AddItem("ParkingBPHAllRun2018B", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/ParkingBPHRun2018B_converted.root", "tree"); 
    filemanager.AddItem("MCSignalMMultipleTau", "/eos/home-m/mhuwiler/DoctoralThesis/Analysis/data/firstAllTau.root", "ntuplizer/tree"); 



	filemanager.OpenAllItems(); 

	gStyle->SetOptStat(0); 


	auto dataframe = RDataFrame(*filemanager.GetItem<TTree*>("MCSignalMMultipleTau")); // tree100k
	

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

	//TFInference TFmodel; //TFEvaluation TFmodel; 

	// Import the python module

	PyObject *dict, *python_class, *object;

	Py_Initialize();
	import_array(); 

	PyObject* myModuleString = PyUnicode_FromString("TFinference");

	PyObject* myModule = PyImport_Import(myModuleString);

	if (myModule == nullptr) 
	{
    	PyErr_Print();
    	std::cerr << "Fails to import the module.\n";
    	return;
  	}
  	Py_DECREF(myModuleString);

  	// dict is a borrowed reference.
  	dict = PyModule_GetDict(myModule);
  	if (dict == nullptr) 
  	{
    	PyErr_Print();
    	std::cerr << "Fails to get the dictionary.\n";
    	return;
  	}
  	Py_DECREF(myModule);

  	// Builds the name of a callable class
  	python_class = PyDict_GetItemString(dict, "TFEvaluation");
  	if (python_class == nullptr) 
  	{
    	PyErr_Print();
    	std::cerr << "Fails to get the Python class.\n";
    	return;
  	}
  	Py_DECREF(dict);

  	// Creates an instance of the class
  	if (PyCallable_Check(python_class)) 
  	{
    	object = PyObject_CallObject(python_class, nullptr);
    	Py_DECREF(python_class);
  	} 
  	else 
  	{
    	std::cout << "Cannot instantiate the Python class" << std::endl;
    	Py_DECREF(python_class);
    	return;
  	}

	//PyObject* args = PyTuple_Pack(1,PyFloat_FromDouble(2.0));

	double a[20];

	PyObject* myResult = PyObject_CallMethod(object, "Add2toNumber", "(d)", a); 

	double result = PyFloat_AsDouble(myResult);

	std::cout << "Result from python: " << result << std::endl; 

	std::vector<double> datavec = {1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2, 2.3}; 

	double *data = datavec.data(); 

    std::cout << "Before making class" << std::endl; 

	PythonInterface pyEvaluation("TestPyTFEval"); 

    std::cout << "After making class" << std::endl; 

	auto response = pyEvaluation.EvaluateArray(datavec);

    std::cout << "After evaluation" << std::endl; 

    for (auto element : response) 
    {
        std::cout << element << ", "; 
    }
    std::cout << std::endl; 

	//ClearPython();  

	double *mydata = data; 

	std::cout << "New function ran" << std::endl; 

	npy_intp dims[1]; 
	dims[0] = 20; 

	PyObject* npdata = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, mydata); 

	PyObject* args = PyTuple_New (1);

	PyTuple_SetItem(args, 1, npdata);

	PyObject* pFunc = PyDict_GetItemString(dict, "DisplayArray"); 

	//PyObject* myResult = PyObject_CallMethod(object, "DisplayArray", "NPY_DOUBLE", npdata); 

	PyObject* newResult = PyObject_CallObject(pFunc, args); 

	std::cout << "After first call" << std::endl; 

	auto TFresponse = [&object](std::vector<float> pt, std::vector<float> eta, std::vector<float> phi, std::vector<float> q, std::vector<float> DOCA2D, std::vector<float> DOCA2DErr, std::vector<float> DOCA3D, std::vector<float> DOCA3DErr, std::vector<float> dzToPV, std::vector<float> dzToClosest, std::vector<float> isAssociate, std::vector<float> assocQualityToPV, std::vector<int> genmatch) 
	{
		std::vector<float> response; 
		std::vector<double> initial = {1., 2., 3.}; 
		TArrayD array(initial.size(), initial.data()); 

		int size = pt.size(); 

		double data[20] = {1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2, 2.3, 3.4, 4.5, 5.6, 6.7, 7.8, 8.9, 9.1, 1.2, 2.3}; 

		npy_intp dims[1]; 
		dims[0] = 20; 

		PyObject* npdata = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data); 

		PyObject* myResult = PyObject_CallMethod(object, "DisplayArray", "NPY_DOUBLE", npdata); 

		//PyObject* pypt= TPython::CPPInstance_FromVoidPtr(&pt, "std::vector< std::vector<float> >"); // Declaring to python what type of object it is
		std::cout << "After assignment" << std::endl; 
		//PyObject* pymain = PyImport_ImportModule("main");
		//PyModule_AddObject(pymain, "pt", pypt);

		//TPython::Prompt(); 

		//std::vector<std::vector<float> > input = {eta, phi, pt, q, DOCA2D, DOCA2DErr, DOCA3D, DOCA3DErr, dzToPV, dzToClosest, isAssociate, assocQualityToPV}; 

		//TFmodel.OtherTest();
		//TFmodel.Eval(size, eta.data(), phi.data(), pt.data(), q.data(), DOCA2D.data(), DOCA2DErr.data(), DOCA3D.data(), DOCA3DErr.data(), dzToPV.data(), dzToClosest.data(), isAssociate.data(), assocQualityToPV.data())
		//TFmodel.Eval(pt);  
		//TFmodel.Test(array); 
		//TFmodel.Condition(-1); 

		response.push_back(-999.); 
		return response; 
	};

	auto withWeight = dataframe.Define("TFscore", TFresponse, {"track_pt", "track_eta", "track_phi", "track_charge", "track_doca2D", "track_doca2Derror", "track_doca", "track_docaerror", "track_dzToPV", "track_dzToClosestVertex", "track_isAssociatedToPV", "track_pvAssociationQuality", "track_isgenmatched"}); 

	withWeight.Snapshot("ntuplizer/tree", "SignalOfficialMC50M_withTHweight.root"); 

	//Pause(5); 

	//PauseUntilEnter(); //system("pause"); 

	filemanager.CloseAll(); 


}

