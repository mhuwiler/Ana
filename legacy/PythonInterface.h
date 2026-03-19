#ifndef PythonInterface_h
#define PythonInterface_h
#include "Python.h"
#include <numpy/arrayobject.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION


class PythonInterface 
{
	public: 
	PyObject *pName, *pModule, *pDict, *pFunc, *pArgs, *python_class, *object, *result, *initFunc;

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
            PyErr_Print(); // TODO: throw an error 
            std::cout << "Cannot instantiate the Python class" << std::endl;
            return;
        }

        PyRun_SimpleString("import sys"); 

        PyRun_SimpleString("print(sys.path)"); 

        PyRun_SimpleString("sys.path = ['/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6', '/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/lib-dynload', '/Users/mhuwiler/Library/Python/3.6/lib/python/site-packages', '/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages']"); 

        PyRun_SimpleString("print(sys.path)"); 

        PyErr_Print();
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

        /*for (int i=0; i<data.size(); i++) 
        {
            std::cout << *(response + i) << ", "; 
            returnvec.push_back(*(response + i)); 
            //response++; 
        }
        std::cout << std::endl; */

        //PyObject* myResult = PyObject_CallMethod(object, "Add2toNumber", "(d)", a); 

	    Py_DECREF (py_array);                             
	    Py_DECREF (pFunc);


	    return returnvec;
	}

    std::vector<float> Evaluate(const vector<float>& data)
    {
        std::vector<double> vec = castVector(data); 
        double *ptr = vec.data(); //const_cast<float*>(data.data());
        npy_intp dims[1] = { static_cast<npy_intp>(vec.size()) };
        PyObject *py_array;

        //std::cout << "Inside Evaluate function " << std::endl; 

        

        py_array = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, ptr);
        

        pArgs = PyTuple_New (1);
        PyTuple_SetItem (pArgs, 0, py_array);

        PyErr_Print();

        //std::cout << "Before function export " << std::endl; 

        pFunc = PyObject_GetAttrString (object, (char*)"Eval"); 

        PyErr_Print();

        //std::cout << "After function export " << std::endl; 

        if (PyCallable_Check (pFunc))
        {
            result = PyObject_CallObject(pFunc, pArgs);
        } 
        else
        {
            cout << "Function is not callable !" << endl;
        }

        double *response = static_cast<double*>(PyArray_DATA((PyArrayObject*)result)); 

        std::vector<float> returnvec; 
        returnvec.reserve(20); 
        std::vector<float> probSig; 
        probSig.reserve(20); 
        std::vector<float> probBkg; 
        probBkg.reserve(20); 

        for (int i=0; i<3*20; i++) 
        {
            //std::cout << *(response + i) << ", "; 
            float value = static_cast<float>(*(response + i)); 
            returnvec.push_back(value); 
            if ((i%3) == 1) 
            {
                probSig.push_back(value); 
            }
            if ((i%3) == 2) 
            {
                probBkg.push_back(value); 
            }
            //response++; 
        }
        //std::cout << std::endl; 

        //PyObject* myResult = PyObject_CallMethod(object, "Add2toNumber", "(d)", a); 

        Py_DECREF (py_array);                             
        Py_DECREF (pFunc);


        return probSig;
    }

    void Initialise(std::string modelname,int givenbatchsize, int givennumpoints) 
    {
        initFunc = PyObject_GetAttrString (object, (char*)"Initialise");

        PyObject *model = PyUnicode_FromStringAndSize(modelname.data(), modelname.size()); 
        PyObject *batchsize = PyLong_FromLong(static_cast<long>(givenbatchsize)); 
        PyObject *numpoints = PyLong_FromLong(static_cast<long>(givennumpoints)); 

        PyObject *localArgs = PyTuple_New (3);
        PyTuple_SetItem (localArgs, 0, model);
        PyTuple_SetItem (localArgs, 1, batchsize);
        PyTuple_SetItem (localArgs, 2, numpoints);

        if (PyCallable_Check (initFunc))
        {
            result = PyObject_CallObject(initFunc, localArgs);
        } 
        else
        {
            cout << "Function is not callable !" << endl;
        }
    }

    std::vector<double> castVector(std::vector<float> vec) 
    {
        std::vector<double> result; 
        result.reserve(vec.size()); 
        for (auto element : vec)
        {
            result.push_back(static_cast<double>(element)); 
        }
        return result; 
    }

};

#endif

