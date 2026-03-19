#include <Python.h>
#include <stdio.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>

float Array [] = {1.2, 3.4, 5.6, 7.8};

// taken from: https://stackoverflow.com/questions/30388170/sending-a-c-array-to-python-and-back-extending-c-with-numpy

//int main (int argc, char *argv[])
int TestPyInclude()
{
    float *ptr = Array;
    PyObject *pName, *pModule, *pDict, *pFunc, *pArgs;
    npy_intp dims[1] = { 4 };
    PyObject *py_array;

    setenv("PYTHONPATH",".",1);
    Py_Initialize ();
    pName = PyUnicode_FromString ("TestPyInclude");

    pModule = PyImport_Import(pName);

    pDict = PyModule_GetDict(pModule);

    import_array ();                                   

    py_array = PyArray_SimpleNewFromData(1, dims, NPY_FLOAT, ptr);
    

    pArgs = PyTuple_New (1);
    PyTuple_SetItem (pArgs, 0, py_array);

    pFunc = PyDict_GetItemString (pDict, (char*)"pyArray"); 

    if (PyCallable_Check (pFunc))
    {
        PyObject_CallObject(pFunc, pArgs);
    } else
    {
        cout << "Function is not callable !" << endl;
    }

    Py_DECREF(pName);
    Py_DECREF (py_array);                             
    Py_DECREF (pModule);
    Py_DECREF (pDict);
    Py_DECREF (pFunc);

    Py_Finalize ();                                    

    return 0;
}