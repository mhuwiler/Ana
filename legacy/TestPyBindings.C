#include "ROOT/RDataFrame.hxx"
#include "TPython.h"
#include "TPyArg.h"
#include "TPyClassGenerator.h"


//class MyPyClass; 

void TestPyBindings() 
{
	// Load the python class
	TPython::LoadMacro("MyPyClass.py"); 
	//TPython::Import("./MyPyClass"); 
	
	MyPyClass pyClass; 


}

