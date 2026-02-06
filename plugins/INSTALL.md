# Installing the plugins for ROOT 


## Downloading the framework 

To install the plugins for ROOT, the easiest (and only) way up to now is to clone the repository from https://gitlab.cern.ch/mhuwiler/plugins/tree/master

## Setting up the plugins to be used in interactive mode

The current way to use the plugins in any interactive session is a bit hacky for now. All pllugins are included (the .C file which itself includes the .h file) in the file libFunctions.C, which furthermore contains some additional functions that are not yet a plugin by themselves. So it is sufficient to include libFunctions.C in any script that wants to use one of the plugins. Iy you want to use the plugins in the interactive ROOT prompt, the libFunctions.C file needs to be loaded: 

	.L path/to/plugins/libFunctions.C

If you want to make ROOT compile the plugins into a library, just type: 

	.L path/to/plugins/libFunctions.C+

This can also be made automatic by creating a file named rootlogon.C containing the code: 

	//rootlogon.C
	{
        gROOT->LoadMacro("path/to/plugins/libFunctions.C+");
	}

And placing it either in your $HOME folder if you want the plugins to be loaded for every version of ROOT you have installed, regardless of the folder you are in, in $ROOTSYS/etc as system.rootlogon.C if you want them to be loaded for the version of ROOT in $ROOTSYS only, or place it in the project folder you want the plugins to be loaded when ROOT is started from this folder. 


## Using the classes in built C++ code 

Once the plugins folder is downloaded, each class can also be used independently in any C++ code based on ROOT, using the standart procedure of including the header files and building the source files. 

## Using the ROOT plugins in python 

Since the plugins can be loaded into the ROOT interpreter, they can also be used in a straightforward way in PyROOT. One just needs to load the libFunctions.C file within the python script or interpreter: 

	import ROOT
	ROOT.gROOT.LoadMacro("path/to/plugins/libFunctions.C+")

And then use them as any other ROOT class: 

	filemanager = ROOT.FileManager()

The plugins are of course also very suitable for use in Jupyter notebooks. The import process is the same as for python. 

