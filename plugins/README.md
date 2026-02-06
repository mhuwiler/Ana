# Plugins for ROOT 

This reporitory hosts classes and modules designed to enhance fucntionalities of the ROOT data analysis framework and address some of its shortcomings. The ROOT framework offers some very good and efficient core functionalities, yet their use by analysers is sometimes not straightforward and can require a significant amount of code. Classes from this repository try to offer an easier and smoother user interface, e.g. for plotting graphs, handling files, etc., making use of the powerful ROOT functionalities in the background. They can easily be integrated into the ROOT environment by sourcing the plugins.C script (adding it to rootlogon.c), and then be called within the interpreter/scripts in the same way as ROOT classes, both in C++ and Python. They can also be integrated as standalone classes into any project, by picking the desired header and source files. 

The documentation of each submodule is found in the corresponding folder. 

[Install instructions](INSTALL.md)

## Repository organisation 

The repository is organised into subfolders, each of them containing all the code of a sumbodule that can be used independently. To use the all plugins at once, including the script plugins.C is all what is needed. The script plugins.C contains some static functions that are not yet ready to be a standalone submodule, and sources in addition all other submodules. Including this script into some code gives thus access to all plugins at once. Using a couple of submodules only is of course possible, by sourcing the corresponding files from within the subfolders. 

## Coding conventions 

For this repo, the coding conventions are inspired by (1). As these modules are plugins for ROOT, the naming conventions of ROOT are followed, and superseed the naming conventions of the previously cited source, with the difference that class names do not need to start with a letter T. Foor indentation, spacing and position of operand, bracket, parenthesis etc. the Allman Style conventions are entirely followed and superseed the conventions found in the ROOT framework. 

The main points of the coding conventions are highlighted below : 

	- Each class name starts with a capital letter 
	- Methods start with a capital letter 
	- CapitalEachWord (CamelCase) 
	- The opening brackets are put to the next line: 
		if (something) 		
		{
			some code; 
		}
	- One space before and after operators, e.g. int x = y
	- Colons are followed, but not preceeded by a space 
	- Getter and Setter methods are declared inline in the header file 
	- One class per file 
	- Same name for file as for class


###souces: 

1. LLVM coding conventions: https://llvm.org/docs/CodingStandards.html#the-high-level-issues
2. ROOT coding convetions: https://root.cern.ch/coding-conventions#Naming_conventions
3. Microsoft coding conventions: https://docs.microsoft.com/en-us/dotnet/csharp/programming-guide/inside-a-program/coding-conventions
4. Allman Style: https://en.wikipedia.org/wiki/Indentation_style#Allman_style




## Unit testing

Since this code may be largely used, unit testing is required for each module. The test should probe all implemented features, and ideally cover all possible states of the program. A nice guide to writing unit tests can be found here: https://www.toptal.com/qa/how-to-write-testable-code-and-why-it-matters All code that is committed or merged to the master should pass the unit tests. 



