#include <string>
#include <iostream>
#include <map>
#include <vector>

#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TGraph2D.h"



class FileManager {

	public : 
	bool AddItem(const TString internalRef, const TString fileName, const TString treeName); 
	

	bool UpdateItem(const TString internalRef, const TString fileName, const TString treeName); 
	

	TObject* GetItem(TString reference, bool autoOpen=false); 

	template<class T> 
	T GetItem(TString reference, bool autoOpen=false); 

	std::string GetObject(const TString& internalRef) const; 

	std::string GetFile(const TString& internalRef) const; 


	bool OpenAllItems(TString permission="READ"); 

	bool OpenItem(TString internalRef, TString permission = "READ"); 

	bool CloseItem(TString internalRef); 

	bool CloseAll(); 

	bool Clear(); 


	void ListOpenFiles(); 

	void ListCollection(); 

	// For debuging purposes 
	void ListObjects(); 

	void ListFiles(); 



	private : 
	bool OpenFile(std::map<TString, std::pair<TString, TString> >::iterator itCollection, TString permission = "READ"); 

	// Should we use unordered maps? see: https://stackoverflow.com/questions/2196995/is-there-any-advantage-of-using-map-over-unordered-map-in-case-of-trivial-keys
	std::map<TString, std::pair<TString, TString> > fFileCollection; 	// Map of Items, with the file name, object path within the file and the reference as key
	std::map<TString, TFile*> fOpenedFiles; 	// Map of pointers to file with filenames as keys. 
	std::map<TString, TObject*> fObjectCollection; 	// Map fo objects opened in the FileManager 
	std::map<TString, TString> fFileUsageBookKeeping; // Map of FileNames in use by internalRef (key). 



}; 
// add a tag for the files 
// get pointer to tree or TObject
// manage opening and closing all or specified file
