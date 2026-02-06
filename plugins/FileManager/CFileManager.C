#include <iostream>

#include "TChain.h"
#include "TROOT.h"
#include "TSystem.h"
#include "CFileManager.h"



bool FileManager::AddItem(const TString internalRef, const TString fileName, const TString treeName)
{ 
	std::map<TString, std::pair<TString, TString> >::iterator it = fFileCollection.find(internalRef); 
	if (it != fFileCollection.end()) 
	{
		std::cerr << "The internal reference " << it->first << " already exists! It points to object: " << it->second.second << " in file: " << it->second.first  << ". " << std::endl 
		<< "The internal reference for an object needs to be unique! Please provide another one. The object: " << treeName << " in file: " << fileName << " was not added to the FileManager and WILL NOT be opened. " << std::endl; 
		return false; 
	}
	else {
		fFileCollection[internalRef] = std::pair<TString, TString>(fileName, treeName); 
		return true; 
	}
}

bool FileManager::UpdateItem(const TString internalRef, const TString fileName, const TString treeName) 
{
	std::map<TString, std::pair<TString, TString> >::iterator it = fFileCollection.find(internalRef); 
	if (it == fFileCollection.end()) std::cout << "There was no object with internal reference: " << internalRef << ", it will be added as a new object. "; 
	
	// Check if the file is open and close it if so 
	bool updated = true; 
	std::map<TString, TString>::iterator itOpenFiles = fFileUsageBookKeeping.find(internalRef); 
	if (itOpenFiles != fFileUsageBookKeeping.end())
	{
		updated = CloseItem(internalRef); 		
	}
	
	fFileCollection[internalRef] = std::pair<TString, TString>(fileName, treeName); 
	return updated; 
}

TObject* FileManager::GetItem(TString reference, bool autoOpen) 
{
	std::map<TString, TObject*>::iterator it = fObjectCollection.find(reference); 
	if (it != fObjectCollection.end()) 
	{
		return it->second;
	}
	else 
	{
		if (autoOpen) 
		{
			OpenItem(reference); 
			return GetItem(reference, false); 
		}
		else 
		{
			std::cerr << "No object with internal reference: " << reference << std::endl << "Please open it first using OpenFile(internalRef), or provide true as second argment to this function: GetObject(internalRef, true). " << std::endl; 
			return nullptr; 
		}
	}
}

template <class T>
T FileManager::GetItem(TString reference, bool autoOpen) 
{
	return static_cast<T>(GetItem(reference, autoOpen)); 
}

std::string FileManager::GetObject(const TString& internalRef) const 
{
	if (fFileCollection.find(internalRef) != fFileCollection.end()) 
	{
		return fFileCollection.at(internalRef).second.Data(); 
	}
	else 
	{
		std::cerr << "No object with internall reference: " << internalRef << ", returning an empty string. " << std::endl; 
		return ""; 
	}
	
}

std::string FileManager::GetFile(const TString& internalRef) const 
{
	if (fFileCollection.find(internalRef) == fFileCollection.end()) return ""; 
	return fFileCollection.at(internalRef).first.Data(); 
}

bool FileManager::OpenAllItems(TString permission)
{
	bool succeded = true; 
	for (std::map<TString, std::pair<TString, TString> >::iterator itCollection=fFileCollection.begin(); itCollection!=fFileCollection.end(); ++itCollection) 
	{
		if(!OpenFile(itCollection, permission)) 
		{
			succeded = false; 
		} 
	}
	return succeded; 
}

bool FileManager::OpenItem(TString internalRef, TString permission) 
{
	std::map<TString, std::pair<TString, TString> >::iterator itCollection = fFileCollection.find(internalRef); 
	if (itCollection != fFileCollection.end()) 
	{
		return OpenFile(itCollection, permission); 
		
	}
	else 
	{
		std::cerr << "The internal reference " << itCollection->first << " does not exists in the file collection! The file could NOT be opened. " << std::endl << "Please make sure you added it previously using : AddFile(fileName, objectName, internalReference)" << std::endl; 
		return false; 
	}
}

bool FileManager::OpenFile(std::map<TString, std::pair<TString, TString> >::iterator itCollection, TString permission) 
{
	// Test if the permission specified does make sense 
	std::vector<TString> possibleCases = {"READ", "RECREATE", "UPDATE", "NEW", "CREATE", "read", "write", "update", "new", "create"}; 
	// Make an enum maybe 
	TString formattedPermission(permission); 
	formattedPermission.ToUpper(); 
	if (std::find(possibleCases.begin(), possibleCases.end(), formattedPermission) == possibleCases.end()) 
	{
		std::cerr << "ERROR: The permission specified does not exist.\n Pleas specify one of the following: READ, RECREATE, UPDATE, NEW, CREATE " << std::endl; 
	}

	TString currentFileName = itCollection->second.first; 
	gSystem->ExpandPathName(currentFileName); 
	TFile *file = nullptr; 

	// Make this a shared method between this and OpenAllItems 
	std::map<TString, TFile*>::iterator it = fOpenedFiles.find(currentFileName); 
	if (it != fOpenedFiles.end()) 
	{
		// The file is already open
		file = it->second; 
	}
	else 
	{
		// Check for file existence 
		file = TFile::Open(currentFileName, permission); 
		if (file->IsOpen()) 
		{
			fOpenedFiles[currentFileName] = file; 
		}
		else 
		{
			std::cerr << "ERROR: The file could not be opened! " << std::endl; 
			return false; 
		}
	}

	// Adding the file in the bookkeeping under this internal ref 
	fFileUsageBookKeeping[itCollection->first] = currentFileName; 
	TString objName = itCollection->second.second; 
	TObject* obj = file->Get(objName); 
	std::pair<std::map<TString, TObject*>::iterator,bool> alreadythere;
	alreadythere = fObjectCollection.insert (std::pair<TString, TObject*>(itCollection->first, obj));
	if (alreadythere.second==false) 
	{
		std::cerr << "ERROR: The key " << alreadythere.first->first << " already axists with value" <<  alreadythere.first->second << ". " << std::endl << "The value " << obj << " was thus not set. " << std::endl; 
	//return false; 
	}
	return true; 
}

bool FileManager::CloseItem(TString internalRef) 
{
	// Delete ref from fObjectCollection
	// Get filename form fFileUsageBookKeeping
	// Delete ref from fFileUaageBookKeeping
	// Search ing FileUsageBookkeeping  for another ref using same file
	// if not close file
	fObjectCollection.erase(internalRef); 
	std::map<TString, TString>::iterator it = fFileUsageBookKeeping.find(internalRef); 
	TString currentFileName = ""; 
	if (it != fFileUsageBookKeeping.end())
	{
		currentFileName = it->second; 
		
	}
	else 
	{
		std::cerr << "ERROR: There is no opened file in the bookkeping for item: " << internalRef << std::endl; 
		return false; 
	}
	fFileUsageBookKeeping.erase(internalRef); 
	bool isStillUsed = false; 
	for (it = fFileUsageBookKeeping.begin(); it != fFileUsageBookKeeping.end(); ++it )
	{
		if (it->second == currentFileName) isStillUsed = true; 
	}
	if (!isStillUsed)
	{
		std::map<TString, TFile*>::iterator it = fOpenedFiles.find(currentFileName); 
		if (it != fOpenedFiles.end()) 
		{
			delete it->second; 
			fOpenedFiles.erase(it); 
		}
		else 
		{
			std::cerr << "ERROR: There is no opened file found in the bookkeeping for file name: " << currentFileName << std::endl; 
			return false; 
		}
	}
   	return true; 
}

bool FileManager::CloseAll() {
	// Close all files properly 
	std::map<TString, TFile*>::iterator it = fOpenedFiles.begin(); 
	for (it = fOpenedFiles.begin(); it != fOpenedFiles.end(); ++it )
	{
		delete it->second; 
	}
	// Errase all maps regardless of what is open. 
	fObjectCollection.clear(); 
	fOpenedFiles.clear(); 
	fFileUsageBookKeeping.clear(); 
	return true; 
}

bool FileManager::Clear() {
	
	bool closed = CloseAll();
	fFileCollection.clear(); 
	return closed; 
}

// Methods for debug
void FileManager::ListOpenFiles() 
{
	for (std::map<TString, TString>::iterator it = fFileUsageBookKeeping.begin(); it != fFileUsageBookKeeping.end(); ++it) 
	{
		cout << it->first << " " << it->second << std::endl; 
	}
}

void FileManager::ListCollection() 
{
	for (std::map<TString, std::pair<TString, TString> >::iterator it=fFileCollection.begin(); it!=fFileCollection.end(); ++it) 
	{
		cout << it->first << " " << it->second.second << " " << it->second.first << std::endl; 
	}
}

// For early stage debug 
void FileManager::ListObjects() 
{
	for (std::map<TString, TObject*>::iterator it = fObjectCollection.begin(); it != fObjectCollection.end(); ++it) 
	{
		cout << it->first << " " << it->second << std::endl; 
	}
}

void FileManager::ListFiles() 
{
	for (std::map<TString, TFile*>::iterator it = fOpenedFiles.begin(); it != fOpenedFiles.end(); ++it) 
	{
		cout << it->first << " " << it->second << std::endl; 
	}
}




// add a tag for the files 
// get pointer to tree or TObject
// manage opening and closing all or specified file
// add a method List() that lists all objects contained in the FileManager (and their status)
// better map (non ordered) 
// implement clear single item
// implement functions with vector 
// update needs to close the file and reopen it 
