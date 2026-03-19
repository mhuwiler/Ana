#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TRandom3.h>
#include <TSystem.h>
#include <TPRegexp.h>
#include <TMVA/Reader.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>


//namespace pt = property_tree; 


// Local set of params and thresholds 
constexpr Int_t maxNumberOfBackup = 50; 
constexpr Int_t autoSaveInterval = -100000; 


struct TMVAWeightfileVariables 
{
    std::vector<std::pair<TString, TString> > variables; 

    void load(const std::string& filename) 
    {
        boost::property_tree::ptree propertyTree; 

        boost::property_tree::read_xml(filename, propertyTree); 

        std::cout << "File name: " << filename << std::endl; 


        for (auto variable : propertyTree.get_child("MethodSetup.Variables")) 
        {
            if (variable.first == "<xmlattr>") continue; 
            // Loading the xml attributes Expression (variable name) and Type into the vector of pairs 
            variables.push_back(std::make_pair<TString, TString>(variable.second.get("<xmlattr>.Label", "").data(), variable.second.get("<xmlattr>.Type", "kNoType").data())); 
        }


    }
};


bool FileNameAlreadyExists(const TString name, bool useXRD=false) 
{
    bool answer = false; //TODO: decide appropriate default 

    if (not useXRD) 
    {
        TString fullName = name; 
        gSystem->ExpandPathName(fullName); 

        FileStat_t dummy; 
        answer = !gSystem->GetPathInfo(fullName, dummy); 
    }
    else 
    {
        // Implement XRootD handling 
    }

    return answer; 
}

bool CreateBackupCopy(const TString filename) 
{
    TString backupName = TString(filename).ReplaceAll(".root", "_bkp.root"); 
    gSystem->ExpandPathName(backupName); 
    // Check if the backup name already exists
    Int_t count = 2; 
    while (FileNameAlreadyExists(backupName) and (count < maxNumberOfBackup)) 
    {
        TPRegexp pattern("_bkp.*.root"); 
        backupName(pattern) = TString::Format("_bkp%i.root", count); 
        //backupName.ReplaceAll("_bkp*.root", TString::Format("_bkp%i.root", count)); 
        count++; 
    }
    TString modifiableSourceName = filename; 
    gSystem->ExpandPathName(modifiableSourceName); 
    if (count >= maxNumberOfBackup) std::cout << "(CreateBackupCopy) ERROR: Maximum number of backups (" << maxNumberOfBackup 
                                        << ") attained, you need to delete the previous ones, or enable overwrite mode! " << std::endl; 


    Int_t result = gSystem->CopyFile(modifiableSourceName, backupName, false); 
    assert(result != -2); // The file should not already exist 
    if (result == -1) std::cerr << "(CreateBackupCopy) ERROR: The file to backup could not be opened: " << filename << std::endl; 
    if (count > (maxNumberOfBackup - 10)) std::cout << "(CreateBackupCopy) Warning: Created more than " << maxNumberOfBackup - 10 << " backups, please consider deleting them... " << std::endl; 

    return (result==0); 
}

TObjArray GetDirInventory(const TDirectory& path) 
{
    TObjArray listOfObjects; // TODO: change to vector<std::pair> with path and object 

    TIter next(path.GetListOfKeys());
    TKey *key=0; 

    while((key = (TKey*)next())) 
    {
        TObject *currentItem = const_cast<TDirectory&>(path).Get(key->GetName()); // TODO: fix this dirty hack
        if (key->ReadObj()->InheritsFrom(TDirectory::Class())) //key->IsFolder() Trees are also folders ... 
        {
            TObjArray subFolderContent = GetDirInventory(*static_cast<TDirectory*>(currentItem)); // TODO: try changing to TKey 
            // Fill back the current collection with the result
            for (Int_t i=0; i<subFolderContent.GetEntries(); i++) 
            {
                listOfObjects.Add(subFolderContent.At(i)); 
            }
        }
        else 
        {
            listOfObjects.Add(new TObjString(TString(key->GetMotherDir()->GetPath())+"/"+TString(key->GetName()))); 
        }
    }
    return listOfObjects; 
}

TString RemovePrefix(const TString& initialString, TString prefix) 
{
    TString finalString = initialString; 
    Size_t length = prefix.Length(); 
    if (initialString.BeginsWith(prefix))
    {
        finalString = finalString(length, initialString.Length() - length); 
    }
    return finalString; 
}

TString RemoveSuffix(TString someString, const TString& suffix) 
{
    if (someString.EndsWith(suffix)) 
    {
        Size_t length = someString.Length() - suffix.Length(); 
        someString = someString(0, length); 
    }
    return someString; 
}

std::string TStringToStdString(const TString& initialString) // TODO: remove 
{
    Size_t len = initialString.Length(); 
    std::string finalString(initialString); 
    std::cout << "Initial String: " << initialString << " of size: " << len << std::endl; 
    std::cout << "Converted file name: " << finalString << std::endl; 
    return finalString; 
}

void FillTrainTestTrees(const Float_t fraction, TTree *inputTree, TDirectory *destination, const Float_t testFraction = -1.) 
{
    assert(fraction < 1.); 
    assert(fraction > 0.); 
    Int_t factor = static_cast<Int_t>(round(1./fraction)); //static_cast<Float_t>(inputTree->GetEntries())*fraction
    std::cout << "Taking every " << factor << " event out of " << inputTree->GetEntries() << std::endl; 
    TString name = inputTree->GetName(); 

    destination->cd(); 
   
	// Creating the train and test trees for signal and background. 
	TTree *train = inputTree->CloneTree(0);
    train->SetName(name+"_train");
	TTree *test = inputTree->CloneTree(0);
    TTree *validation = inputTree->CloneTree(0);
    
    test->SetName(name+"_test");
    validation->SetName(name+"_validation");

	if (autoSaveInterval) 
    {
        train->SetAutoSave(autoSaveInterval); 
        test->SetAutoSave(autoSaveInterval); 
        validation->SetAutoSave(autoSaveInterval); 
    }
    
    
    // Event loop
	for (unsigned int i=0; i<inputTree->GetEntries(); i++ ) 
	{
        if ( i%1000000 == 0 ) cout << "Processing entry " << i << " ... " << std::endl;
        
        //inputTree->GetEntry(i);
        if ( i%factor == 0 )
		{
            inputTree->GetEntry(i);
			train->Fill();
        }
        
        else if ( (i%factor)-1 == 0 )
        {
            inputTree->GetEntry(i);
            test->Fill();
        }
        
        else
        {
            inputTree->GetEntry(i);
            validation->Fill();
        }
        
	}



}

void AddTrainTestTag(const Int_t factor, TTree *tree, const Int_t rangeMax = 10, const Int_t rangeMin = 0, const std::string varName = "trainTestTag", const Int_t randomSeed = 5000) 
{
    // Adding a random number for each entry, between 0 and rangeMax
    
    Int_t myVar;
    TBranch *b_myVar= tree->Branch(varName.c_str(), &myVar, (varName+"/I").c_str()); 

    TRandom3 randomGenerator(randomSeed);   // Set random seed to 5000 for reproducibility 

    Long64_t nentries = tree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) 
    {
        tree->GetEntry(i);

        myVar= static_cast<Int_t>(randomGenerator.Integer(rangeMax));

        b_myVar->Fill();      
   }

}

void AddTrainTestTag(const Int_t factor, const TTree *inputTree, TDirectory *directory, const Int_t rangeMax = 10, const Int_t rangeMin=0, const std::string varName="trainTestTag", const Int_t randomSeed=5000) 
{
    // Adding a random number for each entry, between 0 and rangeMax

    TRandom3 randomGenerator(randomSeed);   // Set random seed to 5000 for reproducibility 

    TTree *readTree = const_cast<TTree*>(inputTree); //TODO: Find a way around this dirty hack 

    //readTree->SetBranchStatus("*", 1); 

    directory->cd(); 

    TTree *outputTree = readTree->CloneTree(0); 

    Int_t myVar;
    TBranch *b_myVar= outputTree->Branch(varName.c_str(), &myVar, (varName+"/I").c_str()); 

    Long64_t nentries = readTree->GetEntries();
    for (Long64_t i=0; i<nentries; i++) 
    {
        readTree->GetEntry(i);

        myVar= static_cast<Int_t>(randomGenerator.Integer(rangeMax));

        //b_myVar->Fill();  
        outputTree->Fill();     
   }

}

void AddWeightFile(const TString& weightFile, TTree *tree, TString branchName) 
{
    tree->SetBranchStatus(branchName, 0); 

    std::unique_ptr<TMVA::Reader> reader = std::make_unique<TMVA::Reader>("!Color:Silent"); 

    // Parse actually the xml to get the variable names 
    TMVAWeightfileVariables weightFileVars ; 

    weightFileVars.load(TStringToStdString(weightFile)); 

    Int_t numI = 0, numF = 0; 

    for (auto entry : weightFileVars.variables) 
    {
        std::cout << "Variable: " << entry.first << " of type: " << entry.second << std::endl; 
        TString type = entry.second; 

        if (type == "F") // The variable is of float type 
        {
            numF++; 
        }
        else if (type == "I") // The variable is of int type 
        {
            numI++; 
        }
        else // The type is not (yet?) supported 
        {
            std::cerr << "(AddWeighFile) ERROR: Unsupported type'" << type << "' for variable: " << entry.first << std::endl; 
        }
    }

    // Using C type arrays for performance 
    Float_t floatVars[numF]; 
    Int_t intVars[numI]; 
    Float_t intToFloat[numI]; // TMVA reader only accepts floats 

    TString floatNames[numF]; 
    TString intNames[numI]; 

    Int_t countF = 0, countI = 0; 
    std::vector<TString> unknownTypes; 

    std::vector<float> *vfloatVar = 0; 

    for (auto entry : weightFileVars.variables) 
    {
        TString varName = entry.first; 
        TString type = entry.second; 

        if (type == "F") // The variable is of float type 
        {
            assert(countF < numF); 
            bool result = tree->SetBranchAddress(varName, &floatVars[countF]); 
            std::cout << "Set branch result: " << result << std::endl; 
            if (result) result = tree->SetBranchAddress(varName, &vfloatVar); 
            std::cout << "Set branch result: " << result << std::endl; 
            floatNames[countF] = varName; 
            countF++; 
        }
        else if (type == "I") // The variable is of int type 
        {
            assert(countI < numI); 
            tree->SetBranchAddress(varName, &intVars[countI]); 
            intNames[countI] = varName; 
            countI++; 
        }
        else // The type is not (yet?) supported 
        {
            unknownTypes.push_back(type); 
        }
    }

    if (!unknownTypes.empty()) 
    {
        std::cout << "(AddWeighFile) ERROR: Variables with the following unsupported types were read from the weight file: "; 
        //std::copy(unknownTypes.begin(), unknownTypes.end(), std::ostream_iterator<TString>(std::cout, ", ")); 
        for (auto it = unknownTypes.begin(); it < unknownTypes.end(); it++) 
        {
            std::cout << it->Data(); 
            if (!(std::next(it) == unknownTypes.end())) std::cout << ", "; 
        }
        std::cout << endl; 
    } 

    for (int i=0; i < numF; i++) 
    {
        reader->AddVariable(floatNames[i], &floatVars[i]); 
    }
    for (int i=0; i < numI; i++) 
    {
        reader->AddVariable(intNames[i], &intToFloat[i]); 
    }

    TTree *newTree = tree->CloneTree(); 

    Float_t weightVar = 0.; 
    TBranch *branch = newTree->Branch(branchName, &weightVar, branchName+"/F"); 

    for (Long64_t i=0; i<tree->GetEntries(); i++)
    {
        tree->GetEntry(i); 

        // Cast int variables to float 
        for (int idx=0; idx < numI; idx++) 
        {
            intToFloat[idx] = static_cast<Float_t>(floatVars[idx]); 
        }

        weightVar = reader->EvaluateMVA("BDT::BDT"); // TODO: get from weight file? 

        newTree->Fill(); 
    }


}

void TrainTestSplit(const Float_t factor, const TString infile, const TString trees, TString copyMethod, TString outfilename="", TString outtrees="") 
{
    // First get the list of comma separated trees to treat
    TObjArray *treeNames = trees.Tokenize(","); 

    // Create a backup copy of the file 
    CreateBackupCopy(infile); 

    TFile *inFile = TFile::Open(infile.Data(), "READ"); 

    if (!inFile->IsOpen()) 
    {
        std::cout << "(TrainTestSplit) ERROR: The input file " << infile << " could not be opened. " << std::endl;  
        return; 
    }

    // Getting the object content of the file
    auto list = GetDirInventory(*static_cast<TDirectory*>(inFile)); 

    if (copyMethod == "newDiff") 
    {
        // Here we create another outputfile 
        TFile *outfile = TFile::Open(TString(infile).ReplaceAll(".root", "_split.root"), "RECREATE"); 

        for (Int_t i=0; i<list.GetEntries(); i++) 
        {
            const TString currentObject = static_cast<TObjString*>(list.At(i))->GetString(); 

            const TObjArray* paths = currentObject.Tokenize(":");  

            assert(paths->GetEntries() == 2); 

            TString relPath = static_cast<TObjString*>(paths->At(1))->GetString(); 

            relPath = RemovePrefix(relPath, "/"); 

            std::cout << "Rel path: " << relPath << std::endl; 

            TObject *found = treeNames->FindObject(relPath); 
            if (found) 
            {
                TString path = static_cast<TObjString*>(found)->GetString(); 
                std::cout << "Found path: " << path << std::endl; 

                TTree *inTree = static_cast<TTree*>(inFile->Get(path)); 

                if (not inTree) 
                {
                    std::cout << "(TrainTestSplit) ERROR: No input tree! Make sure the object exists: " << path << std::endl; 
                    continue; // TODO: Make sure this should not throw an error 
                }

                TString directory = static_cast<TObjString*>(found)->GetString(); 

                TObjArray *folders = directory.Tokenize("/"); 

                TString treeName = static_cast<TObjString*>(folders->At(folders->GetLast()))->GetString(); 

                std::cout << "Tree name: " << treeName << std::endl; 

                delete folders; 

                directory = RemoveSuffix(directory, "/"+treeName); 


                TDirectory *outDir = outfile->mkdir(directory); 
                if (!outDir) 
                {
                    outDir = static_cast<TDirectory*>(outfile->Get(directory)); 
                }
                outDir->cd(); 
                std::cout << "Made directories" << std::endl; 

                //TTree *newTree = inTree->CloneTree(0); 

                FillTrainTestTrees(factor, inTree, outDir); 

                //if (newTree) 
                //{
                //    newTree->Write();    // Check for empty tree (whole data cut away)
                //}
                //else 
                //{
                //    std::cout << "(TrainTestSplit) ERROR: no tree to write. " << std::cout; 
                //}
            }
            //const TString treeName = static_cast<TObjString*>(treeNames->At(i))->GetString(); 
            //TTree *inTree = static_cast<TTree*>(inFile->Get(treeName)); 

            



            // Recreate the tree structure

            delete paths; 
        }

        outfile->Write(); 

        outfile->Close(); 
    }

    if (copyMethod == "newTag") 
    {
        // Here we create another outputfile 
        TFile *outfile = TFile::Open(TString(infile).ReplaceAll(".root", "_tag.root"), "RECREATE"); 

        for (Int_t i=0; i<list.GetEntries(); i++) 
        {
            const TString currentObject = static_cast<TObjString*>(list.At(i))->GetString(); 

            const TObjArray* paths = currentObject.Tokenize(":");  

            assert(paths->GetEntries() == 2); 

            TString relPath = static_cast<TObjString*>(paths->At(1))->GetString(); 

            relPath = RemovePrefix(relPath, "/"); 

            std::cout << "Rel path: " << relPath << std::endl; 

            TObject *found = treeNames->FindObject(relPath); 
            if (found) 
            {
                TString path = static_cast<TObjString*>(found)->GetString(); 
                std::cout << "Found path: " << path << std::endl; 

                TTree *inTree = static_cast<TTree*>(inFile->Get(path)); 

                if (not inTree) 
                {
                    std::cout << "(TrainTestSplit) ERROR: No input tree! Make sure the object exists: " << path << std::endl; 
                    continue; // TODO: Make sure this should not throw an error 
                }

                TString directory = static_cast<TObjString*>(found)->GetString(); 

                TObjArray *folders = directory.Tokenize("/"); 

                TString treeName = static_cast<TObjString*>(folders->At(folders->GetLast()))->GetString(); 

                std::cout << "Tree name: " << treeName << std::endl; 

                delete folders; 

                directory = RemoveSuffix(directory, "/"+treeName); 


                TDirectory *outDir = outfile->mkdir(directory); 
                if (!outDir) 
                {
                    outDir = static_cast<TDirectory*>(outfile->Get(directory)); 
                }
                outDir->cd(); 
                std::cout << "Made directories" << std::endl; 

                //TTree *newTree = inTree->CloneTree(0); 

                AddTrainTestTag(factor, inTree, outDir); 

                //if (newTree) 
                //{
                //    newTree->Write();    // Check for empty tree (whole data cut away)
                //}
                //else 
                //{
                //    std::cout << "(TrainTestSplit) ERROR: no tree to write. " << std::cout; 
                //}
            }
            //const TString treeName = static_cast<TObjString*>(treeNames->At(i))->GetString(); 
            //TTree *inTree = static_cast<TTree*>(inFile->Get(treeName)); 

            



            // Recreate the tree structure

            delete paths; 
        }

        outfile->Write(); 

        outfile->Close(); 
    }

    if (copyMethod == "inplaceCopy") 
    {
        // Here we create another outputfile 
        TFile *outfile = TFile::Open(infile, "UPDATE"); 

        for (Int_t i=0; i<list.GetEntries(); i++) 
        {
            const TString currentObject = static_cast<TObjString*>(list.At(i))->GetString(); 

            const TObjArray* paths = currentObject.Tokenize(":");  

            assert(paths->GetEntries() == 2); 

            TString relPath = static_cast<TObjString*>(paths->At(1))->GetString(); 

            relPath = RemovePrefix(relPath, "/"); 

            std::cout << "Rel path: " << relPath << std::endl; 

            TObject *found = treeNames->FindObject(relPath); 
            if (found) 
            {
                TString path = static_cast<TObjString*>(found)->GetString(); 
                std::cout << "Found path: " << path << std::endl; 

                TTree *inTree = static_cast<TTree*>(inFile->Get(path)); 

                if (not inTree) 
                {
                    std::cout << "(TrainTestSplit) ERROR: No input tree! Make sure the object exists: " << path << std::endl; 
                    continue; // TODO: Make sure this should not throw an error 
                }

                TString directory = RemoveSuffix(path, TString::Format("/%s", inTree->GetName())); 

                outfile->cd(directory); 

                TDirectory *outDir = outfile->CurrentDirectory(); 

                //TTree *newTree = inTree->CloneTree(0); 

                AddTrainTestTag(factor, inTree, outDir); 

                //if (newTree) 
                //{
                //    newTree->Write();    // Check for empty tree (whole data cut away)
                //}
                //else 
                //{
                //    std::cout << "(TrainTestSplit) ERROR: no tree to write. " << std::cout; 
                //}
            }
            //const TString treeName = static_cast<TObjString*>(treeNames->At(i))->GetString(); 
            //TTree *inTree = static_cast<TTree*>(inFile->Get(treeName)); 

            



            // Recreate the tree structure

            delete paths; 
        }

        outfile->Write(); 

        outfile->Close(); 
    }

    if (copyMethod == "inplaceTag") 
    {
        // Here we create another outputfile 
        TFile *outfile = TFile::Open(infile, "UPDATE"); 

        for (Int_t i=0; i<list.GetEntries(); i++) 
        {
            const TString currentObject = static_cast<TObjString*>(list.At(i))->GetString(); 

            const TObjArray* paths = currentObject.Tokenize(":");  

            assert(paths->GetEntries() == 2); 

            TString relPath = static_cast<TObjString*>(paths->At(1))->GetString(); 

            relPath = RemovePrefix(relPath, "/"); 

            std::cout << "Rel path: " << relPath << std::endl; 

            TObject *found = treeNames->FindObject(relPath); 
            if (found) 
            {
                TString path = static_cast<TObjString*>(found)->GetString(); 
                std::cout << "Found path: " << path << std::endl; 

                TTree *tree = static_cast<TTree*>(inFile->Get(path)); 

                if (not tree) 
                {
                    std::cout << "(TrainTestSplit) ERROR: No input tree! Make sure the object exists: " << path << std::endl; 
                    continue; // TODO: Make sure this should not throw an error 
                }

                outfile->cd(tree->GetDirectory()->GetPath()); 

                AddTrainTestTag(factor, tree); 

                tree->Write("", TObject::kOverwrite); // The hardcore way 

                //if (newTree) 
                //{
                //    newTree->Write();    // Check for empty tree (whole data cut away)
                //}
                //else 
                //{
                //    std::cout << "(TrainTestSplit) ERROR: no tree to write. " << std::cout; 
                //}
            }
            //const TString treeName = static_cast<TObjString*>(treeNames->At(i))->GetString(); 
            //TTree *inTree = static_cast<TTree*>(inFile->Get(treeName)); 

            



            // Recreate the tree structure

            delete paths; 
        }

        outfile->Close(); 
    }

    // Just for the sake of testing 
    if (copyMethod == "addTMVAweight") 
    {
        TFile *outfile = TFile::Open(infile, "UPDATE"); 

        for (Int_t i=0; i<list.GetEntries(); i++) 
        {
            const TString currentObject = static_cast<TObjString*>(list.At(i))->GetString(); 

            const TObjArray* paths = currentObject.Tokenize(":");  

            assert(paths->GetEntries() == 2); 

            TString relPath = static_cast<TObjString*>(paths->At(1))->GetString(); 

            relPath = RemovePrefix(relPath, "/"); 

            std::cout << "Rel path: " << relPath << std::endl; 

            TObject *found = treeNames->FindObject(relPath); 
            if (found) 
            {
                TString path = static_cast<TObjString*>(found)->GetString(); 
                std::cout << "Found path: " << path << std::endl; 

                TTree *tree = static_cast<TTree*>(inFile->Get(path)); 

                if (not tree) 
                {
                    std::cout << "(TrainTestSplit) ERROR: No input tree! Make sure the object exists: " << path << std::endl; 
                    continue; // TODO: Make sure this should not throw an error 
                }

                outfile->cd(tree->GetDirectory()->GetPath()); 

                AddWeightFile("/Users/mhuwiler/cernbox/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/ClassificationB2DsDs_firstguessB2DsDs1_BDT.weights.xml", tree, "tmvaWeights"); ///Users/mhuwiler/cernbox/DoctoralThesis/Analysis/scripts/MVA/trainingBayesian/electron_id/model_optimized/weights.xml
            }

            delete paths; 
        }

        outfile->Close(); 

    }

    delete treeNames; 

    delete inFile; 
}
