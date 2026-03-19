# FileManager

FileManager is a class that eases file handling, especially in case of many different objects opened from a variety of files. It makes sure the files are opened only once, and kept open until all objects from within a given file are dereferenced. It is a complementary approach to TChain, which is used for cases where objects (i.e. trees and histograms) are split among multiple files. FileManager proves most useful in case of smaller objects (e.g. histograms, graphs, etc...) opended in several files  with a heteogeneous structure. It was designed in a way that proves especially robust when several objects in use are located in the same file, and used at different places in the code.  

FileManager contains a Collection of Items, which are a tuple of the following three entities: a path of a Root TObject within a file, the file name (path), and a unique reference used to interact with the object through FileManager. The unique reference is the only thing needed to access and interact with the object once it is added to the FileManager. The Items can be added at any time to the collection: 

	FileManager.AddItem(TString uniqueReference, TString filePath, TString object) 

The files are not opened (and one can thus not access to the object) until: 

	FileManager.OpenItem(uniqueReference) 

is called. Equivalently, 

	FileManager.OpenAllItems() 

opens all files existing in the collection at the time of calling. 

It is also possible to specify the permission (e.g. READ, RECREATE, UPDATE, ...)

	FileManager.OpenItem(uniqueReference, "READ") or FileManager.OpenItem(uniqueReference, "RECREATE") 

and

	FileManager.OpenAllItems("READ") 

will open all items with the same permission. 

All methods return in addition to the output and error messages a bool, that is true if then operation succeded, or false if there were some errors. It can be ignored, but it was implemented for upstream code to have an easy feedback on the success of the instruction. 


An Item (pointer to TObject) can be accessed using: 

	TObject* obj = FileManager.GetItem(uniqueReference)

this command assumes the Item was previously opened, and will return nullptr if it is not the case. There is however a way to open the file on the fly: 

	TObject* obj = FileManager.GetItem(uniqueReference, true) 

where the second parameter autoOpen (false by default) will check if the file is open and open it if needed when set to true. 


Items can be closed as following: 

	FileManager.CloseItem(uniqueReference)

This will dereference the TObject (it will no longer be accessible trough GetItem()) and close the file if no other item still uses it. 

As for opening, this can be done for all Items at once: 

	FileManager.CloseAll(); 

The two previous methods do however NOT remove the Items from the Collection, and they can still be opened on request (OpenItem() or OpenAll()) and used again. 

To completely remove Items from the Collection, one has to do: 

	FileManager.Clear()

This removes all Items from the collection in addition to closing all files. 
