from __future__ import division, print_function
import ROOT


def extractVariables(vec, variables): 
	#vec = ROOT.vector("TString")()
	#vec = ROOT.MVAEval.variables
	#vec = []
	ROOT.gInterpreter.Declare("VariableList* castToList(TObject* obj) { return static_cast<VariableList*>(obj); }")
	myvec = ROOT.castToList(vec)
	for item in variables : 
		myvec.Add(ROOT.TString(item[0]+"/"+item[1]))
		#print(item[0])
	return myvec

def getVariables(vec): 
	from anaMVA.BDTvariables import features_save
	return extractVariables(vec, features_save)


if __name__ == "__main__": 

	ROOT.gROOT.LoadMacro("VariableList.h")

	from anaMVA.BDTvariables import features_save

	#extractVariables(features_save)

	#vec = ROOT.vector("TString")()
	vec = ROOT.VariableList()

	result = getVariables(vec)

	print(result)

	for i in range(result.variables.size()): 
		print(result.variables.at(i))


