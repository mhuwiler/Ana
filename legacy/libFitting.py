#!/usr/bin/env python
from __future__ import division, print_function

import ROOT



def WriteDatacard(frames, yields, variables, regions, dataname, datacardname="datacard.txt", workspacefile="workspace.root", workspacename="w"): 
	from datetime import datetime
	from ROOT import RooRealVar, RooArgSet, RooWorkspace, RooDataHist
	MC = list(frames.keys())
	MC.remove(dataname)
	variable = "b_tau_rhomass1"

	with open(datacardname, "w") as datacard: 
		datacard.write("# Datacard generated automatically with {}{} on {}.\n".format(os.getcwd(), __file__, datetime.today().strftime("%d.%m.%y %H:%M:%S")))
		datacard.write("# Simple fit \n\n")
		datacard.write("imax {}\n".format(len(regions)))
		datacard.write("jmax {}\n".format(len(MC)))
		datacard.write("kmax {}\n".format(0)) # For now no systematics
		datacard.write("\n"+"-"*50+"\n")

		datacard.write("# Shapes and RooFit workspace\n")
		file = ROOT.TFile.Open(workspacefile, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		# Creating the variable on which we fit
		var = RooRealVar(variable, variable, 0., 100.)
		fitspace = RooArgSet(var)
		for region in regions: 
			# Writing the data shapes for each region
			name = "data_obs_{}".format(region)
			datacard.write("shapes data_obs {} {} {}\n".format(region, workspacefile, workspacename+":{}".format(name)))
			examplehist = Ana.binning[variable]
			hist = frames[dataname][region].Histo1D(examplehist, variable).GetPtr()
			hist.SetName(name)
			histogram = RooDataHist(name, name, fitspace, hist)
			hist.Write() # Also saving the ROOT hist
			getattr(workspace, "import")(histogram)
			# Writing the MC shapes 
			for item in MC: 
				histname = item+"_"+region
				datacard.write("shapes {} {} {} {}\n".format(item, region, workspacefile, workspacename+":"+histname))
				hist = frames[item][region].Histo1D(examplehist, variable).GetPtr()
				#hist.Scale(yields[item][region].n/hist.Integral())
				hist.SetName(histname)
				roohist = ROOT.RooDataHist(histname, histname, fitspace, hist)
				hist.Write()
				getattr(workspace, "import")(roohist)

		workspace.Write()
		file.Write()
		file.Close()
		datacard.write("\n"+"-"*50+"\n")
		
		datacard.write("# Observed events (data)\n")
		regionstring = "bin "
		for item in regions: 
			regionstring += (item+" ")
		regionstring+="\n"
		datacard.write(regionstring)

		observationstring = "observation "
		for item in regions: 
			print(frames[dataname][item].Count().GetValue())
			observationstring += ("-1 ") # "{} ".format(frames["dataD2"][item].Count().GetValue()) # TODO: fix
		datacard.write(observationstring+"\n")
		# We want to leave a few components floating 
		datacard.write("\n"+"-"*50+"\n")

		datacard.write("# Expected events (MC/model)\n")
		binstring = "bin "
		labelstring = "process "
		indexstring = "process "
		expectedstring = "rate "
		count = 0
		for region in regions: 
			for item in MC: 
				binstring += "{} ".format(region)
				labelstring += "{} ".format(item)
				factor = 1
				if "Sig" in item: 
					factor = -1 # make signal negative
				indexstring += "{} ".format(factor*count)
				expectedstring += "{} ".format(yields[item][region])
				count += 1
		datacard.write(binstring+"\n")
		datacard.write(labelstring+"\n")
		datacard.write(indexstring+"\n")
		datacard.write(expectedstring+"\n")
		#datacard.write("\n"+"-"*50+"\n")
		#datacard.write("lumi     lnN    1.10       1.0 		1.0\n")
		datacard.write("\n"+"-"*50+"\n")
		# Writing out the constraints
		for region in regions: 
			for item in MC: 
				datacard.write("{}_{}_norm rateParam {} {} {} [{},{}]\n".format(item, region, region, item, yields[item][region].n, 0, yields[item][region].n*5.))
		datacard.write("\n"+"-"*50+"\n")
		#for item in variables: 
			#datacard.write("{} flatParam\n".format(item.GetName()))


def WriteWorkspace(frames, yields, variables, regions, dataname, filename="workspace.root"): 
	from ROOT import RooRealVar, RooDataHist, RooArgSet
	from ROOT import Ana

	MC = list(frames.keys())
	print(MC)
	print(dataname)
	MC.remove(dataname)

	for variable in variables: 
		file = ROOT.TFile.Open(filename, "RECREATE")
		# Creating the variable on which we fit
		for region in regions: 
			#directory = file.mkdir(region)
			#directory.cd()
			examplehist = Ana.binning[variable]

			# Writing the data shapes for each region
			name = "data_obs"
	
			hist = frames[dataname][region].Histo1D(examplehist, variable).GetPtr()
			hist.SetName(name)
			hist.Write() # Also saving the ROOT hist
			# Writing the MC shapes 
			for item in MC: 
				histname = item #+"_"+region #+"_"+variable
				hist = frames[item][region].Histo1D(examplehist, variable).GetPtr()
				norm = -1.
				try: 
					norm = yields[item][region].n
					hist.Scale(norm/hist.Integral())
				except:
					pass
				hist.SetName(histname)
				hist.Write()

			# Combinatorial background
			from anaPrepareRegions import GetABCDcomponent
			import anaConfig
			ABCDnorm = -1.
			with open("ABCDnorm.txt", "r") as normfile: 
				yieldtxt = normfile.readlines()
				assert(len(yieldtxt) == 1)
				ABCDnorm = float(yieldtxt[0].rstrip())
				print(ABCDnorm)
			comb = GetABCDcomponent(anaConfig.data, "(b_tau_sumdnn>2.)", "b_B_nmu<1&&b_B_ne<1&&b_B_nh<1", examplehist, variable, [item for item in anaConfig.fitsamples if ((item != anaConfig.Sig) and (item != "data_obs") and (item != "ABCD"))], "ShapefileABCD_fit.root")
			ABCDnorm = 0.9*ABCDnorm
			print("Norm: {}".format(ABCDnorm))
			comb.Scale(ABCDnorm/comb.Integral())
			comb.SetName("ABCD")
			comb.Write()

		file.cd()

		file.Write()
		file.Close()


def WriteWorkspaceDataset(frames, yields, variables, regions, dataname, filename="workspace.root"): 
	from ROOT import RooRealVar, RooDataHist, RooArgSet, RooDataSetHelper
	from ROOT import Ana
	#import anaConfig

	workspacename = "w"

	MC = frames.keys()
	print(MC)
	print(dataname)
	MC.remove(dataname)

	for region in regions: 
		file = ROOT.TFile.Open(filename.replace(".root", "_"+region+".root"), "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)

		roovariables = []
		variablenames = ()
		roovars = RooArgSet()
		for variable in variables: 
			
			# Creating the variable on which we fit
			binning = Ana.binning[variable]
			var = RooRealVar(variable, variable, binning.fXLow, binning.fXUp)
			roovariables.append(var)
			variablenames += (variable, )
			roovars.add(var)
			ROOT.SetOwnership(var, 0)


		print(variablenames)
		print(roovars)


		# https://root.cern/doc/master/rf408__RDataFrameToRooFit_8py.html
		roodatasethelper = RooDataSetHelper("data_obs", "data_obs", roovars)
		dataset = frames[dataname][regions[0]].Book(ROOT.std.move(roodatasethelper), variablenames)
		getattr(workspace, "import")(dataset.GetValue())
		for item in MC: 
			roodatasethelper = RooDataSetHelper(item, item, roovars)
			roodataset = frames[item][regions[0]].Book(ROOT.std.move(roodatasethelper), variablenames)
			getattr(workspace, "import")(roodataset.GetValue())


		workspace.Write()

		file.Write()
		file.Close()


def WriteWorkspaceRooType(frames, yields, variables, regions, dataname, filename="workspace.root", workspacename="w"): 
	from ROOT import RooRealVar, RooDataHist, RooArgSet
	from ROOT import Ana

	MC = frames.keys()
	print(MC)
	print(dataname)
	#MC.remove(dataname)
	#MC.pop(dataname)

	for variable in variables: 
		file = ROOT.TFile.Open(filename, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		# Creating the variable on which we fit
		var = RooRealVar(variable, variable, 0., 100.)
		fitspace = RooArgSet(var)
		for region in regions: 
			examplehist = Ana.binning[variable]

			# Writing the data shapes for each region
			name = "data_obs_{}".format(region, variable)
	
			hist = frames[dataname][region].Histo1D(examplehist, variable).GetPtr()
			hist.SetName(name)
			histogram = RooDataHist(name, name, fitspace, hist)
			hist.Write() # Also saving the ROOT hist
			getattr(workspace, "import")(histogram)
			# Writing the MC shapes 
			for item in MC: 
				histname = item+"_"+region #+"_"+variable
				hist = frames[item][region].Histo1D(examplehist, variable).GetPtr()
				norm = -1.
				try: 
					norm = yields[item][region].n
					hist.Scale(norm/hist.Integral())
				except:
					pass
				hist.SetName(histname)
				roohist = ROOT.RooDataHist(histname, histname, fitspace, hist)
				hist.Write()
				getattr(workspace, "import")(roohist)

			# Combinatorial background
			from anaPrepareRegions import GetABCDcomponent
			import anaConfig
			comb = GetABCDcomponent(anaConfig.data, "(b_tau_sumdnn>2.)", "b_B_nmu<1&&b_B_ne<1&&b_B_nh<1", examplehist, variable)
			comb.SetName("ABCD")
			comb.Write()
			roohist = ROOT.RooDataHist(histname, "ABCD", fitspace, comb)
			getattr(workspace, "import")(roohist)

		workspace.Write()
		file.Write()
		file.Close()


def WriteWorkspaceWithSyst(frames, yields, variables, regions, dataname, filename="workspace.root", workspacename="w"): 
	from ROOT import RooRealVar, RooDataHist, RooArgSet
	from ROOT import Ana

	MC = frames.keys()
	print(MC)
	print(dataname)
	MC.remove(dataname)

	for variable in variables: 
		file = ROOT.TFile.Open(filename, "RECREATE")
		workspace = ROOT.RooWorkspace(workspacename)
		# Creating the variable on which we fit
		var = RooRealVar(variable, variable, 0., 100.)
		fitspace = RooArgSet(var)
		for region in regions: 
			examplehist = Ana.binning[variable]

			# Writing the data shapes for each region
			name = "data_obs_{}".format(region, variable)
	
			hist = frames[dataname][region].Histo1D(examplehist, variable).GetPtr()
			hist.SetName(name)
			histogram = RooDataHist(name, name, fitspace, hist)
			hist.Write() # Also saving the ROOT hist
			getattr(workspace, "import")(histogram)
			# Writing the MC shapes 
			for item in MC: 
				histname = item+"_"+region #+"_"+variable
				hist = frames[item][region].Histo1D(examplehist, variable).GetPtr()
				try: 
					norm = yields[item][region].n
					hist.Scale(norm/hist.Integral())
				except:
					pass
				#hist.Scale(yields[item][region].n/hist.Integral())
				hist.SetName(histname)
				roohist = ROOT.RooDataHist(histname, histname, fitspace, hist)
				hist.Write()
				getattr(workspace, "import")(roohist)

		workspace.Write()
		file.Write()
		file.Close()


def WriteDatacardSimple(frames, yields, variables, regions, fitvariable, dataname, datacardname="datacard.txt", workspacefile="workspace.root", workspacename="w"): 
	import os, copy
	from datetime import datetime
	from ROOT import RooRealVar, RooArgSet, RooWorkspace, RooDataHist
	MC = copy.deepcopy(list(frames.keys()))
	MC.remove(dataname)
	print(MC)

	with open(datacardname, "w") as datacard: 
		datacard.write("# Datacard generated automatically with {}{} on {}.\n".format(os.getcwd(), __file__, datetime.today().strftime("%d.%m.%y %H:%M:%S")))
		datacard.write("# Simple fit \n\n")
		datacard.write("imax {}\n".format(len(regions)))
		datacard.write("jmax {}\n".format(len(MC)-1))
		datacard.write("kmax {}\n".format(0)) # For now no systematics
		datacard.write("\n"+"-"*50+"\n")

		datacard.write("# Shapes and RooFit workspace\n")
		for region in regions: 
			# Writing the data shapes for each region
			name = "data_obs" #"data_obs_{}".format(region)
			datacard.write("shapes data_obs {} {} {}\n".format(region, workspacefile, region+"/{}".format(name)))
			# Writing the MC shapes 
			for item in MC: 
				histname = item #+"_"+region
				datacard.write("shapes {} {} {} {}\n".format(item, region, workspacefile, region+"/"+histname))
		datacard.write("\n"+"-"*50+"\n")
		
		datacard.write("# Observed events (data)\n")
		regionstring = "bin "
		for item in regions: 
			regionstring += (item+" ")
		regionstring+="\n"
		datacard.write(regionstring)

		observationstring = "observation "
		for item in regions: 
			print(frames[dataname][item].Count().GetValue())
			observationstring += ("-1 ") # "{} ".format(frames["dataD2"][item].Count().GetValue()) # TODO: fix
		datacard.write(observationstring+"\n")
		# We want to leave a few components floating 
		datacard.write("\n"+"-"*50+"\n")

		datacard.write("# Expected events (MC/model)\n")
		binstring = "bin "
		labelstring = "process "
		indexstring = "process "
		expectedstring = "rate "
		count = 1
		for region in regions: 
			for item in MC: 
				if (("data" in item) and not ("WS" in item)): 
					from uncertainties import ufloat
					yields[item][region] = ufloat(1., 0.)
				binstring += "{} ".format(region)
				labelstring += "{} ".format(item)
				factor = 1
				if "Sig" in item: 
					factor = -1 # make signal negative
				indexstring += "{} ".format(factor*count)
				# Hack for normalising WS
				if ((not "data" in item) or ("WS" in item)): 
					expectedstring += "{} ".format(yields[item][region].n)
				else: 
					expectedstring += "1. "
				count += 1
		datacard.write(binstring+"\n")
		datacard.write(labelstring+"\n")
		datacard.write(indexstring+"\n")
		datacard.write(expectedstring+"\n")
		#datacard.write("\n"+"-"*50+"\n")
		#datacard.write("lumi     lnN    1.10       1.0 		1.0\n")
		datacard.write("\n"+"-"*50+"\n")
		# Writing out the constraints
		for region in regions: 
			for item in MC: 
				datacard.write("{}_{}_norm rateParam {} {} {} [{},{}]\n".format(item, region, region, item, yields[item][region].n, 0, yields[item][region].n*5.))
		datacard.write("\n"+"-"*50+"\n")
		#for item in variables: 
			#datacard.write("{} flatParam\n".format(item.GetName()))

