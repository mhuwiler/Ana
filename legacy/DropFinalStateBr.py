#!/usr/bin/env python
from __future__ import division, print_function

import json
from argparse import ArgumentParser
from uncertainties import ufloat
from uncertainties.umath import * 


from libEfficiencies import DumpEffs




# Efficiency calculations 



if __name__ == "__main__": 
	from argparse import ArgumentParser

	parser = ArgumentParser(description="GetEfficiency")
	#parser.add_argument("out", action="store", type=str, default="data/etc/FinalStateBranchingFractionsGenerated.json", help="Name of file")
	#parser.add_argument("-N", "--version", dest="iteration", action="store", type=int, default=0, help="Which iteration of inference")
	parser.add_argument("-t", "--table", dest="table", action="store", type=str, default="/Users/mhuwiler/cernbox/DoctoralThesis/Analysis/Presentations/Presentation_23_9_19/finalstatebr.tex", help="Latex fragment with summary table")
	parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="data/etc/FinalStateBranchingFractionsGenerated.json", help="Path for json with yield info")
	
	options = parser.parse_args()


	finalstatebr = {
		"Sig": ufloat(0.0931, 0.0005), # Using pi pi pi nu mode
	}

	finalstatebr["B0toDstarDs"] = ufloat(0.062, 0.006) # Using KKpi pi0 mode (with largest unc.)
	finalstatebr["B0toDstarDsstar"] = finalstatebr["B0toDstarDs"]*ufloat(0.935, 0.007) # Using Ds gamma mode (both modes same unc.)
	finalstatebr["B0toDstarDs1"] = finalstatebr["B0toDstarDsstar"]*ufloat(0.48, 0.11) # Using Ds pi0 mode (largest frac. and unc.)
	finalstatebr["B0toDstarDs0star"] = finalstatebr["B0toDstarDs"]*ufloat(1., 0.2) # Using Ds pi0 mode (only one determined)
	# finalstatebr["B0toDstarD"] = [1.41e-3, 1.6e-4], 
	finalstatebr["B0toDstarD0K"] = ufloat(1., 1.) 
	# finalstatebr["B0toDstarD0Kstar"] = [2.47e-3, 0.21e-3], 
	# finalstatebr["ButoDstarDK"] = [2.61e-3, 2.4e-4], 
	# finalstatebr["ButoDstarD0K"] = [1.3e-2, 1.2e-3], 
	# finalstatebr["ButoDstarpi"] = [1.76e-2, 7e-3], 
	# finalstatebr["BkgDstara1"] = [1.3e-2, 2.7e-3], 
	# finalstatebr["BkgDstar3pi"] = [7.21e-3, 2.9e-4],
	# finalstatebr["B0toDstara1"] = [1.3e-2, 2.7e-3],
	# finalstatebr["B0toDstar3pi"] = [7.21e-3, 2.9e-4],
	# finalstatebr["BkgBuDXc"] = [0.0758, 0.01516], 
	# finalstatebr["BkgB0DD"] = [0.0501, 0.01002]

	DumpEffs(finalstatebr, options.out)


