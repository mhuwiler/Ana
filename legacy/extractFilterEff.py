#!/usr/bin/env python
from __future__ import division, print_function

import pandas as pd
import os
import time
from uncertainties import ufloat
from uncertainties.umath import * 
from argparse import ArgumentParser
#from libEfficiencies import getEff, DumpEffs, ReadEffs
from collections import OrderedDict
import json


def FormatCustom(filecontent): 
	output = filecontent.replace("{","{\n\t").replace("}", "\n}").replace("], ", "], \n\t")
	return output


def DumpEffs(effs, path): 
	effsForWrite = OrderedDict()
	for item, content in effs.iteritems(): 
			eff = effs[item]
			effsForWrite[item] = (eff.n, eff.s)
	with open(path, "w") as file: 
		print(effsForWrite)
		forwrite = FormatCustom(json.dumps(effsForWrite, ensure_ascii=False, encoding="utf8", sort_keys=False)) #indent=4, 
		print(forwrite)
		file.write(forwrite)

def ReadEffs(path): 
	effs = OrderedDict()
	with open(path, "r") as file: 
		effsFromFile = json.load(file, encoding="utf8", object_pairs_hook=OrderedDict)
		for item, content in effsFromFile.iteritems(): 
				eff = effsFromFile[item]
				assert(len(eff)==2)
				effs[item] = ufloat(eff[0], eff[1])
	return effs


def getEff(n, N): 
	eff = float(n)/float(N)
	#print eff
	err = sqrt(eff*(1.-eff)/float(N))
	#print err
	#return eff, err
	return ufloat(eff, err)


def getFilterEffFromSummary(file): 
	frame = pd.read_csv(file, delimiter=" ", header=None)

	#print(frame)

	numbers = pd.to_numeric(frame.iloc[:,6].str.strip('()'))

	initial = pd.to_numeric(frame.iloc[:,8].str.strip('()'))

	#print(numbers)


	N = initial.sum()

	n = numbers.sum()


	eff = getEff(n, N)

	return eff


def updateEfficiency(eff, sample, file, force=False): 
	filterEffs = ReadEffs(file)
	print (filterEffs)

	if ((sample in filterEffs) and (not force)): 
		raise KeyError("The efficiency for {} already exists in the dictionary. Use option -f to overwrite.".format(sample))
	else: 
		filterEffs[sample] = eff
	DumpEffs(filterEffs, file)


def ensureEOSDirectory(directory): 
	if (not os.path.isdir(directory)): 
		os.system("kinit mhuwiler@CERN.CH; eosfusebind -g krb5 $HOME/krb5cc_$UID")
		time.sleep(2)

	return os.path.isdir(directory)


if __name__ == "__main__": 

	parser = ArgumentParser(description="GetEfficiency")
	parser.add_argument("-i", "--file", dest="file", action="store", type=str, default="./summaryfiltereffs.txt", help="File from which to extract filter efficiency")
	parser.add_argument("-s", "--sample", dest="sample", action="store", type=str, default="", help="File from which to extract filter efficiency")
	parser.add_argument("-u", "--update", dest="update", action="store_true", default=False, help="Update filter efficiency in efficiency file ")
	parser.add_argument("-d", "--database", dest="database", action="store", type=str, default="/eos/home-m/mhuwiler/DoctoralThesis/Analysis/scripts/Ana/data/etc/FilterEfficiencies.json", help="Database file where to store the efficiency")
	parser.add_argument("-o", "--out", dest="out", action="store", type=str, default="/eos/home-m/mhuwiler/DoctoralThesis/Analysis/normalisation", help="Folder where to output computed eff")
	parser.add_argument("-f", "--force", dest="force", action="store_true", default=False, help="Force updating efficiency if already present")
	parser.add_argument("-b", "--batch", dest="batch", action="store_true",default=False, help="Run in batch mode")

	options = parser.parse_args()


	eff = getFilterEffFromSummary(options.file)

	if (options.update): 
		assert(options.sample != ""), "ERROR: provided no sample name for updating the efficiency."
		updateEfficiency(eff, options.sample, options.database, options.force)


	ensureEOSDirectory(options.out)
	singleEff = { options.sample: eff }
	DumpEffs(singleEff, options.out+"/{}FilterEff.json".format(options.sample))


	print(eff)


	print("Efficiency for {}:\n\n\t\"{}\": [{}, {}],\n\n".format(options.sample, options.sample, eff.n, eff.s))


	