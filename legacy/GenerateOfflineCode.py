#!/usr/bin/env python


from __future__ import division, print_function

import os
import sys
import time
from datetime import datetime
from argparse import ArgumentParser
from xml.dom.minidom import parse


def stripLine(line, pattern): 
	line = line.replace(pattern, "")
	t = time.strptime(line, "%H:%M:%S")
	dt = datetime.fromtimestamp(time.mktime(t))
	timevalue = dt.hour + dt.minute/60.
	return timevalue



def PromptYesNo(answerasbool=False): 
	# Inspired from Fabrice Couderc 
	rep = ''
	while not rep in [ 'yes', 'no' ]:
		rep = raw_input( "(type 'yes' or 'no'): " ).lower()
	if (answerasbool): 
		if (rep == 'yes'): 
			return True
		else: 
			return False
	return rep

def HoldUntilKeyPress(otherThanEnter=False): 
	# Inspired from Fabrice Couderc 
	rep = ''
	if otherThanEnter: 
		while rep in [ "" ]:
			rep = raw_input( "Press a key to continue... " ).lower()
	else:
		rep = raw_input( "Press a key to continue... " ).lower()
	
	return

if __name__ == "__main__":

	parser = ArgumentParser(description="CheckFiles")
	parser.add_argument("--debug", dest="debug", action="store_true", default=False, help="Turn on debug output")
	parser.add_argument("--taubranchname", dest="taubranchname", action="store", type=str, default="b_tau", help="Name of tau object branch")
	parser.add_argument("--objname", dest="objname", action="store", type=str, default="tau", help="Name of tau object in UpdateTauDNN.C")
	parser.add_argument("--match", dest="match", action="store_true", default=False, help="Do gen matching")

	defaultvalues = {"f": "-999.", "i": "-999", "b": "false", "d": "-999.", "s": "\"\""  } # since code generation given as string

	candidates = ["tau", "B"]

	
	options = parser.parse_args()

	with open("Variables.xml") as variables: 
		document = parse(variables)

		print(document.encoding)

		count = 0

		with open(os.path.expandvars("tauvariables.gcf"), "w") as tauvariables: 
			with open("taubranchwriting.gcf", "w") as taubranchwriting: 
				with open("tauvariableargs.gcf", "w") as tauaffectationargs: 
					with open("taubranchreading.gcf", "w") as taubranchreading: 
						with open("tauvariableaffectation.gcf", "w") as tauaffectation: 
							with open("taubranchcreation.gcf", "w") as taubranchcreation: 
								with open("stringbranches.gcf", "w") as stringbranches: 
									stringbranches.write("std::vector<std::string> stringbranches = {")
									notFirst = False

									for candidate in candidates: 
										element = document.getElementById(candidate)
										
										objectvars = element.getElementsByTagName("Variable")
										lines = []
										for var in objectvars: 
											if (options.debug): print("Attributes: name={} type={}".format(var.getAttribute("name"), var.getAttribute("type")))
											name = var.getAttribute("name")
											if (candidate != "tau"): 
												varname = "{}_{}".format(candidate, name)
											else: 
												varname = name
											branchname = "{}_{}".format(candidate, name)
											assert(name), "Variable has no field \"name\""
											vartype = var.getAttribute("type")
											assert(vartype), "Variable {} no field \"type\"".format(name)
											vartype = vartype.lower()
											assert(vartype in defaultvalues), "Variable {} has unsupported type: {}".format(name)
											branch = var.getAttribute("branch")
											if (not branch): 
												branch = "b_{}".format(branchname)
												if (options.debug): print("Using variable name as branch name: {}".format(branch))
											vbranch = var.getAttribute("vbranch")
											if (not vbranch): 
												vbranch = "v_"+branchname
												if (options.debug): print("Using variable default name as vbranch name: {}".format(vbranch))
											default = var.getAttribute("default")
											if (not default): 
												default = defaultvalues[vartype]
												if (options.debug): print("Using variable default name as default name: {}".format(default))

											# Start writing out the generated code
											cpptype = ""
											branchtypesuffix = ""
											branchtype = ""
											if (vartype == "f"): 
												cpptype = "float"
												branchtypesuffix = "/F"
											elif (vartype == "i"): 
												cpptype = "int"
												branchtypesuffix = "/I"
											elif (vartype == "b"): 
												cpptype = "int"
												branchtypesuffix = "/I" # Use here O for bool branches (if int branch OK, implicit conversion: https://stackoverflow.com/questions/5369770/bool-to-int-conversion)
												branchtype = "int"
											elif (vartype == "d"):
												cpptype = "double"
												branchtypesuffix = "/F" # If really want to store double put D here 
												branchtype = "float"
											elif (vartype == "s"):
												cpptype = "std::string"
												branchtypesuffix = "" 
												if (notFirst): stringbranches.write(", ")
												notFirst = True
												stringbranches.write("\"{}\"".format(vbranch))
											else: 
												raise AttributeError, "No known type: {}".format(vartype)
											assert(cpptype)
											#assert(branchtypesuffix)

											other = "init_{}".format(varname)

											funcName = "Write{}".format(varname.title())

											tauvariables.write("{} {} = {};\n".format(cpptype, varname, default))

											taubranchwriting.write("static {} {}(const Tau& tau)\n{{\n\t return {}.{};\n}};\n\n".format(cpptype, funcName, options.objname, varname))

											if (count > 0): 
												tauaffectationargs.write(", ")
												taubranchreading.write(", ")

											tauaffectationargs.write("ROOT::VecOps::RVec<{}> {}".format(cpptype, other))

											taubranchreading.write("\"{}\"".format(vbranch))

											tauaffectation.write("{}.{} = {}.at(i);\n".format(options.objname, varname, other))

											taubranchcreation.write(".Define(\"{}\", Tau::{}, {{\"{}\"}})".format(branch, funcName, options.taubranchname))

	#										branchinclude.write("{} {};\n".format(branchtype, branch))
	#										branchinclude.write("std::vector<{}> {};\n".format(branchtype, vbranch))

	#										lines.append("tree_->Branch(\"{}\", &{}, \"{}\");\n".format(branch, branch, branch+branchtypesuffix))
	#										bracnhassignement.write("tree_->Branch(\"{}\", &{});\n".format(vbranch, vbranch))

	#										branchcleaning.write("{}.clear();\n".format(vbranch))
	#										branchcleaning.write("{} = {};\n".format(branch, default)) # Some variables need to be reset
										
											count+=1

									stringbranches.write("};")

					for line in lines: 
						bracnhassignement.write(line)


	print("\nSucessfully generated object variables!\n\n")

		


		


