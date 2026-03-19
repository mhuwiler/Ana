#!/bin/bash


if [[ $ENVDEF == "MyMacOS" ]]; then
	NUMPYLIBRARY=/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/numpy/core/include
elif [[ $ENVDEF == "T3PSI" ]]; then
	NUMPYLIBRARY=/t3home/mhuwiler/.local/lib/python3.6/site-packages/numpy/core/include
	export PATH=/work/mhuwiler/software/Analysis/root626/pythonlib/Python-3.6.8/install/include/python3.6m:$PATH
else 
	echo "Unknown environment "$ENVDEF
fi

if [[ $# == 0 ]]; then
	echo "ERROR: No argument provided. You need to at least provide the ref of the bare file."
	return 5
elif [[ $# > 3 ]]; then
	echo "WARNING: Too many arguments provided. Arguments after 2nd position will be ignored. "
#elif [[ $# == 2 ]]; 
fi 


INITIALREF=${1}
VERSION=""
FINALREF="${1}_mva" # Adding suffix to the ref
#Remove suffix from string
#FINALREF=${INITIALREF%".root"}
STARTEVT=0
STOPEVT=0
#if [[ $# == 2 ]]; then
#	FINALREF=${2}
#fi
if [[ $# -ge 2 ]]; then
	VERSION=${2}
fi
if [[ $# -ge 3 ]]; then
	STOPEVT=${3}
fi
if [[ $# -ge 4 ]]; then
	STARTEVT=${4}
fi

echo $INITIALREF
echo $FINALREF

echo $STARTEVT
echo $STOPEVT

olddir=${PWD}
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd ${SCRIPT_DIR}


#root -e 'gInterpreter->AddIncludePath("/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/numpy/core/include");' TestApplyTFweightnew.C
root -q -e 'gInterpreter->AddIncludePath("'$NUMPYLIBRARY'"); gInterpreter->LoadMacro("../../Tau.h+");' 'ApplyXGBOweight.C("'$INITIALREF'", "'$FINALREF'", "'$VERSION'", '$STOPEVT', '$STARTEVT')'
#root -e 'gInterpreter->AddIncludePath("/opt/local/Library/Frameworks/Python.framework/Versions/3.6/lib/python3.6/site-packages/numpy/core/include");' ApplyTFweight.C
#root -e 'gInterpreter->AddIncludePath("/t3home/mhuwiler/.local/lib/python3.6/site-packages/numpy/core/include/");' ApplyTFweight.C

cd ${olddir}
