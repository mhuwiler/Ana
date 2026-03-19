#!/bin/bash

ARG=${1}
VERSION=""
STEP=0
STOPEVT=0
WRONGSIGN=0

if [[ $# -ge 2 ]]; then
	VERSION=${2}
fi
if [[ $# -ge 3 ]]; then
	STEP=${3}
fi
if [[ $# -ge 4 ]]; then
	WRONGSIGN=${4}
fi
if [[ $# -ge 5 ]]; then
	STOPEVT=${5}
fi

function processFull()
{

	if [ $5 -le 0 ]; then
		. ApplyTFweight.sh $1 $2 $3

		RETURNCODE=$?
		if [ $RETURNCODE -ne 0 ]; then
			echo "Abort due to error. "
			return $RETURNCODE
		fi
	fi

	if [ $5 -le 1 ]; then
		. SelectTauCandidate.sh $1 $2 $4

		RETURNCODE=$?
		if [ $RETURNCODE -ne 0 ]; then
			echo "Abort due to error. "
			return $RETURNCODE
		fi
	fi

	if [ $5 -le 2 ]; then
		. AddAnaBDT.sh $1 $2
		#sleep 2
		#./anaMVA/inference/ApplyXGBOweight.sh $1 $2

		RETURNCODE=$?
		if [ $RETURNCODE -ne 0 ]; then
			echo "Abort due to error. "
			return $RETURNCODE
		fi
	fi

	if [ $5 -ge 3 ]; then
		echo "ERROR: Only 3 steps. Starting from step $4 will have no effect. "
	fi

}

# Main 

processFull $ARG $VERSION $STOPEVT $WRONGSIGN $STEP 

