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

	if [ $5 -le 1 ]; then
		. SelectBCandidate.sh $1 $2 $4

		RETURNCODE=$?
		if [ $RETURNCODE -ne 0 ]; then
			echo "Abort due to error. "
			return $RETURNCODE
		fi
	fi

	if [ $5 -ge 2 ]; then
		echo "ERROR: Only 1 steps. Starting from step $4 will have no effect. "
	fi

}

# Main 

processFull $ARG $VERSION $STOPEVT $WRONGSIGN $STEP 

