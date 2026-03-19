#!/bin/bash

INITIALREF=${1}_tf
FINALREF="${1}_DNN"
REF=${1}
VERSION=""
WS=0

if [[ $# -ge 2 ]]; then
	VERSION=${2}
fi

if [[ $# -ge 3 ]]; then
	WS=${3}
fi

echo $INITIALREF
echo $FINALREF

. clean_tau.sh
#root -q -e ".L Tau.h+"  'UpdateTauDNN.C("'$REF'", "'$VERSION'", '$WS')'
root -q 'UpdateTauDNN.C("'$REF'", "'$VERSION'", '$WS')'
