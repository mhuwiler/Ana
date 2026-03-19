#!/bin/bash

REF=${1}
VERSION=""
WS=0

if [[ $# -ge 2 ]]; then
	VERSION=${2}
fi

if [[ $# -ge 3 ]]; then
	WS=${3}
fi

. clean_tau.sh
root -q -e ".L Tau.h+"  'ProcessingDenominatorOffline.C("'$REF'", "'$VERSION'", '$WS')'
