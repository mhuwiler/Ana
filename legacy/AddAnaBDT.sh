#!/bin/bash

INITIALREF=${1}_DNN
FINALREF=${1}
VERSION=""

if [[ $# -ge 2 ]]; then
	VERSION=${2}
fi

root -q -x 'AddMVAVariableSimple.C("'$INITIALREF'", "'$FINALREF'", "'$VERSION'", "anaMVA/NewSelection/model_optimized/weights.xml")'
