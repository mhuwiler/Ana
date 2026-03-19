#!/bin/bash

INITIALREF=${1}_SR
#FINALREF=${1}_AR
VERSION=${2}
WEIGHTFILE=${3}


# if [[ $# -ge 2 ]]; then
# 	VERSION=${2}
# fi

root -q -x 'AddFinalMVASimple.C("'$INITIALREF'", "'$VERSION'", "'$WEIGHTFILE'")'
