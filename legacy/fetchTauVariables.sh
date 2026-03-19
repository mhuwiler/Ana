#!/bin/bash

DATESTRING=$(date '+%Y_%m_%d')

cp Variables.xml Variables_${DATESTRING}.xml
scp t3psi:/work/mhuwiler/software/Analysis/production/CMSSW_10_6_35_patch1/src/EXOVVNtuplizerRunII/Ntuplizer/Variables.xml . 
python GenerateOfflineCode.py

