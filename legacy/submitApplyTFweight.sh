#! /bin/bash
# This is my batch script
#SBATCH --mem=2400
#SBATCH --account=t3
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=standard

# now start our executable
RANDOMSTRING=$(echo $RANDOM | md5sum | head -c 20; echo;)

INITIALREF=${1}
FINALREF="${INITIALREF}_DNN" #${2}
LOGFILE=/work/mhuwiler/software/Analysis/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/submit/$FINALREF.log
TEMP=./test/${RANDOMSTRING}/ #/scratch/mhuwiler/Analysis/production/${RANDOMSTRING}/
DESTINATION=/work/mhuwiler/software/Analysis/CMSSW_11_1_0/src/EXOVVNtuplizerRunII/Ntuplizer/submit/
MAXEVENTS=1000 # 0 means all 

NUMPYLIBRARY=/t3home/mhuwiler/.local/lib/python3.6/site-packages/numpy/core/include
export PATH=/work/mhuwiler/software/Analysis/root626/pythonlib/Python-3.6.8/install/include/python3.6m:$PATH

# Usage: 
# sbatch submithadd.sh ProductionName /pnfs/lcg.cscs.ch/cms/trivcat//store/user/mhuwiler/DATA/BtoDstarTauNu/TauSelection/FirstRun/ParkingBPH2/ParkingBPH2_Run2018B-UL2018_MiniAODv2-v1ULFirstWithFix/220906_152308/

echo Starting job

mkdir -p $TEMP

STARTTIME=$(date +%s.%N)
. /work/mhuwiler/software/Analysis/root626/setuproot.sh
#cmsRun config_MCfull.py 2>&1 | tee runtest.log
root -e 'gInterpreter->AddIncludePath("'$NUMPYLIBRARY'");' 'ApplyTFweight.C("'$INITIALREF'", "'$FINALREF'", '$MAXEVENTS', "'$TEMP'")' 2>&1 | tee $LOGFILE

echo "Working node: "$(hostname) >> $LOGFILE
ENDTIME=$(date +%s.%N)
echo "Execution time:" $(date -u -d "0 $ENDTIME sec - $STARTTIME sec" +"%H:%M:%S") >> $LOGFILE

cp ${TEMP}/*.root ${DESTINATION}
cp ${LOGFILE} ${DESTINATION}
rm -r ${TEMP}

