#!/bin/bash

# Please run this script from the rootNtupleAnalyzerV2 directory by:  
# ./scripts/writeCommandsToRunOnMoreCutFiles.sh

# This scripts creates the whole sets of commands needed to run the analysis on multiple cut files.
# The commands will be written in a text file commandsToRunOnMoreCutFiles.txt in the current directory, 
# to be used by doing cut&paste to a terminal.

# Cut Files should first be created by a script ../rootNtupleMacrosV2/config/eejj/make_eejj_cutFiles.py
# This script will then use those cut files to create the commands needed to run on them.

#### INPUTS HERE ####
#------------
# one file per line
files="/afs/cern.ch/work/s/scooper/private/cmssw/745/LQRootTupleMiniAOD745/src/Leptoquarks/macros/rootNtupleMacrosV2/config2015/ReducedSkims/cutTable_lq1_skim_SingleEle_loose.txt"
#------------
OUTDIRPATH=$LQDATA  # a subdir will be created for each cut file 
SUBDIR="RunII/lq_skim_2015/RootNtuple-v1-1-0-SingleEleLoose_eejjSigAndSomeBackground"
         # output sub-directory (i.e. output will be in OUTDIRPATH/SUBDIR)
         # it is suggested to specify the luminosity in the name of the directory
#------------
FULLEOSDIR="/eos/cms/store/group/phys_exotica/leptonsPlusJets/RootNtuple_skim/RunII/RootNtuple-v1-1-0-SingleEleLoose_eejjSigAndSomeBackground"

#------------
CODENAME=analysisClass_lq1_skim
#------------
#------------
INPUTLIST=config/inputListAllCurrent.txt
#------------
NJOBS=20 #number of jobs for each dataset - for PhysicsDST
WAIT=0 #seconds of delay between submission of different datasets
#------------
QUEUE=8nh #bsub queue  
#------------
#### END OF INPUTS ####

COMMANDFILE=commandsToRunOnMoreCutFiles_MakeSingleEleReducedSkim_nominal_eejjSigAndSomeBackground_`hostname -s |perl -pi -e 's|lxplus[0-9]*|lxplus|'`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_MakeSingleEleReducedSkim_EESdown_LQVector_BetaOneYM500_`hostname -s |perl -pi -e 's|lxplus[0-9]*|lxplus|'`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_MakeSingleEleReducedSkim_EESup_LQVector_BetaOneYM500_`hostname -s |perl -pi -e 's|lxplus[0-9]*|lxplus|'`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_MakeSingleEleReducedSkim_EER_LQVector_BetaOneYM500_`hostname -s |perl -pi -e 's|lxplus[0-9]*|lxplus|'`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_MakeSingleEleReducedSkim_EESdown_LQVector_`hostname -s |perl -pi -e 's|lxplus[0-9]*|lxplus|'`.txt
# etc.
echo "" > $COMMANDFILE

for file in $files
do
suffix=`basename $file`
suffix=${suffix%\.*}
cat >> $COMMANDFILE <<EOF

####################################################
#### launch, check and combine cmds for $suffix ####

  ./scripts/launchAnalysis_batch_ForSkimToEOS.py \
    -i $INPUTLIST \
    -n rootTupleTree/tree \
    -c $file \
    -o $OUTDIRPATH/$SUBDIR/output_$suffix  \
    -j $NJOBS \
    -q $QUEUE \
    -w $WAIT \
    -d $FULLEOSDIR \
    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/launch_${suffix}.log


#### THEN

  ./scripts/check_combine_output_batch_ForSkimToEOS.py \
    -i $INPUTLIST \
    -c $CODENAME \
    -d $OUTDIRPATH/$SUBDIR/output_$suffix \
    -o $OUTDIRPATH/$SUBDIR/output_$suffix \
    -q $QUEUE \
    -s $FULLEOSDIR \
    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/checkcombine_${suffix}.log

EOF
done


echo "The set of commands to run on the cut files:" 
for file in $files
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
