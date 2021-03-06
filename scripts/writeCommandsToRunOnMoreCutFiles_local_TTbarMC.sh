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
# analysis
files="/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleMacrosV2/config2015/TTbarBkg/cutTable_lq_ttbar_emujj.txt
/afs/cern.ch/user/s/scooper/work/private/cmssw/8011/TestRootNTuplizerRecipe/src/Leptoquarks/analyzer/rootNtupleMacrosV2/config2015/TTbarBkg/cutTable_lq_ttbar_eejj.txt
"

#------------
OUTDIRPATH=$LQDATA  # a subdir will be created for each cut file 
SUBDIR=2016ttbar/mar17_reeEmu_ttbarMC
#SUBDIR=2016ttbar/mar2_reeEmu_noTrig_ttbarMC
#SUBDIR=2016ttbar/mar2_reeEmu_ttbarMC
#SUBDIR=2016ttbar/feb13_updateMuonSF_reeEmu_ttbarMC
#SUBDIR=2016ttbar/feb6_gsfEtaCheck_reeEmu_ttbarMC_bugfix
#SUBDIR=2016ttbar/feb6_gsfEtaCheck_reeEmu_ttbarMC_addFinalSels
#SUBDIR=2016ttbar/feb2_gsfEtaCheck_reeEmu_ttbarMC
#SUBDIR=2016ttbar/jan31_gsfEtaCheck_reeEmu_ttbarMC
#SUBDIR=2016ttbar/jan24_gsfEtaCheck_reeEmu_ttbarMC
#SUBDIR=2016ttbar/nov19_ptEE_reeEmu_ttbarMC
#SUBDIR=2016ttbar/oct2_ptEE_ele27wptightOREle115ORPhoton175_reeEmu_ttbarMC
#SUBDIR=2016ttbar/jul4_ele27wptightOREle115ORPhoton175_reeEmu_ttbarMC
#SUBDIR=2016ttbar/may29_ele27wptightOREle115ORPhoton175_reeEmu_ttbarMC
#SUBDIR=2016ttbar/may22_ele27wptightOREle115_reeEmu_ttbarMC
#SUBDIR=2016ttbar/apr11_ele27wptightOREle115_reeEmu_ttbarMC
#SUBDIR=2016ttbar/mar30_reeEmu_ttbarMC
#SUBDIR=2016ttbar/jan21_reeEmu_ttbarMC
# output sub-directory (i.e. output will be in OUTDIRPATH/SUBDIR)
# it is suggested to specify the luminosity in the name of the directory
#------------
ILUM=35867 # [was 36455] ntupleV235 2016B-H rereco runs # integrated luminosity in pb-1 to be used for rescaling/merging MC samples
#ILUM=12900 # ICHEP2016
#ILUM=6910 # ICHEP2016 minus early runs
FACTOR=1000 # numbers in final tables (but *not* in plots) will be multiplied by this scale factor (to see well the decimal digits)
#------------
EXE=mainEEJJttbar
CODENAME=analysisClass_lq_ttbarEst
#------------
INPUTLIST=config/TTbarSk_mar16_v237_local_comb/inputList_ttbarMC.txt
#INPUTLIST=config/TTBarSkim_feb1_SEleL_v237_eoscms_comb/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_jan25_SEleL_v237_eoscms_comb/inputListAllCurrent.txt
#INPUTLIST=config/TTBarSkim_jan23_SEleL_v237_eoscms_comb/inputListAllCurrent.txt
#INPUTLIST=config/TTBarSkim_nov16_tuplev236_comb_eoscms/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_oct2_SEleL_reminiAOD_v236_eoscms/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_may21_SEleL_reminiAOD_v236_eoscms/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_SEleL_reminiAOD_v236/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_SEleL_rereco_v235/inputList_TTBarMC.txt
#INPUTLIST=config/TTBarSkim_SEleL_spring16mc_rereco_v233_heep7/inputList_TTBarMC.txt
#INPUTLIST=config/TestCombinationMay4_TTbarSkim/inputListTTbarMC.txt
#INPUTLIST=config/TestCombinationMay4/inputListAllCurrent.txt
#------------
XSECTION=config/xsection_13TeV_2015.txt #specify cross section file
#XSECTION=versionsOfAnalysis_eejj/1jun_ttbarRescale/xsection_13TeV_2015_TTbarRescale.txt # with ttbar rescaling
#XSECTION=config/xsection_13TeV_2015_Zrescale.txt #first try at Z rescale
#------------
SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_ttbarMC.txt
#SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_eejj.txt
#SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_eejj_QCD.txt
#------------
#NCORES=8 #Number of processor cores to be used to run the job
NCORES=16 #Number of processor cores to be used to run the job
#------------

#### END OF INPUTS ####

COMMANDFILE=commandsToRunOnMoreCutFiles_ttbarMc_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_optimization_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_optimization_QCD_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_QCD_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_trigEff_reHLT_AK5_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_AK5_local_`hostname -s`.txt
#COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_AK4CHS_local_`hostname -s`.txt
echo "" > $COMMANDFILE

for file in $files
do
suffix=`basename $file`
suffix=${suffix%\.*}
cat >> $COMMANDFILE <<EOF

####################################################
#### launch, check and combine cmds for $suffix ####

time python scripts/launchAnalysis.py \
    -e $EXE \
    -k $CODENAME \
    -i $INPUTLIST \
    -n rootTupleTree/tree \
    -c $file \
    -o $OUTDIRPATH/$SUBDIR/output_$suffix  \
    -p $NCORES \
    -k $CODENAME \
    >& launch_${suffix}_`hostname -s`.log

mv launch_${suffix}_`hostname -s`.log $OUTDIRPATH/$SUBDIR/output_$suffix/

time  ./scripts/combinePlots.py \
    -i $INPUTLIST \
    -c $CODENAME \
    -d $OUTDIRPATH/$SUBDIR/output_$suffix \
    -l  `echo "$ILUM*$FACTOR" | bc` \
    -x $XSECTION  \
    -o $OUTDIRPATH/$SUBDIR/output_$suffix \
    -s $SAMPLELISTFORMERGING \
    -f $OUTDIRPATH/$SUBDIR/output_$suffix/launch_${suffix}_`hostname -s`.log \
    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log

EOF
done


echo "The set of commands to run on the cut files:" 
for file in $files
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
