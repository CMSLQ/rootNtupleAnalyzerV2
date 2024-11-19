#!/bin/bash

# Please run this script from the rootNtupleAnalyzerV2 directory by:  
# ./scripts/writeCommandsToRunOnMoreCutFiles.sh

# This scripts creates the whole sets of commands needed to run the analysis on multiple cut files.
# The commands will be written in a text file commandsToRunOnMoreCutFiles.txt in the current directory, 
# to be used by doing cut&paste to a terminal.

# Cut Files should first be created by a script ../rootNtupleMacrosV2/config/eejj/make_eejj_cutFiles.py
# This script will then use those cut files to create the commands needed to run on them.

die () {
    echo >&2 "$@"
    exit 1
}

[ $# -gt 0 ] || die "must specify year as argument; usage: writeCommandsToRunOnMoreCutFiles_batch_eejj_QCD.sh year [doOptimization]"
YEAR=$1

if [[ "$#" -gt 1 ]]; then
  OPT=1
else
  OPT=0
fi

#### INPUTS HERE ####
#------------
#ANANAME=qcd_eejj_11nov2024_presel
#ANANAME=qcd_eejj_13nov2024_bdt_LQToDEle
ANANAME=qcd_eejj_15nov2024_bdt_LQToBEle
#------------
#inputlist2016pre_1FR=config/inputListsRSKQCD_heep_UL16preVFP_24oct24/inputListAllCurrent.txt
#inputlist2016pre_2FR=config/inputListsRSKQCD_heep_UL16preVFP_24oct24/inputListAllCurrent.txt
#inputlist2016post_1FR=config/inputListsRSKQCD_heep_UL16postVFP_24oct24/inputListAllCurrent.txt
#inputlist2016post_2FR=config/inputListsRSKQCD_heep_UL16postVFP_24oct24/inputListAllCurrent.txt
#inputlist2017_1FR=config/inputListsRSKQCD_heep_UL17_24oct24/inputListAllCurrent.txt
#inputlist2017_2FR=config/inputListsRSKQCD_heep_UL17_24oct24/inputListAllCurrent.txt
#inputlist2018_1FR=config/inputListsRSKQCD_heep_UL18_24oct24/inputListAllCurrent.txt
#inputlist2018_2FR=config/inputListsRSKQCD_heep_UL18_24oct24/inputListAllCurrent.txt
#
#inputlist2016pre_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL16preVFP_28oct2024/inputListAllCurrent.txt
#inputlist2016pre_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL16preVFP_28oct2024/inputListAllCurrent.txt
#inputlist2016post_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL16postVFP_28oct2024/inputListAllCurrent.txt
#inputlist2016post_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL16postVFP_28oct2024/inputListAllCurrent.txt
#inputlist2017_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL17_28oct2024/inputListAllCurrent.txt
#inputlist2017_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL17_28oct2024/inputListAllCurrent.txt
#inputlist2018_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL18_28oct2024/inputListAllCurrent.txt
#inputlist2018_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL18_28oct2024/inputListAllCurrent.txt
#
inputlist2016pre_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL16preVFP_11nov2024/inputListAllCurrent.txt
inputlist2016pre_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL16preVFP_11nov2024/inputListAllCurrent.txt
inputlist2016post_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL16postVFP_11nov2024/inputListAllCurrent.txt
inputlist2016post_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL16postVFP_11nov2024/inputListAllCurrent.txt
inputlist2017_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL17_11nov2024/inputListAllCurrent.txt
inputlist2017_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL17_11nov2024/inputListAllCurrent.txt
inputlist2018_1FR=config/inputListsAnalysisQCDPreselectionSkim_1FR_UL18_11nov2024/inputListAllCurrent.txt
inputlist2018_2FR=config/inputListsAnalysisQCDPreselectionSkim_2FR_UL18_11nov2024/inputListAllCurrent.txt
#------------
#xsection2016pre=config/xsection_13TeV_2022.txt
#xsection2016post=config/xsection_13TeV_2022.txt
#xsection2017=config/xsection_13TeV_2022.txt
#xsection2018=config/xsection_13TeV_2022.txt
#
#xsection2016pre=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016preVFP/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_28oct2024_dyjAMCatNLO.txt
#xsection2016post=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016postVFP/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_28oct2024_dyjAMCatNLO.txt
#xsection2017=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2017/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_28oct2024_dyjAMCatNLO.txt
#xsection2018=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2018/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_28oct2024_dyjAMCatNLO.txt
#
xsection2016pre=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016preVFP/eejj_11nov2024_presel_trigSFMod/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_dyjAMCatNLO.txt
xsection2016post=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016postVFP/eejj_11nov2024_presel_trigSFMod/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_dyjAMCatNLO.txt
xsection2017=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2017/eejj_11nov2024_presel_trigSFMod/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_dyjAMCatNLO.txt
xsection2018=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2018/eejj_11nov2024_presel_trigSFMod/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_dyjAMCatNLO.txt
#------------
CODENAME=analysisClass_lq_eejj_QCD
#------------
OUTDIRPATH=$LQDATA  # a subdir will be created for each cut file 
# cut files
#cutFileOpt="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/Optimization/cutTable_lq_eejj_QCD_opt.txt"
# preselection only
#cutFileAna2016pre=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_QCD_singleFR_preselOnly.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_QCD_doubleFR_preselOnly.txt"
#)
#cutFileAna2016post=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_QCD_singleFR_preselOnly.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_QCD_doubleFR_preselOnly.txt"
#)
#cutFileAna2017=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/cutTable_lq_eejj_QCD_singleFR_preselOnly.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/cutTable_lq_eejj_QCD_doubleFR_preselOnly.txt"
#)
#cutFileAna2018=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/cutTable_lq_eejj_QCD_singleFR_preselOnly.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/cutTable_lq_eejj_QCD_doubleFR_preselOnly.txt"
#)

# final selections: LQToDEle
#
#cutFileAna2016pre=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/HTLO-amcatnlo/cutTable_lq_eejj_QCD_singleFR.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/HTLO-amcatnlo/cutTable_lq_eejj_QCD_doubleFR.txt"
#)
#cutFileAna2016post=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/HTLO-amcatnlo/cutTable_lq_eejj_QCD_singleFR.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/HTLO-amcatnlo/cutTable_lq_eejj_QCD_doubleFR.txt"
#)
#cutFileAna2017=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_QCD_singleFR.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_QCD_doubleFR.txt"
#)
#cutFileAna2018=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_QCD_singleFR.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_QCD_doubleFR.txt"
#)

# final selections: LQToBEle
#
cutFileAna2016pre=(
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/LQToBEle/cutTable_lq_eejj_QCD_singleFR.txt"
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/LQToBEle/cutTable_lq_eejj_QCD_doubleFR.txt"
)
cutFileAna2016post=(
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/LQToBEle/cutTable_lq_eejj_QCD_singleFR.txt"
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/LQToBEle/cutTable_lq_eejj_QCD_doubleFR.txt"
)
cutFileAna2017=(
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/LQToBEle/cutTable_lq_eejj_QCD_singleFR.txt"
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/LQToBEle/cutTable_lq_eejj_QCD_doubleFR.txt"
)
cutFileAna2018=(
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/LQToBEle/cutTable_lq_eejj_QCD_singleFR.txt"
  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/LQToBEle/cutTable_lq_eejj_QCD_doubleFR.txt"
)


#cutFileAna2016pre=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_QCD_singleFR_dedicatedMasses.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_QCD_doubleFR_dedicatedMasses.txt"
#)
#cutFileAna2016post=(
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_QCD_singleFR_dedicatedMasses.txt"
#  "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_QCD_doubleFR_dedicatedMasses.txt"
#)

# ilumi
ilumi2016pre=19497.897120
ilumi2016post=16812.151722
ilumi2017=41477.877399
ilumi2018=59827.449483
QUEUEOPT=testmatch
#QUEUEANA=workday
QUEUEANA=tomorrow
#------------
if [ "$YEAR" = "2016preVFP" ]; then
  echo "Doing 2016preVFP!"
  ILUM=$ilumi2016pre
  INPUTLIST1=$inputlist2016pre_1FR
  INPUTLIST2=$inputlist2016pre_2FR
  cutFileAna=( "${cutFileAna2016pre[@]}" )
  XSECTION=$xsection2016pre
elif [ "$YEAR" = "2016postVFP" ]; then
  echo "Doing 2016postVFP!"
  ILUM=$ilumi2016post
  INPUTLIST1=$inputlist2016post_1FR
  INPUTLIST2=$inputlist2016post_2FR
  cutFileAna=( "${cutFileAna2016post[@]}" )
  XSECTION=$xsection2016post
elif [ "$YEAR" = "2017" ]; then
  echo "Doing 2017!"
  ILUM=$ilumi2017
  INPUTLIST1=$inputlist2017_1FR
  INPUTLIST2=$inputlist2017_2FR
  cutFileAna=( "${cutFileAna2017[@]}" )
  XSECTION=$xsection2017
elif [ "$YEAR" = "2018" ]; then
  echo "Doing 2018!"
  ILUM=$ilumi2018
  INPUTLIST1=$inputlist2018_1FR
  INPUTLIST2=$inputlist2018_2FR
  cutFileAna=( "${cutFileAna2018[@]}" )
  XSECTION=$xsection2018
else
  die "year argument ${YEAR} not one of 2016preVFP, 2016postVFP, 2017, 2018"
fi
#------------
if [ "$OPT" = "1" ]; then
  DIRSTR="opt"
  files=( "${cutFileOpt[@]}" )
  queue=$QUEUEOPT
else
  DIRSTR="analysis"
  files=( "${cutFileAna[@]}" )
  queue=$QUEUEANA
fi
SUBDIR=ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
EOSDIR=/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_QCD_${ANANAME}_${YEAR}_${DIRSTR}_batch_$(hostname -s).txt
SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_QCD_dataDriven_${YEAR}.yaml
#------------
FACTOR=1000 # numbers in final tables (but *not* in plots) will be multiplied by this scale factor (to see well the decimal digits)
#------------
EXE=build/mainQCD
#------------
MAXJOBS=50
#------------
if [ -z "$XSECTION" ]; then
  echo "No valid xsection specified"
  exit 1
fi
#------------
#### END OF INPUTS ####

echo "" > $COMMANDFILE

cat >> $COMMANDFILE <<EOF
####################################################
#### launch, check and combine cmds for QCD 1FR+2FR####
EOF

for i in "${!files[@]}"
do
  file=${files[$i]}
  suffix=$(basename $file)
  j=$((i+1))
suffix=${suffix%\.*}
inputlist="INPUTLIST$j"
cat >> $COMMANDFILE <<EOF
python scripts/launchAnalysis_batch_ForSkimToEOS.py -i ${!inputlist} -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree -e $EXE
EOF
done

cat >> $COMMANDFILE <<EOF

####################################################
#### check
EOF

for file in "${files[@]}"
do
  suffix=$(basename $file)
suffix=${suffix%\.*}
cat >> $COMMANDFILE <<EOF
./scripts/checkJobs.sh $OUTDIRPATH/$SUBDIR/condor_$suffix $EOSDIR/$suffix/condor
EOF
done

cat >> $COMMANDFILE <<EOF

####################################################
#### combine
EOF

for i in "${!files[@]}"
do
  file=${files[$i]}
  suffix=$(basename $file)
  j=$((i+1))
suffix=${suffix%\.*}
inputlist="INPUTLIST$j"
cat >> $COMMANDFILE <<EOF
mkdir -p $OUTDIRPATH/$SUBDIR/output_$suffix \
&& time  ./scripts/combinePlots.py \
    -i ${!inputlist} \
    -c $CODENAME \
    -d $EOSDIR/$suffix/condor \
    -l  $(echo "$ILUM*$FACTOR" | bc) \
    -x $XSECTION  \
    -o $EOSDIR/output_$suffix \
    -s $SAMPLELISTFORMERGING \
    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log
EOF
done


# assumes singleFB cutfile is always the first one!
suffix1=$(basename ${files[0]})
suffix1=${suffix1%\.*}
suffix2=$(basename ${files[1]})
suffix2=${suffix2%\.*}
cat >> $COMMANDFILE <<EOF
python scripts/makeQCDYield.py -s $EOSDIR/output_$suffix1 \
    -d $EOSDIR/output_$suffix2 \
    -o $EOSDIR/output_qcdSubtractedYield
EOF

echo "The set of commands to run on the cut files:" 
for file in "${files[@]}"
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
cat $COMMANDFILE
echo "Submitting jobs"
for i in "${!files[@]}"
do
  file=${files[$i]}
  suffix=$(basename $file)
  j=$((i+1))
suffix=${suffix%\.*}
inputlist="INPUTLIST$j"
python scripts/launchAnalysis_batch_ForSkimToEOS.py -i ${!inputlist} -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree -e $EXE
done
echo "Reminder: The set of commands to run on the cut files:" 
for file in "${files[@]}"
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
cat $COMMANDFILE
