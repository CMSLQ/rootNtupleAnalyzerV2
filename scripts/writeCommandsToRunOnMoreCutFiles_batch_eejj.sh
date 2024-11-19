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

[ $# -gt 0 ] || die "must specify year as argument; usage: writeCommandsToRunOnMoreCutFiles_batch_eejj.sh year [doOptimization]"
YEAR=$1

if [[ "$#" -gt 1 ]]; then
  OPT=1
else
  OPT=0
fi

#### INPUTS HERE ####
#------------
#ANANAME=eejj_13nov2024_bdt_LQToDEle
#ANANAME=eejj_11nov2024_presel
ANANAME=eejj_15nov2024_bdt_LQToBEle
#------------
QCDANANAME=qcd_${ANANAME}
#------------
#inputlist2016pre=config/inputListsRSKDoubleEle_heep_UL16preVFP_24oct24/inputListAllCurrent.txt
#inputlist2016post=config/inputListsRSKDoubleEle_heep_UL16postVFP_24oct24/inputListAllCurrent.txt
#inputlist2017=config/inputListsRSKDoubleEle_heep_UL17_24oct24/inputListAllCurrent.txt
#inputlist2018=config/inputListsRSKDoubleEle_heep_UL18_24oct24/inputListAllCurrent.txt
#
#inputlist2016pre=config/inputListsAnalysisPreselSkim_UL16preVFP_28oct2024/inputListAllCurrent.txt
#inputlist2016post=config/inputListsAnalysisPreselSkim_UL16postVFP_28oct2024/inputListAllCurrent.txt
#inputlist2017=config/inputListsAnalysisPreselSkim_UL17_28oct2024/inputListAllCurrent.txt
#inputlist2018=config/inputListsAnalysisPreselSkim_UL18_28oct2024/inputListAllCurrent.txt
#
inputlist2016pre=config/inputListsAnalysisPreselSkim_UL16preVFP_11nov2024/inputListAllCurrent.txt
inputlist2016post=config/inputListsAnalysisPreselSkim_UL16postVFP_11nov2024/inputListAllCurrent.txt
inputlist2017=config/inputListsAnalysisPreselSkim_UL17_11nov2024/inputListAllCurrent.txt
inputlist2018=config/inputListsAnalysisPreselSkim_UL18_11nov2024/inputListAllCurrent.txt
#------------
#xsection2016pre=$LQANA/config/xsection_13TeV_2022.txt
#xsection2016post=$LQANA/config/xsection_13TeV_2022.txt
#xsection2017=$LQANA/config/xsection_13TeV_2022.txt
#xsection2018=$LQANA/config/xsection_13TeV_2022.txt
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
CODENAME=analysisClass_lq_eejj
#CODENAME=analysisClass_lq_eejj_oneBTag
#------------
OUTDIRPATH=$LQDATA  # a subdir will be created for each cut file 
# cut files - preselection only
#cutFileAna2016pre="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_preselOnly.txt"
#cutFileAna2016post="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_preselOnly.txt"
#cutFileAna2017="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/cutTable_lq_eejj_preselOnly.txt"
#cutFileAna2018="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/cutTable_lq_eejj_preselOnly.txt"
# OLD final selections
#cutFileAna2016pre="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/cutTable_lq_eejj_BDT.txt"
#cutFileAna2016post="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/cutTable_lq_eejj_BDT.txt"
#cutFileAna2017="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/cutTable_lq_eejj_BDT.txt"
#cutFileAna2018="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/cutTable_lq_eejj_BDT.txt"
# LQToDEle final selections
#cutFileAna2016pre="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/HTLO-amcatnlo/cutTable_lq_eejj_BDT.txt"
#cutFileAna2016post="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/HTLO-amcatnlo/cutTable_lq_eejj_BDT.txt"
#cutFileAna2017="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_BDT.txt"
#cutFileAna2018="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/HTLO-amcatnlo/cutTable_lq_eejj_BDT.txt"
# LQToBEle final selections
cutFileAna2016pre="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/preVFP/LQToBEle/cutTable_lq_eejj_BDT.txt"
cutFileAna2016post="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/postVFP/LQToBEle/cutTable_lq_eejj_BDT.txt"
cutFileAna2017="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2017/Analysis/LQToBEle/cutTable_lq_eejj_BDT.txt"
cutFileAna2018="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2018/Analysis/LQToBEle/cutTable_lq_eejj_BDT.txt"
#------------
# ilumi
ilumi2016pre=19497.897120
ilumi2016post=16812.151722
ilumi2017=41477.877399
ilumi2018=59827.449483
excludeCombining=""
#------------
QUEUEANA=tomorrow
#QUEUEANA=workday
#------------
if [ "$YEAR" = "2016preVFP" ]; then
  echo "Doing 2016preVFP!"
  ILUM=$ilumi2016pre
  INPUTLIST=$LQANA/$inputlist2016pre
  cutFileAna=$cutFileAna2016pre
  XSECTION=$xsection2016pre
elif [ "$YEAR" = "2016postVFP" ]; then
  echo "Doing 2016postVFP!"
  ILUM=$ilumi2016post
  INPUTLIST=$LQANA/$inputlist2016post
  cutFileAna=$cutFileAna2016post
  XSECTION=$xsection2016post
elif [ "$YEAR" = "2017" ]; then
  echo "Doing 2017!"
  ILUM=$ilumi2017
  INPUTLIST=$LQANA/$inputlist2017
  cutFileAna=$cutFileAna2017
  XSECTION=$xsection2017
elif [ "$YEAR" = "2018" ]; then
  echo "Doing 2018!"
  ILUM=$ilumi2018
  INPUTLIST=$LQANA/$inputlist2018
  cutFileAna=$cutFileAna2018
  XSECTION=$xsection2018
else
  die "year argument ${YEAR} not one of 2016preVFP, 2016postVFP, 2017, 2018"
fi
#------------
DIRSTR="analysis"
files=$cutFileAna
queue=$QUEUEANA
SUBDIR=ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
#EOSDIR=/eos/user/s/scooper/LQ/NanoV7/${YEAR}/$DIRSTR/$ANANAME
EOSDIR=/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
QCDFILE=/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/${DIRSTR}/${YEAR}/${QCDANANAME}/output_qcdSubtractedYield/qcdSubtracted_plots.root
COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_${ANANAME}_${YEAR}_${DIRSTR}_batch_$(hostname -s).txt
SAMPLELISTFORMERGING=$LQANA/config/sampleListForMerging_13TeV_eejj_${YEAR}.yaml
#------------
FACTOR=1000 # numbers in final tables (but *not* in plots) will be multiplied by this scale factor (to see well the decimal digits)
#------------
EXE=build/main
#------------
MAXJOBS=5
#------------
if [ -z "$XSECTION" ]; then
  echo "No valid xsection specified"
  exit 1
fi
#------------
#### END OF INPUTS ####

echo "" > $COMMANDFILE

for file in $files
do
  suffix=$(basename $file)
suffix=${suffix%\.*}
####################################################################################################
# preselection only
####################################################################################################
#cat >> $COMMANDFILE <<EOF
#####################################################
##### launch, check and combine cmds for $suffix ####
## 1 job per dataset
#python scripts/launchAnalysis_batch_ForSkimToEOS.py -i $INPUTLIST -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree
#
#./scripts/checkJobs.sh $OUTDIRPATH/$SUBDIR/condor_$suffix $EOSDIR/$suffix/condor
#
#mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME} && cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME} \
#&& $LQANA/scripts/combinePlotsBatchDriver.py \
#    -i $INPUTLIST \
#    -c $CODENAME \
#    -d $EOSDIR/$suffix/condor/ \
#    -l  $(echo "$ILUM*$FACTOR" | bc) \
#    -x $XSECTION  \
#    -o $EOSDIR/output_$suffix/unscaled \
#    -s $SAMPLELISTFORMERGING
#
#find combinePlotsCondor/error/ -size +0 -type f -print0 -exec echo " had errors; content follows" \; -exec cat {} \;
#
#python $LQANA/scripts/combineAndSortTables.py $EOSDIR/output_$suffix/unscaled/ $SAMPLELISTFORMERGING \
#&& python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py ${QCDFILE} $EOSDIR/output_$suffix/unscaled ${YEAR} > makePlots.log \
#&& python $LQMACRO/plotting2016/calc_DYJetsAndTTBarRescale_And_xsecFile.py ${QCDFILE} $EOSDIR/output_$suffix/unscaled ${YEAR} > rescale.log \
#&& mkdir -p dyj_MGHTLO && cd dyj_MGHTLO \
#&& python $LQMACRO/plotting2016/calc_DYJetsAndTTBarRescale_And_xsecFile.py ${QCDFILE} $EOSDIR/output_$suffix/unscaled ${YEAR} ZJet_madgraphLO_HT > rescale.log \
#
#mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled && cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled \
#&& $LQANA/scripts/combinePlotsBatchDriver.py \
#    -i $INPUTLIST \
#    -c $CODENAME \
#    -d $EOSDIR/$suffix/condor \
#    -l $(echo "$ILUM*$FACTOR" | bc) \
#    -x $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets.txt \
#    -o $EOSDIR/output_$suffix \
#    -s $SAMPLELISTFORMERGING
#
#find combinePlotsCondor/error/ -size +0 -type f -print0 -exec echo " had errors; content follows" \; -exec cat {} \;
#
#python $LQANA/scripts/combineAndSortTables.py $EOSDIR/output_$suffix/ $SAMPLELISTFORMERGING \
#&& python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py ${QCDFILE} $EOSDIR/output_$suffix ${YEAR} > makePlots.log \
#&& python $LQANA/scripts/makeCutflow.py combinedTable.dat > cutflows.txt
#EOF
####################################################################################################
# final selections only
####################################################################################################
cat >> $COMMANDFILE <<EOF
####################################################
#### launch, check and combine cmds for $suffix ####
# 1 job per dataset
python scripts/launchAnalysis_batch_ForSkimToEOS.py -i $INPUTLIST -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree

./scripts/checkJobs.sh $OUTDIRPATH/$SUBDIR/condor_$suffix $EOSDIR/$suffix/condor

mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled && cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled \
&& $LQANA/scripts/combinePlotsBatchDriver.py \
    -i $INPUTLIST \
    -c $CODENAME \
    -d $EOSDIR/$suffix/condor \
    -l  $(echo "$ILUM*$FACTOR" | bc) \
    -x $XSECTION  \
    -o $EOSDIR/output_$suffix \
    -s $SAMPLELISTFORMERGING

find combinePlotsCondor/error/ -size +0 -type f -print0 -exec echo " had errors; content follows" \; -exec cat {} \;

python $LQANA/scripts/combineAndSortTables.py $EOSDIR/output_$suffix/ $SAMPLELISTFORMERGING \
&& python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py ${QCDFILE} $EOSDIR/output_$suffix ${YEAR} > makePlots.log \
&& python $LQANA/scripts/makeCutflow.py combinedTable.dat > cutflows.txt \
&& python $LQANA/scripts/makeDatacard.py ${QCDFILE} $EOSDIR/output_$suffix ${YEAR} > makeDatacard.log
EOF

done


echo "The set of commands to run on the cut files:" 
for file in $files
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
cat $COMMANDFILE
echo "Submitting jobs"
python scripts/launchAnalysis_batch_ForSkimToEOS.py -i $INPUTLIST -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree
echo "Reminder: The set of commands to run on the cut files:" 
for file in $files
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
cat $COMMANDFILE
