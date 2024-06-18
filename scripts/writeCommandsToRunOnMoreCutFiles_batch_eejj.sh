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
#ANANAME=eejj_8mar2024_dedicatedMassBDTs
#ANANAME=eejj_7may2024_preselectionOnly
#ANANAME=eejj_8may2024_dedicatedMassBDTsRedoAprSkims
ANANAME=eejj_24may2024_dedicatedMassBDTs_LQToBEle
#------------
#
#inputlist2016pre=config/inputListsAnalysisPreselSkim_UL16preVFP_4apr2024/inputListAllCurrent.txt
#inputlist2016post=config/inputListsAnalysisPreselSkim_UL16postVFP_4apr2024/inputListAllCurrent.txt
inputlist2016pre=config/inputListsAnalysisPreselSkim_UL16preVFP_7may2024/inputListAllCurrent.txt
inputlist2016post=config/inputListsAnalysisPreselSkim_UL16postVFP_7may2024/inputListAllCurrent.txt
inputlist2017=config/inputListsAnalysisPreselSkim_UL17_7may2024/inputListAllCurrent.txt
inputlist2018=config/inputListsAnalysisPreselSkim_UL18_7may2024/inputListAllCurrent.txt
#------------
#xsection2016pre=config/xsection_13TeV_2022.txt
#xsection2016post=config/xsection_13TeV_2022.txt
#xsection2017=config/xsection_13TeV_2022.txt
#xsection2018=config/xsection_13TeV_2022.txt
xsection2016pre=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016preVFP/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_may7_qcd_dyjAMCatNLO.txt
xsection2016post=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2016postVFP/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_may7_qcd_dyjAMCatNLO.txt
xsection2017=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2017/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_may7_qcd_dyjAMCatNLO.txt
xsection2018=/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/ultralegacy/rescaledCrossSections/2018/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_may7_qcd_dyjAMCatNLO.txt
#------------
CODENAME=analysisClass_lq_eejj
#CODENAME=analysisClass_lq_eejj_oneBTag
#------------
OUTDIRPATH=$LQDATA  # a subdir will be created for each cut file 
## cut files - preselection only
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
#
#cutFileAna2016="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/${YEAR#2016}/cutTable_lq_eejj_looserPresel.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/cutTable_lq_eejj_loosePresel.txt"
#cutFileAna2016="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/${YEAR#2016}/cutTable_lq_eejj_loosePresel.txt"
#cutFileAna2016="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/${YEAR#2016}/cutTable_lq_eejj_loosePresel_noJetRequirements.txt"
#cutFileAna2016="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/${YEAR#2016}/cutTable_lq_eejj_loosePresel_noJetRequirements_noPtEE.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR%%p*}/Analysis/cutTable_lq_eejj.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/Analysis/cutTable_lq_eejj.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/cutTable_lq_eejj_BDT1400.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/cutTable_lq_eejj_BDT1500.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/cutTable_lq_eejj_BDT1600.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/cutTable_lq_eejj_BDT1700.txt"
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/cutTable_lq_eejj_BDT_parametrized.txt"
#
#cutFileAna="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Analysis/cutTable_lq_eejj_MasymTest.txt"
cutFileOpt="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/Optimization/cutTable_lq_eejj_BDT_opt.txt"
#cutFileOpt="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/Optimization/cutTable_lq_eejj_opt.txt"
#cutFileOpt="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/Optimization/cutTable_lq_eejj_oneBTag_opt.txt"
#cutFileOpt="/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/Optimization/cutTable_lq_eejj_twoBTags_opt.txt"
# ilumi
ilumi2016pre=19497.897120
ilumi2016post=16812.151722
ilumi2017=41477.877399
ilumi2018=59827.449483
excludeCombining=""
#------------
QUEUEOPT=testmatch
QUEUEANA=tomorrow
#QUEUEANA=workday
#------------
if [ "$YEAR" = "2016preVFP" ]; then
  echo "Doing 2016preVFP!"
  ILUM=$ilumi2016pre
  INPUTLIST=$inputlist2016pre
  cutFileAna=$cutFileAna2016pre
  XSECTION=$xsection2016pre
elif [ "$YEAR" = "2016postVFP" ]; then
  echo "Doing 2016postVFP!"
  ILUM=$ilumi2016post
  INPUTLIST=$inputlist2016post
  cutFileAna=$cutFileAna2016post
  XSECTION=$xsection2016post
elif [ "$YEAR" = "2017" ]; then
  echo "Doing 2017!"
  ILUM=$ilumi2017
  INPUTLIST=$inputlist2017
  cutFileAna=$cutFileAna2017
  XSECTION=$xsection2017
elif [ "$YEAR" = "2018" ]; then
  echo "Doing 2018!"
  ILUM=$ilumi2018
  INPUTLIST=$inputlist2018
  cutFileAna=$cutFileAna2018
  XSECTION=$xsection2018
else
  die "year argument ${YEAR} not one of 2016preVFP, 2016postVFP, 2017, 2018"
fi
#------------
if [ "$OPT" = "1" ]; then
  DIRSTR="opt"
  files=$cutFileOpt
  queue=$QUEUEOPT
elif [[ "$YEAR" == *2016* ]]; then
  DIRSTR="analysis"
  files=$cutFileAna
  queue=$QUEUEANA
else
  DIRSTR="analysis"
  files=$cutFileAna
  queue=$QUEUEANA
fi
SUBDIR=ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
#EOSDIR=/eos/user/s/scooper/LQ/NanoV7/${YEAR}/$DIRSTR/$ANANAME
EOSDIR=/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/${DIRSTR}/${YEAR}/$ANANAME
COMMANDFILE=commandsToRunOnMoreCutFiles_eejj_${YEAR}_${DIRSTR}_batch_$(hostname -s).txt
SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_eejj_${YEAR}.yaml
#------------
FACTOR=1000 # numbers in final tables (but *not* in plots) will be multiplied by this scale factor (to see well the decimal digits)
#------------
EXE=build/main
#------------
MAXJOBS=10
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
#cat >> $COMMANDFILE <<EOF
#####################################################
##### launch, check and combine cmds for $suffix ####
## 1 job per dataset
#python scripts/launchAnalysis_batch_ForSkimToEOS.py -i $INPUTLIST -o $OUTDIRPATH/$SUBDIR/condor_$suffix -c $file -q $queue -d $EOSDIR/$suffix -j ${MAXJOBS} -n rootTupleTree/tree
#
#./scripts/checkJobs.sh $OUTDIRPATH/$SUBDIR/condor_$suffix $EOSDIR/$suffix/condor
#
#mkdir -p $OUTDIRPATH/$SUBDIR/output_$suffix \
#&& time  ./scripts/combinePlots.py \
#    -i $INPUTLIST \
#    -c $CODENAME \
#    -d $EOSDIR/$suffix/condor \
#    -l  $(echo "$ILUM*$FACTOR" | bc) \
#    -x $XSECTION  \
#    -o $EOSDIR/output_$suffix \
#    -s $SAMPLELISTFORMERGING \
#    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log \
#&& mv -v $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}_unscaled.log \
#&& mv -v $EOSDIR/output_$suffix/${CODENAME}_plots.root $EOSDIR/output_$suffix/${CODENAME}_plots_unscaled.root \
#&& mv -v $EOSDIR/output_$suffix/${CODENAME}_tables.dat $EOSDIR/output_$suffix/${CODENAME}_tables_unscaled.dat \
#&& mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME} && cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME} && python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots_unscaled.root ${YEAR} > makePlots.log \
#&& python $LQMACRO/plotting2016/calc_DYJetsAndTTBarRescale_And_xsecFile.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots_unscaled.root ${YEAR} > rescale.log \
#&& mkdir -p dyj_MGHTLO && cd dyj_MGHTLO \
#&& python $LQMACRO/plotting2016/calc_DYJetsAndTTBarRescale_And_xsecFile.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots_unscaled.root ${YEAR} ZJet_madgraphLO_HT > rescale.log \
#&& cd $LQANA \
#&& mkdir -p $OUTDIRPATH/$SUBDIR/output_$suffix \
#&& time  ./scripts/combinePlots.py \
#    -i $INPUTLIST \
#    -c $CODENAME \
#    -d $EOSDIR/$suffix/condor \
#    -l $(echo "$ILUM*$FACTOR" | bc) \
#    -x $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets.txt \
#    -o $EOSDIR/output_$suffix \
#    -s $SAMPLELISTFORMERGING \
#    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log \
#&& mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled && cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled && python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots.root ${YEAR} > makePlots.log \
#&& python $LQANA/scripts/makeCutflow.py $EOSDIR/output_$suffix/${CODENAME}_tables.dat > cutflows.txt
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

mkdir -p $OUTDIRPATH/$SUBDIR/output_$suffix \
&& time  ./scripts/combinePlots.py \
    -i $INPUTLIST \
    -c $CODENAME \
    -d $EOSDIR/$suffix/condor \
    -l  $(echo "$ILUM*$FACTOR" | bc) \
    -x $XSECTION  \
    -o $EOSDIR/output_$suffix \
    -s $SAMPLELISTFORMERGING \
    | tee $OUTDIRPATH/$SUBDIR/output_$suffix/combineTablesAndPlots_${suffix}.log \
&& mkdir -p $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled \
&& cd $LQANA/versionsOfAnalysis/${YEAR}/eejj/${ANANAME}/scaled \
&& python $LQMACRO/plotting2016/makeStackHistoTemplateV2_eejj.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots.root ${YEAR} > makePlots.log \
&& python $LQANA/scripts/makeCutflow.py $EOSDIR/output_$suffix/${CODENAME}_tables.dat > cutflows.txt \
&& python $LQANA/scripts/makeDatacard.py dummy.root $EOSDIR/output_$suffix/${CODENAME}_plots.root ${YEAR} > makeDatacard.log
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
