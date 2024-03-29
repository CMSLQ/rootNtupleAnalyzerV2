#!/bin/bash

YEAR=$1

if [ "$YEAR" = "2016preVFP" ]; then
  echo "Doing 2016preVFP!"
  SKIMNAME=rskSingleEleL_1dec2022
  INPUTLIST=config/inputListsNanoSkim_UL2016_preVFP_22nov2022_LQOnly/inputListAllCurrent.txt
  #INPUTLIST=config/inputListsNanoSkim_UL2016_preVFP_4oct2022/inputListAllCurrent.txt
  CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/preVFP/cutTable_lq1_noSkim_loose.txt
  #CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/preVFP/cutTable_lq1_skim_SingleEle_loose.txt
  #CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/preVFP/cutTable_lq1_skim_SingleEle_heep.txt
elif [ "$YEAR" = "2016postVFP" ]; then
  echo "Doing 2016postVFP!"
  SKIMNAME=rskSingleEleL_1dec2022
  INPUTLIST=config/inputListsNanoSkim_UL2016_postVFP_22nov2022_LQOnly/inputListAllCurrent.txt
  #INPUTLIST=config/inputListsNanoSkim_UL2016_postVFP_4oct2022/inputListAllCurrent.txt
  CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/postVFP/cutTable_lq1_noSkim_loose.txt
  #CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/postVFP/cutTable_lq1_skim_SingleEle_loose.txt
  #CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config2016/ReducedSkims/postVFP/cutTable_lq1_skim_SingleEle_heep.txt
elif [ "$YEAR" = "2017" ]; then
  SKIMNAME=rskSingleEleL_notDoneYet
  INPUTLIST=config/UL17NanoV9_nanoSkim_1jul2022/inputListAllCurrent.txt
  CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/ReducedSkims/cutTable_lq1_skim_SingleEle_loose.txt
elif [ "$YEAR" = "2018" ]; then
  SKIMNAME=rskSingleEleL_notDoneYet
  #INPUTLIST=config/nanoV7_2018_postProc/inputListAllCurrent.txt
  INPUTLIST=config/nanoV7_2018_postProc/inputList_jun10.txt
  CUTFILE=/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleMacrosV2/config${YEAR}/ReducedSkims/cutTable_lq1_skim_SingleEle_loose.txt
else
  echo "ERROR: did not understand given year of '$YEAR' which is not one of 2016preVFP, 2016postVFP, 2017, 2018"
  echo "Usage: $0 [2016preVFP | 2016postVFP | 2017 | 2018]"
  exit -1
fi

#EOSDIR=/eos/user/s/scooper/LQ/NanoV7/skims/${YEAR}/$SKIMNAME
#EOSDIR=/eos/cms/store/user/scooper/LQ/NanoV7/skims/${YEAR}/$SKIMNAME
EOSDIR=/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/skims/${YEAR}/$SKIMNAME
OUTPUTDIR=/afs/cern.ch/user/s/scooper/work/private/data/Leptoquarks/ultralegacy/skims/${YEAR}/$SKIMNAME

#python scripts/launchAnalysis_batch_ForSkimToEOS.py -j 20 -q tomorrow -i $INPUTLIST -o $OUTPUTDIR -n Events -c $CUTFILE -d $EOSDIR -r
#python scripts/launchAnalysis_batch_ForSkimToEOS.py -j 50 -q tomorrow -i $INPUTLIST -o $OUTPUTDIR -n rootTupleTree/tree -c $CUTFILE -d $EOSDIR -r
#python scripts/launchAnalysis_batch_ForSkimToEOS.py -j -1 -q workday -i $INPUTLIST -o $OUTPUTDIR -n rootTupleTree/tree -c $CUTFILE -d $EOSDIR -r
python scripts/launchAnalysis_batch_ForSkimToEOS.py -j 50 -q tomorrow -i $INPUTLIST -o $OUTPUTDIR -n rootTupleTree/tree -c $CUTFILE -d $EOSDIR -r
# below for nanoSkims
#python3 scripts/launchAnalysis_batch_ForSkimToEOS.py -j 50 -q tomorrow -i $INPUTLIST -o $OUTPUTDIR -n Evemts -c $CUTFILE -d $EOSDIR -s
