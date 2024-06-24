#!/bin/bash

declare -a years=("2016preVFP" "2016postVFP" "2017" "2018")
#declare -a years=("2016preVFP" "2016postVFP")

declare -A inputLists
inputLists[2016preVFP]="config/inputListsNanoAOD_nanoV9_UL2016_preVFP/inputListAllCurrent.txt"
inputLists[2016postVFP]="config/inputListsNanoAOD_nanoV9_UL2016_postVFP/inputListAllCurrent.txt"
inputLists[2017]="config/inputListsNanoAOD_nanoV9_UL2017/inputListAllCurrent.txt"
inputLists[2018]="config/inputListsNanoAOD_nanoV9_UL2018/inputListAllCurrent.txt"
declare -A cutFilesQCD
cutFilesQCD[2016preVFP]="${LQMACRO}/config2016/ReducedSkims/preVFP/cutTable_lq1_skim_QCD_heep.txt"
cutFilesQCD[2016postVFP]="${LQMACRO}/config2016/ReducedSkims/postVFP/cutTable_lq1_skim_QCD_heep.txt"
cutFilesQCD[2017]="${LQMACRO}/config2017/ReducedSkims/cutTable_lq1_skim_QCD_heep.txt"
cutFilesQCD[2018]="${LQMACRO}/config2018/ReducedSkims/cutTable_lq1_skim_QCD_heep.txt"
declare -A cutFilesDoubleEle
cutFilesDoubleEle[2016preVFP]="${LQMACRO}/config2016/ReducedSkims/preVFP/cutTable_lq1_skim_DoubleEle_heep.txt"
cutFilesDoubleEle[2016postVFP]="${LQMACRO}/config2016/ReducedSkims/postVFP/cutTable_lq1_skim_DoubleEle_heep.txt"
cutFilesDoubleEle[2017]="${LQMACRO}/config2017/ReducedSkims/cutTable_lq1_skim_DoubleEle_heep.txt"
cutFilesDoubleEle[2018]="${LQMACRO}/config2018/ReducedSkims/cutTable_lq1_skim_DoubleEle_heep.txt"
skimNameSuffix=$(date +'%d%b%y'| tr A-Z a-z)
#skimNameSuffix=08feb24
nMaxJobsDEle=100
nMaxJobsQCD=200
eosStorageDir="/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/skims"
localStorageDir="${LQDATA}/ultralegacy/skims"

#
# QCD
#
for year in "${years[@]}" ; do
    inputList=${inputLists[$year]}
    cutFile=${cutFilesQCD[$year]}
    echo "python3 scripts/launchAnalysis_batch_ForSkimToEOS.py -j ${nMaxJobsQCD} -i ${inputList} -o ${localStorageDir}/${year}/rskQCD_heep_${skimNameSuffix} -n Events -e build/main -d ${eosStorageDir}/${year}/rskQCD_heep_${skimNameSuffix} -c ${cutFile} -r"
    python3 scripts/launchAnalysis_batch_ForSkimToEOS.py -j ${nMaxJobsQCD} -i ${inputList} -o ${localStorageDir}/${year}/rskQCD_heep_${skimNameSuffix} -n Events -e build/main -d ${eosStorageDir}/${year}/rskQCD_heep_${skimNameSuffix} -c ${cutFile} -r
    retVal=$?
    if [ $retVal -ne 0 ]; then
        echo "Submission command for QCD RSK for year=${year} and cutFile=${cutFile} exited with return code ${retVal}. Stopping here."
        exit $retVal
    fi
done

#
# DoubleEle
#
for year in "${years[@]}" ; do
    inputList=${inputLists[$year]}
    cutFile=${cutFilesDoubleEle[$year]}
    echo "python3 scripts/launchAnalysis_batch_ForSkimToEOS.py -j ${nMaxJobsDEle} -i ${inputList} -o ${localStorageDir}/${year}/rskDoubleEle_heep_${skimNameSuffix} -n Events -e build/main -d ${eosStorageDir}/${year}/rskDoubleEle_heep_${skimNameSuffix} -c ${cutFile} -r"
    python3 scripts/launchAnalysis_batch_ForSkimToEOS.py -j ${nMaxJobsDEle} -i ${inputList} -o ${localStorageDir}/${year}/rskDoubleEle_heep_${skimNameSuffix} -n Events -e build/main -d ${eosStorageDir}/${year}/rskDoubleEle_heep_${skimNameSuffix} -c ${cutFile} -r
    retVal=$?
    if [ $retVal -ne 0 ]; then
        echo "Submission command for DoubleEle RSK for year=${year} and cutFile=${cutFile} exited with return code ${retVal}. Stopping here."
        exit $retVal
    fi
done
