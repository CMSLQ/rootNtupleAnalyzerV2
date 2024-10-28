#!/bin/bash

# year = 2016preVFP/2016postVFP/2017/2018
# skimtype = QCD/DoubleEle
# date = 24oct24

year=$1
skimtype=$2
date=$3

inputListBase="config/inputListsNanoAOD_nanoV9_reNanoST_"
if [ ${year} == "2016preVFP" ]; then
    inputListBase=$inputListBase"UL2016_preVFP/"
elif [ ${year} == "2016postVFP" ]; then
    inputListBase=$inputListBase"UL2016_postVFP/"
elif [ ${year} == "2017" ]; then
    inputListBase=$inputListBase"UL2017/"
elif [ ${year} == "2018" ]; then
    inputListBase=$inputListBase"UL2018/"
else
    echo "Did not understand year given as ${year}; must be one of 2016preVFP/2016postVFP/2017/2018"
    exit -1
fi

inputList=$inputListBase"inputListAllCurrent.txt"
echo "using inputList=${inputList}"

python3 scripts/checkJobs.py ${inputList} /home/scooper/data/LQ/ultralegacy/skims/${year}/rsk${skimtype}_heep_${date}/ /eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/skims/${year}/rsk${skimtype}_heep_${date}/skim/
