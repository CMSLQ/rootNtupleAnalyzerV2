#!/bin/bash

# year = 2016preVFP/2016postVFP/2017/2018
# skimtype = QCD/DoubleEle
# date = 24oct24
# prevDate = 13sep24

year=$1
skimtype=$2
date=$3
prevDate=$4

nLinesToDiff=94

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

dirBase="/home/scooper/data/LQ/ultralegacy/skims/${year}/rsk${skimtype}_heep_${date}/"
prevDirBase="/home/scooper/data/LQ/ultralegacy/skims/${year}/rsk${skimtype}_heep_${prevDate}/"

for dataset in `cat ${inputList}`; do
    f=${dataset##*/}
    dName="${f%.*}"
    echo "validating dat files..."
    datFile1match=${dirBase}"analysisClass_lq1_skim___${dName}/output/analysisClass_lq1_skim___${dName}*.dat"
    prevFileDir=${prevDirBase}"analysisClass_lq1_skim___${dName}/output/"
    for datFile in `/usr/bin/ls ${datFile1match}`; do
        datFileBase=${datFile##*/}
        prevDatFile=${prevFileDir}${datFileBase}
        # echo "diff -q ${datFile} ${prevDatFile}"
        # diff -q ${datFile} ${prevDatFile}
        # diff -q <(head -n ${nLinesToDiff} ${datFile}) <(head -n 1 ${prevDatFile})
        if ! diff -q <(head -n ${nLinesToDiff} ${datFile}) <(head -n ${nLinesToDiff} ${prevDatFile}); then
            echo "WARNING: first ${nLinesToDiff} lines of file ${datFile} differs with previous"
            echo -e "\tcheck with vimdiff ${datFile} ${prevDatFile}"
        fi
    done
done

# for dataset in `cat ${inputList}`; do
#     f=${dataset##*/}
#     dName="${f%.*}"
#     echo "validating root files..."
#     rootFile1match=${dirBase}"analysisClass_lq1_skim___${dName}/output/analysisClass_lq1_skim___${dName}*.root"
#     prevFileDir=${prevDirBase}"analysisClass_lq1_skim___${dName}/output/"
#     for rootFile in `/usr/bin/ls ${rootFile1match}`; do
#         rootFileBase=${rootFile##*/}
#         prevrootFile=${prevFileDir}${rootFileBase}
#         python3 scripts/validateRootPlotFiles.py ${rootFile} ${prevrootFile}
#     done
# done
