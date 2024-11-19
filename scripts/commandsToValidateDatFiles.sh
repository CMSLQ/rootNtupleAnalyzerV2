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

inputListBase="config/inputListsRSKDoubleEle_heep_"
if [ ${year} == "2016preVFP" ]; then
    inputListBase=$inputListBase"UL16preVFP_24oct24/"
elif [ ${year} == "2016postVFP" ]; then
    inputListBase=$inputListBase"UL16postVFP_24oct24/"
elif [ ${year} == "2017" ]; then
    inputListBase=$inputListBase"UL17_24oct24/"
elif [ ${year} == "2018" ]; then
    inputListBase=$inputListBase"UL18_24oct24/"
else
    echo "Did not understand year given as ${year}; must be one of 2016preVFP/2016postVFP/2017/2018"
    exit -1
fi

inputList=$inputListBase"inputListAllCurrent.txt"
echo "using inputList=${inputList}"

# LQToDEle_M-1000_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8_APV/LQToDEle_M-1000_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8_APV_0.dat"
dirBase="/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/${year}/eejj_${date}_bdt_LQToDEle/cutTable_lq_eejj_preselOnly/condor/"
# prevDirBase="/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/${year}/eejj_${prevDate}_bdt_LQToDEle/cutTable_lq_eejj_preselOnly/condor/"
prevDirBase="/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/${year}/eejj_${prevDate}_presel/cutTable_lq_eejj_preselOnly/condor/"

for dataset in `cat ${inputList}`; do
    f=${dataset##*/}
    dName="${f%.*}"
    echo "validating dat files..."
    datFile1match=${dirBase}"${dName}/${dName}*.dat"
    prevFileDir=${prevDirBase}"${dName}/"
    for datFile in `/usr/bin/ls ${datFile1match}`; do
        datFileBase=${datFile##*/}
        prevDatFile=${prevFileDir}${datFileBase}
        # echo "diff -q ${datFile} ${prevDatFile}"
        # diff -q ${datFile} ${prevDatFile}
        # diff -q <(head -n ${nLinesToDiff} ${datFile}) <(head -n 1 ${prevDatFile})
        # if ! diff -q <(head -n ${nLinesToDiff} ${datFile}) <(head -n ${nLinesToDiff} ${prevDatFile}); then
        if ! diff -q ${datFile} ${prevDatFile}; then
            echo "WARNING: file ${datFile} differs with previous"
            echo -e "\tcheck with vimdiff ${datFile} ${prevDatFile}"
        fi
    done
done

# for dataset in `cat ${inputList}`; do
#     f=${dataset##*/}
#     dName="${f%.*}"
#     echo "validating root files..."
#     rootFile1match=${dirBase}"${dName}/output/analysisClass_lq1_skim___${dName}*.root"
#     prevFileDir=${prevDirBase}"${dName}/output/"
#     for rootFile in `/usr/bin/ls ${rootFile1match}`; do
#         rootFileBase=${rootFile##*/}
#         prevrootFile=${prevFileDir}${rootFileBase}
#         python3 scripts/validateRootPlotFiles.py ${rootFile} ${prevrootFile}
#     done
# done
