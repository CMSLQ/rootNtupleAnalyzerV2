#!/bin/bash

#Script loops over all specified years and runs the specified steps of the BDT training/optimization/plots, then copies the dataset directory to EOS if one was created.
makeInputLists=false
makeTrainingTrees=false
doTraining=false
doOptimization=true
doBDTPlots=false
years='2016preVFP 2016postVFP 2017 2018' # '2016preVFP 2016postVFP 2017 2018'
skimDate='7feb'
source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh

if [ "$makeInputLists" = true ]; then
for year in $years; do
	echo "Make input lists for $year"
	echo "Preselection"
	python scripts/createInputLists.py -i /eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/$year/eejj_7feb2025_mjjTest_presel/cutTable_lq_eejj_preselOnly/skim -o config/myDatasets/BDT/$year/$skimDate/trainingTreeInputs/preselOnly
	echo "Single FR"
	python scripts/createInputLists.py -i /eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/$year/qcd_eejj_7feb2025_mjjTest_presel/cutTable_lq_eejj_QCD_singleFR_preselOnly/skim -o config/myDatasets/BDT/$year/$skimDate/trainingTreeInputs/singleFR
	echo "Double FR"
	python scripts/createInputLists.py -i /eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/$year/qcd_eejj_7feb2025_mjjTest_presel/cutTable_lq_eejj_QCD_doubleFR_preselOnly/skim -o config/myDatasets/BDT/$year/$skimDate/trainingTreeInputs/doubleFR
done
fi

if [ "$makeTrainingTrees" = true ]; then
for year in $years; do
	echo "Make training trees for $year"
	python scripts/makeBDTTrainingTrees.py $year
done
fi

years='2016preVFP 2016postVFP 2017 2018'
destDir=$LQDATAEOS/BDT_7FebSkim/LQToBEle/minNB0p5
#destDir=testingBDTs
#rm -r dataset
#cp -r $LQDATAEOS/BDT_7FebSkim/LQToDEle/dataset .

if [ ! -d $destDir ]; then
    echo "EOS directory does not exist. Making one"
    mkdir $destDir
fi

#if [ "$doTraining" = true ]; then
#for year in $years; do
#	echo "Make input lists for $year"
#	./scripts/createInputListsBDT.sh config/myDatasets/BDT/$year/7FebSkim/tmvaInputs $LQDATAEOS/BDTTrainingTrees/LQToDEle/$year/7FebSkims
#	./scripts/createInputListsBDT.sh config/myDatasets/BDT/$year/7FebSkim/tmvaInputsLQToBEle $LQDATAEOS/BDTTrainingTrees/LQToBEle/$year/7FebSkims
#done
#fi

#if [ "$doTraining" = true ]; then
#echo "removing existing directory dataset"
#rm -r dataset
#fi 

if [ "$doTraining" = true ]; then 
echo "train BDT"
python scripts/tmvaBDT.py -t -d $destDir | tee bdtTraining.log
cp bdtTraining.log $destDir
echo "*************************************************************************"
fi

if [ "$doOptimization" = true ]; then
echo "optimize BDT"
python scripts/tmvaBDT.py -o -d $destDir | tee bdtOptimization.log
mv bdtOptimization.log $destDir
echo "*************************************************************************"
fi 

if [ "$doBDTPlots" = true ]; then
echo "do ROC and BDT plots"
python scripts/tmvaBDT.py -r -d $destDir | tee bdtPlots.log
mv bdtPlots.root bdtPlots.log $destDir
echo "*************************************************************************"
fi

if [ "$doTraining" = true ]; then
echo "copy dataset directory to $destDir"
cp -r dataset $destDir
fi
