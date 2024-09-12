#!/bin/bash

#Script loops over all specified years and runs the specified steps of the BDT training/optimization/plots, then copies the dataset directory to EOS if one was created.
doTraining=false
doOptimization=true
doBDTPlots=true
#years='2018' # '2016preVFP 2016postVFP 2017 2018'

#for year in $years; do

#destDir=$LQDATAEOS/BDT_19AugSkim/LQToDEle/$year
destDir=$LQDATAEOS/testBDTAllYears/combinedTest
#rm -r dataset
#cp -r $LQDATAEOS/BDT_7maySkim_10julxsec/LQToDEle/$year/dataset .

if [ ! -d $destDir ]; then
    echo "EOS directory does not exist. Making one"
    mkdir $destDir
fi

if [ "$doTraining" = true ]; then
echo "removing existing directory dataset"
rm -r dataset
fi 

if [ "$doTraining" = true ]; then 
echo "train BDT"
python scripts/tmvaBDT.py -t -d $destDir | tee bdtTraining.log
cp bdtTraining.log $destDir
echo "*************************************************************************"
fi

if [ "$doOptimization" = true ]; then
echo "optimize BDT"
python scripts/tmvaBDT.py -o -d $destDir | tee bdtOptimization-200Bin-minNB1.log
mv bdtOptimization-200Bin-minNB1.log $destDir
echo "*************************************************************************"
fi 

if [ "$doBDTPlots" = true ]; then
echo "do ROC and BDT plots"
python scripts/tmvaBDT.py -r -d $destDir | tee bdtPlots-200Bin.log
mv bdtPlots-200Bin.root bdtPlots-200Bin.log $destDir
echo "*************************************************************************"
fi

#if [ "$doTraining" = true ]; then
#echo "copy dataset directory to $destDir"
#cp -r dataset $destDir
#fi
