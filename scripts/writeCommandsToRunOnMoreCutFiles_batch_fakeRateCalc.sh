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

[ $# -gt 0 ] || die "must specify year as argument; usage: writeCommandsToRunOnMoreCutFiles_batch_fakeRateCalc.sh year"
YEAR=$1

#### INPUTS HERE ####
#------------
ANANAME=fakeRateCalcFinal
#------------
files2016pre="$LQMACRO/config2016/QCDFakeRate/cutTable_lq_QCD_FakeRateCalculation_preVFP.txt"
files2016post="$LQMACRO/config2016/QCDFakeRate/cutTable_lq_QCD_FakeRateCalculation_postVFP.txt"
files2017="$LQMACRO/config2017/QCDFakeRate/cutTable_lq_QCD_FakeRateCalculation.txt"
files2018="$LQMACRO/config2018/QCDFakeRate/cutTable_lq_QCD_FakeRateCalculation.txt"
#------------
QUEUE=workday
#------------
OUTDIRPATH=$LQDATAAFS  # a subdir will be created for each cut file 
excludeCombining=""
#I found it useful to define two separate $LQDATA paths, one for afs and one for eos. - Emma
# output sub-directory (i.e. output will be in OUTDIRPATH/SUBDIR)
# it is suggested to specify the luminosity in the name of the directory
#------------
ilumi2016pre=19501.601622
ilumi2016post=16812.151722
ilumi2017=41540 #FIXME: this number is just from the Egamma twiki
ilumi2018=59830
CODENAME=analysisClass_lq_QCD_FakeRateCalculation
#------------
inputlist2016pre=config/myDatasets/2016HEEPpreVFP/inputListAllCurrent.txt
inputlist2016post=config/myDatasets/2016HEEPpostVFP/inputListAllCurrent.txt
inputlist2017=config/myDatasets/2017HEEP/inputListAllCurrent.txt
inputlist2018=config/myDatasets/2018HEEP/inputListAllCurrent.txt
#------------
if [ "$YEAR" = "2016preVFP" ]; then
  echo "Doing 2016preVFP!"
  ILUM=$ilumi2016pre
  INPUTLIST=$inputlist2016pre
  files=$files2016pre
elif [ "$YEAR" = "2016postVFP" ]; then
  echo "Doing 2016postVFP!"
  ILUM=$ilumi2016post
  INPUTLIST=$inputlist2016post
  files=$files2016post
elif [ "$YEAR" = "2017" ]; then
  ILUM=$ilumi2017
  INPUTLIST=$inputlist2017
  files=$files2017
  excludeCombining="-e LQToBEle*"
elif [ "$YEAR" = "2018" ]; then
  ILUM=$ilumi2018
  INPUTLIST=$inputlist2018
  files=$files2018
fi
SUBDIR=qcdFakeRateCalc/${YEAR}
EOSDIR=$LQDATAEOS/$ANANAME/${YEAR}
COMMANDFILE=commandsToRunOnMoreCutFiles_fakeRateCalc_${YEAR}_batch_`hostname -s`.txt
SAMPLELISTFORMERGING=config/sampleListForMerging_13TeV_QCD_calc_${YEAR}.yaml
#------------
FACTOR=1000 # numbers in final tables (but *not* in plots) will be multiplied by this scale factor (to see well the decimal digits)
#------------
EXE=main
#------------
XSECTION=config/xsection_13TeV_2022.txt #specify cross section file
#------------
#### END OF INPUTS ####

echo "" > $COMMANDFILE

for file in $files
do
suffix=`basename $file`
suffix=${suffix%\.*}
cat >> $COMMANDFILE <<EOF

####################################################
#### launch, check and combine cmds for $suffix ####
# 1 job per file
python scripts/launchAnalysis_batch_ForSkimToEOS.py -i $INPUTLIST -o $OUTDIRPATH/$SUBDIR/condor -c $file -q $QUEUE -d $EOSDIR -j 1 -n rootTupleTree/tree

./scripts/checkJobs.sh $OUTDIRPATH/$SUBDIR/condor $OUTDIRPATH/$SUBDIR/condor

time  ./scripts/combinePlots.py \
    -i $INPUTLIST \
    -c $CODENAME \
    -d $EOSDIR/condor \
    -l  `echo "$ILUM*$FACTOR" | bc` \
    -x $XSECTION  \
    -o $EOSDIR/output_$suffix \
    -s $SAMPLELISTFORMERGING \
    | tee $EOSDIR/output_$suffix/combineTablesAndPlots_${suffix}.log
EOF
done


echo "The set of commands to run on the cut files:" 
for file in $files
do
echo "  " $file
done 
echo "has been written to $COMMANDFILE"
