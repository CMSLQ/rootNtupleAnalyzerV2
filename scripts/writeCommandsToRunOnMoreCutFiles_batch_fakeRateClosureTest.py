import sys
import os

year = str(sys.argv[1])
#analysisType = sys.argv[2] #1P1F or 2F

ananame="qcdFRClosureTest_allYears"
if "2016" in year:
    yearInt = 2016
else:
    yearInt = int(year)
#Script loops over 1P1F and 2F and prints the commands for both to the terminal for the specified year. Commands assume you are working from $LQANA
#Edit to use your input lists and specify where the outputs should be saved.
for analysisType in ["1P1F", "2F"]:
    #Output directory on AFS
    outdirpath = os.getenv("LQDATAAFS") + "/{}/qcdFRClosureTest/frClosureTest_{}/{}/condor".format(year,year,analysisType)
    #EOS directory
    eosdir = os.getenv("LQDATAEOS")+ "/{}_17AugSkims/{}/{}".format(ananame,year,analysisType)
    #Input lists
    if analysisType == "1P1F":
        inputList = "config/myDatasets/{}HEEP/inputListAllCurrent.txt".format(year)
    else:
        inputList = "config/myDatasets/{}DataOnly/inputListAllCurrent.txt".format(year)
#Everything from here down you shouldn't need to modify

    cutfilename = os.getenv("LQMACRO")+"/config{}/QCDFakeRate/cutTable_lq_QCD_FakeRateClosureTest_{}.txt".format(yearInt,analysisType)
    if year == "2016preVFP":
        cutfilename = os.getenv("LQMACRO")+"/config{}/QCDFakeRate/preVFP/cutTable_lq_QCD_FakeRateClosureTest_{}.txt".format(yearInt,analysisType)
    if year == "2016postVFP":
        cutfilename = os.getenv("LQMACRO")+"/config{}/QCDFakeRate/postVFP/cutTable_lq_QCD_FakeRateClosureTest_{}.txt".format(yearInt,analysisType)

    queue = "workday"

    if analysisType == "1P1F":
        sampleList = "config/sampleListForMerging_13TeV_QCD_closureTest_{}.yaml".format(year)
    if analysisType == "2F":
        sampleList = "config/sampleListForMerging_13TeV_QCD_closureTest_{}DataOnly.yaml".format(year) 

    lumi = {}
    lumi["2016preVFP"] = 19501601.622
    lumi["2016postVFP"] = 16812151.722
    lumi["2017"] = 41540000
    lumi["2018"] = 59830000

    launchAnalysis = "python scripts/launchAnalysis_batch_ForSkimToEOS.py -i {} -o {} -c {} -q {} -d {} -j 1 -n rootTupleTree/tree".format(inputList, outdirpath, cutfilename, queue, eosdir)

    combinePlots = "time ./scripts/combinePlots.py -i {} -d {}/condor -c {} -l {} -x config/xsection_13TeV_2022.txt -o {}/output_cutTable_lq_QCD_FakeRateClosureTest -s {} | tee {}/combinePlots.log".format(inputList, eosdir, ananame, lumi[year], eosdir, sampleList, eosdir) 
    
    xsectionFileSF = "xsection_13TeV_2022_MeeControlReg_TTbar_MeeControlReg_DYJets_{}.txt".format(year)
    combinePlotsSF = "time ./scripts/combinePlots.py -i {} -d {}/condor -c {} -l {} -x {} -o {}_SF/output_cutTable_lq_QCD_FakeRateClosureTest -s {} | tee {}/combinePlotsSF.log".format(inputList, eosdir, ananame, lumi[year],xsectionFileSF, eosdir, sampleList, eosdir)

    print("\n+++++++++++++++++++++++++++++\n")
    print(launchAnalysis)
    print("\n+++++++++++++++++++++++++++++\n")
    print(combinePlots)
#    if analysisType=="1P1F":
#        print("\n+++++++++++++++++++++++++++++\n")
#        print(combinePlotsSF)
