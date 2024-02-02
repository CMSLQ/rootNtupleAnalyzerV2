#!/usr/bin/env python3
#
# see: https://github.com/lmoneta/tmva-tutorial/blob/master/notebooks/TMVA_Classification.ipynb
#      https://github.com/lmoneta/tmva-tutorial/blob/master/notebooks/TMVA_Reader_ManyMethods_py.ipynb

from array import array
from collections import OrderedDict
import sys
import os
import math
import multiprocessing
import traceback
import copy
import numpy as np
import ctypes
import time
from optparse import OptionParser

from combineCommon import ParseXSectionFile, lookupXSection

import ROOT
from ROOT import TMVA, TFile, TString, TCut, TChain, TFileCollection, gROOT, gDirectory, gInterpreter, TEntryList, TH1D, TProfile, RDataFrame, TCanvas, TLine, kRed, kBlue, kSpring, TGraph, TGraphErrors, TMultiGraph, gPad, RooStats
ROOT.EnableImplicitMT(6)

class Sample:
    def __init__(self, name, subSampleList, subSampleWeights):
        self.name = name
        self.subSamples = subSampleList
        self.subSampleWeights = subSampleWeights


def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)
    sys.stdout.write("\r"+logString.format(jobCount, len(lqMassesToUse)))
    sys.stdout.write("\t"+str(len(result_list))+" jobs done")
    sys.stdout.flush()


def CalcWeight(fullDatasetName, intLumi, sumWeights):
    fullDatasetName = fullDatasetName.replace("_APV", "")
    xsec = float(lookupXSection(fullDatasetName))
    if xsec > 0:
        print("INFO: weight for {} = {}*{}/{}={}".format(fullDatasetName, intLumi, xsec, sumWeights, intLumi*xsec/sumWeights))
        return intLumi*xsec/sumWeights
    else:
        return 1.0  # in the case of data


def CheckRootFiles(fileList, treeName):
    for i in range(0, fileList.GetEntries()):
        thisFile = fileList.At(i).GetCurrentUrl().GetFile()
        tfile = TFile.Open(thisFile)
        tree = tfile.Get(treeName)
        if not tree or tree is None:
            raise RuntimeError("Couldn't read TTree named {} in file {}".format(treeName, thisFile))
        tfile.Close()


def GetSumWeightsInChain(tchain, cut, weight):
    df = RDataFrame(tchain)
    df = df.Filter(cut.GetTitle())
    df = df.Define('eventWeight', eventWeightExpression)
    sumWeights =  df.Sum("eventWeight").GetValue()
    return sumWeights*weight


def LoadChainFromTxtFile(txtFile, console=None):
    txtFile = os.path.expandvars(txtFile)
    if not os.path.isfile(txtFile):
        raise RuntimeError("File {} does not exist!".format(txtFile))
    fc = TFileCollection("dum","",txtFile)
    if fc.GetNFiles() <= 0:
        raise RuntimeError("Got <= 0 files loaded into the TFileCollection for {}!".format(txtFile))
    treeName = "rootTupleTree/tree"
    CheckRootFiles(fc.GetList(), treeName)
    ch = TChain(treeName)
    ch.AddFileInfoList(fc.GetList())
    # ch.SetBranchStatus("*", 0)
    # for branchName in neededBranches:
    #     ch.SetBranchStatus(branchName, 1)
    if ch.GetEntries() <= 0:
        logString = "WARNING: Got <= 0 entries for dataset={}; returning None!".format(txtFile)
        if console is not None:
            console.print(logString)
        print(logString)
        return None
    # else:
    #     logString = "INFO: Got {} entries for dataset={}".format(ch.GetEntries(), txtFile)
    #     if console is not None:
    #         console.print(logString)
    #     print(logString)
    return ch


def LoadDatasets(datasetDict, neededBranches, signal=False, loader=None, year=None, lqMass=None, nLQPoints=1):
    nTotEvents = 0
    nTotSumWeights = 0
    cut = mycuts if signal else mycutb
    if loader is None:
        totalTChain = TChain("rootTupleTree/tree")
    for key, value in datasetDict.items():
        print("Loading tree for dataset={}; signal={}".format(key, signal))
        if isinstance(value, list):
            nSampleTotEvents = 0
            nSampleSumWeights = 0
            for count, txtFile in enumerate(value):
                txtFile = txtFile.format(year, lqMass)
                ch = LoadChainFromTxtFile(txtFile)
                if ch is None:
                    continue
                nSampleTotEvents += ch.GetEntries()
                datasetName = os.path.basename(txtFile).replace(".txt", "")
                if "data" not in key.lower():
                    sumWeights = GetBackgroundSumWeights(datasetName, txtFile)
                    weight = CalcWeight(datasetName, intLumi, sumWeights)
                else:
                    weight = 1.0
                nSampleSumWeights = GetSumWeightsInChain(ch, cut, weight)
                print("Add file={} with weight*1000={} to collection".format(txtFile, weight*1000))
                if loader is not None:
                    if signal:
                        loader.AddSignalTree    ( ch, weight )
                    else:
                        loader.AddBackgroundTree( ch, weight/nLQPoints )
                    print("Loaded tree from file {} with {} events and {} sumWeights.".format(txtFile, ch.GetEntries(), nSampleSumWeights))
                else:
                    totalTChain.Add(ch)
            print("Loaded tree for sample {} with {} entries; sumWeights={}.".format(key, nSampleTotEvents, nSampleSumWeights))
            nTotEvents += nSampleTotEvents
            nTotSumWeights += nSampleSumWeights
        else:
            txtFile = value
            txtFile = txtFile.format(lqMass)
            ch = LoadChainFromTxtFile(txtFile)
            if ch is None:
                continue
            nSampleTotEvents = ch.GetEntries()
            datasetName = os.path.basename(txtFile).replace(".txt", "")
            if "data" not in key.lower():
                sumWeights = GetBackgroundSumWeights(datasetName, txtFile)
                weight = CalcWeight(datasetName, intLumi, sumWeights)
            else:
                weight = 1.0
            nSampleSumWeights = GetSumWeightsInChain(ch, cut, weight)
            print("Add file={} with weight*1000={} to collection".format(txtFile, weight*1000))
            if loader is not None:
                if signal:
                    loader.AddSignalTree    ( ch, weight )
                else:
                    loader.AddBackgroundTree( ch, weight/nLQPoints )
                print("Loaded tree from file {} with {} events and {} sumWeights.".format(txtFile, ch.GetEntries(), nSampleSumWeights))
            else:
                totalTChain.Add(ch)
            print("Loaded tree for sample {} with {} entries; sumWeights={}".format(key, nSampleTotEvents, nSampleSumWeights))
            nTotEvents += nSampleTotEvents
            nTotSumWeights += nSampleSumWeights
    print("Total: loaded tree with {} entries; sumWeights={}".format(nTotEvents, nTotSumWeights))
    sys.stdout.flush()
    if loader is None:
        return totalTChain


def TrainBDT(args):
    lqMassToUse = args[0]
    year = args[1]
    normVars = args[2]
    # just one use the single LQ signal specified by mass just above
    try:
        signalDatasetsDict = {}
        signalDatasetName = signalNameTemplate.format(lqMassToUse)
        signalDatasetsDict[signalDatasetName] = allSignalDatasetsDict[signalDatasetName]
        print(signalDatasetsDict)
        
        outputFile = TFile.Open("TMVA_ClassificationOutput_"+signalDatasetName+".root", "RECREATE")
        
        TMVA.Tools.Instance()
        factory = TMVA.Factory("TMVAClassification_"+signalDatasetName, outputFile, "!V:ROC:!Silent:Color:DrawProgressBar:AnalysisType=Classification")
        
        loader = TMVA.DataLoader("dataset")
        LoadDatasets(backgroundDatasetsDict, neededBranches, signal=False, loader=loader, lqMass=lqMassToUse, year=year)
        LoadDatasets(signalDatasetsDict, neededBranches, signal=True, loader=loader, lqMass=lqMassToUse, year=year)
        
        #  Set individual event weights (the variables must exist in the original TTree)
        # if(analysisYear < 2018 && hasBranch("PrefireWeight") && !isData()) --> prefire weight
        loader.SetBackgroundWeightExpression( eventWeightExpression )
        loader.SetSignalWeightExpression( eventWeightExpression )
        
        # define variables
        #loader.AddVariable( "myvar1 := var1+var2", 'F' )
        #loader.AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' )
        #loader.AddVariable( "var3",                "Variable 3", "units", 'F' )
        #loader.AddVariable( "var4",                "Variable 4", "units", 'F' )
        #
        #loader.AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' )
        #loader.AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' )
        for var in variableList:
            if normVars:
                if normTo in var:
                    continue
                elif any(exp in var for exp in ["Ele", "Jet", "M", "sT"]) and not "PFMET" in var:
                    loader.AddVariable(var+"/"+normTo, "F")
                    continue
            loader.AddVariable(var, "F")
        
        # Tell the factory how to use the training and testing events
        #
        # If no numbers of events are given, half of the events in the tree are used
        # for training, and the other half for testing:
        # loader.PrepareTrainingAndTestTree( mycuts, mycutb, "V:SplitMode=random:NormMode=None:VerboseLevel=Debug" )
        # loader.PrepareTrainingAndTestTree( mycuts, mycutb, "!V:SplitMode=random:NormMode=NumEvents" )
        # loader.PrepareTrainingAndTestTree( mycuts, mycutb, "V:SplitMode=random:NormMode=EqualNumEvents:VerboseLevel=Debug" )
        loader.PrepareTrainingAndTestTree( mycuts, mycutb, "V:SplitMode=random:NormMode=EqualNumEvents" )
        # To also specify the number of testing events, use:
        # loader.PrepareTrainingAndTestTree( mycut,
        #                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" )
        # loader.PrepareTrainingAndTestTree( mycuts, mycutb,
        #                                     "nTrain_Signal=4000:nTrain_Background=2000:SplitMode=Random:NormMode=NumEvents:!V" )
        
        # Boosted Decision Trees
        # factory.BookMethod(loader,TMVA.Types.kBDT, "BDT",
        #         "V:NTrees=400:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20" )
        # factory.BookMethod(loader,TMVA.Types.kBDT, "BDT",
        #         "!V:NTrees=800:MinNodeSize=2.5%:MaxDepth=2:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20" )
        factory.BookMethod(loader, TMVA.Types.kBDT, "BDTG",
                # "!H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=1:BoostType=Grad:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:NegWeightTreatment=Pray" )
                # "!H:!V:NTrees=400:MinNodeSize=2.5%:MaxDepth=1:BoostType=Grad:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:NegWeightTreatment=IgnoreNegWeightsInTraining" )
                # "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=GiniIndex:NTrees=850:MinNodeSize=2.5%:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=20" )  # LQ2
                # "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=Pray:SeparationType=GiniIndex:NTrees=850:MinNodeSize=2.5%:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=20" )  # LQ2 but use neg. weights in training with pray
                "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=Pray:SeparationType=GiniIndex:NTrees=1000:MinNodeSize=5%:Shrinkage=0.01:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:CreateMVAPdfs:NbinsMVAPdf=20" )
        
        factory.TrainAllMethods()
        factory.TestAllMethods()
        factory.EvaluateAllMethods()
        
        c1 = factory.GetROCCurve(loader)
        # c1.Draw()
        c1.Write("rocCurve_lqm"+str(lqMassToUse)+".png")
        
        outputFile.Close()
    except Exception as e:
        print("ERROR: exception in TrainBDT for lqMass={}".format(lqMassToUse))
        traceback.print_exc()
        raise e

    return True


def TrainParametrizedBDT(lqMassList, year):
   outputFile = TFile.Open("TMVA_ClassificationOutput.root", "RECREATE")
   
   TMVA.Tools.Instance()
   factory = TMVA.Factory("TMVAClassification", outputFile, "!V:ROC:!Silent:Color:DrawProgressBar:AnalysisType=Classification")
   
   loader = TMVA.DataLoader("dataset")
   for lqMass in lqMassList:
       LoadDatasets(backgroundDatasetsDict, neededBranches, signal=False, loader=loader, lqMass=lqMass, nLQPoints=len(lqMassList), year=year)
       signalDatasetName = signalNameTemplate.format(lqMass)
       signalDatasetsDict = {}
       signalDatasetsDict[signalDatasetName] = allSignalDatasetsDict[signalDatasetName]
       LoadDatasets(signalDatasetsDict, neededBranches, signal=True, loader=loader, lqMass=lqMass, year=year)
   
   #  Set individual event weights (the variables must exist in the original TTree)
   # if(analysisYear < 2018 && hasBranch("PrefireWeight") && !isData()) --> prefire weight
   loader.SetBackgroundWeightExpression( eventWeightExpression )
   loader.SetSignalWeightExpression( eventWeightExpression )
   
   for var in variableList:
       loader.AddVariable(var, "F")
   
   loader.PrepareTrainingAndTestTree( mycuts, mycutb, "V:SplitMode=random:NormMode=EqualNumEvents" )
   factory.BookMethod(loader, TMVA.Types.kBDT, "BDTG",
           # "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=IgnoreNegWeightsInTraining:SeparationType=GiniIndex:NTrees=850:MinNodeSize=2.5%:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=20" )  # LQ2
           # "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=Pray:SeparationType=GiniIndex:NTrees=850:MinNodeSize=2.5%:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:CreateMVAPdfs:NbinsMVAPdf=20" )  # LQ2 but use neg. weights in training with pray
           "!H:!V:BoostType=Grad:DoBoostMonitor:NegWeightTreatment=Pray:SeparationType=GiniIndex:NTrees=1000:MinNodeSize=5%:Shrinkage=0.01:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:CreateMVAPdfs:NbinsMVAPdf=20" )
   
   factory.TrainAllMethods()
   factory.TestAllMethods()
   factory.EvaluateAllMethods()
   
   c1 = factory.GetROCCurve(loader)
   # c1.Draw()
   c1.Write("rocCurve.png")
   
   outputFile.Close()


def GetTotalEventsHist(lqMassToUse, year, signalDict, signalNameTemplate):
    signalDatasetName = signalNameTemplate.format(lqMassToUse)
    txtFiles = signalDict[signalDatasetName]
    # profName = "EventsPassingCuts_unscaled"
    hist = None
    histName = "savedHists/EventCounter"
    tfiles = []
    for count, txtFile in enumerate(txtFiles):
        txtFile = txtFile.format(year, lqMassToUse)
        with open(os.path.expandvars(txtFile), "r") as theTxtFile:
            for line in theTxtFile:
                line = line.strip()
                tfile = TFile.Open(os.path.expandvars(line))
                tfiles.append(tfile)
                unscaledEvtsHist = tfile.Get(histName)
                if not unscaledEvtsHist or unscaledEvtsHist.ClassName() != "TH1D":
                    unscaledEvtsHist = tfile.Get("EventCounter")
                    if not unscaledEvtsHist or unscaledEvtsHist is None:
                        raise RuntimeError("Did not find object named {} in file '{}'.".format(histName, tfile.GetName()))
                    if unscaledEvtsHist.ClassName() != "TH1D":
                        raise RuntimeError("Expected class TH1D for object named {} but class is '{}' instead.".format(histName, unscaledEvtsHist.ClassName()))
                if hist is None:
                    hist = copy.deepcopy(unscaledEvtsHist)
                else:
                    hist.Add(unscaledEvtsHist)
    for tfile in tfiles:
        tfile.Close()
    return hist


def GetSignalTotalEvents(lqMassToUse, year):
    hist = GetTotalEventsHist(lqMassToUse, year, allSignalDatasetsDict, signalNameTemplate)
    # for TProfiles
    # unscaledTotalEvts = prof.GetBinContent(1)*prof.GetBinEntries(1)
    unscaledTotalEvts = hist.GetBinContent(1)
    return unscaledTotalEvts


def GetSignalSumWeights(lqMassToUse, year):
    hist = GetTotalEventsHist(lqMassToUse, year, allSignalDatasetsDict, signalNameTemplate)
    # for TProfiles
    # unscaledTotalEvts = prof.GetBinContent(1)*prof.GetBinEntries(1)
    sumWeights = hist.GetBinContent(3)
    return sumWeights


def GetBackgroundSumWeights(sampleName, txtFile):
    hist = GetTotalEventsHist(0, 0, {sampleName: [txtFile]}, sampleName)
    binToUse = 3 if "powhegMiNNLO" not in txtFile else 5
    sumWeights = hist.GetBinContent(binToUse)
    if sumWeights == 0:
        raise RuntimeError("got sumWeights=0 for sampleName={}, txtFile={}; used bin={}".format(sampleName, txtFile, binToUse))
    else:
        print("INFO: got sumWeights={} for sampleName={}, txtFile={}; used bin={}".format(sumWeights, sampleName, txtFile, binToUse))
    return sumWeights


def EvaluateFigureOfMerit(nS, nB, efficiency, bkgEnts, figureOfMerit):
    # see: https://twiki.cern.ch/twiki/bin/view/CMS/FigureOfMerit
    # and https://arxiv.org/pdf/physics/0702156.pdf [1]
    try:
        # s/sqrt(s+b)
        # value = nS / ( math.sqrt ( nS + nB ) )
        if bkgEnts != 0:
            tau = bkgEnts / nB
        # switch to asymptotic formula
        # NB: this approximation doesn't work well with nB < ~ 5 events
        if figureOfMerit == "asymptotic":
            value = math.sqrt(2 * ((nS + nB) * math.log(1 + nS / nB) - nS))
        elif figureOfMerit == "punzi":
            # punzi
            a = 2  # nSigmasExclusion
            b = 5  # nSigmasDiscovery
            # value = efficiency / (nSigmas / 2.0 + math.sqrt(nB))
            smin = a**2/8 + 9*b**2/13 + a*math.sqrt(nB) + (b/2)*math.sqrt(b**2 + 4*a*math.sqrt(nB) + 4*nB)
            value = efficiency / smin
        elif figureOfMerit == "zbi":
            value = RooStats.NumberCountingUtils.BinomialWithTauExpZ(nS, nB, tau)
        elif figureOfMerit == "zpl":  # [1], eqn. 25
            nOff = bkgEnts
            nOn = nS + nB
            nTot = nOff + nOn
            value = math.sqrt(2)*math.sqrt(nOn*math.log(nOn*(1+tau)/nTot) + nOff*math.log(nOff*(1+tau)/(nTot*tau)))
        elif figureOfMerit == "ssb":
            value = nS / math.sqrt(nS+nB)
        else:
            raise RuntimeError("Evaluation of '{}' as figure of merit is not implemented".format(figureOfMerit))
    except ZeroDivisionError:
        value = -999
    except ValueError:
        print("WARNING: had a domain error calculating the value with nS=", nS, "and nB=", nB)
        value = -999
    return value


def OptimizeBDTCut(args):
    bdtWeightFileName, lqMassToUse, sharedOptValsDict, sharedOptHistsDict, sharedFOMInfoDict, year = args
    startTime = time.time()
    try:
        signalDatasetsDict = {}
        signalDatasetName = signalNameTemplate.format(lqMassToUse)
        signalDatasetsDict[signalDatasetName] = allSignalDatasetsDict[signalDatasetName]
        print(signalDatasetsDict)
        TMVA.Tools.Instance()
        # reader = TMVA.Reader("!Color:!Silent")
        # varValueDict = {}
        # for var in neededBranches:
        #     varValueDict[var] = array("f", [0])
        # for var in variableList:
        #     reader.AddVariable(var, varValueDict[var])
        name = "BDTG"
        # methodNames = [name]
        # print("name={}, bdtWeightFileName={}".format(name, bdtWeightFileName))
        # reader.BookMVA(name, bdtWeightFileName )

        binsToUse = 1000 # 10000
        hname = "hsig_" + name + "_" + str(lqMassToUse)
        htitle = "Classifier Output on signal for " + name + ", M_{LQ} = " + str(lqMassToUse) + " GeV"
        hsig = TH1D(hname,htitle,binsToUse,-1,1)
        hsigUnweighted = TH1D(hname+"_unweighted", htitle+" (unweighted)", binsToUse, -1, 1)
        #hname = "hbkg_" + name + "_" + lqMassToUse
        #hbkg = ROOT.RDF.TH1DModel(hname,htitle,binsToUse,-1,1)

        numVars = len(variableList)
        if normalizeVars:
            numVars -= 1
        gInterpreter.ProcessLine(('''
        TMVA::Experimental::RReader BDT{}("{}");
        computeBDT{} = TMVA::Experimental::Compute<{}, float>(BDT{});
        ''').format(lqMassToUse, bdtWeightFileName, lqMassToUse, numVars, lqMassToUse))
        sys.stdout.flush()
        sys.stderr.flush()
        # backgrounds
        histTitle = "Classifier Output on {} background for " + name + ", M_{{LQ}} = " + str(lqMassToUse) + " GeV"
        histName = "BDTOutput{}LQM" + str(lqMassToUse)
        bkgTotal = TH1D(histName.format("TotalBackground"), histTitle.format("all"), binsToUse, -1, 1)
        bkgTotalUnweighted = TH1D(histName.format("TotalBackground")+"Unweighted", histTitle.format("all")+" (unweighted)", binsToUse, -1, 1)
        bkgHists = dict()
        bkgHistsUnweightedUnscaled = dict()
        bkgTotIntegralOverCut = 0
        cutValForIntegral = 0.9940
        for sample in backgroundDatasetsDict.keys():
            bkgSampleIntegralOverCut = 0
            bkgSampleIntegral = 0
            bkgSampleIntegralHist = 0
            bkgHists[sample] = TH1D(histName.format(sample), histTitle.format(sample), binsToUse, -1, 1)
            bkgHistsUnweightedUnscaled[sample] = TH1D(histName.format(sample)+"_unweightedUnscaled", histTitle.format(sample)+", unweighted/unscaled", binsToUse, -1, 1)
            for idx, txtFile in enumerate(backgroundDatasetsDict[sample]):
                txtFile = txtFile.format(year, lqMassToUse)
                #tchainBkg = LoadChainFromTxtFile(txtFile.format(lqMassToUse))
                tchainBkg = LoadChainFromTxtFile(txtFile)
                if tchainBkg is None:
                    continue
                df = RDataFrame(tchainBkg)
                df = df.Filter(mycutb.GetTitle())  # will work for expressions valid in C++
                if "LQCandidateMass" in variableList:
                    df = df.Define("massInt", str(lqMassToUse))
                    df = df.Redefine("LQCandidateMass", "Numba::GetMassFloat(massInt)")
                varNames = getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames()
                if normalizeVars:
                    # print("{} varNames: {}".format(len(varNames), varNames))
                    l_varn = ROOT.std.vector['std::string']()
                    for i_expr, expr in enumerate(varNames):
                        varname = 'v_{}'.format(i_expr)
                        l_varn.push_back(varname)
                        df=df.Define(varname, '(float)({})'.format(expr))
                    df = df.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), l_varn)
                else:
                    df = df.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), varNames)
                df = df.Define('BDT', 'BDTv[0]')
                df = df.Define('eventWeight', eventWeightExpression)
                histName = "BDTVal_{}_{}".format(sample, idx)
                hbkg = TH1D(histName, histName, binsToUse, -1, 1)
                histBkg = df.Histo1D(ROOT.RDF.TH1DModel(hbkg), "BDT", "eventWeight")
                hbkgUnweighted =  TH1D(histName+"_unweighted", histName+"_unweighted", binsToUse, -1, 1)
                histBkgUnweighted = df.Histo1D(ROOT.RDF.TH1DModel(hbkgUnweighted), "BDT")
                #bkgWeight = backgroundDatasetsWeightsTimesOneThousand[os.path.basename(txtFile).replace(".txt", "")]/1000.0
                #bkgWeight = FindWeight(os.path.basename(txtFile).replace(".txt", ""), backgroundDatasetsWeightsTimesOneThousand)/1000.0
                datasetName = os.path.basename(txtFile).replace(".txt", "")
                if "data" not in sample.lower():
                    sumWeights = GetBackgroundSumWeights(datasetName, txtFile)
                    bkgWeight = CalcWeight(datasetName, intLumi, sumWeights)
                else:
                    bkgWeight = 1.0
                histBkg = histBkg.GetValue()
                histBkg.Scale(bkgWeight)
                bkgHists[sample].Add(histBkg)
                bkgHistsUnweightedUnscaled[sample].Add(histBkgUnweighted.GetValue())
                bkgTotal.Add(histBkg)
                #bkgTotalUnweighted.Add(histBkgUnweighted.GetPtr())
                bkgTotalUnweighted.Add(histBkgUnweighted.GetValue())
                #h = df.Histo1D(hbkg, "BDT", "eventWeight")
                #h.Draw()
                bkgIntegral = df.Sum("eventWeight").GetValue()*bkgWeight
                # bkgIntegralOverCut = df.Filter("BDT > {}".format(cutValForIntegral)).Sum("eventWeight").GetValue()*bkgWeight
                # bkgEntriesOverCut = df.Filter("BDT > {}".format(cutValForIntegral)).Count().GetValue()
                print("subsample={}, bkgWeight={}".format(txtFile, bkgWeight), flush=True)
                print("subsample={}, entries = {}, integral unweighted = {}, integral weighted = {}".format(txtFile, histBkg.GetEntries(), histBkg.Integral()/bkgWeight, histBkg.Integral()), flush=True)
                print("subsample={}, df entries = {}, df integral unweighted = {}, df integral weighted = {}".format(txtFile, df.Count().GetValue(), df.Sum("eventWeight").GetValue(), df.Sum("eventWeight").GetValue()*bkgWeight), flush=True)
                # print("subsample={}, entries with BDT > {} = {}, integral unweighted = {}, integral weighted = {}".format(txtFile, cutValForIntegral, bkgEntriesOverCut, bkgIntegralOverCut/bkgWeight, bkgIntegralOverCut))
                sys.stdout.flush()
                # print some entries
                # if sample == "QCDFakes_DATA":
                #     cols = ROOT.vector('string')()
                #     cols.push_back("run")
                #     cols.push_back("event")
                #     cols.push_back("ls")
                #     # cols.push_back("eventWeight")
                #     cols.push_back("sT_eejj")
                #     cols.push_back("PFMET_Type1_Pt")
                #     cols.push_back("M_e1e2")
                #     cols.push_back("M_e1j1")
                #     cols.push_back("M_e1j2")
                #     cols.push_back("M_e2j1")
                #     cols.push_back("M_e2j2")
                #     cols.push_back("Ele1_Pt")
                #     cols.push_back("Ele2_Pt")
                #     cols.push_back("MejMin")
                #     cols.push_back("MejMax")
                #     cols.push_back("Meejj")
                #     cols.push_back("BDT")
                #     df2 = df.Filter("run==276831 && event==2249978927 && ls==1240")
                #     display = df2.Display(cols, 5, 20)
                #     # display.Print()
                #     print(display.AsString(), flush=True)
                #
                # if "photon" in txtFile.lower() and bkgIntegralOverCut > 0:
                #     cols = ROOT.vector('string')()
                #     cols.push_back("run")
                #     cols.push_back("ls")
                #     cols.push_back("event")
                #     cols.push_back("BDT")
                #     cols.push_back("eventWeight")
                #     display = df.Filter("BDT > {}".format(cutValForIntegral)).Display(cols, 20)
                #     print("subsample={}, entries over BDT cut unweighted = {}".format(txtFile, display.Print()))
                #     sys.stdout.flush()
                #bkgTotIntegralOverCut += bkgIntegralOverCut
                #bkgSampleIntegralOverCut += bkgIntegralOverCut
                bkgSampleIntegral += bkgIntegral
                bkgSampleIntegralHist += histBkg.Integral()
            print("sample={}, events = {} [df], from hist = {}".format(sample, bkgSampleIntegral, bkgSampleIntegralHist), flush=True)
            #print("sample={}, events over BDT cut = {}".format(sample, bkgSampleIntegralOverCut))
        # print("bkgIntegralOverCut={}".format(bkgIntegralOverCut), flush=True)

        # signal
        tchainSig = LoadDatasets(signalDatasetsDict, neededBranches, signal=True, loader=None, lqMass=lqMassToUse, year=year)
        dfSig = RDataFrame(tchainSig)
        dfSig = dfSig.Filter(mycuts.GetTitle())  # will work for expressions valid in C++
        # dfSig = dfSig.Define('BDTv', ROOT.computeBDT, ROOT.BDT.GetVariableNames())
        # dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames())
        varNames = getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames()
        if normalizeVars:
            # print("{} varNames: {}".format(len(varNames), varNames))
            l_varn = ROOT.std.vector['std::string']()
            for i_expr, expr in enumerate(varNames):
                varname = 'v_{}'.format(i_expr)
                l_varn.push_back(varname)
                dfSig=dfSig.Define(varname, '(float)({})'.format(expr))
            dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), l_varn)
        else:
            dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), varNames)
        dfSig = dfSig.Define('BDT', 'BDTv[0]')
        dfSig = dfSig.Define('eventWeight', eventWeightExpression)
        histSig = dfSig.Histo1D(ROOT.RDF.TH1DModel(hsig), "BDT", "eventWeight")
        histSigUnweighted = dfSig.Histo1D(ROOT.RDF.TH1DModel(hsigUnweighted), "BDT")
        #datasetName = os.path.basename(txtFile).replace(".txt", "")
        sumWeights = GetSignalSumWeights(lqMassToUse, year)
        signalWeight = CalcWeight(signalDatasetName, intLumi, sumWeights)
        #signalWeight = signalDatasetsWeightsTimesOneThousand[signalDatasetName]/1000.0
        print("multiply hist by signal weight ", signalWeight)
        histSig.Scale(signalWeight)

        # print some entries
        # cols = ROOT.vector('string')()
        # cols.push_back("BDT")
        # cols.push_back("eventWeight")
        # d2 = dfSig.Display(cols)
        # d2.Print()
        print("For LQM={}, totalSignal={}, {} raw events".format(lqMassToUse, histSig.Integral(), histSigUnweighted.Integral()))
        print("For LQM={}, totalBackground={}, {} raw events".format(lqMassToUse, bkgTotal.Integral(), bkgTotalUnweighted.Integral()))

        # now optimize
        #totalSignalEventsUnscaled = GetSignalTotalEvents(lqMassToUse)
        #sumWeights = GetSignalSumWeights(lqMassToUse)
        fomValueToCutInfoDict = {}
        fomList = []
        nSList = []
        effList = []
        nBList = []
        cutValList = []
        nSErrList = []
        nBErrList = []
        effErrList = []
        for iBin in range(1, hbkg.GetNbinsX()+1):
            skipFOMCalc = False
            nSErr = ctypes.c_double()
            nS = histSig.IntegralAndError(iBin, hsig.GetNbinsX(), nSErr)
            efficiency = nS/(signalWeight*sumWeights)
            effErr = nSErr.value/(signalWeight*sumWeights)  # assuming sumWeights and signalWeight have zero error
            # nBErr = ctypes.c_double()
            # nB = bkgTotal.IntegralAndError(iBin, hbkg.GetNbinsX(), nBErr)
            nB = 0
            nBErr = 0
            cutVal = histSig.GetBinLowEdge(iBin)
            for sample, hist in bkgHists.items():
                if "qcd" in sample.lower():
                    continue
                nBThisProcessErr = ctypes.c_double()
                nBThisProcess = hist.IntegralAndError(iBin, hist.GetNbinsX(), nBThisProcessErr)
                nBThisProcessErr = nBThisProcessErr.value
                if nBThisProcess < 0:
                    skipFOMCalc = True
                    break
                nB += nBThisProcess
                nBErr += nBThisProcessErr*nBThisProcessErr
            qcd1FRDataErr = ctypes.c_double()
            qcd1FRDataYield = bkgHists["QCDFakes_DATA"].IntegralAndError(iBin, hist.GetNbinsX(), qcd1FRDataErr)
            # print("INFO: Fot cutVal={}, Got qcd1FRDataYield={} from hist with integral={} and entries={}".format(cutVal, qcd1FRDataYield, bkgHists["QCDFakes_DATA"].Integral(), bkgHists["QCDFakes_DATA"].GetEntries()))
            qcd1FRDYJErr = ctypes.c_double()
            qcd1FRDYJYield = bkgHists["QCDFakes_DYJ"].IntegralAndError(iBin, hist.GetNbinsX(), qcd1FRDYJErr)
            qcd2FRDataErr = ctypes.c_double()
            qcd2FRDataYield = bkgHists["QCDFakes_DATA_2FR"].IntegralAndError(iBin, hist.GetNbinsX(), qcd2FRDataErr)
            qcd1FRYield = qcd1FRDataYield+qcd1FRDYJYield
            if qcd1FRYield < 0:
                # print("INFO: Limiting 1 FR QCD yield for cutVal {} to zero; old qcd1FRYield = 1FRData-DYJ = {} + {} = {}".format(
                #     cutVal, qcd1FRDataYield, qcd1FRDYJYield, qcd1FRYield))
                qcd1FRYield = 0
            limit = 0.5
            if abs(qcd2FRDataYield) > limit*qcd1FRYield:
                # print("INFO: Limiting 2 FR QCD yield for cutVal {} to {} from {}; qcd1FRYield = 1FRData-DYJ = {} + {} = {}; qcdYield = qcd2FRYield + qcd1FRYield = {} + {} = {}".format(
                #     cutVal, -1*limit*qcd1FRYield, qcd2FRDataYield, qcd1FRDataYield, qcd1FRDYJYield, qcd1FRYield, -1*limit*qcd1FRYield, qcd1FRYield, -1*limit*qcd1FRYield+qcd1FRYield))
                qcd2FRDataYield = -1*limit*qcd1FRYield
            nB += qcd2FRDataYield+qcd1FRYield
            nBErr += pow(qcd1FRDataErr.value, 2)+pow(qcd1FRDYJErr.value, 2)+pow(qcd2FRDataErr.value, 2)
            nBErr = math.sqrt(nBErr)
            # if nB < 3:
            # if nS < 5:
            #      fomValueToCutInfoDict[iBin] = [-1.0, cutVal, nS, efficiency, nB]
            #      continue
            # if iBin==9963:
            #     nSUnweightedErr = ctypes.c_double()
            #     unweightedNs = histSigUnweighted.IntegralAndError(iBin, hsigUnweighted.GetNbinsX(), nSUnweightedErr)
            #     nBUnweightedErr = ctypes.c_double()
            #     unweightedNb = bkgTotalUnweighted.IntegralAndError(iBin, hbkgUnweighted.GetNbinsX(), nBUnweightedErr)
            #     print("Evaluate figure of merit for nS={}, nB={}, unweightedNs={}, unweightedNb={}".format(nS, nB, unweightedNs, unweightedNb), flush=True)
            #     exit(0)
            # require at least one background event expected
            if nB > 1 and not skipFOMCalc:
                # fom = EvaluateFigureOfMerit(nS, nB if nB > 0.0 else 0.0, efficiency, bkgTotalUnweighted.Integral(iBin, hbkg.GetNbinsX()), "punzi")
                fom = EvaluateFigureOfMerit(nS, nB if nB > 0.0 else 0.0, efficiency, bkgTotalUnweighted.Integral(iBin, hbkg.GetNbinsX()), "asymptotic")
                # fom = EvaluateFigureOfMerit(nS, nB if nB > 0.0 else 0.0, efficiency, bkgTotalUnweighted.Integral(iBin, hbkg.GetNbinsX()), "zpl")
                # fom = EvaluateFigureOfMerit(nS, nB if nB > 0.0 else 0.0, efficiency, bkgTotalUnweighted.Integral(iBin, hbkg.GetNbinsX()), "zbi")
                # fom = EvaluateFigureOfMerit(nS, nB if nB > 0.0 else 0.0, efficiency, bkgTotalUnweighted.Integral(iBin, hbkg.GetNbinsX()), "ssb")
            else:
                fom = 0.0
            fomValueToCutInfoDict[iBin] = [fom, cutVal, nS, efficiency, nB]
            fomList.append(fom)
            nSList.append(nS)
            effList.append(efficiency)
            nBList.append(nB)
            cutValList.append(cutVal)
            nSErrList.append(nSErr.value)
            nBErrList.append(nBErr)
            effErrList.append(effErr)
        # sort by FOM
        sortedDict = OrderedDict(sorted(fomValueToCutInfoDict.items(), key=lambda t: float(t[1][0]), reverse=True))
        # now the max FOM value should be the first entry
        maxVal = next(iter(sortedDict.items()))
        #print("For lqMass={}, max FOM: ibin={} with FOM={}, cutVal={}, nS={}, eff={}, nB={}".format(lqMassToUse, maxVal[0], *maxVal[1]))
        # testVals=list(fomValueToCutInfoDict.items())
        # testVal=testVals[9949]
        # print("test FOM: ibin={} with FOM={}, cutVal={}, nS={}, eff={}, nB={}".format(testVal[0], *testVal[1]))
        # testVal=testVals[9948]
        # print("test FOM: ibin={} with FOM={}, cutVal={}, nS={}, eff={}, nB={}".format(testVal[0], *testVal[1]))
        # testVal=testVals[9950]
        # print("test FOM: ibin={} with FOM={}, cutVal={}, nS={}, eff={}, nB={}".format(testVal[0], *testVal[1]))
        valList = [maxVal[0]]
        valList.extend(maxVal[1])
        sharedOptValsDict[lqMassToUse] = valList
        sharedOptHistsDict[lqMassToUse] = [histSig.GetValue(), bkgTotal, histSigUnweighted.GetValue(), bkgTotalUnweighted]
        sharedFOMInfoDict[lqMassToUse]["FOM"] = fomList
        sharedFOMInfoDict[lqMassToUse]["nS"] = nSList
        sharedFOMInfoDict[lqMassToUse]["eff"] = effList
        sharedFOMInfoDict[lqMassToUse]["nB"] = nBList
        sharedFOMInfoDict[lqMassToUse]["cutVal"] = cutValList
        sharedFOMInfoDict[lqMassToUse]["nSErr"] = nSErrList
        sharedFOMInfoDict[lqMassToUse]["nBErr"] = nBErrList
        sharedFOMInfoDict[lqMassToUse]["effErr"] = effErrList
        # single values
        nSErr = ctypes.c_double()
        sharedFOMInfoDict[lqMassToUse]["nSNoBDTCut"] = histSig.IntegralAndError(1, histSig.GetNbinsX(), nSErr)
        sharedFOMInfoDict[lqMassToUse]["nSErrNoBDTCut"] = nSErr.value
        nBErr = ctypes.c_double()
        sharedFOMInfoDict[lqMassToUse]["nBNoBDTCut"] = bkgTotal.IntegralAndError(1, bkgTotal.GetNbinsX(), nBErr)
        sharedFOMInfoDict[lqMassToUse]["nBErrNoBDTCut"] = nBErr.value
        nSErr = ctypes.c_double()
        sharedFOMInfoDict[lqMassToUse]["nSUnweightedNoBDTCut"] = histSigUnweighted.IntegralAndError(1, histSig.GetNbinsX(), nSErr)
        sharedFOMInfoDict[lqMassToUse]["nSErrUnweightedNoBDTCut"] = nSErr.value
        nBErr = ctypes.c_double()
        sharedFOMInfoDict[lqMassToUse]["nBUnweightedNoBDTCut"] = bkgTotalUnweighted.IntegralAndError(1, bkgTotal.GetNbinsX(), nBErr)
        sharedFOMInfoDict[lqMassToUse]["nBErrUnweightedNoBDTCut"] = nBErr.value
        cutVal = maxVal[1][1]
        print(f"For LQM={lqMassToUse:4}, cutVal={cutVal:4.3f}", flush=True)
        for sample, hist in bkgHists.items():
            cutBin = hist.FindFixBin(cutVal)
            nBErr = ctypes.c_double()
            nB = hist.IntegralAndError(cutBin, hist.GetNbinsX(), nBErr)
            nBErr = nBErr.value
            bkgIntegral = hist.Integral()
            rawEventsHist = bkgHistsUnweightedUnscaled[sample]
            nBEventsErr = ctypes.c_double()
            nBEvents = rawEventsHist.IntegralAndError(cutBin, rawEventsHist.GetNbinsX(), nBEventsErr)
            nBEventsErr = nBEventsErr.value
            print(f"Background yield for optimized BDT cut for background={sample:20}: yield={nB:4.6f}+/-{nBErr:4.6f} [raw events={nBEvents:4.6f}+/-{nBEventsErr:4.6f}], fullYield={bkgIntegral:4.6f}", flush=True)
        cutBin = bkgTotal.FindFixBin(cutVal)
        nBErr = ctypes.c_double()
        nB = bkgTotal.IntegralAndError(cutBin, bkgTotal.GetNbinsX(), nBErr)
        nBErr = nBErr.value
        bkgIntegral = bkgTotal.Integral()
        rawEventsHist = bkgTotalUnweighted
        nBEventsErr = ctypes.c_double()
        nBEvents = rawEventsHist.IntegralAndError(cutBin, rawEventsHist.GetNbinsX(), nBEventsErr)
        nBEventsErr = nBEventsErr.value
        print(f"Background yield for optimized BDT cut for background={'Total':20}: yield={nB:4.6f}+/-{nBErr:4.6f} [raw events={nBEvents:4.6f}+/-{nBEventsErr:4.6f}], fullYield={bkgIntegral:4.6f}", flush=True)
    except Exception as e:
        print("ERROR: exception in OptimizeBDTCut for lqMass={}".format(lqMassToUse))
        traceback.print_exc()
        raise e
    endTime = time.time()
    totTime = endTime - startTime
    print("total time = ", totTime)
    return True


def DoROCAndBDTPlots(args):
    rootFileName, bdtWeightFileName, lqMassToUse, year = args
    startTime = time.time()
    timeMakingHistos = 0
    try:
        rootFile = TFile.Open(rootFileName, "update")
        rootFile.cd()
        rootFile.mkdir("LQM{}".format(lqMassToUse))
        rootFile.cd("LQM{}".format(lqMassToUse))
        signalDatasetsDict = {}
        signalDatasetName = signalNameTemplate.format(lqMassToUse)
        signalDatasetsDict[signalDatasetName] = allSignalDatasetsDict[signalDatasetName]
        # print(signalDatasetsDict)
        TMVA.Tools.Instance()
        numVars = len(variableList)
        if normalizeVars:
            numVars -= 1
        bdt = getattr(ROOT, "BDT{}".format(lqMassToUse), None)
        if bdt is None:
            gInterpreter.ProcessLine(('''
            TMVA::Experimental::RReader BDT{}("{}");
            computeBDT{} = TMVA::Experimental::Compute<{}, float>(BDT{});
            ''').format(lqMassToUse, bdtWeightFileName, lqMassToUse, numVars, lqMassToUse))
        sys.stdout.flush()
        sys.stderr.flush()
        bkgMvaValues = []
        bkgWeights = []
        # backgrounds
        varHistsBkg = {}
        for var in variableList+["BDT"]:
            if var == "LQCandidateMass":
                continue
            varHistsBkg[var] = TH1D(var+"_bkg", var+" bkg", *variableHistInfo[var])            
        for sample in backgroundDatasetsDict.keys():
            for idx, txtFile in enumerate(backgroundDatasetsDict[sample]):
                txtFile = txtFile.format(year, lqMassToUse)
                #tchainBkg = LoadChainFromTxtFile(txtFile.format(lqMassToUse))
                tchainBkg = LoadChainFromTxtFile(txtFile)
                if tchainBkg is None:
                    continue
                df = RDataFrame(tchainBkg)
                df = df.Filter(mycutb.GetTitle())  # will work for expressions valid in C++
                if "LQCandidateMass" in variableList:
                    df = df.Define("massInt", str(lqMassToUse))
                    df = df.Redefine("LQCandidateMass", "Numba::GetMassFloat(massInt)")
                # df = df.Define('BDTv', ROOT.computeBDT, ROOT.BDT.GetVariableNames())
                # df = df.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames())
                varNames = getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames()
                if normalizeVars:
                    l_varn = ROOT.std.vector['std::string']()
                    for i_expr, expr in enumerate(varNames):
                        varname = 'v_{}'.format(i_expr)
                        l_varn.push_back(varname)
                        df=df.Define(varname, '(float)({})'.format(expr))
                    df = df.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), l_varn)
                else:
                    df = df.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), varNames)
                df = df.Define('BDT', 'BDTv[0]')
                df = df.Define('eventWeight', eventWeightExpression)
                # df = df.Filter("BDT > -0.2 && BDT < 0.2")
                datasetName = os.path.basename(txtFile).replace(".txt", "")
                if "data" not in sample.lower():
                    sumWeights = GetBackgroundSumWeights(datasetName, txtFile)
                    bkgWeight = CalcWeight(datasetName, intLumi, sumWeights)
                else:
                    bkgWeight = 1.0
                df = df.Define("datasetWeight", str(bkgWeight))
                df = df.Redefine("eventWeight", "eventWeight*datasetWeight")
                bkgMvaValues.extend(df.Take["float"]("BDT").GetValue())
                bkgWeights.extend(df.Take["double"]("eventWeight").GetValue())
                
                st = time.time()
                histosToRun = []
                hbkg = {}
                histBkg = {}
                for var in variableList+["BDT"]:
                    if var == "LQCandidateMass":
                        continue
                    histName = "{}_{}_{}".format(var, sample, idx)
                    hbkg[var] = TH1D(histName, histName, *variableHistInfo[var])
                    histBkg[var] = df.Histo1D(ROOT.RDF.TH1DModel(hbkg[var]), var, "eventWeight")
                    # histBkg.Scale(bkgWeight)
                    histosToRun.append(histBkg[var])
                ROOT.RDF.RunGraphs(histosToRun)
                for var in variableList+["BDT"]:
                    if var == "LQCandidateMass":
                        continue
                    varHistsBkg[var].Add(histBkg[var].GetPtr())
                et = time.time()
                timeMakingHistos += et - st    
        rootFile.cd("LQM{}".format(lqMassToUse))

        st = time.time()
        for varHist in varHistsBkg.values():
            varHist.Write()
        et = time.time()
        timeMakingHistos+= et - st
        # signal
        tchainSig = LoadDatasets(signalDatasetsDict, neededBranches, signal=True, loader=None, lqMass=lqMassToUse, year=year)
        dfSig = RDataFrame(tchainSig)
        dfSig = dfSig.Filter(mycuts.GetTitle())  # will work for expressions valid in C++
        # dfSig = dfSig.Define('BDTv', ROOT.computeBDT, ROOT.BDT.GetVariableNames())
        # dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), getattr(ROOT, "BDT{}".format(lqMassToUse)).GetVariableNames())
        if normalizeVars:
            l_varn = ROOT.std.vector['std::string']()
            for i_expr, expr in enumerate(varNames):
                varname = 'v_{}'.format(i_expr)
                l_varn.push_back(varname)
                dfSig=dfSig.Define(varname, '(float)({})'.format(expr))
            dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), l_varn)
        else:
            dfSig = dfSig.Define('BDTv', getattr(ROOT, "computeBDT{}".format(lqMassToUse)), varNames)
        dfSig = dfSig.Define('BDT', 'BDTv[0]')
        dfSig = dfSig.Define('eventWeight', eventWeightExpression)
        # dfSig = dfSig.Filter("BDT > -0.2 && BDT < 0.2")
        sumWeights = GetSignalSumWeights(lqMassToUse, year)
        signalWeight = CalcWeight(signalDatasetName, intLumi, sumWeights)
        dfSig = dfSig.Define("datasetWeight", str(signalWeight))
        dfSig = dfSig.Redefine("eventWeight", "eventWeight*datasetWeight")
        # vars
        rootFile.cd("LQM{}".format(lqMassToUse))
        st = time.time()
        for var in variableList+["BDT"]:
            if var == "LQCandidateMass":
                continue
            hsig = TH1D(var+"_LQ{}".format(lqMassToUse), var+" LQ{}".format(lqMassToUse), *variableHistInfo[var])
            histSig = dfSig.Histo1D(ROOT.RDF.TH1DModel(hsig), var, "eventWeight")
            histSig.Scale(signalWeight)
            histSig.Write()
        et = time.time()
        timeMakingHistos += et - st
        # ROC
        startROCcode = time.time()
        rocCurve = TMVA.ROCCurve(dfSig.Take["float"]("BDT").GetValue(), bkgMvaValues, dfSig.Take["double"]("eventWeight").GetValue(), bkgWeights)
        rocGraph = rocCurve.GetROCCurve(500)
        c = TCanvas("rocLQ"+str(lqMassToUse))
        c.cd()
        rocGraph.Draw("ap")
        c.Update()
        if not os.path.isdir("dataset/plots"):
            os.mkdir("dataset/plots")
        c.Print("dataset/plots/"+"rocLQ"+str(lqMassToUse)+".png")
        rootFile.cd("LQM{}".format(lqMassToUse))
        rocGraph.Write()
        rocAUC = rocCurve.GetROCIntegral(500)
        print("For lqMass={}, got ROC AUC = {}".format(lqMassToUse, rocAUC))

        rootFile.Close()
    except Exception as e:
        print("ERROR: exception in DoROCAndBDTPlots for lqMass={}".format(lqMassToUse))
        traceback.print_exc()
        raise e
    endTime = time.time()
    ROCTime = endTime - startROCcode
    totTime = endTime - startTime
    print("total time = ", totTime)
    print("time spent making and writing histograms = ", timeMakingHistos)
    print("time to make ROC curve = ", ROCTime)
    return True

@ROOT.Numba.Declare(["int"], "float")
def GetMassFloat(mass):
    return float(mass)


def PrintBDTCuts(optValsDict, parametrized):
    for mass, valList in optValsDict.items():
        valListFormatted = ["{:0.4f}".format(i) for i in valList]
        print("For lqMass={}, max FOM: ibin={} with FOM={}, cutVal={}, nS={}, eff={}, nB={}".format(mass, *valListFormatted))
    sortedDict = OrderedDict(sorted(optValsDict.items(), key=lambda t: float(t[1][2])))
    for mass, valList in sortedDict.items():
        print("#"+114*"-")
        print("# LQ M {} optimization".format(mass))
        print("#"+114*"-")
        if parametrized:
            print("BDTOutput_LQ{}                         {}                +inf 		-		-	2	100 -1 1    TMVACut:BDT,BDTWeightFileName,trainingSelection,LQCandidateMass={}".format(mass, round(valList[2], 4), mass))
        else:
            print("BDTOutput_LQ{}                         {}                +inf 		-		-	2	100 -1 1    TMVACut:BDT,BDTWeightFileLQM{},trainingSelection".format(mass, round(valList[2], 4), mass))


def WriteOptimizationHists(rootFileName, optHistsDict, optValsDict, fomInfoDict):
    rebinFactor = 100
    signalColor = kRed
    backgroundColor = kBlue
    rootFile = TFile(rootFileName, "recreate")
    rootFile.cd()
    lqMasses = [float(key) for key in optHistsDict.keys()]
    optCutVals = [float(optValsDict[lqMass][2]) for lqMass in lqMasses]
    optCutValsVsLQMassGraph = TGraph(len(lqMasses), np.array(lqMasses), np.array(optCutVals))
    optCutValsVsLQMassGraph.SetName("optCutValVsLQMass")
    optCutValsVsLQMassGraph.SetTitle("Opt. BDT cut vs. M_{LQ}; M_{LQ} [GeV]")
    optCutValsVsLQMassGraph.SetMarkerStyle(22)
    optCutValsVsLQMassGraph.SetMarkerColor(kSpring-1)
    optCutValsVsLQMassGraph.Draw("ap")
    optCutValsVsLQMassGraph.Write()
    optNSVals = [float(optValsDict[lqMass][3]) for lqMass in lqMasses]
    optNSValsVsLQMassGraph = TGraph(len(lqMasses), np.array(lqMasses), np.array(optNSVals))
    optNSValsVsLQMassGraph.SetName("optNSVsLQMass")
    optNSValsVsLQMassGraph.SetTitle("Opt. cut n_{S} vs. M_{LQ}; M_{LQ} [GeV]; Opt. n_{S}")
    optNSValsVsLQMassGraph.SetMarkerStyle(22)
    optNSValsVsLQMassGraph.SetMarkerColor(signalColor)
    optNSValsVsLQMassGraph.Draw("ap")
    optNSValsVsLQMassGraph.Write()
    optEffVals = [float(optValsDict[lqMass][4]) for lqMass in lqMasses]
    optEffValsVsLQMassGraph = TGraph(len(lqMasses), np.array(lqMasses), np.array(optEffVals))
    optEffValsVsLQMassGraph.SetName("optEffVsLQMass")
    optEffValsVsLQMassGraph.SetTitle("Opt. cut eff. vs. M_{LQ}; M_{LQ} [GeV]; Opt. eff.")
    optEffValsVsLQMassGraph.SetMarkerStyle(22)
    optEffValsVsLQMassGraph.SetMarkerColor(kSpring-1)
    optEffValsVsLQMassGraph.Draw("ap")
    optEffValsVsLQMassGraph.Write()
    optNBVals = [float(optValsDict[lqMass][5]) for lqMass in lqMasses]
    optNBValsVsLQMassGraph = TGraph(len(lqMasses), np.array(lqMasses), np.array(optNBVals))
    optNBValsVsLQMassGraph.SetName("optNBVsLQMass")
    optNBValsVsLQMassGraph.SetTitle("Opt. cut n_{B} vs. M_{LQ}; M_{LQ} [GeV]; Opt. n_{B}")
    optNBValsVsLQMassGraph.SetMarkerStyle(22)
    optNBValsVsLQMassGraph.SetMarkerColor(kSpring-1)
    optNBValsVsLQMassGraph.Draw("ap")
    optNBValsVsLQMassGraph.Write()
    nSBeforeBDTCutList = [fomInfoDict[lqMass]["nSNoBDTCut"] for lqMass in lqMasses]
    nSErrBeforeBDTCutList = [fomInfoDict[lqMass]["nSErrNoBDTCut"] for lqMass in lqMasses]
    nSBeforeBDTCutVsLQMassGraph = TGraphErrors(len(lqMasses), np.array(lqMasses), np.array(nSBeforeBDTCutList), np.zeros(len(lqMasses)), np.array(nSErrBeforeBDTCutList))
    nSBeforeBDTCutVsLQMassGraph.SetName("nSBeforeBDTCutVsLQMass")
    nSBeforeBDTCutVsLQMassGraph.SetTitle("n_{S} w/o BDT cut vs. M_{LQ}; M_{LQ} [GeV]; n_{S}")
    nSBeforeBDTCutVsLQMassGraph.SetMarkerStyle(22)
    nSBeforeBDTCutVsLQMassGraph.SetMarkerColor(signalColor)
    nSBeforeBDTCutVsLQMassGraph.Draw("ap")
    nSBeforeBDTCutVsLQMassGraph.Write()
    nBBeforeBDTCutList = [fomInfoDict[lqMass]["nBNoBDTCut"] for lqMass in lqMasses]
    nBErrBeforeBDTCutList = [fomInfoDict[lqMass]["nBErrNoBDTCut"] for lqMass in lqMasses]
    nBBeforeBDTCutVsLQMassGraph = TGraphErrors(len(lqMasses), np.array(lqMasses), np.array(nBBeforeBDTCutList), np.zeros(len(lqMasses)), np.array(nBErrBeforeBDTCutList))
    nBBeforeBDTCutVsLQMassGraph.SetName("nBBeforeBDTCutVsLQMass")
    nBBeforeBDTCutVsLQMassGraph.SetTitle("n_{B} w/o BDT cut vs. M_{LQ}; M_{LQ} [GeV]; n_{B}")
    nBBeforeBDTCutVsLQMassGraph.SetMarkerStyle(22)
    nBBeforeBDTCutVsLQMassGraph.SetMarkerColor(backgroundColor)
    nBBeforeBDTCutVsLQMassGraph.Draw("ap")
    nBBeforeBDTCutVsLQMassGraph.Write()
    for lqMass, optHistsList in optHistsDict.items():
        rootFile.mkdir("LQM"+str(lqMass))
        rootFile.cd("LQM"+str(lqMass))
        optCutVal = float(optValsDict[lqMass][2])
        c = TCanvas("optLQ"+str(lqMass))
        c.cd()
        c.SetLogy()
        signalHist = optHistsList[0]
        signalHist.SetLineWidth(2)
        signalHist.SetLineColor(signalColor)
        signalHist.SetMarkerColor(signalColor)
        signalHist.Write()
        bkgHist = optHistsList[1]
        bkgHist.SetLineWidth(2)
        bkgHist.SetLineColor(backgroundColor)
        bkgHist.Write()
        signalUnweightedHist = optHistsList[2]
        signalUnweightedHist.SetLineWidth(2)
        signalUnweightedHist.SetLineColor(signalColor)
        signalUnweightedHist.SetMarkerColor(signalColor)
        signalUnweightedHist.Write()
        bkgUnweightedHist = optHistsList[3]
        bkgUnweightedHist.SetLineWidth(2)
        bkgUnweightedHist.SetLineColor(backgroundColor)
        bkgUnweightedHist.Write()
        # now make overlay canvas
        signalHist.Rebin(rebinFactor)
        bkgHist.Rebin(rebinFactor)
        bkgHist.SetTitle("Classifier output for M_{LQ} = " + str(lqMass) + " GeV")
        minY = min(signalHist.GetMinimum(), bkgHist.GetMinimum())
        #maxY = 1.1*max(signalHist.GetYaxis().GetXmax(), bkgHist.GetYaxis().GetXmax())
        maxY = 1.1*max(signalHist.GetMaximum(), bkgHist.GetMaximum())
        bkgHist.GetYaxis().SetRangeUser(minY, maxY)
        bkgHist.Draw()
        signalHist.Draw("sames")
        line = TLine(optCutVal, minY, optCutVal, maxY)
        line.SetLineColor(kSpring-5)
        line.SetLineWidth(2)
        line.Draw()
        c.Modified()
        c.Update()
        c.Write()
        #
        fomVals = fomInfoDict[lqMass]["FOM"]
        bdtCutVals = fomInfoDict[lqMass]["cutVal"]
        fomValVsBDTCutGraph = TGraph(len(bdtCutVals), np.array(bdtCutVals), np.array(fomVals))
        fomValVsBDTCutGraph.SetName("fomValVsBDTCutGraphLQM"+str(lqMass))
        fomValVsBDTCutGraph.SetTitle("FOM value vs. BDT cut for M_{LQ} = "+str(lqMass)+"; M_{LQ}; fom.")
        fomValVsBDTCutGraph.SetMarkerStyle(21)
        fomValVsBDTCutGraph.SetMarkerColor(kSpring-5)
        fomValVsBDTCutGraph.Draw("ap")
        fomValVsBDTCutGraph.Write()
        nSigVals = fomInfoDict[lqMass]["nS"]
        nSigErrVals = fomInfoDict[lqMass]["nSErr"]
        nSigVsBDTCutGraph = TGraphErrors(len(bdtCutVals), np.array(bdtCutVals), np.array(nSigVals), np.zeros(len(bdtCutVals)), np.array(nSigErrVals))
        nSigVsBDTCutGraph.SetName("nSigVsBDTCutGraphLQM"+str(lqMass))
        nSigVsBDTCutGraph.SetTitle("N_{S} vs. BDT cut for M_{LQ} = "+str(lqMass)+"; M_{LQ}; N_{S}")
        nSigVsBDTCutGraph.SetMarkerStyle(26)
        nSigVsBDTCutGraph.SetMarkerColor(backgroundColor)
        nSigVsBDTCutGraph.Draw("ap")
        nSigVsBDTCutGraph.Write()
        effVals = fomInfoDict[lqMass]["eff"]
        effErrVals = fomInfoDict[lqMass]["effErr"]
        effVsBDTCutGraph = TGraphErrors(len(bdtCutVals), np.array(bdtCutVals), np.array(effVals), np.zeros(len(bdtCutVals)), np.array(effErrVals))
        effVsBDTCutGraph.SetName("effVsBDTCutGraphLQM"+str(lqMass))
        effVsBDTCutGraph.SetTitle("Sig. eff. vs. BDT cut for M_{LQ} = "+str(lqMass)+"; M_{LQ}; eff.")
        effVsBDTCutGraph.SetMarkerStyle(22)
        effVsBDTCutGraph.SetMarkerColor(signalColor)
        effVsBDTCutGraph.Draw("ap")
        effVsBDTCutGraph.Write()
        nBkgVals = fomInfoDict[lqMass]["nB"]
        nBkgErrVals = fomInfoDict[lqMass]["nBErr"]
        nBkgVsBDTCutGraph = TGraphErrors(len(bdtCutVals), np.array(bdtCutVals), np.array(nBkgVals), np.zeros(len(bdtCutVals)), np.array(nBkgErrVals))
        nBkgVsBDTCutGraph.SetName("nBkgVsBDTCutGraphLQM"+str(lqMass))
        nBkgVsBDTCutGraph.SetTitle("N_{B} vs. BDT cut for M_{LQ} = "+str(lqMass)+"; M_{LQ}; n_{B}")
        nBkgVsBDTCutGraph.SetMarkerStyle(23)
        nBkgVsBDTCutGraph.SetMarkerColor(backgroundColor)
        nBkgVsBDTCutGraph.Draw("ap")
        nBkgVsBDTCutGraph.Write()
    rootFile.Close()


def GetBackgroundDatasets(inputListBkgBase):
    # preselection-skimmed background datasets
    #FIXME: comment out datasets with zero events; should handle this a bit better
    #       this came from the fact that the makeBDTTrainingTrees script doesn't write out files for trees with zero entries
    backgroundDatasetsDict = {
            "ZJet_amcatnlo_ptBinned" : [
                #inclusive stitched yields no events
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                ],
            #"ZJet_powhegminnlo" : [
                #inputListBkgBase+"DYJetsToEE_M-1000to1500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-100to200_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                # inputListBkgBase+"DYJetsToEE_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-1500to2000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-2000toInf_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-200to400_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-400to500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-500to700_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-700to800_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #inputListBkgBase+"DYJetsToEE_M-800to1000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
                #],
            "TTbar_powheg" : [
                inputListBkgBase+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.txt",
                #inputListBkgBase+"TTToHadronic_TuneCP5_13TeV-powheg-pythia8.txt",
                # inputListBkgBase+"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.txt",
                ],
            "DIBOSON_nlo" : [
                #inputListBkgBase+"WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                # inputListBkgBase+"WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"ZZTo4L_TuneCP5_13TeV_powheg_pythia8.txt",
                #inputListBkgBase+"ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8.txt",
                # inputListBkgBase+"WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                #inputListBkgBase+"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8.txt",
                # inputListBkgBase+"WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                ],
            # "TRIBOSON" : [
            #     inputListBkgBase+"WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"WZZ_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"ZZZ_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     ],
            # "TTW" : [
            #     inputListBkgBase+"TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.txt",
            #     inputListBkgBase+"TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.txt",
            #     ],
            # "TTZ" : [
            #     inputListBkgBase+"ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8.txt"
            #     ],
            "SingleTop" : [
                inputListBkgBase+"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                # inputListBkgBase+"ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                # inputListBkgBase+"ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                # inputListBkgBase+"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.txt",
                ],
            # "WJet_amcatnlo_jetBinned" : [
            #     #inputListBkgBase+"WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     #inputListBkgBase+"WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     ],
            #"PhotonJets_Madgraph" : [
            #    #inputListBkgBase+"GJets_HT-40To100_madgraphMLM.txt",
            #    #inputListBkgBase+"GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #    inputListBkgBase+"GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #    inputListBkgBase+"GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #    inputListBkgBase+"GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #    ],
    }
    return backgroundDatasetsDict


def GetQCDDatasetsDict(inputListQCD1FRBase, inputListQCD2FRBase, year):
    # if year == "2016preVFP":
    #     qcdFakes = {
    #             "QCDFakes_DATA" : [
    #                 #inputListQCDBase+"SinglePhoton_Run2016B-ver1_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 ]
    #     }
    # elif year == "2016postVFP":
    #     qcdFakes = {
    #             "QCDFakes_DATA" : [
    #                 inputListQCDBase+"SinglePhoton_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2.txt",
    #                 inputListQCDBase+"SinglePhoton_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1.txt",
    #                 ]
    #     }
    qcdFakes = {
        "QCDFakes_DYJ" : [
            #inputListQCD1FRBase+"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            ],
        #"QCDFakes_DYJ" : [
            #inputListQCD1FRBase+"DYJetsToEE_M-1000to1500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-100to200_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            # inputListQCD1FRBase+"DYJetsToEE_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-1500to2000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-2000toInf_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-200to400_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-400to500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-500to700_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-700to800_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #inputListQCD1FRBase+"DYJetsToEE_M-800to1000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
            #]
    }
    if year == "2016preVFP":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"SinglePhoton_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v4.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
            inputListQCD2FRBase+"SinglePhoton_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v4.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            ]
        qcdFakes["QCDFakes_DYJ"] = [txtFile.replace(".txt", "_APV.txt") for txtFile in qcdFakes["QCDFakes_DYJ"]]
    elif year == "2016postVFP":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"SinglePhoton_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
                    inputListQCD2FRBase+"SinglePhoton_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    ]
    return qcdFakes

    
if __name__ == "__main__":
    ####################################################################################################
    # Configurables
    ####################################################################################################
    variableList = [
        "sT_eejj",
        "PFMET_Type1_Pt",
        # "PFMET_Type1_Phi",
        "M_e1e2",
        "M_e1j1",
        "M_e1j2",
        "M_e2j1",
        "M_e2j2",
        "Ele1_Pt",
        "Ele2_Pt",
#        "Ele1_Eta",
#        "Ele2_Eta",
#        "Ele1_Phi",
#        "Ele2_Phi",
        "Jet1_Pt",
        "Jet2_Pt",
        "Jet3_Pt",
#        "Jet1_Eta",
#        "Jet2_Eta",
#        "Jet3_Eta",
#        "Jet1_Phi",
#        "Jet2_Phi",
#        "Jet3_Phi",
#        "DR_Ele1Jet1",
#        "DR_Ele1Jet2",
#        "DR_Ele2Jet1",
#        "DR_Ele2Jet2",
#        "DR_Jet1Jet2",
#        "Masym",
        "MejMin",
        "MejMax",
        "Meejj",
        "LQCandidateMass"
    ]
    variableHistInfo = {}
    variableHistInfo["sT_eejj"] = [700, 0, 7000]  # nbins, min, max
    variableHistInfo["PFMET_Type1_Pt"] = [200, 0, 2000]
    variableHistInfo["PFMET_Type1_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["M_e1e2"] = [500, 0, 5000]
    variableHistInfo["M_e1j1"] = [500, 0, 5000]
    variableHistInfo["M_e1j2"] = [500, 0, 5000]
    variableHistInfo["M_e2j1"] = [500, 0, 5000]
    variableHistInfo["M_e2j2"] = [500, 0, 5000]
    variableHistInfo["Ele1_Pt"] = [200, 0, 2000]
    variableHistInfo["Ele2_Pt"] = [200, 0, 2000]
    variableHistInfo["Ele1_Eta"] = [100, -5, 5]
    variableHistInfo["Ele2_Eta"] = [100, -5, 5]
    variableHistInfo["Ele1_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["Ele2_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["Jet1_Pt"] = [200, 0, 2000]
    variableHistInfo["Jet2_Pt"] = [200, 0, 2000]
    variableHistInfo["Jet3_Pt"] = [200, 0, 2000]
    variableHistInfo["Jet1_Eta"] = [100, -5, 5]
    variableHistInfo["Jet2_Eta"] = [100, -5, 5]
    variableHistInfo["Jet3_Eta"] = [100, -5, 5]
    variableHistInfo["Jet1_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["Jet2_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["Jet3_Phi"] = [60, -3.1416, 3.1416]
    variableHistInfo["DR_Ele1Jet1"] = [100, 0, 10]
    variableHistInfo["DR_Ele1Jet2"] = [100, 0, 10]
    variableHistInfo["DR_Ele2Jet1"] = [100, 0, 10]
    variableHistInfo["DR_Ele2Jet2"] = [100, 0, 10]
    variableHistInfo["DR_Jet1Jet2"] = [100, 0, 10]
    variableHistInfo["Masym"] = [200, 0, 2000]
    variableHistInfo["MejMin"] = [350, 0, 3500]
    variableHistInfo["MejMax"] = [350, 0, 3500]
    variableHistInfo["Meejj"] = [700, 0, 7000]
    variableHistInfo["BDT"] = [200, -1, 1]
    eventWeightExpression = "EventWeight"
    #eventWeightExpression = "Weight*PrefireWeight*puWeight*FakeRateEffective*MinPrescale*Ele1_RecoSF*Ele2_RecoSF*EventTriggerScaleFactor*ZVtxSF"
    ## EGM Loose ID
    ##eventWeightExpression += "*Ele1_EGMLooseIDSF*Ele2_EGMLooseIDSF"
    ## HEEP
    #eventWeightExpression += "*Ele1_HEEPSF*Ele2_HEEPSF"
    #
    # neededBranches = ["Weight", "PrefireWeight", "puWeight", "FakeRateEffective", "MinPrescale", "Ele1_RecoSF", "Ele2_RecoSF", "EventTriggerScaleFactor", "ZVtxSF"]
    neededBranches = ["EventWeight"]
    neededBranches.extend(["run", "ls", "event"])
    # loose
    #neededBranches.extend(["Ele1_EGMLooseIDSF", "Ele2_EGMLooseIDSF"])
    # HEEP
    #neededBranches.extend(["Ele1_HEEPSF", "Ele2_HEEPSF"])
    neededBranches.extend(variableList)
    # Apply additional cuts on the signal and background samples (can be different)
    mycuts = TCut("M_e1e2 > 220 && sT_eejj > 400 && PassTrigger==1")
    mycutb = mycuts
    # mycuts = TCut("M_e1e2 > 200 && sT_eejj > 400 && PassTrigger==1")
    # mycutb = TCut("M_e1e2 > 200 && sT_eejj > 400 && PassTrigger==1")
    # mycuts = TCut("M_e1e2 > 220 && PassTrigger==1")
    # mycutb = TCut("M_e1e2 > 220 && PassTrigger==1")
        
    result_list = []
    logString = "INFO: running {} parallel jobs for {} separate LQ masses requested..."
    
    
    ####################################################################################################
    # Run
    ####################################################################################################
    usage = "usage: %prog [options] year"
    parser = OptionParser(usage=usage)
    # parser.add_option(
    #     "-y",
    #     "--year",
    #     dest="year",
    #     help="year",
    #     metavar="YEAR",
    # )
    parser.add_option(
        "-t",
        "--train",
        action="store_true",
        dest="train",
        default=False,
        help="train the BDT",
        metavar="TRAIN",
    )
    parser.add_option(
        "-o",
        "--optimize",
        action="store_true",
        dest="optimize",
        default=False,
        help="optimize the BDT",
        metavar="OPTIMIZE",
    )
    parser.add_option(
        "-r",
        "--roc",
        action="store_true",
        dest="roc",
        default=False,
        help="make ROC plots",
        metavar="ROC",
    )
    (options, args) = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
        raise RuntimeError("Must specify year")
    year = args[0]

    gROOT.SetBatch()
    dateStr = "9oct2023"
    inputListBkgBase = os.getenv("LQANA")+"/config/myDatasets/BDT/{}_amcatnloDY/{}/"
    inputListQCD1FRBase = os.getenv("LQANA")+"/config/myDatasets/BDT/{}_amcatnloDY/{}/QCDFakes_1FR/"
    inputListQCD2FRBase = os.getenv("LQANA")+"/config/myDatasets/BDT/{}_amcatnloDY/{}/QCDFakes_DATA_2FR/"
    backgroundDatasetsDict = GetBackgroundDatasets(inputListBkgBase)
    xsectionFiles = dict()
    xsectionTxt = "xsection_13TeV_2022_Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar_Mee_BkgControlRegion_DYJets_2016preVFP.txt"
    xsectionFiles["2016preVFP"] = os.getenv("LQANA")+"/config/" + xsectionTxt
    xsectionFiles["2016postVFP"] = os.getenv("LQANA")+"/versionsOfAnalysis/2016postVFP/eejj/eejj_4oct2023_heep_preselOnly/" + xsectionTxt
    train = options.train
    optimize = options.optimize
    roc = options.roc
    parallelize = True
    parametrized = False
    includeQCD = True
    normalizeVars = False
    # normTo = "Meejj"
    #lqMassesToUse = [1000]#,1100,1200]
    #lqMassesToUse = list(range(2000, 3100, 100))
    #lqMassesToUse = list(range(600, 1100, 100))
    #lqMassesToUse = list(range(800, 1300, 100))
    #lqMassesToUse = list(range(300, 1000, 100))
    #lqMassesToUse = list(range(700, 1100, 100))
    #lqMassesToUse = list(range(300, 1100, 100))
    lqMassesToUse = list(range(300,3100,100))
    #lqMassesToUse = [300]#,400,500,600]
    signalNameTemplate = "LQToDEle_M-{}_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8"
    weightFile = os.path.abspath(os.getcwd())+"/dataset/weights/TMVAClassification_BDTG.weights.xml"
    # weightFile = "dataset/weights/TMVAClassification_"+signalNameTemplate.format(300)+"_APV_BDTG.weights.xml"
    # weightFile = "dataset/weights/TMVAClassification_"+signalNameTemplate.format(1500)+"_APV_BDTG.weights.xml"
    optimizationPlotFile = os.path.abspath(os.getcwd())+"/optimizationPlots.root"
    bdtPlotFile = os.path.abspath(os.getcwd())+"/bdtPlots.root"

    if not parametrized:
        neededBranches.remove("LQCandidateMass")
        variableList.remove("LQCandidateMass")
    xsectionFile = xsectionFiles[year]
    if year == "2016preVFP":
        intLumi = 19497.897120
        for dataset in backgroundDatasetsDict.keys():
            if "QCDFakes_DATA" in dataset:
                continue
            backgroundDatasetsDict[dataset] = [txtFile.replace(".txt", "_APV.txt") for txtFile in backgroundDatasetsDict[dataset]]
        if includeQCD:
            backgroundDatasetsDict.update(GetQCDDatasetsDict(inputListQCD1FRBase, inputListQCD2FRBase, year))
        signalNameTemplate+="_APV"
    elif year == "2016postVFP":
        intLumi = 16812.151722
        if includeQCD:
            backgroundDatasetsDict.update(GetQCDDatasetsDict(inputListQCD1FRBase, inputListQCD2FRBase, year))
    elif year == "2017":
        intLumi = 41477.877399
    elif year == "2018":
        intLumi = 59827.449483
    else:
        raise RuntimeError("Did not understand 'year' parameter whose value is {}; must be one of 2016preVFP, 2016postVFP, 2017, 2018".format(year))

    ParseXSectionFile(xsectionFile)
    inputListSignalBase = inputListBkgBase
    allSignalDatasetsDict = {}
    massList = list(range(300, 3100, 100))
    massList.extend([3500, 4000])
    for mass in massList:
        signalName = signalNameTemplate.format(mass)
        allSignalDatasetsDict[signalName] = [inputListSignalBase+signalName+".txt"]

    if train:
        print("INFO: Begin {} training.".format("parametrized BDT" if parametrized else "BDT"))
        if parametrized:
            TrainParametrizedBDT(lqMassesToUse, year)
        else:
            if parallelize:
                # ncores = multiprocessing.cpu_count()
                ncores = 4  # only use 4 parallel jobs to be nice
                pool = multiprocessing.Pool(ncores)
                jobCount = 0
                for mass in lqMassesToUse:
                    try:
                        pool.apply_async(TrainBDT, [[mass, year, normalizeVars]], callback=log_result)
                        jobCount += 1
                    except KeyboardInterrupt:
                        print("\n\nCtrl-C detected: Bailing.")
                        pool.terminate()
                        exit(-1)
                    except Exception as e:
                        print("ERROR: caught exception in job for LQ mass: {}; exiting".format(mass))
                        traceback.print_exc()
                        exit(-2)
                
                # now close the pool and wait for jobs to finish
                pool.close()
                sys.stdout.write(logString.format(jobCount, len(lqMassesToUse)))
                sys.stdout.write("\t"+str(len(result_list))+" jobs done")
                sys.stdout.flush()
                pool.join()
                # check results?
                if len(result_list) < jobCount:
                    raise RuntimeError("ERROR: {} jobs had errors. Exiting.".format(jobCount-len(result_list)))
            else:
                for mass in lqMassesToUse:
                    TrainBDT([mass, year, normalizeVars])
        print("INFO: Training {} done.".format("parametrized BDT" if parametrized else "BDT"))
    
    if optimize:
        print("INFO: Begin optimization.")
        manager = multiprocessing.Manager()
        dictOptValues = manager.dict()
        dictOptHists = manager.dict()
        dictOptFOMInfo = manager.dict()
        if parallelize:
            # ncores = multiprocessing.cpu_count()
            ncores = 4  # only use 4 parallel jobs to be nice
            pool = multiprocessing.Pool(ncores)
            jobCount = 0
            for mass in lqMassesToUse:
                if not parametrized:
                    signalDatasetName = signalNameTemplate.format(mass)
                    weightFile = "dataset/weights/TMVAClassification_"+signalDatasetName+"_BDTG.weights.xml"
                dictOptFOMInfo[mass] = manager.dict()
                try:
                    pool.apply_async(OptimizeBDTCut, [[weightFile.format(mass), mass, dictOptValues, dictOptHists, dictOptFOMInfo, year]], callback=log_result)
                    jobCount += 1
                except KeyboardInterrupt:
                    print("\n\nCtrl-C detected: Bailing.")
                    pool.terminate()
                    exit(-1)
                except Exception as e:
                    print("ERROR: caught exception in job for LQ mass: {}; exiting".format(mass))
                    traceback.print_exc()
                    exit(-2)
            
            # now close the pool and wait for jobs to finish
            pool.close()
            sys.stdout.write(logString.format(jobCount, len(lqMassesToUse)))
            sys.stdout.write("\t"+str(len(result_list))+" jobs done")
            sys.stdout.flush()
            pool.join()
            # check results?
            if len(result_list) < jobCount:
                raise RuntimeError("ERROR: {} jobs had errors. Exiting.".format(jobCount-len(result_list)))
        else:
            for mass in lqMassesToUse:
                dictOptFOMInfo[mass] = manager.dict()
                OptimizeBDTCut([weightFile.format(mass), mass, dictOptValues, dictOptHists, dictOptFOMInfo, year])
        PrintBDTCuts(dictOptValues, parametrized)
        WriteOptimizationHists(optimizationPlotFile, dictOptHists, dictOptValues, dictOptFOMInfo)

    if roc:
        print("INFO: Do ROC and BDT plots.")
        if os.path.isfile(bdtPlotFile):
            os.unlink(bdtPlotFile)
        # if parallelize:
        #     ncores = 4  # only use 4 parallel jobs to be nice
        #     pool = multiprocessing.Pool(ncores)
        #     jobCount = 0
        #     for mass in lqMassesToUse:
        #         try:
        #             pool.apply_async(DoROCAndBDTPlots, [[bdtPlotFile, weightFile.format(mass), mass, year]], callback=log_result)
        #             jobCount += 1
        #         except KeyboardInterrupt:
        #             print("\n\nCtrl-C detected: Bailing.")
        #             pool.terminate()
        #             exit(-1)
        #         except Exception as e:
        #             print("ERROR: caught exception in job for LQ mass: {}; exiting".format(mass))
        #             traceback.print_exc()
        #             exit(-2)
        #     
        #     # now close the pool and wait for jobs to finish
        #     pool.close()
        #     sys.stdout.write(logString.format(jobCount, len(lqMassesToUse)))
        #     sys.stdout.write("\t"+str(len(result_list))+" jobs done")
        #     sys.stdout.flush()
        #     pool.join()
        #     # check results?
        #     if len(result_list) < jobCount:
        #         raise RuntimeError("ERROR: {} jobs had errors. Exiting.".format(jobCount-len(result_list)))
        # else:
        #     for mass in lqMassesToUse:
        #         DoROCAndBDTPlots([bdtPlotFile, weightFile.format(mass), mass, year])
        for mass in lqMassesToUse:
            if not parametrized:
                signalDatasetName = signalNameTemplate.format(mass)
                weightFile = "dataset/weights/TMVAClassification_"+signalDatasetName+"_BDTG.weights.xml"
            DoROCAndBDTPlots([bdtPlotFile, weightFile, mass, year])
