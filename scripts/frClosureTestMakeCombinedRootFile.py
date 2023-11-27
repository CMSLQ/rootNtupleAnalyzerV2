
from __future__ import print_function

from ROOT import (
    gROOT,
    TFile,
    TH1D,
    TCanvas,
    TColor,
    TLegend
)
import ROOT
import numpy as np
import copy
import ctypes
import sys
import os
import math

#runs on the output of combinePlots for 1P1F and 2F. Copies the histos you need for the closure test into a single ROOT file. Also adds all the MC together to make an MC total histo, and puts that in the output file too.
input_file2F = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/2F/output_cutTable_lq_QCD_FakeRateClosureTest/analysisClass_lq_QCD_FakeRateClosureTest_plots.root"
#input_file1P1F = "$LQDATA/2016/qcdFRClosureTest/frClosureTest_2016pre_Aug23/MCfix/output_cutTable_lq_QCD_FakeRateClosureTest/analysisClass_lq_QCD_FakeRateClosureTest_plots.root"
input_file1P1F = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/1P1F/output_cutTable_lq_QCD_FakeRateClosureTest/analysisClass_lq_QCD_FakeRateClosureTest_plots.root"
output_name = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/FRCTCombined.root"

#run
gROOT.SetBatch(True)

output_file = TFile(output_name,"recreate")
variableNameList = [
    "Pt1stEle_PAS",
    "Me1j1_PAS",
    "Mee_PAS",
    "sT_PAS",
    "MET_PAS",
    "Me2j1_PAS",
    "Pt2ndEle_PAS", 
    "HT",
    "Mt_MET_Ele1_PAS",
    "Mt_MET_Ele2_PAS",
    "Phi1stEle_PAS",
    "Phi2ndEle_PAS",
    "METPhi_PAS",
    "nJet_PAS",
    "Pt1stEle_tight", 
    "Me1j1_tight", 
    "Mee_tight", 
    "sT_tight", 
    "MET_tight", 
    "Me2j1_tight", 
    "Pt2ndEle_tight",
    "HT_tight",
    "Mt_MET_Ele1_tight",
    "Mt_MET_Ele2_tight",
    "Phi1stEle_tight",
    "Phi2ndEle_tight",
    "METPhi_tight",
    "MeeControlReg",
]
histoNameData = "histo1D__QCDFakes_DATA__"
#both eles barrel, one barrel one endcap, both eles endcap
etaRegions = ["_BB","_BE","_EE",""]
fileList = [input_file2F, input_file1P1F]

mcSamples = [
    "ZJet_NNLO_IncStitch",
    "WJet_sherpa", 
    "WJet_amcatnlo_jetBinned",
    "WJet_HTBinned",
    "WJet_HTBinned_IncStitch",
    "WJet_madgraph_Inc",
    "TTBar_powheg", 
    "SingleTop", 
    "GJets", 
    "DIBOSON_nlo",
]#script makes separate MCTotal histo for each WJet sample.
histoNamesMC = []
for name in mcSamples:
    histoNamesMC.append("histo1D__"+name+"__")

for filename in fileList:
    tfile = TFile.Open(filename)
    print("opening "+filename)
    for variableName in variableNameList:
        for region in etaRegions:
            if "control" in variableName.lower() and not region == "":
                continue
            histoData = tfile.Get(histoNameData+variableName+region)
            if not histoData:
                print("ERROR: could not find histo ",histoNameData+variableName+region)
            output_file.cd()
            if "1P1F" in filename:
                histoData.SetName("histo1D__QCDFakes_DATA1P1F__"+variableName+region)
                histoMCTotal = 0
                isFirstMCHist = True
                for name in histoNamesMC:
                    histo = tfile.Get(name + variableName+region)
                    if not histo:
                        print("ERROR: could not find histo ", name + variableName+region)
                    print("writing histo "+histo.GetName())
                    histo.Write()
                    if isFirstMCHist:
                        histoMCTotal_WHT = copy.deepcopy(histo)
                        histoMCTotal_WJetBin = copy.deepcopy(histo)
                        histoMCTotal_WSherpa = copy.deepcopy(histo)
                        histoMCTotal_WHTincStitch = copy.deepcopy(histo)
                        histoMCTotal_WHT.SetName("histo1D__MCTotal-WHTBinned__"+variableName+region)
                        histoMCTotal_WJetBin.SetName("histo1D__MCTotal-WJetBinned__"+variableName+region)
                        histoMCTotal_WSherpa.SetName("histo1D__MCTotal-WSherpa__"+variableName+region)
                        histoMCTotal_WHTincStitch.SetName("histo1D__MCTotal-WHTBinnedIncStitch__"+variableName+region)
                        isFirstMCHist = False
                    else:
                        if not "WJet" in name:
                            histoMCTotal_WHT.Add(histo)
                            histoMCTotal_WJetBin.Add(histo)
                            histoMCTotal_WSherpa.Add(histo)
                            histoMCTotal_WHTincStitch.Add(histo)
                        else:
                            if "ht" in name.lower() and not "stitch" in name.lower():
                                histoMCTotal_WHT.Add(histo)
                            elif "jetbinned" in name.lower():
                                histoMCTotal_WJetBin.Add(histo)
                            elif "sherpa" in name.lower():
                                histoMCTotal_WSherpa.Add(histo)
                            elif "ht" in name.lower() and "stitch" in name.lower():
                                histoMCTotal_WHTincStitch.Add(histo)
                            else:
                                print("cannot determine WJet sample from MC name ", name)
                print("writing histo "+histoMCTotal_WHT.GetName())
                if "control" in variableName.lower():
                    histoAllBkgContReg = copy.deepcopy(histoMCTotal_WHTincStitch)
                histoMCTotal_WHT.Write()
                print("writing histo "+histoMCTotal_WJetBin.GetName())
                histoMCTotal_WJetBin.Write()
                print("writing histo "+histoMCTotal_WSherpa.GetName())
                histoMCTotal_WSherpa.Write()
                print("writing histo "+histoMCTotal_WHTincStitch.GetName())
                histoMCTotal_WHTincStitch.Write()
            elif "2F" in filename:
                histoData2F = copy.deepcopy(histoData)
                histoData2F.SetName("histo1D__QCDFakes_DATA2F__"+variableName+region)
                #each bin in this histo is the sum of the squares of the errors in that bin. 
                #to get the error for the bin, I need to take the sqrt
                #data hist and error hist have the same binning.
                shortVarName = variableName.replace("_PAS","")
                histoNameErr = histoNameData+"errFRsq_"
                #print("histoNameErr: "+histoNameErr+shortVarName+region)
                if not "HT" in variableName and not "control" in variableName.lower() and not "nJet" in variableName and not "Phi" in variableName:
                    histoErrSQ2F = tfile.Get(histoNameErr+shortVarName+region)
                    if not histoErrSQ2F:
                        print("ERROR: could not find histo ",histoNameErr+shortVarName+region)
                    histoErrSQ2F.SetName("histo1D__QCDFakes_DATA2F__errFRsq_"+shortVarName+region)
                    histoData2FWErr = copy.deepcopy(histoData2F)
                    histoData2FWErr.SetName("histo1D__QCDFakes_DATA2F_WError__"+variableName+region)
                #if "Pt" in variableName:
                   # print("bin contents for ", variableName)
                    for i in range(histoData2F.GetNbinsX()):
                        errSQ = histoErrSQ2F.GetBinContent(i)
                        err = math.sqrt(errSQ)
                        #if "Pt" in variableName:
                        #   print("bin ",i,": ",histoData2F.GetBinContent(i)," +/- ",err)
                        histoData2FWErr.SetBinError(i, err)
                else:
                    histoData2FWErr = copy.deepcopy(histoData2F)
                    histoData2FWErr.SetName("histo1D__QCDFakes_DATA2F_WError__"+variableName+region)
                print("writing histo "+histoData2F.GetName())
                if "control" in variableName.lower():
                    histoFRContReg = copy.deepcopy(histoData2F)
                histoData2F.Write()
                print("writing histo "+histoData2FWErr.GetName())
                histoData2FWErr.Write()
                print("writing histo "+histoErrSQ2F.GetName())
                histoErrSQ2F.Write()
            else:
                print("cannot determine region (2F or 1P1F) from filename")
                break
histoAllBkgContReg.SetName("histo1D__AllBkgIncFR__MeeControlReg")
histoAllBkgContReg.Add(histoFRContReg)
print("writing histo "+histoAllBkgContReg.GetName())
histoAllBkgContReg.Write()
print("histos written to "+output_name)
