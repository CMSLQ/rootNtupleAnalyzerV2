#!/usr/bin/env python3

# validate systematics calculations for LHEScale and LHEPDF uncertainty calculations

from ROOT import TFile
import sys
import glob
import math
import copy
import re
import numpy as np
from combineCommon import CalculateShapeSystematic, CalculatePDFVariationHessian, ExtractBranchTitles, IsHistEmpty, AddHistoBins, GetPDFVariationType, GetBranchTitle, CalculatePDFSystematic


# 2016preVFP
scaleFactorsForPiece = {}
# scaleFactorsForPiece["DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 0.07665937664707981/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 0.022305129048936105/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 0.025731861831311143/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 0.07403710668708309/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 0.3728395389392373/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 2.837623603324894/1000.
# scaleFactorsForPiece["DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8"] = 2.6702508813590455/1000.
# scaleFactorsForPiece["ST_tW_top_5f_NoFullyHadronicDecays"] = 5.918836820183422/1000.
# scaleFactorsForPiece["ST_tW_antitop_5f_NoFullyHadronicDecays"] = 6.138639777333125/1000.
scaleFactorsForPiece["WWTo2L2Nu_TuneCP5_13TeV-powheg"] = 7.097138092178493/1000.
scaleFactorsForPiece["ZZTo2Q2L"] = 0.8147357027915243/1000.
scaleFactorsForPiece["ZZTo4L_TuneCP5_13TeV_powheg"] = 0.39241351738237223/1000.
scaleFactorsForPiece["ZZTo2L2Nu_TuneCP5_13TeV_powheg"] = 1.1563916821351499/1000.
scaleFactorsForPiece["WZTo2Q2L"] = 0.8331161555423957/1000.
# scaleFactorsForPiece["WZTo3LNu_powheg"] = 18.0562937086093/1000.
scaleFactorsForPiece["WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg"] = 18.0562937086093/1000.
# overall scale factor
for key in scaleFactorsForPiece.keys():
    if "DYJ" in key:
        scaleFactorsForPiece[key] *= 0.9783443485220265
    elif "TTT" in key:
        scaleFactorsForPiece[key] *= 0.8014439346650566


def GetSystematicsHisto(tfile):
    for key in tfile.GetListOfKeys():
        histoName = key.GetName()
        if histoName.endswith("systematics"):
            hist = key.ReadObj()
            return copy.deepcopy(hist)
    raise RuntimeError("Histogram with name ending in 'systematics' not found in the file '{}'".format(file.GetName()))


def GetTMap(tfile):
    histoKeys = tfile.GetListOfKeys()
    histoList = [x.ReadObj() for x in histoKeys]
    return next((x for x in histoList if x.ClassName() == "TMap" and "systematicNameToBranchesMap" in x.GetName()), None)


def GetScaleFactorForPiece(piece):
    return scaleFactorsForPiece[piece]


# analysisFilePath = sys.argv[1]
# pieceNames = sys.argv[2]  # corresponds to pieces to combine, so similar to dataset
# combinedFile = sys.argv[2]  # after combinePlots

# analysisFilePath = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_11nov2024_presel/cutTable_lq_eejj_preselOnly/condor/"
analysisFilePath = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_4nov2024_bdt_LQToDEle/cutTable_lq_eejj_BDT/condor/"
pieceNames = scaleFactorsForPiece.keys()
# combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_11nov2024_presel/output_cutTable_lq_eejj_preselOnly/unscaled/analysisClass_lq_eejj_ZJet_amcatnlo_ptBinned_IncStitch_plots.root"
# combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_11nov2024_presel/output_cutTable_lq_eejj_preselOnly/unscaled/analysisClass_lq_eejj_ZJet_amcatnlo_Inc_plots.root"
# combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_4nov2024_bdt_LQToDEle/output_cutTable_lq_eejj_BDT/analysisClass_lq_eejj_ZJet_amcatnlo_ptBinned_IncStitch_plots.root"
# combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_4nov2024_bdt_LQToDEle/output_cutTable_lq_eejj_BDT/modded2/analysisClass_lq_eejj_SingleTop_plots.root"
# combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_4nov2024_bdt_LQToDEle/output_cutTable_lq_eejj_BDT/modded2/analysisClass_lq_eejj_WZ_nlo_plots.root"
combinedFile = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_4nov2024_bdt_LQToDEle/output_cutTable_lq_eejj_BDT/modded2/analysisClass_lq_eejj_DIBOSON_nlo_plots.root"
sample = "DIBOSON_nlo"

corrLHESysts = False


print("Using combined file:", combinedFile)

allAnaFiles = []
totalSystHist = None
for piece in pieceNames:
    analysisFiles = glob.glob(analysisFilePath+"/{}*/{}*.root".format(piece, piece))
    print("INFO: found files {} for piece={}".format(analysisFiles, piece))
    analysisTFiles = [TFile.Open(f) for f in analysisFiles]
    if len(analysisFiles) <= 0:
        raise RuntimeError("Did not find any files for piece={} in path={}".format(piece, analysisFilePath))
    pieceHist = None
    for f in analysisTFiles:
        systHist = GetSystematicsHisto(f)
        systHist.Scale(GetScaleFactorForPiece(piece))
        # print("DEBUG: systHist={} from file {}: xbin=1 ybin=2, binc={}, bine={}".format(systHist.GetName(), f.GetName(), systHist.GetBinContent(1, 2), systHist.GetBinError(1, 2)))
        if pieceHist is None:
            pieceHist = systHist
        else:
            success = pieceHist.Add(systHist)
            if not success:
                raise RuntimeError("Failed adding systHist from file {} to pieceHist:".format(f.GetName()))
    # labelsToAdd = ["LHEPdfWeightMC_UpComb", "LHEPdfWeightMC_DownComb", "LHEPdfWeightHessian_NominalComb", "LHEScaleWeight_maxComb", "LHEScaleWeight_maxIndex"]
    labelsToAdd = ["LHEScaleWeight_maxIndex"]
    if not any(substring in list(pieceHist.GetYaxis().GetLabels()) for substring in labelsToAdd):
        pieceHist = AddHistoBins(pieceHist, "y", labelsToAdd)
    verbose = False
    sampleTMap = GetTMap(analysisTFiles[0])
    systNameToBranchTitleDict = ExtractBranchTitles(pieceHist, sampleTMap)
    # LHEScaleWeight systs
    scaleSystDeltas, yBins, maxYIndices = CalculateShapeSystematic(pieceHist, piece, systNameToBranchTitleDict, verbose)
    lheScaleMaxIndexBin = pieceHist.GetYaxis().FindFixBin("LHEScaleWeight_maxIndex")
    lheScaleMaxCombBin = pieceHist.GetYaxis().FindFixBin("LHEScaleWeight_maxComb")
    for xBin in range(0, pieceHist.GetNbinsX()+2):
       nominal = pieceHist.GetBinContent(xBin, 1)  # y-bin 1 always being nominal
       pieceHist.SetBinContent(xBin, lheScaleMaxCombBin, scaleSystDeltas[xBin]+nominal)
       pieceHist.SetBinError(xBin, lheScaleMaxCombBin, scaleSystDeltas[xBin])
       pieceHist.SetBinContent(xBin, lheScaleMaxIndexBin, maxYIndices[xBin])
    # LHEPdfWeight systs
    hessianNominalYields = np.zeros(pieceHist.GetNbinsX()+2)
    pdfSystDeltasUp = np.zeros(pieceHist.GetNbinsX()+2)
    pdfSystDeltasDown = np.zeros(pieceHist.GetNbinsX()+2)
    pdfVariationType = ""
    pdfVariationType, pdfName = GetPDFVariationType(GetBranchTitle("LHEPdfWeight", piece, systNameToBranchTitleDict)[0])
    if "mc" in pdfVariationType:
        pdfSystDeltasUp, pdfSystDeltasDown, yBins = CalculatePDFSystematic(pieceHist, piece, systNameToBranchTitleDict)
    elif "hessian" in pdfVariationType:
        hessianNominalYields = [pieceHist.GetBinContent(xBin, 1) for xBin in range(0, pieceHist.GetNbinsX()+2)]
    pdfWeightMCDownCombBin = pieceHist.GetYaxis().FindFixBin("LHEPdfWeightMC_DownComb")
    pdfWeightMCUpCombBin = pieceHist.GetYaxis().FindFixBin("LHEPdfWeightMC_UpComb")
    pdfWeightHessianCombBin = pieceHist.GetYaxis().FindFixBin("LHEPdfWeightHessian_NominalComb")
    for xBin in range(0, pieceHist.GetNbinsX()+2):
        pieceHist.SetBinContent(xBin, pdfWeightMCUpCombBin, pdfSystDeltasUp[xBin])
        pieceHist.SetBinContent(xBin, pdfWeightMCDownCombBin, pdfSystDeltasDown[xBin])
        pieceHist.SetBinContent(xBin, pdfWeightHessianCombBin, hessianNominalYields[xBin])
        pieceHist.SetBinError(xBin, pdfWeightMCUpCombBin, pdfSystDeltasUp[xBin])
        pieceHist.SetBinError(xBin, pdfWeightMCDownCombBin, pdfSystDeltasDown[xBin])
        pieceHist.SetBinError(xBin, pdfWeightHessianCombBin, hessianNominalYields[xBin])
    if "mc" in pdfVariationType:
        pdfWeightLabels = [label.GetString().Data() for label in pieceHist.GetYaxis().GetLabels() if "LHEPdfWeight" in label.GetString().Data()]
        pdfWeightLabels.remove("LHEPdfWeightMC_UpComb")
        pdfWeightLabels.remove("LHEPdfWeightMC_DownComb")
        pdfWeightLabels.remove("LHEPdfWeightHessian_NominalComb")
        for xBin in range(0, pieceHist.GetNbinsX()+2):
            for label in pdfWeightLabels:
                yBin = pieceHist.GetYaxis().FindFixBin(label)
                pieceHist.SetBinContent(xBin, yBin, 0)
                pieceHist.SetBinError(xBin, yBin, 0)
    #
    if totalSystHist is None:
        totalSystHist = pieceHist
    else:
        success = totalSystHist.Add(pieceHist)
        # print("DEBUG: summed totalHistSyst={} from file {}: xbin=1 ybin=2, binc={}, bine={}".format(totalSystHist.GetName(), f.GetName(), totalSystHist.GetBinContent(1, 2), totalSystHist.GetBinError(1, 2)))
        if not success:
            raise RuntimeError("Failed adding pieceHist for piece={}".format(piece))
    allAnaFiles.extend(analysisTFiles)
if corrLHESysts:
    scaleSystDeltas, yBins, maxYIndices = CalculateShapeSystematic(totalSystHist, piece, systNameToBranchTitleDict, verbose)
    print("DEBUG: scaleSystDeltas: calculate from scale weight variations; deltas[1]={} ".format(scaleSystDeltas[1]))
else:
    scaleSystDeltas = [totalSystHist.GetBinError(xBin, totalSystHist.GetYaxis().FindFixBin("LHEScaleWeight_maxComb")) for xBin in range(0, totalSystHist.GetNbinsX()+2)]
lheScaleMaxIndexBin = pieceHist.GetYaxis().FindFixBin("LHEScaleWeight_maxIndex")
scaleUpCombBin = totalSystHist.GetYaxis().FindFixBin("LHEScale_UpComb")
scaleDownCombBin = totalSystHist.GetYaxis().FindFixBin("LHEScale_DownComb")
for xBin in range(0, totalSystHist.GetNbinsX()+2):
    nominal = totalSystHist.GetBinContent(xBin, 1)  # y-bin 1 always being nominal
    totalSystHist.SetBinContent(xBin, scaleUpCombBin, scaleSystDeltas[xBin]+nominal)
    totalSystHist.SetBinError(xBin, scaleUpCombBin, scaleSystDeltas[xBin])  # for plotting, we rely on the sumQuad addition of bin errors
    totalSystHist.SetBinContent(xBin, scaleDownCombBin, scaleSystDeltas[xBin]+nominal)
    totalSystHist.SetBinError(xBin, scaleDownCombBin, scaleSystDeltas[xBin])
    if corrLHESysts:
        totalSystHist.SetBinContent(xBin, lheScaleMaxIndexBin, maxYIndices[xBin])
    else:
        totalSystHist.SetBinContent(xBin, lheScaleMaxIndexBin, 0.0)

sampleTMap = GetTMap(allAnaFiles[0])
# sampleSystHist = next((x for x in histoList if "systematics" == x.GetName().split("__")[-1] ), None)
#systObjects = [obj.GetName() for obj in sampleHistoDict.values() if "syst" in obj.GetName().lower()]
#print("for sample={}. found systObjects: {}".format(sample, systObjects))
systNameToBranchTitleDict = {}
systNameToBranchTitleDict = ExtractBranchTitles(totalSystHist, sampleTMap)
pdfWeightLabels = [label.GetString().Data() for label in totalSystHist.GetYaxis().GetLabels() if "LHEPdfWeight" in label.GetString().Data() and "comb" not in label.GetString().Data().lower()]
pdfWeightBins = [totalSystHist.GetYaxis().FindFixBin(label) for label in pdfWeightLabels]
pdfMCVarLabels = ["LHEPdfWeightMC_UpComb", "LHEPdfWeightMC_DownComb"]
pdfMCVarBins = [totalSystHist.GetYaxis().FindFixBin(label) for label in pdfMCVarLabels]
hessianPDFSystDeltasUp = np.zeros(totalSystHist.GetNbinsX()+2)
scaleSystDeltas = np.zeros(totalSystHist.GetNbinsX()+2)
mcPDFSystDeltasUp = [totalSystHist.GetBinError(xBin, totalSystHist.GetYaxis().FindFixBin(pdfMCVarLabels[0])) for xBin in range(0, totalSystHist.GetNbinsX()+2)]
mcPDFSystDeltasDown = [totalSystHist.GetBinError(xBin, totalSystHist.GetYaxis().FindFixBin(pdfMCVarLabels[1])) for xBin in range(0, totalSystHist.GetNbinsX()+2)]

hasHessianPDFVar = False
for xBin in range(0, totalSystHist.GetNbinsX()+2):
    if totalSystHist.GetBinContent(xBin, totalSystHist.GetYaxis().FindFixBin("LHEPdfWeightHessian_NominalComb")) != 0:
        hasHessianPDFVar = True
        break
if hasHessianPDFVar:
    pdfKeyLabels = pdfWeightLabels.copy()
    pdfKeyLabels.remove("LHEPdfWeight_0")
    hessianPDFSystDeltasUp, pdfSystDeltasDown, yBins = CalculatePDFVariationHessian(totalSystHist, sample, pdfKeyLabels, totalSystHist.GetYaxis().FindFixBin("LHEPdfWeightHessian_NominalComb"))
    # hessianPDFSystDeltasUp, pdfSystDeltasDown, yBins = CalculatePDFVariationHessian(totalSystHist, sample, pdfKeyLabels, 1)
if not IsHistEmpty(totalSystHist):
    pdfUpCombBin = totalSystHist.GetYaxis().FindFixBin("LHEPdf_UpComb")
    pdfDownCombBin = totalSystHist.GetYaxis().FindFixBin("LHEPdf_DownComb")
    for xBin in range(0, totalSystHist.GetNbinsX()+2):
        nominal = totalSystHist.GetBinContent(xBin, 1)  # y-bin 1 always being nominal
        pdfSystDeltaUp = math.sqrt(pow(hessianPDFSystDeltasUp[xBin], 2)+pow(mcPDFSystDeltasUp[xBin], 2))
        pdfSystDeltaDown = math.sqrt(pow(hessianPDFSystDeltasUp[xBin], 2)+pow(mcPDFSystDeltasDown[xBin], 2))
        pdfSystTotUp = pdfSystDeltaUp+nominal  # convention of filling with x' = deltaX + x
        pdfSystTotDown = pdfSystDeltaDown+nominal
        totalSystHist.SetBinContent(xBin, pdfUpCombBin, pdfSystTotUp)
        totalSystHist.SetBinContent(xBin, pdfDownCombBin, pdfSystTotDown)
        totalSystHist.SetBinError(xBin, pdfUpCombBin, pdfSystDeltaUp)
        totalSystHist.SetBinError(xBin, pdfDownCombBin, pdfSystDeltaDown)  # for plotting, we rely on the sumQuad addition of bin errors

if corrLHESysts:
    verbose = False
    scaleSystDeltas, yBins, maxYIndices = CalculateShapeSystematic(totalSystHist, sample, systNameToBranchTitleDict, verbose)
    print("DEBUG: scaleSystDeltas: calculate from scale weight variations; deltas[1]={} ".format(scaleSystDeltas[1]))
    scaleUpCombBin = totalSystHist.GetYaxis().FindFixBin("LHEScale_UpComb")
    scaleDownCombBin = totalSystHist.GetYaxis().FindFixBin("LHEScale_DownComb")
    for xBin in range(0, totalSystHist.GetNbinsX()+2):
        nominal = totalSystHist.GetBinContent(xBin, 1)  # y-bin 1 always being nominal
        totalSystHist.SetBinContent(xBin, scaleUpCombBin, scaleSystDeltas[xBin]+nominal)
        totalSystHist.SetBinError(xBin, scaleUpCombBin, scaleSystDeltas[xBin])  # for plotting, we rely on the sumQuad addition of bin errors
        totalSystHist.SetBinContent(xBin, scaleDownCombBin, scaleSystDeltas[xBin]+nominal)
        totalSystHist.SetBinError(xBin, scaleDownCombBin, scaleSystDeltas[xBin])
    lheScaleMaxCombBin = totalSystHist.GetYaxis().FindFixBin("LHEScaleWeight_maxComb")
    for xBin in range(0, totalSystHist.GetNbinsX()+2):
        nominal = totalSystHist.GetBinContent(xBin, 1)  # y-bin 1 always being nominal
        totalSystHist.SetBinContent(xBin, lheScaleMaxCombBin, scaleSystDeltas[xBin]+nominal)
        totalSystHist.SetBinError(xBin, lheScaleMaxCombBin, scaleSystDeltas[xBin])

# need to ignore the special combination systematics bin
# yBinLabelsToIgnore = ["LHEPdfWeightHessian_NominalComb", "LHEPdf_UpComb", "LHEPdf_DownComb", "LHEScaleWeight_maxComb", "LHEScale_UpComb", "LHEScale_DownComb",]
# yBinLabelsToIgnore = ["LHEPdfWeightHessian_NominalComb", "LHEScaleWeight_maxComb"]
yBinLabelsToIgnore = []
combTFile = TFile.Open(combinedFile)
systHistComb = GetSystematicsHisto(combTFile)
xBin = totalSystHist.GetXaxis().FindFixBin("preselection")
for yBin in range(0, systHistComb.GetNbinsY()+2):
    yBinLabel = systHistComb.GetYaxis().GetBinLabel(yBin)
    if any(yBinLabel in label for label in yBinLabelsToIgnore):
        continue
    if systHistComb.GetYaxis().GetBinLabel(yBin) != totalSystHist.GetYaxis().GetBinLabel(yBin):
        raise RuntimeError("Bin label (y) of comb hist '{}' (bin {}) doesn't match that of the summed hist '{}'; comb hist has {} yBins, summed hist has {}".format(systHistComb.GetYaxis().GetBinLabel(yBin), yBin,
            totalSystHist.GetYaxis().GetBinLabel(yBin), systHistComb.GetNbinsY(), totalSystHist.GetNbinsY()))
    if systHistComb.GetXaxis().GetBinLabel(xBin) != totalSystHist.GetXaxis().GetBinLabel(xBin):
        raise RuntimeError("Bin label (x) of comb hist '{}' doesn't match that of the summed hist '{}'".format(systHistComb.GetXaxis().GetBinLabel(xBin),
            totalSystHist.GetXaxis().GetBinLabel(xBin)))
    if not math.isclose(systHistComb.GetBinContent(xBin, yBin), totalSystHist.GetBinContent(xBin, yBin), rel_tol=1e-5):
        raise RuntimeError("Bin content of comb hist {} doesn't match that of the summed hist {} for xBin={} yBin={} with labels xBin={} yBin={}".format(
            systHistComb.GetBinContent(xBin, yBin), totalSystHist.GetBinContent(xBin, yBin), xBin, yBin, systHistComb.GetXaxis().GetBinLabel(xBin),
            systHistComb.GetYaxis().GetBinLabel(yBin)))
    # we can only check the noninal yield for bin errors, as the bin errors are reset in combinePlots
    if yBin == 1:
        if not math.isclose(systHistComb.GetBinError(xBin, yBin), totalSystHist.GetBinError(xBin, yBin), rel_tol=1e-5):
            raise RuntimeError("Bin error of comb hist {} doesn't match that of the summed hist {} for xBin={} yBin={} with labels xBin={} yBin={}; bin contents = {}, {}".format(
                systHistComb.GetBinError(xBin, yBin), totalSystHist.GetBinError(xBin, yBin), xBin, yBin, systHistComb.GetXaxis().GetBinLabel(xBin),
                systHistComb.GetYaxis().GetBinLabel(yBin),
                systHistComb.GetBinContent(xBin, yBin), totalSystHist.GetBinContent(xBin, yBin)
                ))
    if "comb" in yBinLabel.lower() or yBin == 1 or "scale" in yBinLabel.lower():
        print("INFO: Bin error of comb hist {} matches that of the summed hist {} for xBin={} yBin={} with labels xBin={} yBin={}; bin contents = {}, {}".format(
            systHistComb.GetBinError(xBin, yBin), totalSystHist.GetBinError(xBin, yBin), xBin, yBin, systHistComb.GetXaxis().GetBinLabel(xBin),
            systHistComb.GetYaxis().GetBinLabel(yBin),
            systHistComb.GetBinContent(xBin, yBin), totalSystHist.GetBinContent(xBin, yBin)
            ))

combTFile.Close()
for f in allAnaFiles:
    f.Close()
