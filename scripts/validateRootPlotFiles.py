#/usr/bin/env python3

import sys
import math
import ctypes
import numpy as np
from ROOT import TIter, TFile
import ROOT as r

from combineCommon import ComparePDFBranches


def GetSampleHistoNamesFromTFile(tfileName):
    # histNameToHistDict = {}
    if tfileName.startswith("/eos/cms"):
        tfileName = "root://eoscms/" + tfileName
    elif tfileName.startswith("/eos/user"):
        tfileName = "root://eosuser/" + tfileName
    histNames = []
    tfile = TFile.Open(tfileName)
    for key in tfile.GetListOfKeys():
        histNames.append(key.GetName())
    tfile.Close()
    if len(histNames) < 1:
        raise RuntimeError(
                "GetSampleHistoNamesFromTFile({}) -- failed to read any histos from this file!".format(tfile.GetName()))
    return histNames, tfileName


def CheckIfIdenticalTObjs(h1, h2):
    if h1.GetName() != h2.GetName():
        raise RuntimeError("names of objects do not agree: '{}' and '{}'".format(h1.GetName(), h2.GetName()))
    relTol = 1e-3
    digits = 5
    if h1.InheritsFrom("TH1"):
        if h1.GetNbinsX() != h2.GetNbinsX():
            print("WARN: hist 1 {} has {} x bins while hist 2 {} has {} x bins".format(h1.GetName(), h1.GetNbinsX(), h2.GetName(), h2.GetNbinsX()))
        if h1.GetNbinsY() != h2.GetNbinsY():
            print("WARN: hist 1 has {} y bins while hist 2 {} has {} y bins".format(h1.GetName(), h1.GetNbinsY(), h2.GetName(), h2.GetNbinsY()))
        if h1.GetNbinsZ() != h2.GetNbinsZ():
            print("WARN: hist 1 {} has {} z bins while hist 2 {} has {} z bins".format(h1.GetName(), h1.GetNbinsZ(), h2.GetName(), h2.GetNbinsZ()))
        for globalBin in range(0, h1.GetNcells()+1):
            h1BinContent = np.float32(h1.GetBinContent(globalBin))
            h2BinContent = np.float32(h2.GetBinContent(globalBin))
            if h1BinContent != h2BinContent:
                # compare first 6 digits only (float precision from TH1F)
                # h1bc = str(h1BinContent)[:digits+2] if abs(h1BinContent) < 1 else str(h1BinContent)[:digits+1]
                # h2bc = str(h2BinContent)[:digits+2] if abs(h2BinContent) < 1 else str(h2BinContent)[:digits+1]
                # if h1bc != h2bc:
                # if '{:.6}'.format(h1BinContent) != '{:.6}'.format(h2BinContent):
                # try to compare contents with tolerance
                if not math.isclose(h1BinContent, h2BinContent, rel_tol=relTol):
                    systName = ""
                    if "systematics" in h1.GetName().lower():
                        xbin  = ctypes.c_int()
                        ybin  = ctypes.c_int()
                        zbin  = ctypes.c_int()
                        h1.GetBinXYZ(globalBin, xbin, ybin, zbin)
                        if ybin.value > 1:
                            systName = h1.GetYaxis().GetBinLabel(ybin.value) + ","
                    print("WARN: For histo {}, bin {}, {} content 1 '{}' != content 2 '{}' and exceeds tolerance ({}).".format(h1.GetName(), globalBin, systName, h1BinContent, h2BinContent, relTol))
                    # print("WARN: For histo {}, bin {}, {} content 1 '{}' != content 2 '{}' and exceeds tolerance ({} != {}).".format(h1.GetName(), globalBin, systName, h1BinContent, h2BinContent, h1bc, h2bc))
            if "TProfile" in h1.ClassName():
                h1BinEntries = h1.GetBinEntries(globalBin)
                h2BinEntries = h2.GetBinEntries(globalBin)
                if h1BinEntries != h2BinEntries:
                    raise RuntimeError("For profile {}, bin {}, entries 1 '{}' != entries 2 '{}'.".format(h1.GetName(), globalBin, h1BinEntries, h2BinEntries))
                h1BinSumw2 = h1.GetSumw2.At(globalBin)
                h2BinSumw2 = h2.GetSumw2.At(globalBin)
                if h1BinSumw2 != h2BinSumw2:
                    raise RuntimeError("For profile {}, bin {}, sumw2 1 '{}' != sumw2 2 '{}'.".format(h1.GetName(), globalBin, h1BinSumw2, h2BinSumw2))
            else:
                h1BinError = h1.GetBinError(globalBin)
                h2BinError = h2.GetBinError(globalBin)
                if h1BinError != h2BinError:
                    if not math.isclose(h1BinError, h2BinError, rel_tol=relTol):
                    # if '{:.6}'.format(h1BinError) != '{:.6}'.format(h2BinError):
                        systName = ""
                        if "systematics" in h1.GetName().lower():
                            xbin  = ctypes.c_int()
                            ybin  = ctypes.c_int()
                            zbin  = ctypes.c_int()
                            h1.GetBinXYZ(globalBin, xbin, ybin, zbin)
                            if ybin.value > 1:
                                systName = h1.GetYaxis().GetBinLabel(ybin.value) + ","
                        print("WARN: For histo {}, bin {}, {} error 1 '{}' != error 2 '{}' and exceeds tolerance ({}).".format(h1.GetName(), globalBin, systName, h1BinError, h2BinError, relTol))
    elif h1.InheritsFrom("TMap"):
        CheckSystTMapConsistency(h1, h2)
    else:
        print("Found object name '{}', class '{}'; name2: '{}', class2 '{}' [NOT CHECKING!]".format(h1.GetName(), h1.ClassName(), h2.GetName(), h2.ClassName()))


def CheckSystTMapConsistency(combinedSampleMap, mapToCheck):
    combIter = TIter(combinedSampleMap)
    combKey = combIter.Next()
    while combKey:
        sampleMapObject = mapToCheck.FindObject(combKey.GetName())
        if not sampleMapObject:
            sampleItr = TIter(mapToCheck)
            sampleKey = sampleItr.Next()
            while sampleKey:
                # print("sampleKey: '{}'".format(sampleKey.GetName()))
                sampleKey = sampleItr.Next()
            raise RuntimeError("could not find syst '{}' in sampleMap. systematics TMap in combined sample is inconsistent in input root file".format(
                combKey.GetName()))
        sampleMapVal = sampleMapObject.Value()
        combVal = combinedSampleMap.GetValue(combKey)
        # now check the TLists
        combListIter = TIter(combVal)
        combListItem = combListIter.Next()
        while combListItem:
            if "lha id" in combListItem.GetName().lower():
                return ComparePDFBranches([item for item in combVal], [item for item in sampleMapVal])
            sampleListItem = sampleMapVal.FindObject(combListItem.GetName())
            if not sampleListItem:
                raise RuntimeError("branch title used in combined sample '{}' for syst '{}' not found in sample list: {}".format(
                    combListItem.GetName(), combKey.GetName(), [item for item in sampleMapVal]))
            if sampleListItem.GetName() != combListItem.GetName():
                raise RuntimeError("branch title for combined sample '{}' does not equal branch name for candidate sample '{}' for syst '{}'".format(
                    combListItem.GetName(), sampleListItem.GetName(), combKey.GetName()))
            combListItem = combListIter.Next()
        combKey = combIter.Next()


def PrintCompressionInfo(tfileName):
    if tfileName.startswith("/eos/cms"):
        tfileName = "root://eoscms/" + tfileName
    elif tfileName.startswith("/eos/user"):
        tfileName = "root://eosuser/" + tfileName
    tfile = TFile.Open(tfileName)
    print("INFO: TFile {}:\ncompression info: algorithm {}, level {}, settings {}; compression factor {}".format(tfileName, tfile.GetCompressionAlgorithm(),
        tfile.GetCompressionLevel(), tfile.GetCompressionSettings(), tfile.GetCompressionFactor()))
    tfile.Close()


####################################################################################################
# Run
####################################################################################################
# sampleToSelect = "ALLBKG_powhegTTBar_ZJetPtIncStitch_NLODiboson"
# sampleToSelect = "ZJet_amcatnlo_ptBinned_IncStitch"
sampleToSelect = "TTTo2L2Nu"

if len(sys.argv) < 3:
    raise RuntimeError("Need to specify the two root files to compare as command-line arguments")
file1 = sys.argv[1]
file2 = sys.argv[2]

# PrintCompressionInfo(file1)
# PrintCompressionInfo(file2)

# print("INFO: Get TObj names from file1...", end="", flush=True)
sampleHistoNames1, tfileName1 = GetSampleHistoNamesFromTFile(file1)
# print("done.", flush=True)
# print("INFO: Get TObj names from file2...", end="", flush=True)
sampleHistoNames2, tfileName2 = GetSampleHistoNamesFromTFile(file2)
# print("done.", flush=True)

nObjs1 = len(sampleHistoNames1)
nObjs2 = len(sampleHistoNames2)

if nObjs1 != nObjs2:
    sample1Set = set(sampleHistoNames1)
    sample2Set = set(sampleHistoNames2)
    print("ERROR: {} TObjs found in file1 but {} in file2".format(nObjs1, nObjs2))
    print("Objects in file1 but not file2:", sample1Set - sample2Set)
    print("Objects in file2 but not file1:", sample2Set - sample1Set)
    # raise RuntimeError("Number of objects in files differ")
    print("Proceeding anyway to compare objects")
    
print("INFO: Found {} TObjs from file1 and TObjs {} in file2".format(nObjs1, nObjs2), flush=True)
tfile1 = TFile.Open(tfileName1)
tfile2 = TFile.Open(tfileName2)
# compare only the shorter list of hist names in case of hists missing from one file or another
sampleHistNamesToCompare = sampleHistoNames1 if len(sampleHistoNames1) <= len(sampleHistoNames2) else sampleHistoNames2
if len(sampleToSelect):
    print("INFO: comparing histograms for sample '{}' only.".format(sampleToSelect))
    sampleHistNamesToCompare = [histName for histName in sampleHistNamesToCompare if sampleToSelect in histName]
nHistsToCompare = len(sampleHistNamesToCompare)
# split into chunks
n = 1000
sampleHistNamesToCompare = [sampleHistNamesToCompare[i * n:(i + 1) * n] for i in range((len(sampleHistNamesToCompare) + n - 1) // n )]
print("INFO: Comparing {} objects in {} chunks of {}...".format(nHistsToCompare, len(sampleHistNamesToCompare), n))
for idx, chunk in enumerate(sampleHistNamesToCompare):
    print("INFO: Comparing chunk {}/{}".format(idx+1, len(sampleHistNamesToCompare)))
    histList1 = []
    histList2 = []
    # batch read chunks of histos and compare them
    for histName in chunk:
        tobj = tfile1.Get(histName)
        tobj2 = tfile2.Get(histName)
        r.SetOwnership(tobj, True)
        r.SetOwnership(tobj2, True)
        if tobj.InheritsFrom("TH1"):
            tobj.SetDirectory(0)
            tobj2.SetDirectory(0)
        histList1.append(tobj)
        histList2.append(tobj2)
    for idx, tobj in enumerate(histList1):
        CheckIfIdenticalTObjs(tobj, histList2[idx])
tfile1.Close()
tfile2.Close()
