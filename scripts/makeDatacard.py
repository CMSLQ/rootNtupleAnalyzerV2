#!/usr/bin/env python3

import os
import sys
import math
import string
import re
from prettytable import PrettyTable
from tabulate import tabulate
from pathlib import Path
from decimal import Decimal
from termcolor import colored
from bisect import bisect
import numpy as np
import ROOT as r
import matplotlib.pyplot as plt
import mplhep as hep
import combineCommon as cc

sys.path.append(os.getenv("LQMACRO").rstrip("/") + "/plotting2016/")
from plotSignalEfficiencyTimesAcceptance import plotSignalEfficiencyTimesAcceptance

r.gROOT.SetBatch()


def EntryIsValid(entry):
    if entry == "-":
        return False
    elif entry == "no effect":
        return False
    return True


def GetFullSignalName(signal_name, mass_point):
    verbose = False
    fullSignalName = signal_name.replace("[masspoint]", mass_point)
    signalNameForFile = ""
    if verbose:
        print("GetFullSignalName(): signal_name=", signal_name, "fullSignalName=", fullSignalName)
    if "BetaHalf" in signal_name:
        signalNameForFile = "LQToUE_ENuJJFilter_M-" + mass_point + "_BetaHalf"
    elif "LQ" in signal_name:
        signalNameForFile = signalNameTemplate.format(mass_point)
        fullSignalName = "LQ_M"+str(mass_point)  # RunStatsBasicCLS requires this signal name
    elif "Stop" in signal_name:
        ctau = signal_name[signal_name.find("CTau") + 4:]
        # print 'found ctau=',ctau,'in signal_name:',signal_name,'mass point:',mass_point
        signalNameForFile = "DisplacedSUSY_StopToBL_M-" + mass_point + "_CTau-" + ctau
    else:
        print("WARN: GetFullSignalName(): not sure how to define signalNameForFile for signal_name={}; returning empty string.".format(signal_name))
    return fullSignalName, signalNameForFile


def GetStatErrors(nevts, nUnweighted, theta=1.0):
    # use parameters as in datacard gamma distribution
    # see: https://cds.cern.ch/record/1379837/files/NOTE2011_005.pdf (eq. 18)
    verbose = False
    alpha = 1.0 - 0.682689492
    lower = (
        0 if nUnweighted == 0
        else nevts - r.Math.gamma_quantile(alpha/2, nUnweighted, theta)
    )
    upper = r.Math.gamma_quantile_c(alpha/2, nUnweighted+1, theta) - nevts
    if verbose:
        print("calculate upper gamma quantile_c for nevts={}, nUnweighted={}, theta={}; upper error={}".format(
                nevts, nUnweighted, theta, upper))
    return upper, lower


def GetStatErrorFromDict(statErrDict, mass):
    availMasses = sorted(statErrDict.keys())
    if mass not in availMasses and mass > availMasses[-1]:
        mass = availMasses[-1]
    return statErrDict[mass]


def GetStatErrorsFromDatacard(dictEntry, nevts):
    pdfType = dictEntry[0]
    if pdfType == "lnN" and len(dictEntry) == 2:
        return dictEntry[1], dictEntry[1]
    elif pdfType == "gmN" and len(dictEntry) == 3:
        # print "GetStatErrors({}, {}, {})".format(nevts, dictEntry[1], dictEntry[2])
        return GetStatErrors(nevts, dictEntry[1], dictEntry[2])
    else:
        raise RuntimeError("Could not GetStatErrorsFromDatacard: didn't understand dictEntry={}".format(dictEntry))


# returns dict like systDict[systName][selection] = yield
def GetSystematicsDict(rootFile, sampleName, selections, verbose=False):
    systHistName = "histo2D__{}__systematics".format(sampleName)
    systHist = rootFile.Get(systHistName)
    systDict = {}
    xBinPresel = systHist.GetXaxis().FindFixBin("preselection")
    for selection in selections:
        # expect that selections are either "preselection" or "LQXXXX"
        finalSelection = cc.GetFinalSelection(selection, doEEJJ)
        systDict[selection] = {}
        xBin = systHist.GetXaxis().FindFixBin(finalSelection)
        if xBin < 1:
            raise RuntimeError("Could not find requested selection name '{}' in hist {} in file {}".format(finalSelection, systHistName, rootFile.GetName()))
        lheScaleYBin = None
        for yBin in range(1, systHist.GetNbinsY()+1):
            systName = systHist.GetYaxis().GetBinLabel(yBin)
            preselYield = None
            if "LHEPdfWeight" in systName or "LHEScaleWeight" in systName:
                continue
            elif "LHEPdf" in systName or "LHEScale" in systName:
                systName = systName.replace("_UpComb", "CombUp").replace("_DownComb", "CombDown")
                if "LHEScale" in systName:
                    preselYieldBin = systHist.GetYaxis().FindFixBin("LHEScaleWeight_preselYield")
                    preselYield = systHist.GetBinContent(xBin, preselYieldBin)
                    # print("DEBUG: for systName={}, selection={}, sampleName={}, with xBinPresel={}, preselYieldBin={}; preselYield={}".format(
                    #     systName, selection, sampleName, xBinPresel, preselYieldBin, preselYield))
            systDict[selection][systName] = {}
            systDict[selection][systName]["yield"] = systHist.GetBinContent(xBin, yBin)
            systDict[selection][systName]["preselYield"] = preselYield
            # print("DEBUG: ASSIGNED for systName={}, selection={}, preselYield={}".format(systName, selection, preselYield))
            # if "ZJet" in sampleName and "LQ3000" in selection:
            #     print("INFO: for sampleName={} and selection={}, xBin={} yBin={} ({}), got content={}".format(sampleName, selection, xBin, yBin, systName, systDict[selection][systName]))
    # add special entry for branch titles
    systDict["branchTitles"] = {}
    tmapName = "tmap__{}__systematicNameToBranchesMap".format(sampleName)
    tmap = rootFile.Get(tmapName)
    if not tmap or tmap is None:
        raise RuntimeError("Could not find TMap '{}' in file {}".format(tmapName, rootFile.GetName()))
    for yBin in range(1, systHist.GetNbinsY()+1):
        branchTitleList = []
        systName = systHist.GetYaxis().GetBinLabel(yBin)
        # special handling for LHEPdf/LHEScalesysts
        if "LHEPdfWeight" in systName or "LHEScaleWeight" in systName:
            continue
        elif "LHEPdf" in systName or "LHEScale" in systName:
            systName = systName.replace("_UpComb", "CombUp").replace("_DownComb", "CombDown")
            systDict["branchTitles"][systName] = [systName]
            continue
        mapObject = tmap.FindObject(systName)
        if not mapObject:
            # print("\tINFO: look again for matching TObject for syst {}".format(systName[:systName.rfind("_")]), flush=True)
            # assume it's an array syst, so try to match stripping off the _N part
            mapObject = tmap.FindObject(systName[:systName.rfind("_")])
        if not mapObject:
            tmap.Print()
            raise RuntimeError("Could not find matching TObject in map for syst {} or {}; see map content above".format(systName, systName[:systName.rfind("_")]))
        branchTitleListItr = r.TIter(mapObject.Value())
        branchTitle = branchTitleListItr.Next()
        while branchTitle:
            branchTitleList.append(branchTitle.GetName())
            branchTitle = branchTitleListItr.Next()
        systDict["branchTitles"][systName] = branchTitleList
        # print "INFO: systDict[\"branchTitles\"][{}] = {}".format(systName, branchTitleList)
    #if verbose:
    #    print("sampleName={}: systDict={}".format(sampleName, systDict["LQ300"]))
    # print("DEBUG: sampleName={}: systDict={}".format(sampleName, systDict))
    # reindex by syst name
    systDict = {syst: {sel: systDict[sel][syst] for sel in systDict} for syst in list(systDict[selections[0]].keys())}
    # print("sampleName={}: systDict=".format(sampleName), systDict)
    # if verbose:
    #     print "sampleName={}: systDict=".format(sampleName)
    #     for systName, selYieldDicts in systDict.items():
    #         print "{} : {}".format(systName, selYieldDicts["LQ300"]),
    #     print
    return systDict


def GetSystYield(deltaOverNom, systNomYield):
    delta = deltaOverNom * systNomYield
    systYield = delta + systNomYield
    return systYield


def HasPreselectionSystYield(systName, systDict, selection, up):
    if "shape" in systName.lower():
        systName = "LHEScaleComb"
    elif "LHEPdfWeight" == systName:
        systName = "LHEPdfComb"
    systName += "Up" if up else "Down"
    preselSystYield = None
    if systName in systDict.keys():
        preselSystYield = systDict[systName][selection]["preselYield"]
    if preselSystYield is not None:
        return True, preselSystYield
    else:
        return False, preselSystYield
    # elif "lumi" in systName.lower():
    #     return False, None
    # if "preselYield" in systDict[systName].keys():
    #     preselSystYield = systDict[systName]["preselYield"]
    #     return True, preselSystYield
    # return False, None


def GetSignalSystDeltaOverNominal(signalName, selectionName, verbose=False):
    return GetTotalSystDeltaOverNominal(signalName, selectionName, systematicsNamesSignal, d_signal_systs, verbose)


def GetBackgroundSystDeltaOverNominal(background_name, selectionName, verbose=False):
    return GetTotalSystDeltaOverNominal(background_name, selectionName, systematicsNamesBackground, d_background_systs, verbose)


def GetTotalSystDeltaOverNominal(sampleName, selectionName, systematicsNames, d_systs, verbose=False):
    totalSyst = 0
    for syst in systematicsNames:
        systEntry, deltaOverNominalUp, deltaOverNominalDown, _, _ = GetSystematicEffectAbs(syst, sampleName, selectionName, d_systs)
        deltaOverNominalMax = max(deltaOverNominalUp, deltaOverNominalDown)
        if deltaOverNominalMax > 0:
            totalSyst += deltaOverNominalMax * deltaOverNominalMax
        elif deltaOverNominalMax > 1:
            raise RuntimeError(
                    "deltaOverNominalMax > 1 for sampleName={} syst={} selection={} deltaOverNominalUp={} deltaOverNominalDown={}".format(
                        sampleName, syst, selectionName, deltaOverNominalUp, deltaOverNominalDown))
        # if verbose:
        #     print "GetSystematicEffect({}, {}, {}, {})".format(syst, sampleName, selectionName, d_systs.keys())
        #     print "\t result={}".format(GetSystematicEffect(syst, sampleName, selectionName, d_systs))
        #     print "\t totalSyst={}".format(totalSyst)
    if verbose:
        print("GetTotalSystDeltaOverNominal(): {} -- return sqrt(totalSyst) = sqrt({}) = {}".format(
                sampleName, totalSyst, math.sqrt(totalSyst)))
    return math.sqrt(totalSyst)


# return abs val of syst; also checks for deltaOverNom==1 and does renorm if needed
def GetSystematicEffectAbs(systName, sampleName, selection, fullSystDict, verbose=False):
    # verbose = True
    entry, deltaOverNomUp, deltaOverNomDown, symmetric, systNominal, systSelection = GetSystematicEffect(systName, sampleName, selection, fullSystDict)
    if verbose or "pdf" in systName.lower():
        print("For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}, systNominal={}, systSelection={}".format(
            sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown, systNominal, systSelection))
    if entry == "-":
        # this means the systematic doesn't apply
        return "-", 0, 0, systNominal, systSelection
    nomYield = systNominal
    systNomSelection = systSelection
    # returned deltaOverNom is always systYield - nominal
    # if "norm" not in systName.lower() and "shape" not in systName.lower() and "lumi" not in systName.lower():
    if "norm" not in systName.lower() and selection != "preselection" and selection != "trainingSelection":
        if sampleName in backgroundsToRenormSystAtPresel:
            if verbose:
                print("\tINFO: renormalizing syst={} for background={}".format(systName, sampleName))
                print("\t original entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(entry, deltaOverNomUp, deltaOverNomDown))
            preselNomYield, preselNomSelection, preselRawEvents = GetSystNominalYield(systName, fullSystDict[sampleName], "preselection", sampleName)
            preselEntry, preselDOverNUp, preselDOverNDown, preselSymm, preselSystNominal, preselSystSelection = GetSystematicEffect(systName, sampleName, "preselection", fullSystDict)
            if preselSystSelection != preselNomSelection:
                raise RuntimeError("Something strange happened: the selection used for the preselection systematic '{}' was not the same as the presel syst nominal yield selection '{}'.".format(
                    preselSystSelection, preselNomSelection))
            if preselSystNominal != preselNomYield:
                raise RuntimeError("Something strange happened: the selection used for the preselection systematic yield '{}' was not the same as the presel syst nominal yield '{}'.".format(
                    preselSystNominal, preselNomYield))
            preselSystYieldUp = GetSystYield(preselDOverNUp, preselNomYield)
            hasPreselSystYield, preselSystYield = HasPreselectionSystYield(systName, fullSystDict[sampleName], selection, True)
            if hasPreselSystYield:
                if verbose:
                    print("\tOverriding preselSystYieldUp for sampleName={} systName={} to {}".format(sampleName, systName, preselSystYield))
                preselSystYieldUp = preselSystYield
            nomYield, systNomSelection, systRawEvents = GetSystNominalYield(systName, fullSystDict[sampleName], selection, sampleName)
            systYieldUp = GetSystYield(deltaOverNomUp, nomYield)
            systYieldUpRenorm = systYieldUp * preselNomYield/preselSystYieldUp
            if verbose:
                print("\trenorm systYieldUp to {} = {}*{}/{} = systYieldUp*preselNomYield/preselSystYieldUp".format(
                        systYieldUpRenorm, systYieldUp, preselNomYield, preselSystYieldUp))
            deltaOverNomUp = (systYieldUpRenorm - nomYield) / nomYield
            if verbose:
                print("\t nominal yield at selection={}; deltaX/X={}".format(nomYield, deltaOverNomUp))
            entry = str(1 + math.fabs(deltaOverNomUp))
            if not symmetric:
                deltaDown = deltaOverNomDown * nomYield
                systYieldDown = deltaDown + nomYield
                preselSystYieldDown = GetSystYield(preselDOverNDown, preselNomYield)
                hasPreselSystYield, preselSystYield = HasPreselectionSystYield(systName, fullSystDict[sampleName], selection, False)
                if hasPreselSystYield:
                    if verbose:
                        print("\tOverriding preselSystYieldDown up for sampleName={} systName={} to {}".format(sampleName, systName, preselSystYield))
                    preselSystYieldDown = preselSystYield
                if verbose:
                    print("\trenorm systYieldDown to {}*{}/{} = systYieldDown*preselNomYield/preselSystYieldDown".format(
                            systYieldDown, preselNomYield, preselSystYieldDown))
                systYieldDownRenorm = systYieldDown * preselNomYield/preselSystYieldDown
                deltaOverNomDown = (systYieldDownRenorm - nomYield) / nomYield
                if verbose:
                    print("\trenorm systYieldDown to {} = {}*{}/{} = systYieldDown*preselNomYield/preselSystYieldDown".format(
                            systYieldDownRenorm, systYieldDown, preselNomYield, preselSystYieldDown))
                kUp = systYieldUpRenorm/nomYield
                kDown = systYieldDownRenorm/nomYield
                entry = str(kDown/kUp)
                if verbose:
                    print("\t new kUp={} = {} / {} = systYieldDownRenorm/nomYield, new kDown={} = {} / {}".format(
                            kUp, systYieldUpRenorm, nomYield, kDown, systYieldDownRenorm, nomYield))
            else:
                deltaOverNomDown = deltaOverNomUp
            if verbose:
                print("\t new entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(entry, deltaOverNomUp, deltaOverNomDown))
    if deltaOverNomUp == 0 and deltaOverNomDown == 0:
        # print("WARN: For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown))
        return "no effect", 0, 0, nomYield, systNomSelection
    return entry, math.fabs(deltaOverNomUp), math.fabs(deltaOverNomDown), nomYield, systNomSelection


def DoesSystematicApply(systName, sampleName):
    return systName in d_applicableSystematics[sampleName]


def GetSystematicEffect(systName, sampleName, selection, fullSystDict):
    # print("GetSystematicEffect({}, {}. {})".format(systName, sampleName, selection))
    verbose = False
    systDict = fullSystDict[sampleName]
    try:
        nominal = systDict["nominal"][selection]["yield"]
    except KeyError as e:
        raise RuntimeError("Could not find key 'nominal' in systDict for sampleName={}; systDict.keys()={}".format(sampleName, systDict.keys()))
    except Exception as e:
        raise RuntimeError("Got exception accessing systDict[{}][{}][{}] for sampleName={}; systDict[{}]={}".format("nominal", selection, "yield", sampleName, "nominal", systDict["nominal"]))
    if not DoesSystematicApply(systName, sampleName):
        return "-", 0, 0, True, nominal, selection
    if "shape" in systName.lower():
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = CalculateUpDownSystematic("LHEScaleComb", systDict, selection, sampleName, verbose)
    elif systName == "LHEPdfWeight":
        # print("INFO GetSystematicEffect(): for sample={}, selection={}, systName={}, fullSystDict[{}][LHEPdfCombUp] = {}".format(sampleName, selection, systName, sampleName, fullSystDict[sampleName]["LHEPdfCombUp"][selection]))
        # print("INFO GetSystematicEffect(): for sample={}, selection={}, systName={}, fullSystDict[{}][LHEPdfCombDown] = {}".format(sampleName, selection, systName, sampleName, fullSystDict[sampleName]["LHEPdfCombDown"][selection]))
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = CalculateUpDownSystematic("LHEPdfComb", systDict, selection, sampleName, verbose)
        # if verbose:
        #     components = dictSamples[sampleName]["pieces"]
        #     for sampleName in components:
        #         if sampleName not in dictSamples.keys():
        #             print("INFO: Cannot print PDF systematic output for sample {} as it was not saved.".format(sampleName))
        #         else:
        #             CalculateUpDownSystematic("LHEPdfComb", fullSystDict[sampleName], selection, sampleName, True)  # this call is broken somehow
    elif systName == "Lumi"+str(yearInt):
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = str(1+lumiDeltaXOverX), lumiDeltaXOverX, lumiDeltaXOverX, True, nominal, selection
    elif systName == "LumiCorrelated":
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = str(1+lumiCorrelatedDeltaXOverX), lumiCorrelatedDeltaXOverX, lumiCorrelatedDeltaXOverX, True, nominal, selection
    elif systName == "LumiCorrelated1718":
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = str(1+lumi1718CorrelatedDeltaXOverX), lumi1718CorrelatedDeltaXOverX, lumi1718CorrelatedDeltaXOverX, True, nominal, selection
    elif "norm" in systName.lower():
        if "tt" in systName.lower() or "dy" in systName.lower() or "qcd" in systName.lower():
            entry, deltaNomUp, deltaNomDown, symmetric = CalculateFlatSystematic(systName, systDict, selection)
            newNominal = nominal if nominal > 0 else 1
            newSelection = selection
    else:
        entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = CalculateUpDownSystematic(systName, systDict, selection, sampleName)
    return entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection


def CalculateShiftSystematic(systName, systDict, selection, sampleName):
    nominal, selection, rawEvents = GetSystNominalYield(systName, systDict, selection, sampleName)
    try:
        systYield = systDict[systName][selection]["yield"]
    except KeyError as e:
        if systName not in systDict.keys():
            raise RuntimeError("Could not find syst named {} in systDict keys {}".format(systName, systDict.keys()))
        elif selection not in systDict[systName].keys():
            raise RuntimeError("Could not find selection {} in systDict[{}].keys() {}".format(selection, systName, systDict[systName].keys()))
        else:
            raise RuntimeError("Strange KeyError occurred")
    delta = systYield - nominal
    # print "CalculateShiftSystematic({}, {}, {}): nominal={}, systYield={}, delta={}, delta/nominal={}".format(
    #         systName, selection, sampleName, nominal, systYield, delta, delta/nominal)
    return str(1 + math.fabs(delta)/nominal), delta/nominal, delta/nominal, True


def CalculateUpDownSystematic(systName, systDict, selection, sampleName, verbose=False):
    symmetric = False
    # nominal = systDict["nominal"][selection]
    # print("CalculateUpDownSystematic({}, {}) for sample: {}; nominal yield={}".format(systName, selection, sampleName, nominal))
    nominal, selection, rawEvents = GetSystNominalYield(systName, systDict, selection, sampleName)
    # print("CalculateUpDownSystematic(): nominal={}, selection={} after GetNominalYield({})".format(nominal, selection, systName))
    # verbose = True
    try:
        systYieldUp = systDict[systName+"Up"][selection]["yield"]
        systYieldDown = systDict[systName+"Down"][selection]["yield"]
    except KeyError:
        raise RuntimeError("Could not find Up or Down key for syst={}; keys={}".format(systName, list(systDict.keys())))
    logString = "CalculateUpDownSystematic(): selection={}, sample={} syst={}, systYieldUp={}, systYieldDown={}, systNominal={}, origSystNominal={}".format(
                selection, sampleName, systName, systYieldUp, systYieldDown, nominal, systDict["nominal"][selection])
    kUp = systYieldUp/nominal
    kDown = systYieldDown/nominal
    upDelta = systYieldUp-nominal
    downDelta = systYieldDown-nominal
    if kUp <= 0:
        if verbose:
            print(logString)
            print("CalculateUpDownSystematic(): kUp <= 0 = {}; return symmetric kDown syst '{}', {}, {}".format(kUp, str(1 + math.fabs(downDelta)/nominal), downDelta/nominal, downDelta/nominal))
        return str(1 + math.fabs(downDelta)/nominal), downDelta/nominal, downDelta/nominal, True, nominal, selection
    elif kDown <= 0:
        if verbose:
            print(logString)
            print("CalculateUpDownSystematic(): kDown <= 0 = {}; return symmetric kDown syst '{}', {}, {}".format(kDown, str(1 + math.fabs(upDelta)/nominal), upDelta/nominal, upDelta/nominal))
        return str(1 + math.fabs(upDelta)/nominal), upDelta/nominal, upDelta/nominal, True, nominal, selection
    try:
        entry = kDown/kUp
        deltaNomUp = upDelta/nominal
        deltaNomDown = downDelta/nominal
    except ZeroDivisionError:
        raise RuntimeError("CalculateUpDownSystematic(): caught ZeroDivisionError for a quantity. kDown/kUp = {}/{}, upDelta/nominal={}/{}, downDelta/nominal={}/{}.".format(kDown, kUp, upDelta, nominal, downDelta, nominal)+ " "+logString)
    if deltaNomUp == deltaNomDown:
        entry = 1 + math.fabs(upDelta)/nominal
        symmetric = True
    if verbose:
        print(logString)
        print("CalculateUpDownSystematic(): kDown/kUp = {}/{}; return '{}', {}, {}, {}".format(kDown, kUp, str(entry), deltaNomUp, deltaNomDown, symmetric))
    return str(entry), deltaNomUp, deltaNomDown, symmetric, nominal, selection


def CalculateFlatSystematic(systName, systDict, selection):
    # assumes that the number here is < 1
    # print "d_applicableSystematics["+sampleName+"]=", d_applicableSystematics[sampleName]
    # print "for sample={}, systName={}, systDict.keys()={}".format(sampleName, systName, systDict.keys())
    return str(1 + systDict[systName][selection]["yield"]), systDict[systName][selection]["yield"], systDict[systName][selection]["yield"], True


def GetSystNominalYield(systName, systDict, selection, sampleName):
    verbose = False
    minRawEvents = 2
    requestedSelection = selection
    try:
        nominal = systDict["nominal"][selection]["yield"]
    except KeyError:
        raise RuntimeError("Could not find nominal key for systName={} sampleName={} selection={}; keys={}".format(
            systName, sampleName, selection, list(systDict.keys()))
            )
    if sampleName in list(d_background_rates.keys()):
        unscaledRatesDict =  d_background_unscaledRates[sampleName]
    elif sampleName in list(d_signal_rates.keys()):
        unscaledRatesDict =  d_signal_unscaledRates[sampleName]
    else:
        raise RuntimeError("Could not find sampleName={} in background keys={} or signal keys={}".format(
            sampleName, list(d_background_rates.keys()), list(d_signal_rates.keys())))
    try:
        rawEvents = unscaledRatesDict[selection]
    except KeyError:
        raise RuntimeError("Could not find nominal key for sampleName={} selection={}; keys={}".format(
            sampleName, selection, list(systDict.keys()))
            )
    if nominal <= 0 or rawEvents < minRawEvents:
        # take the last selection for which we have some rate, and use that for systematic evaluation
        lastNonzeroSelection, rate, err, events = GetNearestPositiveSelectionYields(sampleName, selection)
        nominal = rate
        selection = lastNonzeroSelection
        rawEvents = events
    #if sampleName in list(d_background_unscaledRates.keys()):
    #    unscaledRate = d_background_unscaledRates[sampleName][selection]
    #    err = d_background_rateErrs[sampleName][selection]
    #else:
    #    unscaledRate = d_signal_unscaledRates[sampleName][selection]
    #    err = d_signal_rateErrs[sampleName][selection]
    #    if verbose:
    #        print("GetSystNominalYield({}, {}, {}), lastSelection: {}, raw events: {}, rate: {} +/- {}".format(
    #            systName, requestedSelection, sampleName, selection, unscaledRate, nominal, err))
    return nominal, selection, rawEvents


def RoundToN(x, n):
    # if n < 1:
    #    raise ValueError("can't round to less than 1 sig digit!")
    # # number of digits given by n
    # return "%.*e" % (n-1, x)
    if isinstance(x, float):
        return round(x, n)
    else:
        return x


def GetSampleNameFromSubstring(substring, background_names):
    sampleName = [sample for sample in background_names if substring in sample]
    if len(sampleName) != 1:
        raise RuntimeError("Could not find unique {} sample name in background_names with the given keys: {}".format(substring, background_names))
    return sampleName[0]


def GetTableEntryStr(evts, errStatUp="-", errStatDown="-", errSyst=0, latex=False):
    if evts == "-":
        return evts
    # rounding
    evtsR = RoundToN(evts, 2)
    errStatUpR = RoundToN(errStatUp, 2)
    errStatDownR = RoundToN(errStatDown, 2)
    # add additional decimal place if it's zero after rounding
    if evtsR == 0.0:
        evtsR = RoundToN(evts, 3)
    if errStatUpR == 0.0:
        errStatUpR = RoundToN(errStatUp, 3)
    if errStatDownR == 0.0:
        errStatDownR = RoundToN(errStatDown, 3)
    # try again
    if evtsR == 0.0:
        evtsR = RoundToN(evts, 4)
    if errStatUpR == 0.0:
        errStatUpR = RoundToN(errStatUp, 4)
    if errStatDownR == 0.0:
        errStatDownR = RoundToN(errStatDown, 4)
    # handle cases where we don't specify stat or syst
    if errStatUp == "-":
        return str(evtsR)
    elif errSyst == 0:
        if errStatUp == errStatDown:
            if not latex:
                return str(evtsR) + " +/- " + str(errStatUpR)
            else:
                return str(evtsR) + " \\pm " + str(errStatUpR)
        else:
            if not latex:
                return str(evtsR) + " + " + str(errStatUpR) + " - " + str(errStatDownR)
            else:
                return (
                    str(evtsR)
                    + "^{+"
                    + str(errStatUpR)
                    + "}_{-"
                    + str(errStatDownR)
                    + "}"
                )
    else:
        errSystR = RoundToN(errSyst, 2)
        if errStatUp == errStatDown:
            if not latex:
                return str(evtsR) + " +/- " + str(errStatUpR) + " +/- " + str(errSystR)
            else:
                return (
                    str(evtsR) + " \\pm " + str(errStatUpR) + " \\pm " + str(errSystR)
                )
        else:
            return (
                str(evtsR)
                + "^{+"
                + str(errStatUpR)
                + "}_{-"
                + str(errStatDownR)
                + "} \\pm "
                + str(errSystR)
            )


def GetLatexHeaderFromColumnNames(columnNames):
    headers = []
    for col in columnNames:
        if "mlq" in col.lower():
            headers.append("$M_{LQ}$")
        elif "ttbar" in col.lower():
            headers.append(r"$t\bar{t}$")
        else:
            headers.append(col)
    headerLine = " & ".join(headers) + r" \\"
    return headerLine


def GetLastNonzeroSelectionYieldsFromDicts(rateDict, rateErrDict, unscaledRateDict):
    massPointsRev = list(reversed(mass_points))
    idx = 0
    rate = 0
    err = 0
    while rate == 0 or err == 0:
        lastSelectionName = "LQ" + massPointsRev[idx]
        rate = rateDict[lastSelectionName]
        err = rateErrDict[lastSelectionName]
        rawEvents = unscaledRateDict[lastSelectionName]
        idx += 1
        if idx >= len(massPointsRev):
            print("WARN GetLastNonzeroSelectionYieldsFromDicts: could not find final selection point with nonzero yield for this background; returning trainingSelection")
            return "trainingSelection", rateDict["trainingSelection"], rateErrDict["trainingSelection"], unscaledRateDict["trainingSelection"]
    return lastSelectionName, rate, err, rawEvents


def GetLastNonzeroSelectionYields(sampleName):
    if sampleName in list(d_background_rates.keys()):
        return GetLastNonzeroSelectionYieldsFromDicts(
                d_background_rates[sampleName], d_background_rateErrs[sampleName], d_background_unscaledRates[sampleName])
    if sampleName in list(d_signal_rates.keys()):
        return GetLastNonzeroSelectionYieldsFromDicts(
                d_signal_rates[sampleName], d_signal_rateErrs[sampleName], d_signal_unscaledRates[sampleName])
    raise RuntimeError("Could not find sampleName={} in background keys={} or signal keys={}".format(
        sampleName, list(d_background_rates.keys()), list(d_signal_rates.keys())))


def GetNearestPositiveSelectionYieldsFromDicts(sampleName, rateDict, rateErrDict, unscaledRateDict, selection):
    massPointsRev = list(reversed(mass_points))
    # idx = massPointsRev.index(selection.replace("LQ", ""))
    idxInc = mass_points.index(selection.replace("LQ", ""))
    idxInc += 1
    idxDec = massPointsRev.index(selection.replace("LQ", ""))
    idxDec += 1
    rate = 0
    err = 0
    rawEvents = 0
    minRawEvents = 2
    selectionsToCheck = []
    # for index in range(idx, len(mass_points)):
    #     selectionsToCheck.append("LQ" + mass_points[index])
    #     selectionsToCheck.append("LQ" + mass_points[index])
    trimmedMassesInc = mass_points[idxInc:]
    trimmedMassesDec = massPointsRev[idxDec:]
    # selectionsToCheck = list(set(["LQ"+val for pair in zip(massPointsRev, mass_points) for val in pair]))
    # unscaledRate = 0
    # look around given selection
    lastSelectionsChecked = []
    index = 0
    while rate <= 0 or err <= 0 or rawEvents < minRawEvents:
        # lastSelectionName = "LQ" + massPointsRev[idx]
        # lastSelectionName = selectionsToCheck[index]
        if index < len(trimmedMassesInc):
            lastSelectionName = "LQ"+trimmedMassesInc[index]
            rate = rateDict[lastSelectionName]
            err = rateErrDict[lastSelectionName]
            rawEvents = unscaledRateDict[lastSelectionName]
            lastSelectionsChecked.append(lastSelectionName)
            # print("INFO: GetNearestPositiveSelectionYieldsFromDicts: for sample {}, with initial selection={}, check selectionName={}: rate = {} +/- {}".format(sampleName, selection, lastSelectionName, rate, err))
        if rate <= 0 or err <= 0 or rawEvents < minRawEvents:
            if index < len(trimmedMassesDec):
                lastSelectionName = "LQ"+trimmedMassesDec[index]
                rate = rateDict[lastSelectionName]
                err = rateErrDict[lastSelectionName]
                rawEvents = unscaledRateDict[lastSelectionName]
                lastSelectionsChecked.append(lastSelectionName)
                # print("INFO: GetNearestPositiveSelectionYieldsFromDicts: for sample {}, with initial selection={}, check selectionName={}: rate = {} +/- {}".format(sampleName, selection, lastSelectionName, rate, err))
        index += 1
    if len(lastSelectionsChecked) <= 0:
        raise RuntimeError("Could not find nonzero selection for sample={}; rates look like {} and errors look like {}; checked selections={}".format(
            sampleName, rateDict, rateErrDict, lastSelectionsChecked))
    print("INFO: GetNearestPositiveSelectionYieldsFromDicts: for sample {}, with zero nominal rate for selection={}, found lastSelectionName={} with rate = {} +/- {} [{} evts]".format(sampleName, selection, lastSelectionName, rate, err, rawEvents))
    # print "INFO: GetNearestPositiveSelectionYieldsFromDicts: found last selectionName={} with {} unscaled events".format(selectionName, unscaledRate)
    return lastSelectionName, rate, err, rawEvents


def GetNearestPositiveSelectionYields(sampleName, selection):
    # print "INFO: GetNearestPositiveSelectionYields({}, {})".format(sampleName, selection)
    if sampleName in list(d_background_rates.keys()):
        return GetNearestPositiveSelectionYieldsFromDicts(sampleName,
                d_background_rates[sampleName], d_background_rateErrs[sampleName], d_background_unscaledRates[sampleName], selection)
    elif sampleName in list(d_signal_rates.keys()):
        return GetNearestPositiveSelectionYieldsFromDicts(sampleName,
                d_signal_rates[sampleName], d_signal_rateErrs[sampleName], d_signal_unscaledRates[sampleName], selection)
    else:
        raise RuntimeError("Could not find sampleName={} in background keys={} or signal keys={}".format(
            sampleName, list(d_background_rates.keys()), list(d_signal_rates.keys())))


def CheckHistBins(hist):
    for iBin in range(0, hist.GetNbinsX()+2):
        if hist.GetBinContent(iBin) <= 0:
            hist.SetBinContent(iBin, 1e-10)
    return hist


def CreateAndWriteHist(rootFile, mass, sample, content, err, name="yieldFinalSelection", titleStart="Yield at LQ"):
    newHist = r.TH1D("{}_{}_{}".format(name, mass, sample), "{} {} final selection for {}".format(titleStart, mass, sample), 1, 0, 1)
    newHist.SetBinContent(1, content)
    newHist.SetBinError(1, err)
    if sample == "DATA":
        newHist.SetName(newHist.GetName().replace("DATA", "data_obs"))
        # print("Hist {} now has {} bins after rebinning with {} bins".format(hist.GetName(), hist.GetNbinsX(), len(xbins)))
    newHist = CheckHistBins(newHist)
    rootFile.cd()
    newHist.Write()


def CreateAndWriteHistograms(outputRootFilename):
    if not os.path.exists(datacardHistoDir):
        os.makedirs(datacardHistoDir)
    outputRootFilename = datacardHistoDir + "/" + outputRootFilename
    for i_signal_name, signal_name in enumerate(signal_names):
        outputRootFilename = outputRootFilename.replace(".root", "_{}.root".format(signal_name))
        for iSel, selectionName in enumerate(selectionNames):
            if selectionName == "preselection" or selectionName == "trainingSelection":
                continue
            mass_point = selectionName.replace("LQ", "")
            outputRootFile = r.TFile.Open(outputRootFilename.format(mass_point), "recreate")
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, mass_point)
            signalEvts = d_signal_rates[signalNameForFile][selectionName]
            signalEvtErrs = d_signal_rateErrs[signalNameForFile][selectionName]
            CreateAndWriteHist(outputRootFile, mass_point, fullSignalName, signalEvts, signalEvtErrs)
            signalFailEvts = d_signal_failRates[signalNameForFile][selectionName]
            signalFailEvtErrs = d_signal_failRateErrs[signalNameForFile][selectionName]
            CreateAndWriteHist(outputRootFile, mass_point, fullSignalName, signalFailEvts, signalFailEvtErrs, "yieldFailFinalSelection", "Yield failing final selection for LQ")
            for ibkg, background_name in enumerate(background_names):
                thisBkgEvts = d_background_rates[background_name][selectionName]
                thisBkgEvtsErr = d_background_rateErrs[background_name][selectionName]
                CreateAndWriteHist(outputRootFile, mass_point, background_name, thisBkgEvts, thisBkgEvtsErr)
                thisBkgFailEvts = d_background_failRates[background_name][selectionName]
                thisBkgFailEvtsErr = d_background_failRateErrs[background_name][selectionName]
                CreateAndWriteHist(outputRootFile, mass_point, background_name, thisBkgFailEvts, thisBkgFailEvtsErr, "yieldFailFinalSelection", "Yield failing final selection for LQ")
            CreateAndWriteHist(outputRootFile, mass_point, "DATA", d_data_rates["DATA"][selectionName], d_data_rateErrs["DATA"][selectionName])
            CreateAndWriteHist(outputRootFile, mass_point, "DATA", d_data_failRates["DATA"][selectionName], d_data_failRateErrs["DATA"][selectionName], "yieldFailFinalSelection", "Yield failing final selection for LQ")
            outputRootFile.Close()


def FillDicts(rootFilepath, sampleNames, bkgType, verbose=False):
    isData = False if "mc" in bkgType.lower() or "signal" in bkgType.lower() else True
    d_rates = {}
    d_rateErrs = {}
    d_unscaledRates = {}
    d_totalEvents = {}
    d_systematics = {}
    d_failRates = {}
    d_failRateErrs = {}
    d_unscaledFailRates = {}
    # start sample
    for i_sample, sampleName in enumerate(sampleNames):
        thisSampleRootFilepath = rootFilepath
        if ".root" not in rootFilepath:
            thisSampleRootFilepath += "analysisClass_lq_eejj_{}_plots.root".format(sampleName)
        scaledRootFile = r.TFile.Open(thisSampleRootFilepath)
        if not scaledRootFile or scaledRootFile.IsZombie():
            raise RuntimeError("Could not open root file: {}".format(scaledRootFile.GetName()))
        unscaledTotalEvts = cc.GetUnscaledTotalEvents(scaledRootFile, sampleName)
        if unscaledTotalEvts < 0:
            print(colored("WARN: for sample {}, found negative sampleUnscaledTotalEvents: {}; set to zero.".format(sampleName, unscaledTotalEvts), "red"))
            unscaledTotalEvts = 0.0
        ratesDict = {}
        rateErrsDict = {}
        unscaledRatesDict = {}
        failRatesDict = {}
        failRateErrsDict = {}
        unscaledFailRatesDict = {}
        systematicsNominalDict = {}
        # do all selections
        for i_signal_name, signal_name in enumerate(signal_names):
            for i_mass_point, mass_point in enumerate(selectionPoints):
                selectionName = selectionNames[i_mass_point]
                # print '------>Call GetRatesAndErrors for sampleName=',bkgSample
                sampleRate, sampleRateErr, sampleUnscaledRate = cc.GetRatesAndErrors(
                    scaledRootFile,
                    sampleName,
                    selectionName,
                    doEEJJ,
                    isData,
                    bkgType == "TTData"
                )
                sampleFailRate, sampleFailRateErr, sampleFailUnscaledRate = cc.GetFailingRatesAndErrors(
                    scaledRootFile,
                    sampleName,
                    selectionName,
                    trainingSelectionCutName,
                    doEEJJ,
                    isData,
                    bkgType == "TTData"
                )
                # print '------>rate=',rate,'rateErr=',rateErr,'unscaledRate=',unscaledRate
                # if isQCD:
                #  print 'for sample:',bkgSample,'got unscaled entries=',unscaledRate
                # print 'sampleName={}, sampleRate:',sampleRate,'sampleRateErr=',sampleRateErr,'sampleUnscaledRate=',sampleUnscaledRate
                if selectionName == "preselection" and verbose:
                    print("INFO: for sampleName={}, PRESELECTION ------>rate={} rateErr={} unscaledRate={} unscaledTotalEvts={}".format(sampleName, sampleRate, sampleRateErr, sampleUnscaledRate, unscaledTotalEvts))
                    print("INFO: for sampleName={}, PRESELECTION ------>failRate={} failRateErr={} unscaledFailRate={}".format(sampleName, sampleFailRate, sampleFailRateErr, sampleUnscaledFailRate))
                elif selectionName == "LQ1000" and verbose:
                    print("INFO: for sampleName={}, {} ------>rate={} rateErr={} unscaledRate={} unscaledTotalEvts={}".format(sampleName, selectionName, sampleRate, sampleRateErr, sampleUnscaledRate, unscaledTotalEvts))
                    print("INFO: for sampleName={}, {} ------>failRate={} failRateErr={} unscaledFailRate={}".format(sampleName, selectionName, sampleFailRate, sampleFailRateErr, sampleUnscaledFailRate))
                elif selectionName == "LQ1100" and verbose:
                    print("INFO: for sampleName={}, {} ------>rate={} rateErr={} unscaledRate={} unscaledTotalEvts={}".format(sampleName, selectionName, sampleRate, sampleRateErr, sampleUnscaledRate, unscaledTotalEvts))
                    print("INFO: for sampleName={}, {} ------>failRate={} failRateErr={} unscaledFailRate={}".format(sampleName, selectionName, sampleFailRate, sampleFailRateErr, sampleUnscaledFailRate))
                elif selectionName == "LQ1200" and verbose:
                    print("INFO: for sampleName={}, {} ------>rate={} rateErr={} unscaledRate={} unscaledTotalEvts={}".format(sampleName, selectionName, sampleRate, sampleRateErr, sampleUnscaledRate, unscaledTotalEvts))
                    print("INFO: for sampleName={}, {} ------>failRate={} failRateErr={} unscaledFailRate={}".format(sampleName, selectionName, sampleFailRate, sampleFailRateErr, sampleUnscaledFailRate))
                ratesDict[selectionName] = sampleRate
                rateErrsDict[selectionName] = sampleRateErr
                unscaledRatesDict[selectionName] = sampleUnscaledRate
                failRatesDict[selectionName] = sampleFailRate
                failRateErrsDict[selectionName] = sampleFailRateErr
                unscaledFailRatesDict[selectionName] = sampleFailUnscaledRate
                systematicsNominalDict[selectionName] = {"yield": sampleRate, "preselYield": None}
                if ratesDict[selectionName] < 0:
                    print("WARN: for sample {}, selection {}: found negative rate: {}; set to zero. Had {} unscaled events.".format(sampleName, selectionName, sampleRate, sampleUnscaledRate))
                    ratesDict[selectionName] = 0.0
                if unscaledRatesDict[selectionName] < 0:
                    print("WARN: for sample {}, selection {}: found negative unscaled rate: {}; set to zero.".format(sampleName, selectionName, sampleUnscaledRate))
                    unscaledRatesDict[selectionName] = 0.0
                if unscaledFailRatesDict[selectionName] < 0 and selectionName != "preselection":
                    print("WARN: for sample {}, selection {}: found negative unscaled fail rate: {}; set to zero.".format(sampleName, selectionName, sampleFailUnscaledRate))
                    unscaledFailRatesDict[selectionName] = 0.0
        if not isData and doSystematics:
            d_systematics[sampleName] = GetSystematicsDict(scaledRootFile, sampleName, selectionNames)  # , sampleName == "LQToDEle_M-300_pair")
            if "zjet" in sampleName.lower():
                d_systematics[sampleName].update({"DY_Norm": {sel: {"yield": dyNormDeltaXOverX, "preselYield": None} for sel in selectionNames}})
            if "ttbar" in sampleName.lower() or "ttto2l2nu" in sampleName.lower():
                d_systematics[sampleName].update({"TT_Norm": {sel: {"yield": ttBarNormDeltaXOverX, "preselYield": None} for sel in selectionNames}})
        elif sampleName == "QCDFakes_DATA":
            # special handling for QCDFakes_DATA
            d_selsAndSysts = {sel: {"yield": qcdNormDeltaXOverX, "preselYield": None} for sel in selectionNames}
            d_systematics[sampleName] = {"QCD_Norm": d_selsAndSysts}
            d_systematics[sampleName].update({"nominal": systematicsNominalDict})
        # fill full dicts
        d_rates[sampleName] = ratesDict
        d_rateErrs[sampleName] = rateErrsDict
        d_unscaledRates[sampleName] = unscaledRatesDict
        d_totalEvents[sampleName] = unscaledTotalEvts
        d_failRates[sampleName] = failRatesDict
        d_failRateErrs[sampleName] = failRateErrsDict
        d_unscaledFailRates[sampleName] = unscaledFailRatesDict
        scaledRootFile.Close()
    return d_rates, d_rateErrs, d_unscaledRates, d_totalEvents, d_systematics, d_failRates, d_failRateErrs, d_unscaledFailRates


def WriteDatacard(card_file_path):
    thresholdForAutoMCStats = 10
    yearStr = "# " + year + "\n"
    lumiStr = "# " + str(intLumi) + "\n\n"
    card_file = open(card_file_path, "w")
    card_file.write(yearStr)
    card_file.write(lumiStr)
    for i_signal_name, signal_name in enumerate(signal_names):
        doMassPointLoop = True
        #for i_mass_point, mass_point in enumerate(mass_points):
        for iSel, selectionName in enumerate(selectionNames):
            # fullSignalName = signal_name.replace('[masspoint]',mass_point)
            mass_point = selectionName.replace("LQ", "") if "LQ" in selectionName else "0"
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, mass_point)
            # print "consider fullSignalName={}".format(fullSignalName)
            #selectionName = "LQ" + mass_point
            # this will need to be fixed later if needed
            # else:
            #     # figure out mass point from name. currently the only case for this is RPV stop, where they are like 'Stop_M100_CTau100'
            #     mass = int(signal_name.split("_")[1].replace("M", ""))
            #     if mass < 200:
            #         mass = 200
            #     selectionName = "LQ" + str(mass)
            #     # print 'use selection name=',selectionName,'for fullSignalName=',fullSignalName
            #     doMassPointLoop = False
            if selectionName != "preselection" and selectionName != "trainingSelection":
                txt_file_name = datacardHistoDir.rstrip("/") + "/" + fullSignalName + ".txt"
                indivMassPointCard = open(txt_file_name, "w")
                indivMassPointCard.write(yearStr)
                indivMassPointCard.write(lumiStr)
    
                datacardLines = []
                datacardLines.append("# " + txt_file_name + "\n\n")
                datacardLines.append("imax " + str(n_channels) + "\n")
                datacardLines.append("jmax " + str(n_background) + "\n")
                if doSystematics:
                    datacardLines.append("kmax " + str(n_systematics) + "\n\n")
                else:
                    datacardLines.append("kmax 0\n\n")
                datacardLines.append("---------------\n")
                # datacardLines.append("shapes * * "+str(Path(shapeHistos_filePath).resolve())+" yieldFinalSelection_$MASS_$PROCESS\n")
                datacardLines.append("shapes * * "+shapeHistos_filePath+" yieldFinalSelection_$MASS_$PROCESS\n")
                datacardLines.append("---------------\n")
                datacardLines.append("bin bin1\n\n")
    
                total_data = d_data_rates["DATA"][selectionName]
                datacardLines.append("observation " + str(total_data) + "\n\n")
    
                line = "bin "
                for i_channel in range(0, n_background + 1):
                    line = line + "bin1 "
                datacardLines.append(line + "\n")
    
                line = "process " + fullSignalName + " "
                for background_name in background_names:
                    line = line + background_name + " "
                datacardLines.append(line + "\n")
    
                line = "process 0 "
                for background_name in background_names:
                    line = line + "1 "
                datacardLines.append(line + "\n\n")
    
                # rate line
                signalYield = d_signal_rates[signalNameForFile][selectionName]
                line = "rate {} ".format(signalYield)
                totalBackgroundYield = 0
                for background_name in background_names:
                    bkgYield = d_background_rates[background_name][selectionName]
                    line += "{} ".format(bkgYield)
                    totalBackgroundYield += bkgYield
                datacardLines.append(line.strip() + "\n")
                datacardLines.append("------------------------------\n")
                datacardLines.append("* autoMCStats "+str(thresholdForAutoMCStats)+"\n")
    
                # print signal_name, mass_point, total_signal, total_bkg, total_data
                # print signal_name+str(mass_point), total_signal, total_bkg
    
            # recall the form: systDict['PileupUp'/systematicFromHist]['ZJet_amcatnlo_ptBinned'/sampleName]['LQXXXX'/selection] = yield
            # for RPV, select proper signalSystDict based on ctau of signal
            # FIXME
            # if doRPV:
            #     ctau = int(signal_name[signal_name.find("CTau") + 4:])
            #     signalSystDict = signalSystDictByCTau[ctau]
            if doSystematics:
                d_systematicsApplied[selectionName] = {}
                d_systematicsApplied[selectionName][signalNameForFile] = []
                for syst in systematicsNamesBackground:
                    # if int(mass_point) > maxLQSelectionMass:
                    #     selectionNameSyst = "LQ"+str(maxLQSelectionMass)
                    # else:
                    #     selectionNameSyst = selectionName
                    selectionNameSyst = selectionName
                    line = syst + " lnN "
                    if syst in systematicsNamesSignal and selectionName != "preselection" and selectionName != "trainingSelection":
                        if signalNameForFile not in list(d_systNominals.keys()):
                            d_systNominals[signalNameForFile] = {}
                            d_systNominalErrs[signalNameForFile] = {}
                            d_systUpDeltas[signalNameForFile] = {}
                            d_systDownDeltas[signalNameForFile] = {}
                        if syst not in list(d_systNominals[signalNameForFile].keys()):
                            d_systNominals[signalNameForFile][syst] = {}
                            d_systNominalErrs[signalNameForFile][syst] = {}
                            d_systUpDeltas[signalNameForFile][syst] = {}
                            d_systDownDeltas[signalNameForFile][syst] = {}
                        systEntry, deltaOverNominalUp, deltaOverNominalDown, systNomYield, systSelection = GetSystematicEffectAbs(syst, signalNameForFile, selectionNameSyst, d_signal_systs)
                        thisSigEvts = d_signal_rates[signalNameForFile][selectionNameSyst]
                        thisSigEvtErr = d_signal_rateErrs[signalNameForFile][selectionNameSyst]
                        thisSigSystUp = deltaOverNominalUp*systNomYield
                        thisSigSystDown = deltaOverNominalDown*systNomYield
                        d_systNominals[signalNameForFile][syst][selectionNameSyst] = systNomYield
                        d_systNominalErrs[signalNameForFile][syst][selectionNameSyst] = d_signal_rateErrs[signalNameForFile][systSelection]
                        d_systUpDeltas[signalNameForFile][syst][selectionNameSyst] = thisSigSystUp
                        d_systDownDeltas[signalNameForFile][syst][selectionNameSyst] = thisSigSystDown
                        if EntryIsValid(systEntry):
                            d_systematicsApplied[selectionName][signalNameForFile].append(syst)
                            # if syst == "LHEPdfWeight": systEntry = 1.1  # hack signal PDF uncertainty
                            line += str(systEntry) + " "
                        else:
                            line += "- "
                    else:
                        line += "- "
                    for ibkg, background_name in enumerate(background_names):
                        if background_name not in d_systematicsApplied[selectionName]:
                            d_systematicsApplied[selectionName][background_name] = []
                        if background_name not in list(d_systNominals.keys()):
                            d_systNominals[background_name] = {}
                            d_systNominalErrs[background_name] = {}
                            d_systUpDeltas[background_name] = {}
                            d_systDownDeltas[background_name] = {}
                        if syst not in list(d_systNominals[background_name].keys()):
                            d_systNominals[background_name][syst] = {}
                            d_systNominalErrs[background_name][syst] = {}
                            d_systUpDeltas[background_name][syst] = {}
                            d_systDownDeltas[background_name][syst] = {}
                        systEntry, deltaOverNominalUp, deltaOverNominalDown, systNomYield, systSelection = GetSystematicEffectAbs(syst, background_name, selectionNameSyst, d_background_systs)
                        thisBkgEvts = d_background_rates[background_name][selectionNameSyst]
                        thisBkgEvtErr = d_background_rateErrs[background_name][selectionNameSyst]
                        thisBkgSystUp = deltaOverNominalUp*systNomYield
                        thisBkgSystDown = deltaOverNominalDown*systNomYield
                        d_systNominals[background_name][syst][selectionNameSyst] = systNomYield
                        d_systNominalErrs[background_name][syst][selectionNameSyst] = d_background_rateErrs[background_name][systSelection]
                        d_systUpDeltas[background_name][syst][selectionNameSyst] = thisBkgSystUp
                        d_systDownDeltas[background_name][syst][selectionNameSyst] = thisBkgSystDown
                        if EntryIsValid(systEntry):
                            d_systematicsApplied[selectionName][background_name].append(syst)
                            line += str(systEntry) + " "
                        else:
                            line += "- "
                        # if "ZJet" in background_name and "800" in selectionNameSyst and "EES" in syst:
                        try:
                            if float(systEntry) < 0:
                                 #print("INFO: For sample={} selection={} syst={}, systEntry={}, thisBkgEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                 raise RuntimeError("INFO: For sample={} selection={} syst={}, systEntry={}, thisBkgEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                         background_name, selectionNameSyst, syst, systEntry, thisBkgEvts, thisBkgSystUp, thisBkgSystDown))
                        except ValueError:
                            continue
                    # need to always fill the syst dicts, but only write the datacard if we have a BDT selection
                    if selectionName != "preselection" and selectionName != "trainingSelection":
                        datacardLines.append(line + "\n")
    
            if selectionName != "preselection" and selectionName != "trainingSelection":
                # rateParam for signal scaling
                signalScaleParam = 1.0
                datacardLines.append("signalScaleParam rateParam bin1 {} {}".format(fullSignalName, signalScaleParam))
                datacardLines.append("\n\n\n")
                for line in datacardLines:
                    card_file.write(line)
                    indivMassPointCard.write(line.replace("$MASS", mass_point).replace(".root", "_" + signal_name.format(mass_point) + ".root"))
                indivMassPointCard.close()
            if not doMassPointLoop:
                break
    card_file.close()


###################################################################################################
# CONFIGURABLES
###################################################################################################

blinded = True
doSystematics = True
doQCD = True
doEEJJ = True
signalNameTemplate = "LQToDEle_M-{}_pair"
# signalNameTemplate = "LQToBEle_M-{}_pair"
doRPV = False  # to do RPV, set doEEJJ and doRPV to True
# forceGmNNormBkgStatUncert = False
# cc.finalSelectionName = "sT_eejj"  # "min_M_ej"
cc.finalSelectionName = "BDTOutput_"  # has LQ300 (for example) appended to it inside the code
# cc.finalSelectionName = "Meejj"
trainingSelectionCutName = "trainingSelection"
if len(sys.argv) < 4:
    print("ERROR: need to specify qcdPath, dataMCPath, year")
    print("Example: makeDatacard.py root://path/to/qcdSubtracted_plots.root root://path/to/analysisClass_lq_eejj_plots.root 2016preVFP")
    raise RuntimeError("ERROR: need to specify qcdPath, dataMCPath, year")

qcdPath = sys.argv[1]
dataMCPath = sys.argv[2]
year = sys.argv[3]

if not signalNameTemplate.split("_")[0] in dataMCPath or not signalNameTemplate.split("_")[0] in qcdPath:
    raise RuntimeError("signalNameTemplate specified is {}, while it does not appear to match the dataMC or qcdPaths given.".format(signalNameTemplate))

sampleListForMerging = "$LQANA/config/sampleListForMerging_13TeV_eejj_{}.yaml"
sampleListsForMergingQCD = "$LQANA/config/sampleListForMerging_13TeV_QCD_dataDriven_{}.yaml"

if doEEJJ:
    sampleListForMerging = os.path.expandvars(sampleListForMerging.format(year))
    # SIC 6 Jul 2020 remove
    # sampleListForMergingTTBar = (
    #     os.environ["LQANA"] + "/config/sampleListForMerging_13TeV_ttbarBkg_emujj.txt"
    # )
    if doQCD:
        sampleListForMergingQCD = os.path.expandvars(sampleListsForMergingQCD.format(year))
        qcdFilePath = os.path.expandvars(qcdPath)
    filePath = os.path.expandvars(dataMCPath)
    dataMC_filepath = filePath.rstrip("/") + "/"
    # ttbarFilePath = (
    #     os.environ["LQDATA"]
    #     + "/2016ttbar/mar17_emujj_fixMuons/output_cutTable_lq_ttbar_emujj_correctTrig/"
    # )
    # # calculated with fitForStatErrs.py script. mass, stat. uncert.
    # # statErrorsSingleTop = { 800: 1.10, 850: 0.85, 900: 0.67, 950: 0.53, 1000: 0.42, 1050: 0.33 }
    # # here we increased the fit range
    # statErrorsPhotonJets = {
    #     400: 0.23,
    #     450: 0.33,
    #     500: 0.36,
    #     550: 0.35,
    #     600: 0.32,
    #     650: 0.28,
    #     700: 0.24,
    #     750: 0.21,
    #     800: 0.17,
    #     850: 0.14,
    #     900: 0.12,
    #     950: 0.095,
    #     1000: 0.077,
    #     1050: 0.062,
    # }
else:
    sampleListForMerging = (
        os.environ["LQANA"] + "/config/sampleListForMerging_13TeV_enujj.txt"
    )
    sampleListForMergingQCD = (
        os.environ["LQANA"] + "/config/sampleListForMerging_13TeV_QCD_dataDriven.txt"
    )
    filePath = (
        os.environ["LQDATA"]
        + "/2016analysis/enujj_psk_mar16_fixMuons/output_cutTable_lq_enujj_MT/"
    )
    qcdFilePath = (
        os.environ["LQDATA"]
        + "/2016qcd/enujj_mar16_fixMuons/output_cutTable_lq_enujj_MT_QCD/"
    )
    # # calculated with fitForStatErrs.py script. mass, stat. uncert.
    # statErrorsSingleTop = { 650:0.213, 700:0.330, 750:0.399, 800:0.431, 850:0.442, 900:0.438, 950:0.425, 1000:0.407, 1050:0.385, 1100:0.363, 1150:0.342, 1200:0.321 }

# this has the TopPtReweight+updatedSF and the Z+jets St corrections at final selections
# filePath = os.environ["LQDATA"] + '/RunII/eejj_analysis_zJetsStCorrectionFinalSelections_21jul/output_cutTable_lq_eejj/'

if doEEJJ:
    # ttbar_data_filepath = ttbarFilePath + "analysisClass_lq_ttbarEst_plots.root"
    # SIC 6 Jul 2020 remove
    ttbar_data_filepath = ""
else:
    ttbar_data_filepath = ""
    # ttbarFilePath = filePath

do2016preVFP = False
do2016postVFP = False
do2017 = False
do2018 = False
if year == "2016preVFP":
    do2016preVFP = True
    cc.intLumi = 19501.601622000
    yearInt = 2016
elif year == "2016postVFP":
    do2016postVFP = True
    cc.intLumi = 16812.151722000
    yearInt = 2016
elif year == "2017":
    do2017 = True
    cc.intLumi = 41477.877399
    yearInt = 2017
elif year == "2018":
    do2018 = True
    cc.intLumi = 59827.449483
    yearInt = 2018
else:
    print("ERROR: could not find one of 2016preVFP/2016postVFP/2017/2018 as year. Quitting.")
    exit(-1)
intLumi = cc.intLumi
do2016 = do2016preVFP or do2016postVFP
if doRPV:
    mass_points = [
        str(i) for i in range(200, 1250, 100)
    ]  # go from 200-1200 in 100 GeV steps
else:
    # LQ case
    mass_points = [
        str(i) for i in range(300, 3100, 100)
    ]  # go from 300-3000 in 100 GeV steps
    ## mass_points.extend(["3500", "4000"])

if doEEJJ:
    if doRPV:
        signal_names = [
            "Stop_M[masspoint]_CTau1000",
            "Stop_M[masspoint]_CTau100",
            "Stop_M[masspoint]_CTau10",
            "Stop_M[masspoint]_CTau1",
        ]
        # signal_names = [ "Stop_M[masspoint]_CTau10","Stop_M[masspoint]_CTau1"]
        # put in some more signals that don't fit the general pattern
        # signal_names = ['Stop_M100_CTau100','Stop_M125_CTau100','Stop_M150_CTau100','Stop_M175_CTau100','Stop_M200_CTau50'] + signal_names
        signalNameTemplate = "stop"
    else:
        signal_names = [signalNameTemplate]
    # maxLQSelectionMass = 1200  # max background selection point used
    maxLQSelectionMass = 100000 # FIXME remove
    # background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTBarFromDATA", "ZJet_amcatnlo_ptBinned", "WJet_amcatnlo_ptBinned", "DIBOSON","SingleTop"  ]
    background_names = [
        # "ZJet_powhegminnlo",
        "ZJet_amcatnlo_ptBinned_IncStitch",
        # "TTBarFromDATA",
        "TTTo2L2Nu"] + (["QCDFakes_DATA"] if doQCD else []) + [
        "DIBOSON_nlo",
        # "TRIBOSON",
        # "TTX",
        "SingleTop",
        # "WJet_amcatnlo_jetBinned",
        # "PhotonJets_Madgraph",
    ]
    background_fromMC_names = [bkg for bkg in background_names if "data" not in bkg.lower()]
    background_QCDfromData = [bkg for bkg in background_names if "data" in bkg.lower() and "qcd" in bkg.lower()]
    d_systTitles = {
            "EleTrigSF": "Trigger",
            "Pileup": "Pileup",
            "Prefire": "L1EcalPrefiring",
            "LHEPdfWeight": "PDF",
            "Lumi"+str(yearInt): "Lumi",
            "LumiCorrelated": "Lumi correlated",
            "JER": "Jet energy resolution",
            "JES": "Jet energy scale",
            "EleRecoSF": "Electron reconstruction",
            "EleIDSF": "Electron identification",
            "EES": "Electron energy scale",
            "EER": "Electron energy resolution",
            "UnclusteredEne": "Unclustered energy (MET)",
            "TT_Norm": "TTbar normalization",
            "DY_Norm": "DYJ normalization",
            "DY_Shape": "DYJ shape",
            "TT_Shape": "TTbar shape",
            "Diboson_Shape": "Diboson shape",
            "ST_Shape": "SingleTop shape",
            "QCD_Norm": "QCD bkg. normalization",
    }
    if do2017 or do2018:
        d_systTitles["LumiCorrelated1718"] =  "Lumi correlated 2017-2018"
    otherBackgrounds = [
        # "PhotonJets_Madgraph",
        # "WJet_amcatnlo_jetBinned",
        "DIBOSON_nlo",
        # "TRIBOSON",
        # "TTW",
        # "TTZ",
        "SingleTop",
    ]
else:
    signal_names = ["LQ_BetaHalf_M[masspoint]"]
    signalNameTemplate = "betaHalfLQ"
    # background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTbar_amcatnlo_Inc", "ZJet_amcatnlo_ptBinned", "WJet_amcatnlo_ptBinned", "DIBOSON","SingleTop"  ]
    # background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTbar_powheg", "ZJet_amcatnlo_ptBinned", "WJet_amcatnlo_ptBinned", "DIBOSON","SingleTop"  ]
    background_names = [
        "PhotonJets_Madgraph",
        "QCDFakes_DATA",
        "TTbar_powheg",
        "ZJet_amcatnlo_ptBinned",
        "WJet_amcatnlo_ptBinned",
        "DIBOSON_amcatnlo",
        "SingleTop",
    ]
    syst_background_names = [
        "GJets",
        "QCDFakes_DATA",
        "TTbar",
        "DY",
        "WJets",
        "Diboson",
        "Singletop",
    ]
    maxLQSelectionMass = 900  # max background selection point used
    systematicsNamesSignal = [
        "Trigger",
        "Reco",
        "PU",
        # by EXO convention, no PDF systs for signals
        # "PDF",
        "Lumi"+str(yearInt),
        "LumiCorrelated",
        "JER",
        "JEC",
        "HEEP",
        "E_scale",
        "EER",
        "MET",
    ]
    if do2017 or do2018:
        systematicsNamesSignal.append("LumiCorrelated1718")
    otherBackgrounds = [
        "PhotonJets_Madgraph",
        "ZJet_amcatnlo_ptBinned",
        "DIBOSON_amcatnlo",
        "SingleTop",
    ]

minLQselectionBkg = "LQ200"
zjetsSampleName = GetSampleNameFromSubstring("ZJet", background_names)
#ttbarSampleName = GetSampleNameFromSubstring("TTbar", background_names)
ttbarSampleName = GetSampleNameFromSubstring("TTTo2L2Nu", background_names)
dibosonSampleName = GetSampleNameFromSubstring("DIBOSON", background_names)
singleTopSampleName = GetSampleNameFromSubstring("SingleTop", background_names)
systematicsNamesBackground = d_systTitles.keys()
systematicsNamesSignal = [syst for syst in systematicsNamesBackground if "shape" not in syst.lower() and "norm" not in syst.lower()]
systematicsNamesSignal.remove("LHEPdfWeight")

backgroundTitlesDict = {zjetsSampleName: "Z+jets", ttbarSampleName: "TTbar", "QCDFakes_DATA": "QCD(data)", "WJet_amcatnlo_ptBinned": "W+jets", "WJet_amcatnlo_jetBinned": "W+jets",
        dibosonSampleName: "DIBOSON", "TRIBOSON": "TRIBOSON", "TTX": "TTX", "SingleTop": "SingleTop", "PhotonJets_Madgraph": "Gamma+jets"}
backgroundTitles = [backgroundTitlesDict[bkg] for bkg in background_names]
# SIC 6 Jul 2020 remove
# if doEEJJ:
#     ttBarNormDeltaXOverX = 0.01
#     ttbarSampleName = "TTBarFromDATA"
#     ttBarUnscaledRawSampleName = "TTBarUnscaledRawFromDATA"
#     # nonTTBarSampleName='NONTTBARBKG_amcatnloPt_emujj'
#     nonTTBarSampleName = "NONTTBARBKG_amcatnloPt_amcAtNLODiboson_emujj"

selectionPoints = ["preselection", "trainingSelection"] + mass_points
selectionNames = ["LQ"+sel if "selection" not in sel.lower() else sel for sel in selectionPoints]
additionalBkgSystsDict = {}
# QCDNorm is 0.40 [40% norm uncertainty for eejj = uncertaintyPerElectron*2]
# lumi uncertainty from https://twiki.cern.ch/twiki/bin/view/CMS/LumiRecommendationsRun2#Combination_and_correlations
dyNormDeltaXOverX = 0.1
ttBarNormDeltaXOverX = 0.1
if doEEJJ:
    qcdNormDeltaXOverX = 0.4
    if do2016:
        lumiDeltaXOverX = 0.01
        lumiCorrelatedDeltaXOverX = 0.006
        # from: https://arxiv.org/pdf/1506.04072.pdf tables 1-3: 5% is probably good enough for SingleTop PDF syst
        # additionalBkgSystsDict["SingleTop"] = {"LHEPdfWeight": {sel: 0.10 for sel in selectionNames}}
    elif do2017:
        lumiDeltaXOverX = 0.02
        lumiCorrelatedDeltaXOverX = 0.009
        lumi1718CorrelatedDeltaXOverX = 0.006
    elif do2018:
        lumiDeltaXOverX = 0.015
        lumiCorrelatedDeltaXOverX = 0.02
        lumi1718CorrelatedDeltaXOverX = 0.002
else:
    qcdNormDeltaXOverX = 0.2

n_background = len(background_names)
# all bkg systematics, plus stat 'systs' for all bkg plus signal plus 3 backNormSysts
# update July/Sep 2020: no more data-driven TBar, so only 1 extra syst for eejj
if doSystematics:
    if doEEJJ:
        # systematics + stat_errors_background + stat_errors_signal
        # n_systematics = len(systematicsNamesBackground) + n_background + 1
        n_systematics = len(systematicsNamesBackground)
    else:
        n_systematics = (
            len(systematicsNamesBackground) + n_background + 1 + 1
        )  # QCD norm only
else:
    n_systematics = n_background + 1
n_channels = 1

signalNameList = [GetFullSignalName(signalNameTemplate, massPoint)[1] for massPoint in mass_points]
allBkgSysts = [syst for syst in systematicsNamesBackground if "norm" not in syst.lower() and "shape" not in syst.lower()]
d_applicableSystematics = {bkg: list(allBkgSysts) for bkg in background_fromMC_names}
d_applicableSystematics.update({bkg: "QCD_Norm" for bkg in background_QCDfromData})
d_applicableSystematics[zjetsSampleName].append("DY_Norm")
d_applicableSystematics[zjetsSampleName].append("DY_Shape")
backgroundsToRenormSystAtPresel = [zjetsSampleName, ttbarSampleName]
d_applicableSystematics[ttbarSampleName].append("TT_Norm")
d_applicableSystematics[ttbarSampleName].append("TT_Shape")
d_applicableSystematics[dibosonSampleName].append("Diboson_Shape")
d_applicableSystematics[singleTopSampleName].append("ST_Shape")
d_applicableSystematics.update({sig: systematicsNamesSignal for sig in signalNameList})
datacard_filePath = "tmp_card_file_{}.txt".format(signalNameTemplate.split("_")[0])
plots_filePath = "plots.root"
plotsDir = "makeDatacard_plots"
datacardHistoDir = "makeDatacard_separateMassPoints"
tablesDir = "makeDatacard_tables"
systematics_dictFilePath = "systematics_dict.txt"
shapeHistos_filePathBase = "shapeHistos_{}.root"
###################################################################################################
# RUN
###################################################################################################

# ---Check if sampleListForMerging file exist
if os.path.isfile(sampleListForMerging) is False:
    print("ERROR: file " + sampleListForMerging + " not found")
    print("exiting...")
    sys.exit()

# ---Check if sampleListForMergingQCD file exist
if doQCD:
    if os.path.isfile(sampleListForMergingQCD) is False:
        print("ERROR: file " + sampleListForMergingQCD + " not found")
        print("exiting...")
        sys.exit()

print("Launched like:")
for arg in sys.argv:
    print("\t" + arg)
print("Using input root paths/files:")
print("\t Data/MC:", dataMC_filepath)
if doQCD:
    print("\t QCD(data):", qcdFilePath)
else:
    print("\t No QCD used.")
    qcdFilePath = None

dictSamples = cc.GetSamplesToCombineDict(sampleListForMerging)
dictSamples = {k: v for k, v in dictSamples.items() if dictSamples[k]["save"]}

d_systUpDeltas = {}
d_systDownDeltas = {}
d_systNominals = {}  # same as rates, but more conveniently indexed
d_systNominalErrs = {}
d_systematicsApplied = {}
d_datacardStatErrs = {}

# SIC 6 Jul 2020 remove
# if doEEJJ:
#     dictSamplesTTBarRaw = GetSamplesToCombineDict(sampleListForMergingTTBar)
#     # only care about the TTBar parts
#     dictSamplesTTBar = {}
#     # for key in dictSamplesTTBarRaw.iterkeys():
#     #  if 'ttbar' in key.lower():
#     #    if 'ttbarunscaledrawfromdata' in key.lower():
#     #      print 'set dictSamplesTTBar[TTBarFromDATA] =',dictSamplesTTBarRaw[key],'for key:',key
#     #      dictSamplesTTBar['TTBarFromDATA'] = dictSamplesTTBarRaw[key]
#     # NB: we rely on this exact sample name for the total TTBar data-driven sample
#     # dictSamplesTTBar['TTBarFromDATA'] = dictSamplesTTBarRaw['TTBarUnscaledRawFromDATA']
#     # dictSamplesTTBar['TTBarFromDATA'] = ['TTBarFromDATA']
#     dictSamplesTTBar[ttbarSampleName] = [ttbarSampleName]
#     dictSamples.update(dictSamplesTTBar)
#     print "found ttbar samples:", dictSamplesTTBar

# rates/etc.
print("INFO: Filling background [MC] information...")
#d_background_rates, d_background_rateErrs, d_background_unscaledRates, d_background_totalEvents, d_background_systs = FillDicts(dataMC_filepath, background_fromMC_names, "MC")
d_background_rates, d_background_rateErrs, d_background_unscaledRates, d_background_totalEvents, d_background_systs, d_background_failRates, d_background_failRateErrs, d_background_unscaledFailRates = FillDicts(dataMC_filepath, list(dictSamples.keys()), "MC")
if doQCD:
    print("INFO: Filling QCD[data] information...")
    bgFromData_rates, bgFromData_rateErrs, bgFromData_unscaledRates, bgFromData_totalEvents, bgFromData_systs, bgFromData_failRates, bgFromData_failRateErrs, bgFromData_unscaledFailRates = FillDicts(qcdFilePath, background_QCDfromData, "DATA")
    d_background_rates.update(bgFromData_rates)
    d_background_rateErrs.update(bgFromData_rateErrs)
    d_background_unscaledRates.update(bgFromData_unscaledRates)
    d_background_totalEvents.update(bgFromData_totalEvents)
    d_background_systs.update(bgFromData_systs)
    d_background_failRates.update(bgFromData_failRates)
    d_background_failRateErrs.update(bgFromData_failRateErrs)
    d_background_unscaledFailRates.update(bgFromData_unscaledFailRates)
if doSystematics:
    for sampleName in additionalBkgSystsDict.keys():
        d_background_systs[sampleName].update(additionalBkgSystsDict[sampleName])
# above would be similar for TTBarFromDATA
print("INFO: Filling signal information...")
d_signal_rates, d_signal_rateErrs, d_signal_unscaledRates, d_signal_totalEvents, d_signal_systs, d_signal_failRates, d_signal_failRateErrs, d_signal_unscaledFailRates = FillDicts(dataMC_filepath, signalNameList, "signal")
print("INFO: Filling data information...")
d_data_rates, d_data_rateErrs, d_data_unscaledRates, d_data_totalEvents, d_data_systs_garbage, d_data_failRates, d_data_failRateErrs, d_data_unscaledFailRates = FillDicts(dataMC_filepath, ["DATA"], "DATA")
# print one of the systematics for checking
# for syst in backgroundSystDict.keys():
#    print 'Syst is:',syst
#    print 'selection\t\tvalue'
#    for selection in sorted(backgroundSystDict[syst].keys()):
#        print selection+'\t\t'+str(backgroundSystDict[syst][selection])
#    break
# print signalSystDict
# print backgroundSystDict

# selectionNames = ["preselection"]
# selectionNames.extend(["LQ"+str(mass) for mass in mass_points])

print("INFO: Preparing shape histograms...", end=' ')
shapeHistos_filePath = shapeHistos_filePathBase.format(year)
CreateAndWriteHistograms(shapeHistos_filePath)
print("Done")

print("INFO: Preparing datacard...", end=' ')
WriteDatacard(datacard_filePath)
print("Done")

if not os.path.exists(tablesDir):
    os.makedirs(tablesDir)

if doSystematics:
    # tables
    columnNames = ["Systematic", "Signal (%)", "Background (%)"]
    for selectionName in selectionNames:
        table = []
        for syst in systematicsNamesBackground:
            print("INFO: total syst table for selection {}".format(selectionName))
            selectionNameSyst = selectionName
            if selectionName != "preselection" and selectionName != "trainingSelection":
                massPoint = selectionName.replace("LQ", "")
                if int(massPoint) > maxLQSelectionMass:
                    selectionNameSyst = "LQ"+str(maxLQSelectionMass)
            if syst in systematicsNamesSignal and selectionName != "preselection" and selectionName != "trainingSelection":
                for i_signal_name, signal_name in enumerate(signal_names):
                    fullSignalName, signalNameForFile = GetFullSignalName(signal_name, massPoint)
                    thisEntry, sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown, _, _ = GetSystematicEffectAbs(syst, signalNameForFile, selectionNameSyst, d_signal_systs)
                    if thisEntry == "-":
                        thisSigSystPercent = 0
                        continue
                    thisSigSyst = max(sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown)
                    thisSigSystPercent = round(100*thisSigSyst, 1)
            else:
                thisSigSystPercent = "  - "
            totalBkgSyst = 0
            totalBkgNominal = 0
            for background_name in background_names:
                # compute effect of this systematic on the total background yield
                thisBkgNominalEvts = d_background_rates[background_name][selectionName]
                totalBkgNominal += d_systNominals[background_name][syst][selectionNameSyst]
                # print("INFO: for background {}, systematic {}, selection {}, bkgYieldNominalAtSyst={}, bkgEvts={}".format(background_name, syst, selectionNameSyst, d_systNominals[background_name][syst][selectionNameSyst], thisBkgNominalEvts))
                thisEntry, bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown, _, _ = GetSystematicEffectAbs(syst, background_name, selectionNameSyst, d_background_systs)
                if thisEntry != "-":
                    thisBkgSyst = max(bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown)
                    thisBkgDelta = thisBkgSyst*d_systNominals[background_name][syst][selectionNameSyst]
                    totalBkgSyst += thisBkgDelta*thisBkgDelta
                    # print("INFO: for background {}, systematic {}, selection {}, thisBkgSyst={}, bkgYieldNominalAtSyst={}, thisBkgDelta={}.".format(background_name, syst, selectionNameSyst, thisBkgSyst, d_systNominals[background_name][syst][selectionNameSyst], thisBkgDelta))
                else:
                    print("INFO: for background {}, systematic {}, selection {}, syst not applied as thisEntry='{}'.".format(background_name, syst, selectionNameSyst, thisEntry))
            if totalBkgNominal > 0:
                totalBkgSystPercent = 100*(math.sqrt(totalBkgSyst))/totalBkgNominal
                print("INFO: for systematic {}, selection {}, totalBkgSyst={}, totalBkgNominal={}.".format(syst, selectionNameSyst, math.sqrt(totalBkgSyst), totalBkgNominal))
            else:
                totalBkgSystPercent = -1
            table.append([d_systTitles[syst], thisSigSystPercent, totalBkgSystPercent])
        print("Selection: {}".format(selectionName))
        print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".1f"))
        print(tabulate(table, headers=columnNames, tablefmt="latex", floatfmt=".1f"))
        print()
        with open(tablesDir+"/"+selectionName+".txt", "w") as table_file:
            table_file.write(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".1f"))
        with open(tablesDir+"/"+selectionName+".tex", "w") as table_file:
            table_file.write(tabulate(table, headers=columnNames, tablefmt="latex", floatfmt=".1f"))

    # print info on systematics used
    if "LQ300" in d_systematicsApplied.keys():
        print("For selection LQ300")
        lq300systsDict = d_systematicsApplied["LQ300"]
        print("{0:40}\t{1}".format("sampleName", "systematics applied"))
        for sampleName in sorted(lq300systsDict.keys()):
            if "LQ" in sampleName:
                continue
            print("{0:40}\t{1}".format(sampleName, lq300systsDict[sampleName]))
        for sampleName in sorted(lq300systsDict.keys()):
            if "LQ" not in sampleName:
                continue
            print("{0:40}\t{1}".format(sampleName, lq300systsDict[sampleName]))
        print()

    # compute total syst for all backgrounds at each selection
    d_totalSystDeltaXOverX = {}
    for selectionName in selectionNames:
        if selectionName not in list(d_totalSystDeltaXOverX.keys()):
            d_totalSystDeltaXOverX[selectionName] = {}
        otherBkgTotalDeltaXOverXSqr = 0
        for background_name in background_names:
            thisBkgTotalSystOverNomSqr = 0
            for syst in systematicsNamesBackground:
                selectionNameSyst = selectionName
                if selectionName != "preselection" and selectionName != "trainingSelection":
                    massPoint = selectionName.replace("LQ", "")
                    if int(massPoint) > maxLQSelectionMass:
                        selectionNameSyst = "LQ"+str(maxLQSelectionMass)
                thisEntry, bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown, _, _ = GetSystematicEffectAbs(syst, background_name, selectionNameSyst, d_background_systs)
                if thisEntry != "-":
                    thisBkgSystOverNom = max(bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown)
                    print("INFO: selection={}, background_name={}, syst={}, deltaX/X={}".format(selectionName, background_name, syst, thisBkgSystOverNom))
                    thisBkgTotalSystOverNomSqr += thisBkgSystOverNom*thisBkgSystOverNom
            print("INFO: selection={}, background_name={}, total deltaX/X={}".format(selectionName, background_name, math.sqrt(thisBkgTotalSystOverNomSqr)))
            d_totalSystDeltaXOverX[selectionName][background_name] = math.sqrt(thisBkgTotalSystOverNomSqr)
            if background_name in otherBackgrounds:
                otherBkgTotalDeltaXOverXSqr += thisBkgSystOverNom*thisBkgSystOverNom
        d_totalSystDeltaXOverX[selectionName]["otherBackgrounds"] = math.sqrt(otherBkgTotalDeltaXOverXSqr)
        for i_signal_name, signal_name in enumerate(signal_names):
            if selectionName == "preselection" and selectionName != "trainingSelection":
                d_totalSystDeltaXOverX[selectionName][signalNameForFile] = 0.0
                continue
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, massPoint)
            thisSigTotalSystOverNomSqr = 0
            for syst in systematicsNamesSignal:
                thisEntry, sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown, _, _ = GetSystematicEffectAbs(syst, signalNameForFile, selectionNameSyst, d_signal_systs)
                if thisEntry != "-":
                    thisSigSystOverNom = max(sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown)
                    print("INFO: selection={}, signal_name={}, syst={}, deltaX/X={}".format(selectionName, signal_name, syst, thisSigSystOverNom))
                    thisSigTotalSystOverNomSqr += thisSigSystOverNom*thisSigSystOverNom
            print("INFO: selection={}, signal_name={}, total deltaX/X={}".format(selectionName, signal_name, math.sqrt(thisSigTotalSystOverNomSqr)))
            d_totalSystDeltaXOverX[selectionName][signalNameForFile] = math.sqrt(thisSigTotalSystOverNomSqr)

    with open(systematics_dictFilePath, "w") as systematics_dictFile:
        systematics_dictFile.write(str(d_totalSystDeltaXOverX))
    print("systematics dict written to {}".format(systematics_dictFilePath))

    # make the plots
    if not os.path.exists(plotsDir):
        os.makedirs(plotsDir)
    cc.SetDistinctiveTColorPalette()
    plots_tfile = r.TFile(plots_filePath, "recreate")
    systDir = plots_tfile.mkdir("systematics")
    allSystCanvasNames = []
    d_totalBkgSysts = {}
    totalBkgSystsHeaders = []
    includeYieldLine = True
    for sampleName in list(d_systNominals.keys()):
        yieldRow = ["yield"]
        mcRow = ["mcEvts"]
        yieldTable = []
        table = []
        systDeltaTable = []
        tableHeaders = ["systematic"]
        yieldTableHeaders = ["yield"]
        systDeltaTableHeaders = ["delta"]
        systStacks = []
        systStacks.append(r.THStack())
        systStacks[-1].SetName("allSystsStack_"+sampleName)
        systematicsByNameAndMass = {}
        systematicsSelectionsByNameAndMass = {}
        for iSyst, syst in enumerate(d_systNominals[sampleName].keys()):
            if not DoesSystematicApply(syst, sampleName):
                continue
            if not syst in d_totalBkgSysts.keys():
                d_totalBkgSysts[syst] = {}
            if not "totalYield" in d_totalBkgSysts.keys():
                d_totalBkgSysts["totalYield"] = {}
            tableRow = [d_systTitles[syst]]
            systYieldRow = [d_systTitles[syst]]
            systematicsByNameAndMass[d_systTitles[syst]] = []
            systematicsSelectionsByNameAndMass[d_systTitles[syst]] = []
            sortedSelectionDict = sorted(
                    iter(d_systNominals[sampleName][syst].items()), key=lambda x: float(x[0].replace("LQ", "").replace("preselection", "0").replace("trainingSelection", "1")))
            systDir.cd()
            systDir.mkdir(syst, syst, True).cd()
            if "LQ" in sampleName:
                signalMass = int(sampleName.split("-")[-1].split("_")[0])
                systHists = [r.TH1D(syst+"_uncertVsMassHist_{}".format(sampleName), syst, 2, signalMass-100, signalMass+100)]
                massRanges = [signalMass-100, signalMass+100]
            else:
                # split into 3 ranges
                splitInto = 3
                massPoints = [int(selection.replace("LQ", "")) for selection, _ in sortedSelectionDict if "selection" not in selection.lower()]
                minMass = min(massPoints)
                maxMass = max(massPoints)
                subRange = math.ceil((maxMass+100-minMass) / (1000*splitInto))*1000
                bins = int(subRange / 100)  # bins 100 GeV wide (one per mass point)
                systHists = []
                massRanges = []
                if len(systStacks) < len(systHists):
                    systStacks[0].SetName(systStacks[0].GetName().replace("Stack_", "Stack_MassRange1_"))
                for idx in range(0, splitInto):
                    minRange = idx*subRange+minMass-50
                    maxRange = (idx+1)*subRange+minMass-50
                    systHists.append(r.TH1D(syst+"_uncertVsMassHist{}_{}".format(idx, sampleName), syst, bins, minRange, maxRange))
                    massRanges.append(minRange)
                    if len(systStacks) < len(systHists) and idx > 0:
                        systStacks.append(r.THStack())
                        systStacks[-1].SetName("allSystsStack_MassRange{}_".format(idx+1)+sampleName)
            massList = []
            nominals = []
            upVariation = []
            downVariation = []
            thisSystTableHeaders = ["% "+sName[0] for sName in sortedSelectionDict]
            if len(tableHeaders) == 1:
                tableHeaders.extend(thisSystTableHeaders)
                yieldTableHeaders.extend([sName[0] for sName in sortedSelectionDict])
                systDeltaTableHeaders.extend([sName[0] for sName in sortedSelectionDict])
            elif len(tableHeaders) - 1 != len(thisSystTableHeaders):
                raise RuntimeError("For sample {}, number of selections {} for this syst {} is larger than the other systs [{}]".format(
                    sampleName, syst, len(thisSystTableHeaders), len(tableHeaders) - 1))
            for selection, value in sortedSelectionDict:
                if selection not in d_totalBkgSysts[syst].keys():
                    d_totalBkgSysts[syst][selection] = 0
                if selection not in d_totalBkgSysts["totalYield"].keys():
                    d_totalBkgSysts["totalYield"][selection] = 0
                #if selection == "preselection":
                #    continue
                # massList.append(mass)
                # nominals.append(0.0)  # nominals.append(float(value))
                fillVal = "N/A"
                valid = False
                if value > 0:
                    upVariation.append(100*float(d_systUpDeltas[sampleName][syst][selection])/value)
                    downVariation.append(100*float(d_systDownDeltas[sampleName][syst][selection])/value)
                    fillVal = max(100*float(d_systUpDeltas[sampleName][syst][selection])/value,
                                  100*float(d_systDownDeltas[sampleName][syst][selection])/value)
                    systYieldVal = max(float(d_systUpDeltas[sampleName][syst][selection]),
                                  float(d_systDownDeltas[sampleName][syst][selection]))
                    # print("DEBUG: for sample {}, syst {}, fill value (% deltaX/X) = {}".format(sampleName, syst, fillVal), flush=True)
                    # tableRow.append(fillVal)
                    if selection != "preselection" and selection != "trainingSelection":
                        mass = float(selection.replace("LQ", ""))
                        idxToFill = bisect(massRanges, mass) - 1
                        systHists[idxToFill].Fill(mass, fillVal)
                        valid = True
                if valid:
                    systematicsByNameAndMass[d_systTitles[syst]].append(fillVal if valid else 0)
                    systematicsSelectionsByNameAndMass[d_systTitles[syst]].append(mass)
                tableRow.append(fillVal)
                # systYieldRow.append(systYieldVal+value)
                systYieldRow.append(systYieldVal)
                if sampleName in background_names:
                    toAdd = max(float(d_systUpDeltas[sampleName][syst][selection]), float(d_systDownDeltas[sampleName][syst][selection]))
                    d_totalBkgSysts[syst][selection] += toAdd*toAdd
                    d_totalBkgSysts["totalYield"][selection] += value
                # else:
                #     lastSelection, rate, err, events = GetLastNonzeroSelectionYields(background_name)
                #     value = rate
                #     upVariation.append(100*float(d_systUpDeltas[sampleName][syst][lastSelection])/value)
                #     downVariation.append(100*float(d_systDownDeltas[sampleName][syst][lastSelection])/value)
                #     fillVal = max(100*float(d_systUpDeltas[sampleName][syst][lastSelection])/value,
                #                   100*float(d_systDownDeltas[sampleName][syst][lastSelection])/value)
                #     tableRow.append(fillVal)
                #     if selection != "preselection":
                #         mass = float(selection.replace("LQ", ""))
                #         systHist.Fill(mass, fillVal)
                # if "ZJet" in sampleName and "800" in selection and "EES" in syst:
                #     print "INFO: For sample={} selection={} syst={}, nominal={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                #             sampleName, selection, syst, value, d_systUpDeltas[sampleName][syst][selection],
                #             d_systDownDeltas[sampleName][syst][selection])
                if len(yieldRow) < len(tableHeaders):
                    # err = d_systNominalErrs[sampleName][syst][selection]
                    # err = format(Decimal(format(err, '#.4g')), 'f')
                    err = format(Decimal(format(d_background_rateErrs[sampleName][selection], '#.4g')), 'f')
                    thisBkgRawEvts = d_background_unscaledRates[sampleName][selection]
                    # val = format(Decimal(format(value, '#.4g')), 'f')
                    val = format(Decimal(format(d_background_rates[sampleName][selection], '#.4g')), 'f')
                    # yieldRow.append("{:.4f} +/- {:.4f} [MC evts: {:.1f}]".format(value, err, thisBkgRawEvts))
                    yieldRow.append("{:.2f} +/- {:.2f}".format(value, float(err)))
                    mcRow.append("{:.1f}".format(thisBkgRawEvts))
                    # if "TTbar" in sampleName:
                    #     # print("systNominals[{}][{}][{}] = {}; d_background_rates[{}][{}] = {}".format(sampleName, syst, selection, value, sampleName, selection, val))
                    #     print("systNominals[{}][{}][{}] = {}; d_background_rates[{}][{}] = {}".format(sampleName, syst, selection, value, sampleName, selection, d_background_rates[sampleName][selection]))
            # if "Z" in sampleName:
            #     if "EES" in syst:
            #         print "INFO for sampleName={} syst={}".format(sampleName, syst)
            #         print "massList:", massList
            #         print "nominals:", nominals
            #         print "downVariation:", downVariation
            #         print "upVariation:", upVariation
            for idx, systHist in enumerate(systHists):
                systHist.GetXaxis().SetTitle("M_{LQ} [GeV]")
                systHist.GetYaxis().SetTitle(syst+" Uncertainty [%]")
                systHist.SetLineWidth(2)
                systHist.SetBarWidth(0.08)
                systHist.Write()
                # print("Sample Name: {}; systHist name: {}; idx={}; nStacks={}; nSystHists={}".format(sampleName, systHist.GetName(), idx, len(systStacks), len(systHists)))
                systStacks[idx].Add(systHist, "hist")
            table.append(tableRow)
            systDeltaTable.append(systYieldRow)
        yieldTable.append(yieldRow)
        yieldTable.append(mcRow)
        print("Systematics tables for sample Name: {}".format(sampleName))
        print(tabulate(table, headers=tableHeaders, tablefmt="fancy_grid", floatfmt=".2f"))
        print(tabulate(systDeltaTable, headers=systDeltaTableHeaders, tablefmt="fancy_grid", floatfmt=".2f"))
        if includeYieldLine:
            print(tabulate(yieldTable, headers=yieldTableHeaders, tablefmt="fancy_grid", floatfmt=".4f"))
        print()

        if len(totalBkgSystsHeaders) <= 0:
            totalBkgSystsHeaders = tableHeaders
        systDir.cd()
        for idx, systStack in enumerate(systStacks):
            canvas = r.TCanvas("c", "c", 1200, 600)
            canvas.SetName("allSystsCanvas{}_".format(idx+1)+sampleName)
            allSystCanvasNames.append(canvas.GetName())
            canvas.cd()
            canvas.SetGridy()
            canvas.Draw()
            r.gStyle.SetPaintTextFormat(".0f%")
            systStack.Paint()
            systStack.Draw("pfc plc nostackb text")
            systStack.GetHistogram().GetXaxis().SetTitle("M_{LQ} [GeV]")
            systStack.GetHistogram().GetYaxis().SetTitle("Max. Syst. Uncertainty [%]")
            systStack.SetMaximum(100)
            systStack.SetTitle(sampleName)
            canvas.BuildLegend(0.15, 0.55, 0.5, 0.85, "", "l")
            canvas.Write()
            # this broke for some reason. it just hangs.
            # if len(systStacks) > 1:
            #     # canvas.SaveAs(plotsDir+"/systematics_massRange{}_{}.pdf".format(idx+1, sampleName))
            #     # canvas.SaveAs(plotsDir+"/systematics_massRange{}_{}.png".format(idx+1, sampleName))
            #     # print("\tDEBUG: canvas.SaveAs({})".format(plotsDir+"/systematics_massRange{}_{}".format(idx+1, sampleName)), flush=True)
            #     # canvas.SaveAs(plotsDir+"/systematics_massRange{}_{}.png".format(idx+1, sampleName))
            # else:
            #     # print("\tDEBUG: canvas.SaveAs({})".format(plotsDir+"/systematics_{}.pdf".format(sampleName)), flush=True)
            #     # canvas.SaveAs(plotsDir+"/systematics_{}.pdf".format(sampleName))
            systStack.Write()

        # matplotlib bar chart for saving to png, since the root version above broke
        # split the systematics into two groups, as otherwise there are too many bars
        half = len(systematicsByNameAndMass) // 2
        systsGroup1 = {k: v for i, (k, v) in enumerate(systematicsByNameAndMass.items()) if i < half}
        systsGroup2 = {k: v for i, (k, v) in enumerate(systematicsByNameAndMass.items()) if i >= half}
        systGroups = [systsGroup1, systsGroup2]
        if "PDF" in systematicsByNameAndMass.keys():
            systsGroupPDF = {"PDF": systematicsByNameAndMass["PDF"]}
            systGroups.append(systsGroupPDF)
        width = 6
        multiplier = 0
        for idx, systDict in enumerate(systGroups):
            if len(systDict.keys()) < 1:
                continue
            fig, ax = plt.subplots(layout='constrained')
            fig.set_size_inches(12, 5)
            plt.style.use(hep.style.CMS)
            plt.grid(axis = 'y', linestyle = ':')
            maxSyst = 10  # approx. minimum y-axis range (gets multiplied by 1.1 later)
            for systName, effect in systDict.items():
                offset = width * multiplier
                x = np.array(systematicsSelectionsByNameAndMass[systName])
                rects = ax.bar(x + offset, effect, width, label=systName)
                ax.bar_label(rects, padding=3, fmt=lambda x: '{:.1f}%'.format(x) if x > 0 else '', rotation=90, fontsize=6)
                multiplier += 1
                maxEffect = max(effect)
                if maxEffect > maxSyst:
                    maxSyst = maxEffect
            # Add some text for labels, title and custom x-axis tick labels, etc.
            ax = plt.gca()
            ax.set_ylabel('Max. syst. uncert. [%]')
            ax.set_xlabel('$\\mathrm{M}_{LQ}$ [GeV]')
            ax.set_title(sampleName, fontsize=14)
            ax.set_yticks(np.arange(start=0, stop=110, step=10))
            # ax.set_xlim([300, 600])
            # ax.set_ylim([0, 100])
            ax.set_ylim([0, int(round(1.1*maxSyst))])
            ax.legend(loc='upper left', ncols=2, fontsize=10, framealpha=1, frameon=True, facecolor='white')
            plt.savefig(plotsDir+"/systematics_{}_group{}.png".format(sampleName, idx), dpi=100)
            plt.close()

        massList = []
        nominals = []
        deltaXOverX = []
        for selection in list(d_totalSystDeltaXOverX.keys()):
            if sampleName not in list(d_totalSystDeltaXOverX[selection].keys()):
                continue  # skip other signals besides the one corresponding to this selection
            if selection == "preselection":
                massList.append(0.0)
            elif selectionName != "trainingSelection":
                massList.append(1.0)
            else:
                massList.append(float(selection.replace("LQ", "")))
            nominals.append(0.0)
            deltaXOverX.append(100*d_totalSystDeltaXOverX[selection][sampleName])
        systDir.cd()
        systDir.mkdir("totalDeltaXOverX", "totalDeltaXOverX", True).cd()
        plots_tfile.cd("systematics/totalDeltaXOverX")
        # totalSystGraph = r.TGraphErrors(len(massList), np.array(massList), np.array(nominals), np.array([0.0]*len(massList)), np.array(deltaXOverX))
        # totalSystGraph.SetNameTitle("uncertVsMass_{}".format(sampleName))
        # totalSystGraph.GetXaxis().SetTitle("M_{LQ} [GeV]")
        # totalSystGraph.GetYaxis().SetTitle("Rel. Uncertainty [%]")
        # totalSystGraph.SetFillColor(9)
        # totalSystGraph.SetFillStyle(3003)
        # totalSystGraph.Write()
    plots_tfile.Close()
    print("plots written to: {}".format(plots_filePath))
    print()
    print("Total background")
    totalBkgSystTable = []
    for syst in systematicsNamesBackground:
        if syst == "totalYield":
            continue
        tableRow = [syst]
        for selection in d_totalBkgSysts[syst].keys():
            # total = d_totalBkgSysts["totalYield"][selection]
            total = sum([d_systNominals[sampleName][syst][selection] for sampleName in background_names])
            try:
                tableRow.append(100*math.sqrt(d_totalBkgSysts[syst][selection])/total)
            except ZeroDivisionError as e:
                print("Caught a ZeroDivisionError! total is {}. Examine background yields for syst={}, selection={}".format(total, syst, selection))
                for sampleName in background_names:
                    print("sampleName={}, systNominal (bkgYield) = {}".format(sampleName, d_systNominals[sampleName][syst][selection]))
                raise e
        totalBkgSystTable.append(tableRow)
        if "norm" in syst.lower():
            print("DEBUG: QCDFakes_DATA, syst={}, d_systNominals[{}][{}], total sum = {}".format(
                syst, "QCDFakes_DATA", syst, d_systNominals["QCDFakes_DATA"][syst], total))
    print(tabulate(totalBkgSystTable, headers=totalBkgSystsHeaders, tablefmt="fancy_grid", floatfmt=".2f"))

# make final selection tables
if doEEJJ:
    columnNames = [
        "MLQ",
        "signal"]
    columnNames.extend(backgroundTitles)
    columnNames.extend([
        "Total BG (stat) (syst)",
        "Data",
    ])
else:
    columnNames = [
        "MLQ",
        "signal",
        "W+jets",
        "ttbar"] + (["QCD(data)"] if doQCD else []) + ["DIBOSON",
        "SingleTop",
        "Z+Jets",
        "PhotonJets",
        "Total BG",
        "Data",
    ]
latexRowsAN = []
latexRowsPaper = []
t = PrettyTable(columnNames)
t.float_format = "4.3"
# selectionNames = ["preselection"]
# selectionNames.extend(["LQ"+str(mass) for mass in mass_points])
for i_signal_name, signal_name in enumerate(signal_names):
    for i_sel_name, selectionName in enumerate(selectionNames):
        # signal events
        thisSigEvts = "-"
        thisSigEvtsErrUp = "-"
        thisSigEvtsErrDown = "-"
        if selectionName != "preselection" and selectionName != "trainingSelection":
            massPoint = selectionName.replace("LQ", "")
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, massPoint)
            thisSigEvts = d_signal_rates[signalNameForFile][selectionName]
            thisSigEvtsErr = d_signal_rateErrs[signalNameForFile][selectionName]
            # thisSigEvtsErrUp, thisSigEvtsErrDown = GetStatErrorsFromDatacard(d_datacardStatErrs[selectionName][signalNameForFile], thisSigEvts)
            thisSigEvtsErrUp = thisSigEvtsErr
            thisSigEvtsErrDown = thisSigEvtsErr
        # print 'INFO: thisSignal=',fullSignalName,'selection=',selectionName
        # print 'd_data_rates[data]['+selectionName+']'
        if blinded:
            thisDataEvents = -1
        else:
            thisDataEvents = d_data_rates["DATA"][selectionName]
        backgroundEvts = {}
        backgroundEvtsErrUp = {}
        backgroundEvtsErrDown = {}
        totalBackground = 0.0
        totalBackgroundErrStatUp = 0.0
        totalBackgroundErrStatDown = 0.0
        totalBackgroundErrSyst = 0.0
        otherBackground = 0.0
        otherBackgroundErrStatUp = 0.0
        otherBackgroundErrStatDown = 0.0
        for i_background_name, background_name in enumerate(background_names):
            thisBkgEvts = d_background_rates[background_name][selectionName]
            thisBkgEvtsErr = d_background_rateErrs[background_name][selectionName]
            # thisBkgEvtsErrUp, thisBkgEvtsErrDown = GetStatErrorsFromDatacard(d_datacardStatErrs[selectionName][background_name], thisBkgEvts)
            # print "GetStatErrorsFromDatacard for selection={}, background={}, thisBkgEvts={} + {} - {}".format(selectionName, background_name, thisBkgEvts, thisBkgEvtsErrUp, thisBkgEvtsErrDown)
            thisBkgEvtsErrUp = thisBkgEvtsErr
            thisBkgEvtsErrDown = thisBkgEvtsErr
            # thisBkgTotalEntries = d_background_unscaledRates[background_name][selectionName]
            thisBkgEffEntries = thisBkgEvts**2/thisBkgEvtsErr**2 if thisBkgEvtsErr != 0 else 0
            thisBkgSyst = 0
            if doSystematics:
                thisBkgSyst = GetBackgroundSystDeltaOverNominal(background_name, selectionName)
            thisBkgSystErr = thisBkgEvts * thisBkgSyst
            totalBackgroundErrSyst += thisBkgSystErr * thisBkgSystErr
            totalBackground += thisBkgEvts
            totalBackgroundErrStatUp += thisBkgEvtsErrUp * thisBkgEvtsErrUp
            totalBackgroundErrStatDown += thisBkgEvtsErrDown * thisBkgEvtsErrDown
            if background_name in otherBackgrounds:
                otherBackground += thisBkgEvts
                otherBackgroundErrStatUp += thisBkgEvtsErrUp * thisBkgEvtsErrUp
                otherBackgroundErrStatDown += thisBkgEvtsErrDown * thisBkgEvtsErrDown
            if thisBkgEvts < 0:
                print("WARNING: Found", thisBkgEvts, "events for", background_name, "; setting to zero")
                thisBkgEvts = 0.0
            backgroundEvts[background_name] = thisBkgEvts
            backgroundEvtsErrUp[background_name] = thisBkgEvtsErrUp
            backgroundEvtsErrDown[background_name] = thisBkgEvtsErrDown
        totalBackgroundErrStatUp = math.sqrt(totalBackgroundErrStatUp)
        totalBackgroundErrStatDown = math.sqrt(totalBackgroundErrStatDown)
        totalBackgroundErrSyst = math.sqrt(totalBackgroundErrSyst)
        otherBackgroundErrStatUp = math.sqrt(otherBackgroundErrStatUp)
        otherBackgroundErrStatDown = math.sqrt(otherBackgroundErrStatDown)

        row = [selectionName]
        row = [
            selectionName.replace("LQ", ""),
            GetTableEntryStr(
                thisSigEvts, thisSigEvtsErrUp, thisSigEvtsErrDown
            ),  # assumes we always have > 0 signal events
        ]
        for background_name in background_names:
            row.append(
                GetTableEntryStr(
                    backgroundEvts[background_name],
                    backgroundEvtsErrUp[background_name],
                    backgroundEvtsErrDown[background_name],
                )
            )
        row.append(
                GetTableEntryStr(
                    totalBackground,
                    totalBackgroundErrStatUp,
                    totalBackgroundErrStatDown,
                    totalBackgroundErrSyst,
                    )
                )
        row.append(GetTableEntryStr(thisDataEvents))
        t.add_row(row)
        # latex tables
        latexRow = [
            selectionName.replace("LQ", ""),
            GetTableEntryStr(
                thisSigEvts, thisSigEvtsErrUp, thisSigEvtsErrDown, latex=True
            ),  # assumes we always have > 0 signal events
        ]
        latexRowPaper = list(latexRow)
        for background_name in background_names:
            thisEntry = GetTableEntryStr(
                    backgroundEvts[background_name],
                    backgroundEvtsErrUp[background_name],
                    backgroundEvtsErrDown[background_name],
                    latex=True,
                    )
            if background_name not in otherBackgrounds:
                latexRowPaper.append(thisEntry)
            latexRow.append(thisEntry)
        latexRowPaper.append(
                GetTableEntryStr(
                    otherBackground,
                    otherBackgroundErrStatUp,
                    otherBackgroundErrStatDown,
                    latex=True,
                    )
                )
        totalBkgEntry = GetTableEntryStr(
                totalBackground,
                totalBackgroundErrStatUp,
                totalBackgroundErrStatDown,
                totalBackgroundErrSyst,
                True,
                )
        dataEntry = GetTableEntryStr(thisDataEvents, latex=True)
        latexRow.append(totalBkgEntry)
        latexRow.append(dataEntry)
        latexRowPaper.append(totalBkgEntry)
        latexRowPaper.append(dataEntry)

        latexRow = [
            "$" + entry + "$" if "LQ" not in entry and "pres" not in entry and "train" not in entry else entry
            for entry in latexRow
        ]
        latexRowPaper = [
            "$" + entry + "$" if "LQ" not in entry and "pres" not in entry and "train" not in entry else entry
            for entry in latexRowPaper
        ]
        for i, rowEntry in enumerate(latexRow):
            if i < len(latexRow) - 1:
                # rowEntry+=' & '
                latexRow[i] += " & "
            else:
                # rowEntry+=' \\\\ '
                latexRow[i] += " \\\\ "
        latexRowsAN.append("".join(latexRow))
        #
        for i, rowEntry in enumerate(latexRowPaper):
            if i < len(latexRowPaper) - 1:
                latexRowPaper[i] += " & "
            else:
                latexRowPaper[i] += " \\\\ "
        latexRowsPaper.append("".join(latexRowPaper))
        if selectionName == "preselection" or selectionName == "trainingSelection":
            latexRowsAN.append("\\hline")
            latexRowsPaper.append("\\hline")
print(t)
# table_txt = t.get_string()
# table_file.write(table_txt+"\n\n")
 
with open(tablesDir + "/eventYieldsAN.tex", "w") as table_file:
    print()
    print("Latex table: AN")
    print()
    # latex table -- AN
    prelims = [r"\setlength\tabcolsep{2pt}"]
    prelims.append(r"\resizebox{\textwidth}{!}{")
    prelims.append(r"\begin{tabular}{| l | c | c | c | c | c | c | c | c | c | c | c | c | c |}")
    prelims.append(r"\hline")
    for line in prelims:
        print(line)
        table_file.write(line+"\n")
    headers = GetLatexHeaderFromColumnNames(columnNames)
    print(headers)
    table_file.write(headers+"\n")
    for line in latexRowsAN:
        print(line)
        table_file.write(line+"\n")
    ending = r"\end{tabular}}"  # extra } to end resizebox
    print(ending)
    table_file.write(ending+"\n")
    print()
    # table_file.write("\n")

with open(tablesDir + "/eventYieldsPaper.tex", "w") as table_file:
    print()
    print("Latex table: Paper")
    print()
    table_file.write("Latex table: Paper\n")
    # latex table -- Paper
    for line in latexRowsPaper:
        print (line)
        table_file.write(line+"\n")
    print()
    # table_file.write("\n")

# acc * eff
plotSignalEfficiencyTimesAcceptance(dataMC_filepath, signalNameTemplate, [int(m) for m in mass_points], year)
print("datacard written to {}".format(datacard_filePath))
