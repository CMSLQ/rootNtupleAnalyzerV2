#!/usr/bin/env python3

import os
import sys
import math
import string
import re
import copy
import json
from prettytable import PrettyTable
from tabulate import tabulate
from pathlib import Path
from decimal import Decimal
from termcolor import colored
from bisect import bisect
from collections import OrderedDict
import numpy as np
import ROOT as r
import matplotlib.pyplot as plt
import mplhep as hep
import combineCommon as cc

sys.path.append(os.getenv("LQMACRO").rstrip("/") + "/plotting2016/")
from plotSignalEfficiencyTimesAcceptance import plotSignalEfficiencyTimesAcceptance

r.gROOT.SetBatch()


class SampleInfo:
    def __init__(self, name, rates, rateErrs, unscaledRates, totalEvents, failRates, failRateErrs, unscaledFailRates, systematics):
        self.sampleName = name
        self.rates = rates
        self.rateErrs = rateErrs
        self.unscaledRates = unscaledRates
        self.totalEvents = totalEvents
        self.failRates = failRates
        self.failRateErrs = failRateErrs
        self.unscaledFailRates = unscaledFailRates
        self.systematics = systematics
        self.systematicsApplied = set()  # we assume the same systematics for every selection
        # self.systNominals = {}
        # self.systNominalErrs = {}
        # self.systUpDeltas = {}
        # self.systDownDeltas = {}
        # self.minRawEvents = 3  # min number of events necessary to evaluate systematics
        self.minRawEvents = 20  # min number of events necessary to evaluate systematics

    def __add__(self, other):
        new = SampleInfo(self.sampleName, self.rates, self.rateErrs, self.unscaledRates, self.totalEvents, self.failRates, self.failRateErrs, self.unscaledFailRates, self.systematics)
        new += other
        return new

    def __iadd__(self, other):
        for selection in self.rates.keys():
            self.rates[selection] += other.rates[selection]
            self.rateErrs[selection] = math.sqrt(self.rateErrs[selection]**2 + other.rateErrs[selection]**2)
            self.unscaledRates[selection] += other.unscaledRates[selection]
            self.failRates[selection] += other.failRates[selection]
            self.failRateErrs[selection] = math.sqrt(self.failRateErrs[selection]**2 + other.failRateErrs[selection]**2)
            self.unscaledFailRates[selection] += other.unscaledFailRates[selection]
        self.systematics = self.CheckAndAddSystematics(other)
        for syst in other.systematics.keys():
            self.systematicsApplied.add(syst)
        self.totalEvents += other.totalEvents
        return self

    def CheckAndAddSystematics(self, other):
        systDictSelf = self.systematics
        systDictOther = other.systematics
        systDictSelfOrig = copy.deepcopy(systDictSelf)
        # if sorted(systDictSelf.keys()) != sorted(systDictOther.keys()):
        #     systsSample1 = [x for x in sorted(systDictSelf.keys()) if x not in sorted(systDictOther.keys())]
        #     systsSample2 = [x for x in sorted(systDictOther.keys()) if x not in sorted(systDictSelf.keys())]
        #     raise RuntimeError("systematics dicts for sample={} and sample={} have different systematics: \n'{}'\npresent in sample1 and not in 2 while \n'{}'\npresent in sample2 and not in 1. Cannot add the systematics.".format(
        #         self.sampleName, other.sampleName, systsSample1, systsSample2))
        for syst in set(systDictSelf.keys()) & set(systDictOther.keys()):
            for selection in systDictSelf[syst].keys():
                if selection not in list(systDictOther[syst].keys()):
                    raise RuntimeError("cannot find selection '{}' in sample={}; systematics dicts for sample={} and sample={} have different selections for syst={}: '{}' vs. '{}'. Cannot add them.".format(
                        selection, other.sampleName, self.sampleName, other.sampleName, syst, list(systDictSelf[syst].keys()), list(systDictOther[syst].keys())))
                if "branchTitles" not in selection:
                    # if "nominal" in syst.lower() and "selection" in selection.lower():
                    # if "JES" in syst and ("LQ2700" in selection or "LQ2800" in selection):
                    #     print("SIC DEBUG: CheckAndAddSystematics() - sample={}, other={}, systDictSelf[{}][{}]={}".format(self.sampleName, other.sampleName, syst, selection, systDictSelf[syst][selection]))
                    #     print("SIC DEBUG: CheckAndAddSystematics() - sample={}, other={}, systDictOther[{}][{}]={}".format(self.sampleName, other.sampleName, syst, selection, systDictOther[syst][selection]))
                    systDictSelf[syst][selection]["yield"] += systDictOther[syst][selection]["yield"]
                    if systDictSelf[syst][selection]["preselYield"] is None and systDictOther[syst][selection]["preselYield"] is not None:
                        systDictSelf[syst][selection]["preselYield"] = systDictOther[syst][selection]["preselYield"]
                    elif systDictSelf[syst][selection]["preselYield"] is not None and systDictOther[syst][selection]["preselYield"] is not None:
                        systDictSelf[syst][selection]["preselYield"] += systDictOther[syst][selection]["preselYield"]
                        # if "lumicorrelated" in syst.lower() and selection == "preselection":
                        #     nominal = systDictSelf["nominal"][selection]["yield"]
                        #     print("SIC DEBUG: CheckAndAddSystematics() - AFTER ADDING - sample={}, other={}, systDictSelf[{}][{}]={}, entry/nominal - 1 (if nominal != 0)".format(self.sampleName, other.sampleName, syst, selection, systDictSelf[syst][selection], systDictSelf[syst][selection]["yield"]/nominal-1 if nominal != 0 else 0))
                else:
                    # if the branch titles are different, then we just glom them together here
                    # again, for composite samples, this might make no sense
                    if "lhepdf" in syst.lower() and len(systDictOther[syst][selection]) > 1:
                        raise RuntimeError("Not sure how to handle systDictOther for sampleName={} branchTitles={} which have more than 1 entry".format(other.sampleName, systDictOther[syst][selection]))
                    elif len(systDictOther[syst][selection]) == 0:
                        continue
                    if systDictOther[syst][selection][0] not in systDictSelf[syst][selection]:
                        systDictSelf[syst][selection].extend(systDictOther[syst][selection])
                # print("SIC DEBUG: sample={}, systDictSelf[{}][{}]={}".format(self.sampleName, syst, selection, systDictSelf[syst][selection]))
        # handle systs only appearing in sample 2
        systsSample2 = [x for x in sorted(systDictOther.keys()) if x not in systDictSelf.keys()]
        for syst in systsSample2:
            systDictSelf[syst] = systDictOther[syst]
            # need to add nominal of self onto this
            for sel in systDictSelf[syst].keys():
                if "branchTitles" not in sel:
                    systDictSelf[syst][sel]["yield"] += systDictSelfOrig["nominal"][sel]["yield"]
        # handle systs only appearing in sample 1
        systsSample1 = [x for x in sorted(systDictSelf.keys()) if x not in systDictOther.keys()]
        for syst in systsSample1:
            # need to add nominal of other onto this
            for sel in systDictSelf[syst].keys():
                if "branchTitles" not in sel:
                    systDictSelf[syst][sel]["yield"] += systDictOther["nominal"][sel]["yield"]
        return systDictSelf

    def GetRateAndErr(self, selectionName):
        return self.rates[selectionName], self.rateErrs[selectionName]

    def GetFailRateAndErr(self, selectionName):
        return self.failRates[selectionName], self.failRateErrs[selectionName]

    def GetUnscaledRate(self, selectionName):
        return self.unscaledRates[selectionName]

    def UpdateSystsApplied(self, syst):
        self.systematicsApplied.add(syst)

    def GetNearestPositiveSelectionYield(self, selection, minRawEvents=-1):
        if minRawEvents == -1:
            minRawEvents = self.minRawEvents
        massPointsRev = list(reversed(mass_points))
        idxInc = mass_points.index(selection.replace("LQ", ""))
        idxInc += 1
        idxDec = massPointsRev.index(selection.replace("LQ", ""))
        idxDec += 1
        rate = self.rates[selection]
        err = self.rateErrs[selection]
        rawEvents = self.unscaledRates[selection]
        selectionsToCheck = []
        trimmedMassesInc = mass_points[idxInc:]
        trimmedMassesInc = ["LQ"+mass for mass in trimmedMassesInc]
        trimmedMassesDec = massPointsRev[idxDec:]
        trimmedMassesDec = ["LQ"+mass for mass in trimmedMassesDec] + ["trainingSelection", "preselection"]
        # look around given selection
        lastSelectionsChecked = []
        index = 0
        while (rate <= 0 or err <= 0 or rawEvents < minRawEvents) and (index < len(trimmedMassesInc) or index < len(trimmedMassesDec)):
            # lastSelectionName = "LQ" + massPointsRev[idx]
            # lastSelectionName = selectionsToCheck[index]
            if index < len(trimmedMassesInc):
                lastSelectionName = trimmedMassesInc[index]
                rate = self.rates[lastSelectionName]
                err = self.rateErrs[lastSelectionName]
                rawEvents = self.unscaledRates[lastSelectionName]
                lastSelectionsChecked.append(lastSelectionName)
                # print("INFO: GetNearestPositiveSelectionYieldsFromDicts: for sample {}, with initial selection={}, check selectionName={}: rate = {} +/- {}".format(sampleName, selection, lastSelectionName, rate, err))
            if rate <= 0 or err <= 0 or rawEvents < minRawEvents:
                if index < len(trimmedMassesDec):
                    lastSelectionName = trimmedMassesDec[index]
                    rate = self.rates[lastSelectionName]
                    err = self.rateErrs[lastSelectionName]
                    rawEvents = self.unscaledRates[lastSelectionName]
                    lastSelectionsChecked.append(lastSelectionName)
                    # print("INFO: GetNearestPositiveSelectionYieldsFromDicts: for sample {}, with initial selection={}, check selectionName={}: rate = {} +/- {}".format(sampleName, selection, lastSelectionName, rate, err))
            index += 1
        if rate <= 0 or err <= 0 or rawEvents < minRawEvents:
            raise RuntimeError("Could not find valid selection for sample={}; rates look like {} and errors look like {}; checked selections={}".format(
                self.sampleName, self.rates, self.rateErrs, lastSelectionsChecked))
        # print("INFO: GetNearestPositiveSelectionYield: for sample {}, with zero nominal rate for selection={}, found lastSelectionName={} with rate = {} +/- {} [{} evts]".format(self.sampleName, selection, lastSelectionName, rate, err, rawEvents))
        print("INFO: GetNearestPositiveSelectionYield: for sample {}, with rate {} +/- {} [ {} evts] for selection={}, found lastSelectionName={} with rate = {} +/- {} [{} evts]".format(self.sampleName, self.GetRateAndErr(selection)[0], self.GetRateAndErr(selection)[1], self.GetUnscaledRate(selection), selection, lastSelectionName, rate, err, rawEvents))
        # print "INFO: GetNearestPositiveSelectionYieldsFromDicts: found last selectionName={} with {} unscaled events".format(selectionName, unscaledRate)
        return lastSelectionName, rate, err, rawEvents

    # return abs val of syst; also checks for deltaOverNom==1 and does renorm if needed
    def GetSystematicEffectAbs(self, year, systName, selection, applicableSysts, verbose=False):
        # verbose = True
        entry, deltaOverNomUp, deltaOverNomDown, symmetric, systNominal, systSelection = self.GetSystematicEffect(year, systName, selection, applicableSysts[year])
        if verbose:
            print("INFO GetSystematicEffectAbs(): For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}, systNominal={}, systSelection={}".format(
                self.sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown, systNominal, systSelection))
        if entry == "-":
            # this means the systematic doesn't apply
            # print("INFO: GetSystematicEffectAbs() systematic doesn't apply: For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}, systNominal={}, systSelection={}".format(
            #     sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown, systNominal, systSelection))
            return "-", 0, 0, systNominal, systSelection, None
        nomYield = systNominal
        systNomSelection = systSelection
        # returned deltaOverNom is always systYield - nominal
        # if "norm" not in systName.lower() and "shape" not in systName.lower() and "lumi" not in systName.lower():
        preselRatioSystUp = None
        preselRatioSystDown = None
        if "norm" not in systName.lower() and selection != "preselection" and selection != "trainingSelection":
            if self.sampleName in backgroundsToRenormSystAtPresel:
                if verbose:
                    print("\tINFO: renormalizing syst={} for background={}".format(systName, self.sampleName))
                    print("\t original entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(entry, deltaOverNomUp, deltaOverNomDown))
                preselNomYield, preselNomSelection, preselRawEvents = self.GetSystNominalYield(systName, "preselection")
                preselEntry, preselDOverNUp, preselDOverNDown, preselSymm, preselSystNominal, preselSystSelection = self.GetSystematicEffect(year, systName, "preselection", applicableSysts[year])
                if preselSystSelection != preselNomSelection:
                    raise RuntimeError("Something strange happened: the selection used for the preselection systematic '{}' was not the same as the presel syst nominal yield selection '{}'.".format(
                        preselSystSelection, preselNomSelection))
                if not math.isclose(preselSystNominal, preselNomYield, rel_tol=1e-7):
                    raise RuntimeError("Something strange happened: the selection used for the preselection systematic yield '{}' was not the same as the presel syst nominal yield '{}'.".format(
                        preselSystNominal, preselNomYield))
                preselSystYieldUp = GetSystYield(preselDOverNUp, preselNomYield)
                hasPreselSystYield, preselSystYield = self.HasPreselectionSystYield(systName, selection, True)
                if hasPreselSystYield:
                    if verbose:
                        print("\tOverriding preselSystYieldUp for sampleName={} systName={} to {}".format(self.sampleName, systName, preselSystYield))
                    preselSystYieldUp = preselSystYield
                nomYield, systNomSelection, systRawEvents = self.GetSystNominalYield(systName, selection)
                systYieldUp = GetSystYield(deltaOverNomUp, nomYield)
                preselRatioSystUp = preselNomYield/preselSystYieldUp
                systYieldUpRenorm = systYieldUp * preselRatioSystUp
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
                    hasPreselSystYield, preselSystYield = self.HasPreselectionSystYield(systName, selection, False)
                    if hasPreselSystYield:
                        if verbose:
                            print("\tOverriding preselSystYieldDown up for sampleName={} systName={} to {}".format(self.sampleName, systName, preselSystYield))
                        preselSystYieldDown = preselSystYield
                    if verbose:
                        print("\trenorm systYieldDown to {}*{}/{} = systYieldDown*preselNomYield/preselSystYieldDown".format(
                                systYieldDown, preselNomYield, preselSystYieldDown))
                    preselRatioSystDown = preselNomYield/preselSystYieldDown
                    systYieldDownRenorm = systYieldDown * preselRatioSystDown
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
                    preselRatioSystDown = preselRatioSystUp
                if verbose:
                    print("\t new entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(entry, deltaOverNomUp, deltaOverNomDown))
        if (deltaOverNomUp == 0 and deltaOverNomDown == 0) or math.isclose(float(entry), 1):
            # print("WARN: For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}".format(sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown))
            return "no effect", 0, 0, nomYield, systNomSelection, None
        tolerance = 0.5
        if math.fabs(deltaOverNomUp) > tolerance or math.fabs(deltaOverNomDown) > tolerance:
            if not "qcdfakes_data" in self.sampleName.lower():
                rawEvents = self.unscaledRates[systSelection]
                print("WARN: For sample={}, selection={}, syst={}: deltaOverNomUp or deltaOverNomDown is larger than {}%! entry={}, deltaOverNomUp={}, deltaOverNomDown={}, systNominal={}, systSelection={}, rawEvents={}".format(
                    self.sampleName, selection, systName, 100*tolerance, entry, deltaOverNomUp, deltaOverNomDown, systNominal, systSelection, rawEvents))
        # if "zjet" in self.sampleName.lower():
        #     print("DEBUGINFO GetSystematicEffectAbs(): For sample={}, selection={}, syst={}: entry={}, deltaOverNomUp={}, deltaOverNomDown={}, systNominal={}, systSelection={}, preselRatioSystUp={}, preselRatioSystDown={}".format(
        #         self.sampleName, selection, systName, entry, deltaOverNomUp, deltaOverNomDown, systNominal, systSelection, preselRatioSystUp, preselRatioSystDown))
        return entry, math.fabs(deltaOverNomUp), math.fabs(deltaOverNomDown), nomYield, systNomSelection, [preselRatioSystUp, preselRatioSystDown]

    def GetSystematicEffect(self, year, systName, selection, applicableSystematics):
        # print("GetSystematicEffect(systName={}, selection={}. d_applicableSystematics={})".format(systName, selection, d_applicableSystematics))
        verbose = False
        systDict = self.systematics
        try:
            nominal = systDict["nominal"][selection]["yield"]
        except KeyError as e:
            raise RuntimeError("Could not find key 'nominal' in systDict for sampleName={}; systDict.keys()={}".format(sampleName, systDict.keys()))
        except Exception as e:
            raise RuntimeError("Got exception accessing systDict[{}][{}][{}] for sampleName={}; systDict[{}]={}".format("nominal", selection, "yield", sampleName, "nominal", systDict["nominal"]))
        if not DoesSystematicApply(systName, self.sampleName, applicableSystematics):
            return "-", 0, 0, True, nominal, selection
        if "shape" in systName.lower():
            entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = self.CalculateUpDownSystematic("LHEScaleComb", selection, verbose)
        elif systName == "LHEPdfWeight":
            # print("INFO GetSystematicEffect(): for sample={}, selection={}, systName={}, fullSystDict[{}][LHEPdfCombUp] = {}".format(sampleName, selection, systName, sampleName, fullSystDict[sampleName]["LHEPdfCombUp"][selection]))
            # print("INFO GetSystematicEffect(): for sample={}, selection={}, systName={}, fullSystDict[{}][LHEPdfCombDown] = {}".format(sampleName, selection, systName, sampleName, fullSystDict[sampleName]["LHEPdfCombDown"][selection]))
            entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = self.CalculateUpDownSystematic("LHEPdfComb", selection, verbose)
            # if verbose:
            #     components = dictSamples[sampleName]["pieces"]
            #     for sampleName in components:
            #         if sampleName not in dictSamples.keys():
            #             print("INFO: Cannot print PDF systematic output for sample {} as it was not saved.".format(sampleName))
            #         else:
            #             CalculateUpDownSystematic("LHEPdfComb", fullSystDict[sampleName], selection, sampleName, True)  # this call is broken somehow
        elif "lumi" in systName.lower():
            entry, deltaNomUp, deltaNomDown, symmetric = self.CalculateFlatSystematic(systName, selection)
            newNominal = nominal if nominal > 0 else 1
            # newNominal = nominal
            newSelection = selection
        elif "norm" in systName.lower():
            if "tt" in systName.lower() or "dy" in systName.lower() or "qcd" in systName.lower():
                entry, deltaNomUp, deltaNomDown, symmetric = self.CalculateFlatSystematic(systName, selection)
                newNominal = nominal if nominal > 0 else 1  # XXX SICF FIXME WHY?!?!
                newSelection = selection
        else:
            entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection = self.CalculateUpDownSystematic(systName, selection)
        return entry, deltaNomUp, deltaNomDown, symmetric, newNominal, newSelection

    def GetTotalSystDeltaOverNominal(self, year, selectionName, systematicsNames, applicableSysts, verbose=False):
        totalSyst = 0
        for syst in systematicsNames:
            systEntry, deltaOverNominalUp, deltaOverNominalDown, _, _, _ = self.GetSystematicEffectAbs(year, syst, selectionName, applicableSysts)
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
    
    # def CalculateShiftSystematic(systName, selection):
    #     nominal, selection, rawEvents = GetSystNominalYield(systName, selection)
    #     try:
    #         systYield = systDict[systName][selection]["yield"]
    #     except KeyError as e:
    #         if systName not in systDict.keys():
    #             raise RuntimeError("Could not find syst named {} in systDict keys {}".format(systName, systDict.keys()))
    #         elif selection not in systDict[systName].keys():
    #             raise RuntimeError("Could not find selection {} in systDict[{}].keys() {}".format(selection, systName, systDict[systName].keys()))
    #         else:
    #             raise RuntimeError("Strange KeyError occurred")
    #     delta = systYield - nominal
    #     # print "CalculateShiftSystematic({}, {}, {}): nominal={}, systYield={}, delta={}, delta/nominal={}".format(
    #     #         systName, selection, sampleName, nominal, systYield, delta, delta/nominal)
    #     return str(1 + math.fabs(delta)/nominal), delta/nominal, delta/nominal, True
    
    def CalculateUpDownSystematic(self, systName, selection, verbose=False):
        symmetric = False
        systDict = self.systematics
        nominal = systDict["nominal"][selection]
        # print("CalculateUpDownSystematic({}, {}) for sample: {}; nominal yield={}".format(systName, selection, self.sampleName, nominal))
        nominal, selection, rawEvents = self.GetSystNominalYield(systName, selection)
        # print("CalculateUpDownSystematic(): nominal={}, selection={} after GetNominalYield({})".format(nominal, selection, systName))
        # if "OTHERBKG" in self.sampleName:
        #     verbose = True
        try:
            systYieldUp = systDict[systName+"Up"][selection]["yield"]
            systYieldDown = systDict[systName+"Down"][selection]["yield"]
        except KeyError:
            raise RuntimeError("Could not find Up or Down key for syst={}; keys={}".format(systName, list(systDict.keys())))
        logString = "CalculateUpDownSystematic(): selection={}, sample={} syst={}, systYieldUp={}, systYieldDown={}, systNominal={}, origSystNominal={}".format(
                    selection, self.sampleName, systName, systYieldUp, systYieldDown, nominal, systDict["nominal"][selection])
        kUp = systYieldUp/nominal if nominal != 0 else 0
        kDown = systYieldDown/nominal if nominal != 0 else 0
        upDelta = systYieldUp-nominal if nominal != 0 else 0
        downDelta = systYieldDown-nominal if nominal != 0 else 0
        entry = kDown/kUp if kUp != 0 else 0
        deltaNomUp = upDelta/nominal if nominal != 0 else 0
        deltaNomDown = downDelta/nominal if nominal != 0 else 0
        if nominal == 0 and not IsComponentBackground(self.sampleName):
            raise RuntimeError("CalculateUpDownSystematic(): for sample={}, ZeroDivisionError for a quantity, as nominal={}. kDown/kUp = {}/{}, upDelta/nominal={}/{}, downDelta/nominal={}/{}.".format(self.sampleName, nominal, kDown, kUp, upDelta, nominal, downDelta, nominal)+ " "+logString)
        if kUp <= 0:
            if verbose:
                print(logString)
                print("CalculateUpDownSystematic(): kUp <= 0 = {}; return symmetric kDown syst '{}', {}, {}".format(kUp, str(1 + math.fabs(downDelta)/nominal), downDelta/nominal, downDelta/nominal))
            return str(1 + math.fabs(deltaNomDown)), deltaNomDown, deltaNomDown, True, nominal, selection
        elif kDown <= 0:
            if verbose:
                print(logString)
                print("CalculateUpDownSystematic(): kDown <= 0 = {}; return symmetric kDown syst '{}', {}, {}".format(kDown, str(1 + math.fabs(upDelta)/nominal), upDelta/nominal, upDelta/nominal))
            return str(1 + math.fabs(deltaNomUp)), deltaNomUp, deltaNomUp, True, nominal, selection
        if deltaNomUp == deltaNomDown:
            entry = 1 + math.fabs(upDelta)/nominal
            symmetric = True
        if verbose:
            print(logString)
            print("CalculateUpDownSystematic(): kDown/kUp = {}/{}; return '{}', {}, {}, {}".format(kDown, kUp, str(entry), deltaNomUp, deltaNomDown, symmetric))
        return str(entry), deltaNomUp, deltaNomDown, symmetric, nominal, selection
    
    def CalculateFlatSystematic(self, systName, selection):
        # assumes that the number here is < 1
        systDict = self.systematics
        nominal, selection, rawEvents = self.GetSystNominalYield(systName, selection)
        # print("INFO2: CalculateFlatSystematic() - for sample={}, systName={}, selection={}, syst nominal = {}, systDict[systName][selection]={}".format(self.sampleName, systName, selection, nominal, systDict[systName][selection]))
        deltaOverNom = systDict[systName][selection]["yield"]/nominal - 1 if nominal != 0 else 0
        # if "lumicorrelated" in systName.lower():
        #     print("CalculateFlatSystematic() - for sample={}, systName={}, selection={}, yield={}, nom={}, yield/nom -1 = {}".format(self.sampleName, systName, selection, systDict[systName][selection]["yield"], nominal, deltaOverNom))
        return str(1 + deltaOverNom), deltaOverNom, deltaOverNom, True
    
    def GetSystNominalYield(self, systName, selection):
        verbose = False
        minRawEvents = self.minRawEvents
        requestedSelection = selection
        systDict = self.systematics
        nominal = systDict["nominal"][selection]["yield"]
        rawEvents = self.GetUnscaledRate(selection)
        # if self.sampleName in dictSampleComponents.keys():
        #     allSelections = set()
        #     for componentBkg in dictSampleComponents[self.sampleName]:
        #         compSampleInfo = d_backgroundSampleInfos[componentBkg]
        #         nominal = compSampleInfo.systematics["nominal"][selection]["yield"]
        #         rawEvents = compSampleInfo.GetUnscaledRate(selection)
        #         if nominal <= 0 or rawEvents < minRawEvents:
        #             selections = FindAllSelectionsPassingThresholds(compSampleInfo, selectionNames, minRawEvents)
        #             print("selections='{}'".format(selections))
        #             allSelections.update(selections)
        #         else:
        #             allSelections.update([selection])
        #     if "trainingSelection" in allSelections:
        #         selection = "trainingSelection"
        #     elif "preselection" in allSelections:
        #         selection = "preselection"
        #     else:
        #         print("allSelections='{}'".format(allSelections))
        #         selection = "LQ" + str(max([int(sel[2:]) for sel in list(allSelections)]))
        #     # nominal = systDict["nominal"][selection]["yield"]
        #     nominal, err = self.GetRateAndErr(selection)
        #     rawEvents = self.GetUnscaledRate(selection)
        #     print("GetSystNominalYield({}, {}), for sample={} with components - found best lastSelection: {}, raw events: {}, rate: {} +/- {}".format(
        #         systName, requestedSelection, self.sampleName, selection, rawEvents, nominal, err))
        # elif not IsComponentBackground(self.sampleName):
        #     nominal, selection, rawEvents = self.FindNearestSelectionPassingThresholds(selection, nominal, rawEvents, minRawEvents)

        # if "OTHERBKG" in self.sampleName:
        #     print("GetSystNominalYield({}, {}, {}), before finding nearest valid selection -- selection: {}, raw events: {}, rate: {} +/- {}, systNominal={}".format(
        #         systName, requestedSelection, self.sampleName, selection, rawEvents, self.GetRateAndErr(selection)[0], self.GetRateAndErr(selection)[1], nominal))
        if not IsComponentBackground(self.sampleName):
            nominal, selection, rawEvents = self.FindNearestSelectionPassingThresholds(selection, nominal, rawEvents, minRawEvents)

        #if sampleName in list(d_background_unscaledRates.keys()):
        #    unscaledRate = d_background_unscaledRates[sampleName][selection]
        #    err = d_background_rateErrs[sampleName][selection]
        #else:
        #    unscaledRate = d_signal_unscaledRates[sampleName][selection]
        #    err = d_signal_rateErrs[sampleName][selection]
        systNominal = systDict["nominal"][selection]["yield"]
        # if "OTHERBKG" in self.sampleName:
        #     print("GetSystNominalYield({}, {}, {}), lastSelection: {}, raw events: {}, rate: {} +/- {}; systNominal={}".format(
        #         systName, requestedSelection, self.sampleName, selection, rawEvents, nominal, self.GetRateAndErr(selection)[1], systNominal))
        return systNominal, selection, rawEvents

    def FindNearestSelectionPassingThresholds(self, selection, nominal, rawEvents, minRawEvents):
        if nominal <= 0 or rawEvents < minRawEvents:
            # take the last selection for which we have some rate, and use that for systematic evaluation
            lastNonzeroSelection, rate, err, events = self.GetNearestPositiveSelectionYield(selection, minRawEvents)
            nominal = rate
            selection = lastNonzeroSelection
            rawEvents = events
        # print("SICDEBUG: FindNearestSelectionPassingThresholds(): return nominal={}, selection={}, rawEvents={}".format(nominal, selection, rawEvents))
        return nominal, selection, rawEvents

    def HasPreselectionSystYield(self, systName, selection, up=None, verbose=False):
        if "shape" in systName.lower():
            systName = "LHEScaleComb"
        elif "LHEPdfWeight" == systName:
            systName = "LHEPdfComb"
        systName += "Up" if up else "Down" if not up else ""
        preselSystYield = None
        systDict = self.systematics
        if systName in systDict.keys():
            if verbose:
                print("INFO: sample {} hasPreselSystYield for systName={}, selection={}, up={}; systDict[{}]={}".format(self.sampleName, systName, selection, up, systName, systDict[systName]))
            preselSystYield = systDict[systName][selection]["preselYield"]
        elif verbose:
            print("INFO: sample {} hasPreselSystYield [systName not in systDict keys] for systName={}, selection={}, up={}; systDict.keys={}".format(self.sampleName, systName, selection, up, systName, systDict.keys()))
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

    def GetLastNonzeroSelectionYields(self):
        massPointsRev = list(reversed(mass_points))
        idx = 0
        rate = 0
        err = 0
        while rate == 0 or err == 0:
            lastSelectionName = "LQ" + massPointsRev[idx]
            rate, err = self.GetRateAndErr(lastSelectionName)
            rawEvents = self.GetUnscaledRate(lastSelectionName)
            idx += 1
            if idx >= len(massPointsRev):
                print("WARN GetLastNonzeroSelectionYieldsFromDicts: could not find final selection point with nonzero yield for this background; returning trainingSelection")
                return "trainingSelection", self.GetRateAndErr("trainingSelection"), self.GetUnscaledRate("trainingSelection")
        return lastSelectionName, rate, err, rawEvents
    

def FindAllSelectionsPassingThresholds(sampleInfo, selectionsToCheck, minRawEvents):
    systDict = sampleInfo.systematics
    selections = []
    for selection in selectionsToCheck:
        nominal = systDict["nominal"][selection]["yield"]
        rawEvents = sampleInfo.GetUnscaledRate(selection)
        if nominal > 0 and rawEvents >= minRawEvents:
            selections.append(selection)
    return selections


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


def GetSingleTopYieldAndUncertainty(year, mass):
    func = r.TF1("func", "expo(0)", 0, 5000)
    constTerm = 3.90306
    constTermErr = 1.25629
    slopeTerm = -0.00491907
    slopeTermErr = 0.00170827
    # simple lumi scaling from full Run II fit
    if year == "2016preVFP":
        scaleFactor = 0.1416840016416827
    elif year == "2016postVFP":
        scaleFactor = 0.1221676839055999
    elif year == "2017":
        scaleFactor = 0.3014043829098545
    elif year == "2018":
        scaleFactor = 0.4347439315428627
    elif year == "all":
        scaleFactor = 1.0
    else:
        raise RuntimeError("Couldn't understand year={}, so couldn't scale the fit function.".format(year))
    expConstant = math.log(scaleFactor*math.exp(constTerm))
    func.SetParameters(expConstant, slopeTerm)
    predictedYield = func.Eval(mass+50)  # 100 being the bin width, evaluate in the middle of the bin
    uncertaintiesOnYield = []
    for i in range(1, 9):
        if i == 1:
            toAddConst = constTermErr
            toAddSlope = 0
        elif i == 2:
            toAddConst = -1*constTermErr
            toAddSlope = 0
        elif i == 3:
            toAddConst = 0
            toAddSlope = slopeTermErr
        elif i == 4:
            toAddConst = 0
            toAddSlope = -1*slopeTermErr
        else:
            break
        # elif i == 5:
        #     toAddConst = constTermErr
        #     toAddSlope = slopeTermErr
        # elif i == 6:
        #     toAddConst = constTermErr
        #     toAddSlope = -1*slopeTermErr
        # elif i == 7:
        #     toAddConst = -1*constTermErr
        #     toAddSlope = slopeTermErr
        # elif i == 8:
        #     toAddConst = -1*constTermErr
        #     toAddSlope = -1*slopeTermErr
        expConstant = math.log(scaleFactor*math.exp(constTerm+toAddConst))
        newFunc = func.Clone("newFunc"+str(i))
        newFunc.SetParameters(expConstant, slopeTerm+toAddSlope)
        uncertaintiesOnYield.append(math.fabs(predictedYield-newFunc.Eval(mass+50)))
    maxUncertainty = max(uncertaintiesOnYield)
    maxUncertainty = min(maxUncertainty, predictedYield)  # cap uncertainty at 100%
    return predictedYield, maxUncertainty


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
def GetSystematicsDict(rootFile, sampleName, selections, year, verbose=False):
    systHistName = "histo2D__{}__systematics".format(sampleName)
    systHist = rootFile.Get(systHistName)
    systDict = {}
    systDict["branchTitles"] = {}
    xBinPresel = systHist.GetXaxis().FindFixBin("preselection")
    for selection in selections:
        # expect that selections are either "preselection" or "LQXXXX"
        finalSelection = cc.GetFinalSelection(selection)
        systDict[selection] = {}
        xBin = systHist.GetXaxis().FindFixBin(finalSelection)
        if xBin < 1:
            raise RuntimeError("Could not find requested selection name '{}' in hist {} in file {}".format(finalSelection, systHistName, rootFile.GetName()))
        lheScaleYBin = None
        for yBin in range(1, systHist.GetNbinsY()+1):
            systName = systHist.GetYaxis().GetBinLabel(yBin)
            preselYield = None
            # if "LHEPdfWeight" in systName or "LHEScaleWeight" in systName:
            #     continue
            # if "LHEPdf" in systName or "LHEScale" in systName:
            #     systName = systName.replace("_UpComb", "CombUp").replace("_DownComb", "CombDown")
            #     if "LHEScale" in systName:
            #         preselYieldBin = systHist.GetYaxis().FindFixBin("LHEScaleWeight_preselYield")
            #         preselYield = systHist.GetBinContent(xBin, preselYieldBin)
                    # print("DEBUG: for systName={}, selection={}, sampleName={}, with xBinPresel={}, preselYieldBin={}; preselYield={}".format(
                    #     systName, selection, sampleName, xBinPresel, preselYieldBin, preselYield))
            if "LHEPdf" in systName or "LHEScale" in systName:
                if "comb" in systName.lower():
                    systName = systName.replace("_UpComb", "CombUp").replace("_DownComb", "CombDown")
                if "LHEScale" in systName:
                    preselYieldBin = systHist.GetYaxis().FindFixBin("LHEScaleWeight_preselYield")
                    preselYield = systHist.GetBinContent(xBin, preselYieldBin)
            systDict[selection][systName] = {}
            systDict[selection][systName]["yield"] = systHist.GetBinContent(xBin, yBin)
            systDict[selection][systName]["preselYield"] = preselYield
            # print("DEBUG: ASSIGNED for systName={}, selection={}, preselYield={}".format(systName, selection, preselYield))
            # if "ZJet" in sampleName and "LQ3000" in selection:
            #     print("INFO: for sampleName={} and selection={}, xBin={} yBin={} ({}), got content={}".format(sampleName, selection, xBin, yBin, systName, systDict[selection][systName]))
        # add flat systematics here
        for systName, effect in d_flatSystematics[year].items():
            systDict[selection][systName] = {}
            systDict[selection][systName]["yield"] = (effect+1) * systDict[selection]["nominal"]["yield"]  # effect being deltaX/X = (X'-X)/X = X'/X - 1
            # print("INFO: flat syst: year={}, systName={}, effect={}, nominal={}, x'={}".format(year, systName, effect, systDict[selection]["nominal"]["yield"], (effect+1) * systDict[selection]["nominal"]["yield"]))
            systDict[selection][systName]["preselYield"] = None
            systDict["branchTitles"][systName] = []
    # add special entry for branch titles
    tmapName = "tmap__{}__systematicNameToBranchesMap".format(sampleName)
    tmap = rootFile.Get(tmapName)
    if not tmap or tmap is None:
        raise RuntimeError("Could not find TMap '{}' in file {}".format(tmapName, rootFile.GetName()))
    for yBin in range(1, systHist.GetNbinsY()+1):
        branchTitleList = []
        systName = systHist.GetYaxis().GetBinLabel(yBin)
        # for branch titles, because of how we combine the samples, the TMap is just taken from the first file in the first sample
        #   so it's not necessarily correct for composite samples!
        systNameForLookup = systName
        if "lhe" in systName.lower():
            if not "comb" in systName.lower():
                if "LHEPdfWeight" in systName:
                    systNameForLookup = "LHEPdfWeight"
                elif "LHEScaleWeight" in systName:
                    systNameForLookup = "LHEScaleWeight"
            else:
                systName = systName.replace("_UpComb", "CombUp").replace("_DownComb", "CombDown")
                systDict["branchTitles"][systName] = [systName]
                continue
        mapObject = tmap.FindObject(systNameForLookup)
        if not mapObject:
            # print("\tINFO: look again for matching TObject for syst {}".format(systName[:systName.rfind("_")]), flush=True)
            # assume it's an array syst, so try to match stripping off the _N part
            mapObject = tmap.FindObject(systNameForLookup[:systNameForLookup.rfind("_")])
        if not mapObject:
            tmap.Print()
            raise RuntimeError("Could not find matching TObject in map for syst {} using lookup {} or {}; see map content above".format(systName, systNameForLookup, systNameForLookup[:systNameForLookup.rfind("_")]))
        branchTitleListItr = r.TIter(mapObject.Value())
        branchTitle = branchTitleListItr.Next()
        while branchTitle:
            branchTitleList.append(branchTitle.GetName())
            branchTitle = branchTitleListItr.Next()
        systDict["branchTitles"][systName] = branchTitleList
        # print("INFO: sampleName={}, systDict[\"branchTitles\"][{}] = {}".format(sampleName, systName, branchTitleList))

    #if verbose:
    # print("sampleName={}: systDict['LQ300']={}".format(sampleName, systDict["LQ300"]))

    # print("DEBUG: sampleName={}: systDict.keys()={}".format(sampleName, systDict.keys()))
    # print("DEBUG: sampleName={}: systDict[systDict.keys()[0]].keys()={}".format(sampleName, systDict[list(systDict.keys())[0]].keys()))
    # print("DEBUG: sampleName={}: list(systDict[selections[0]].keys())={}".format(sampleName, list(systDict[selections[0]].keys())))
    # for syst in list(systDict[selections[0]].keys()):
    #     print("DEBUG: syst={}".format(syst))
    #     for sel in systDict:
    #         print("\tDEBUG: syst={}, sel={}, systDict[sel].keys()={}".format(syst, sel, systDict[sel].keys()))
    #         print("\t\tDEBUG2:  systDict[sel][syst]={}".format(systDict[sel][syst]))

        
    
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


def DoesSystematicApply(systName, sampleName, applicableSysts):
    # if "LQ" in sampleName:
    #     print("DoesSystematicApply? {} for systName={}, sampleName={}, applicableSysts={}".format(systName in applicableSysts[sampleName], systName, sampleName, applicableSysts))
    try:
        return systName in applicableSysts[sampleName]
    except KeyError as e:
        if not sampleName in applicableSysts.keys():
            raise RuntimeError("Could not find sampleName '{}' in applicable systs dict; keys = {}".format(sampleName, applicableSysts.keys()))
        raise RuntimeError("Could not find systName '{}' in keys for sampleName = {}; keys = {}".format(systName, sampleName, applicableSysts[sampleName].keys()))


def DoesSystematicApplyAnyYear(systName, sampleName, d_applicableSystematicsAllYears):
    for year, applicableSysts in d_applicableSystematicsAllYears.items():
        if DoesSystematicApply(systName, sampleName, applicableSysts):
            return True
    return False


def IsComponentBackground(sampleName):
    if sampleName in dictSamplesContainingSample.keys():
        return True
    return False


def IsBackground(sampleName):
    if sampleName in background_names:
        return True
    return False


def GetSampleNameFromSubstring(substring, background_names):
    sampleName = [sample for sample in background_names if substring in sample]
    if len(sampleName) != 1:
        raise RuntimeError("Could not find unique {} sample name in background_names with the given keys: {}".format(substring, background_names))
    return sampleName[0]


def RoundToN(x, n, returnFloat=False):
    # if n < 1:
    #    raise ValueError("can't round to less than 1 sig digit!")
    # # number of digits given by n
    # return "%.*e" % (n-1, x)
    if isinstance(x, float):
        rounded = round(x, n)
        # return "{0:.{1}f}".format(rounded, n)
        if returnFloat:
            return rounded
        # return str(rounded)
        return "{0:.{1}f}".format(rounded, n)
    else:
        return str(x)


def GetNDecimalDigits(num):
    d = Decimal(str(num))
    positiveExponent = abs(d.as_tuple().exponent)
    return positiveExponent


def FormatToNDigits(num, n):
    if isinstance(num, str) or isinstance(num, int):
        return str(num)
    rounded = RoundToN(num, n, returnFloat=True)
    return "{0:.{1}f}".format(rounded, n)


def RoundToNSigFigs(num, n=1):
    if isinstance(num, str):
        return num, -1
    nDecDigits = GetNDecimalDigits(num)
    d = Decimal(str(num))
    nDigits = len(d.normalize().as_tuple().digits)
    digitsBeforeDecimal = nDigits - nDecDigits
    d = round(d, n-digitsBeforeDecimal)
    # dStr = str(np.format_float_positional(float(d), trim='-'))
    dStr = str(int(d))
    if d < 10:
        dStr = "{0:.{1}f}".format(float(d), n-digitsBeforeDecimal)
    # print("DEBUG: RoundToNSigFigs for num={} --> {}; n-digitsBeforeDecimal = {}-{} = {}".format(num, dStr, n, digitsBeforeDecimal, n-digitsBeforeDecimal))
    # but now, digits before decimal could have changed, e.g., 0.99 --> 1.0
    if "." in dStr:
        d = Decimal(dStr)
        nDecDigitsUpdated = GetNDecimalDigits(d)
        dUpdated = Decimal(dStr)
        nDigitsUpdated = len(dUpdated.as_tuple().digits) if "." in dStr else len(d.normalize().as_tuple().digits)
        digitsBeforeDecimalUpdated = nDigitsUpdated - nDecDigitsUpdated
        if digitsBeforeDecimalUpdated > digitsBeforeDecimal:
            dStr = str(int(d))
        # print("DEBUG: RoundToNSigFigs for num={} --> {}; [updated] n-digitsBeforeDecimal = {}-{} = {}".format(num, dStr, n, digitsBeforeDecimal, n-digitsBeforeDecimal))
    return dStr, n-digitsBeforeDecimal


def GetTableEntryStr(evts, errStatUp="-", errStatDown="-", errSyst=0, addDecimalsUntilNonzero=False, latex=False):
    if evts == "-":
        return evts
    # print("DEBUG: RoundToNSigFigs for errStatUp={}".format(errStatUp))
    errStatUpR, digitsAwayFromDecimal = RoundToNSigFigs(errStatUp)
    # print("DEBUG: AFTER RoundToNSigFigs for errStatUp={}: got errStatUpR={}, digitsAwayFromDecimal={}".format(errStatUp, errStatUpR, digitsAwayFromDecimal))
    errStatDownR = str(float(round(Decimal(errStatDown), digitsAwayFromDecimal))) if not isinstance(errStatDown, str) else errStatDown
    # print("DEBUG: Now for errStatDown={}: got errStatDownR={}, digitsAwayFromDecimal={}".format(errStatDown, errStatDownR, digitsAwayFromDecimal))
    evtsR = str(float(round(Decimal(evts), digitsAwayFromDecimal)))
    if float(evtsR) > 1 and evtsR.endswith(".0") and len(evtsR) > len(errStatUpR):
        evtsR = evtsR[:-2]
    elif errStatUpR != "-" and len(evtsR) != len(errStatUpR):
        # print("DEBUG: reformat evtsR to {} digits: evtsR={} --> {}".format(GetNDecimalDigits(errStatUpR), evtsR, "{0:.{1}f}".format(float(evtsR), GetNDecimalDigits(errStatUpR))))
        evtsR = "{0:.{1}f}".format(float(evtsR), GetNDecimalDigits(errStatUpR))
    # print("DEBUG: Now for evts={}: got evtsR={}, digitsAwayFromDecimal={}".format(evts, evtsR, digitsAwayFromDecimal))
    # # rounding
    # errStatUpR = RoundToN(errStatUp, 2)
    # if GetNDecimalDigits(errStatUpR) > 1:
    #     errStatDownR = FormatToNDigits(errStatDown, 2)
    #     evtsR = FormatToNDigits(evts, 2)
    #     errStatUpR = FormatToNDigits(errStatUp, 2)
    # else:
    #     errStatDownR = FormatToNDigits(errStatDown, 1)
    #     evtsR = FormatToNDigits(evts, 1)
    #     errStatUpR = FormatToNDigits(errStatUp, 1)

    # suppress
    suppressed = False
    if errStatUpR != "-" and float(errStatUpR) < 0.01:
        errStatUpR = "0.0"
        suppressed = True
    if errStatDownR != "-" and float(errStatDownR) < 0.01:
        errStatDownR = "0.0"
        suppressed = True
    if suppressed:
        evtsR, _ = RoundToNSigFigs(evts, n=2)
        if float(evtsR) < 0.01:
            evtsR, _ = RoundToNSigFigs(evts, n=1)
    # add additional decimal place if it's zero after rounding
    if addDecimalsUntilNonzero and float(evtsR) == 0.0:
        # print("DEBUG: adding decimals until nonzero for evts={}; evtsR={}; errStatUp={}, errStatDown={}".format(evts, evtsR, errStatUp, errStatDown))
        nDigits = 2
        while float(evtsR) == 0.0 and nDigits < 5:
            nDigits += 1
            evtsR = RoundToN(evts, nDigits)
            # print("\tDEBUG: rounded to n={}, evtsR={}".format(nDigits, evtsR))
        # errStatUpR = RoundToN(errStatUp, nDigits)
        # errStatDownR = RoundToN(errStatDown, nDigits)
        # print("\tDEBUG: for evts={}; ended up with evtsR={}; errStatUpR={}, errStatDownR={}".format(evts, evtsR, errStatUpR, errStatDownR))
    # handle cases where we don't specify stat or syst
    if errStatUp == "-":
        return evtsR
    elif errSyst == 0:
        if errStatUp == errStatDown:
            if not latex:
                return evtsR + " +/- " + errStatUpR if float(errStatUpR) != 0 else evtsR
            else:
                return evtsR + " \\pm " + errStatUpR if float(errStatUpR) != 0 else evtsR
        else:
            if not latex:
                return evtsR + " + " + errStatUpR + " - " + errStatDownR
            else:
                return (
                    evtsR
                    + "^{+"
                    + errStatUpR
                    + "}_{-"
                    + errStatDownR
                    + "}"
                )
    else:
        # errSystR = FormatToNDigits(errSyst, 2)
        errSystR = str(float(round(Decimal(errSyst), digitsAwayFromDecimal))) if not isinstance(errSyst, str) else errSyst
        # suppress or remove trailing zero
        if errSystR != "-" and float(errSystR) < 0.01:
            errSystR = "0.0"
        elif float(errSystR) > 1 and errSystR.endswith(".0") and len(errSystR) > len(evtsR):
            errSystR = errSystR[:-2]
        elif errStatUpR != "-" and len(errSystR) != len(errStatUpR):
            # print("DEBUG: reformat errSystR to {} digits: errSystR={} --> {}".format(GetNDecimalDigits(errStatUpR), errSystR, "{0:.{1}f}".format(float(errSystR), GetNDecimalDigits(errStatUpR))))
            errSystR = "{0:.{1}f}".format(float(errSystR), GetNDecimalDigits(errStatUpR))
        if errStatUp == errStatDown:
            if not latex:
                return evtsR + " +/- " + errStatUpR + " +/- " + errSystR
            else:
                return (
                    evtsR + " \\pm " + errStatUpR + " \\pm " + errSystR
                )
        else:
            return (
                evtsR
                + "^{+"
                + errStatUpR
                + "}_{-"
                + errStatDownR
                + "} \\pm "
                + errSystR
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


def SumSampleInfoOverYears(sampleInfos, years):
    yearsDone = []
    for year in years:
        for sample, sampleInfo in sampleInfos[year].items():
            try:
                if not sample in sampleInfos.keys():
                    sampleInfos[sample] = sampleInfo
                else:
                    sampleInfos[sample] += sampleInfo
            except RuntimeError as e:
                raise RuntimeError("Caught RuntimeError when trying to add sample {} from year {} to sample {} from previously-summed years {}: '{}'".format(sampleInfo.sampleName, year, sampleInfos[sample].sampleName, yearsDone, e))
        yearsDone.append(year)
    firstSample = list(sampleInfos.keys())[1]


def CheckHistBins(hist):
    for iBin in range(0, hist.GetNbinsX()+2):
        if hist.GetBinContent(iBin) <= 0:
            hist.SetBinContent(iBin, 1e-10)
    return hist


def CreateAndWriteBackgroundVsMLQHist(rootFiles, sample, bkgEvtsDict, bkgErrsDict):
    newHist = r.TH1D("yieldVsMLQ_{}".format(sample), "yieldVsMLQ for {}".format(sample), len(mass_points_int)-1, min(mass_points_int), max(mass_points_int))
    for selectionName in bkgEvtsDict.keys():
        mass_point = int(selectionName.replace("LQ", ""))
        content = bkgEvtsDict[selectionName]
        err = bkgErrsDict[selectionName]
        newHist.SetBinContent(newHist.FindBin(mass_point), content)
        newHist.SetBinError(newHist.FindBin(mass_point), err)
    # newHist = CheckHistBins(newHist)
    for f in rootFiles:
        f.cd()
        newHist.Write()


def CreateAndWriteHist(rootFiles, mass, sample, content, err, name="yieldFinalSelection", titleStart="Yield at LQ"):
    newHist = r.TH1D("{}_{}_{}".format(name, mass, sample), "{} {} final selection for {}".format(titleStart, mass, sample), 1, 0, 1)
    newHist.SetBinContent(1, content)
    newHist.SetBinError(1, err)
    if sample == "DATA":
        newHist.SetName(newHist.GetName().replace("DATA", "data_obs"))
        # print("Hist {} now has {} bins after rebinning with {} bins".format(hist.GetName(), hist.GetNbinsX(), len(xbins)))
    newHist = CheckHistBins(newHist)
    for f in rootFiles:
        f.cd()
        newHist.Write()


def CreateAndWriteHistograms(year):
    outputRootFilename = shapeHistos_filePathBase.format(year)
    if not os.path.exists(datacardHistoDir):
        os.makedirs(datacardHistoDir)
    allOutputRootFile = r.TFile.Open(outputRootFilename, "recreate")
    outputRootFilename = datacardHistoDir + "/" + outputRootFilename
    bkgEvtsDict = {}
    bkgErrsDict = {}
    for i_signal_name, signal_name in enumerate(signal_names):
        outputRootFilename = outputRootFilename.replace(".root", "_{}.root".format(signal_name))
        for iSel, selectionName in enumerate(selectionNames):
            if selectionName == "preselection" or selectionName == "trainingSelection":
                continue
            mass_point = selectionName.replace("LQ", "")
            outputRootFile = r.TFile.Open(outputRootFilename.format(mass_point), "recreate")
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, mass_point)
            signalEvts, signalEvtErrs = d_signalSampleInfos[signalNameForFile].GetRateAndErr(selectionName)
            CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, fullSignalName, signalEvts, signalEvtErrs)
            signalFailEvts, signalFailEvtErrs = d_signalSampleInfos[signalNameForFile].GetFailRateAndErr(selectionName)
            CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, fullSignalName, signalFailEvts, signalFailEvtErrs, "yieldFailFinalSelection", "Yield failing final selection for LQ")
            for ibkg, background_name in enumerate(background_names):
                thisBkgEvts, thisBkgEvtsErr = d_backgroundSampleInfos[background_name].GetRateAndErr(selectionName)
                if thisBkgEvts <= 0:
                    print("INFO: CreateAndWriteHistograms() - Got thisBkgEvts={} for background_name={}, selectionName={}".format(thisBkgEvts, background_name, selectionName))
                    lastSel, rate, err, rawEvents = d_backgroundSampleInfos[background_name].GetNearestPositiveSelectionYield(selectionName)
                    # take upper Poisson 68% CL limit of 1.8410216450092634 and scale it by the factor obtained from the nearest nonzero selection
                    thisBkgEvtsErr = 1.8410216450092634 * rate/rawEvents
                CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, background_name, thisBkgEvts, thisBkgEvtsErr)
                thisBkgFailEvts, thisBkgFailEvtsErr = d_backgroundSampleInfos[background_name].GetFailRateAndErr(selectionName)
                CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, background_name, thisBkgFailEvts, thisBkgFailEvtsErr, "yieldFailFinalSelection", "Yield failing final selection for LQ")
                if background_name not in bkgEvtsDict.keys():
                    bkgEvtsDict[background_name] = {}
                    bkgErrsDict[background_name] = {}
                bkgEvtsDict[background_name][selectionName] = thisBkgEvts
                bkgErrsDict[background_name][selectionName] = thisBkgEvtsErr
            dataRate, dataRateErr = d_dataSampleInfos['DATA'].GetRateAndErr(selectionName)
            dataFailRate, dataFailRateErr = d_dataSampleInfos['DATA'].GetFailRateAndErr(selectionName)
            CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, "DATA", dataRate, dataRateErr)
            CreateAndWriteHist([outputRootFile, allOutputRootFile], mass_point, "DATA", dataFailRate, dataFailRateErr, "yieldFailFinalSelection", "Yield failing final selection for LQ")
            outputRootFile.Close()
    for background_name in background_names:
        CreateAndWriteBackgroundVsMLQHist([allOutputRootFile], background_name, bkgEvtsDict[background_name], bkgErrsDict[background_name])
    allOutputRootFile.Close()


def GetPlotFilename(sampleName, dictSavedSamples, dictSamplesContainingSample):
    if sampleName in dictSavedSamples.keys():
        return "analysisClass_lq_eejj_{}_plots.root".format(sampleName)
    else:
        # find saved sample which has this as a piece
        return "analysisClass_lq_eejj_{}_plots.root".format(dictSamplesContainingSample[sampleName])
    raise RuntimeError("Could not find sample '{}' included as a piece in any of the saved samples '{}'".format(sampleName, dictSavedSamples.keys()))


def FillDicts(rootFilepath, sampleNames, bkgType, dictSavedSamples, dictSamplesContainingSample, year, verbose=False):
    isData = False if "mc" in bkgType.lower() or "signal" in bkgType.lower() else True
    sampleInfos = {}
    # start sample
    for i_sample, sampleName in enumerate(sampleNames):
        thisSampleRootFilepath = rootFilepath
        if ".root" not in rootFilepath:
            thisSampleRootFilepath += GetPlotFilename(sampleName, dictSavedSamples, dictSamplesContainingSample)
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
        systematics = {}
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
                    isData,
                    bkgType == "TTData"
                )
                sampleFailRate, sampleFailRateErr, sampleFailUnscaledRate = cc.GetFailingRatesAndErrors(
                    scaledRootFile,
                    sampleName,
                    selectionName,
                    trainingSelectionCutName,
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
                if sampleName == singleTopSampleName and overrideSingleTopYields and sampleUnscaledRate < overrideThreshold:
                    # print("SICDEBUG: override SingleTop yield for sampleName={}, selectionName={}, unscaledRate={}: {} +/- {}".format(sampleName, selectionName, sampleUnscaledRate, sampleRate, sampleRateErr), end='')
                    sampleRate, sampleRateErr = GetSingleTopYieldAndUncertainty(year, int(mass_point))
                    # print(" --> {} +/- {}".format(sampleRate, sampleRateErr))
                ratesDict[selectionName] = sampleRate
                rateErrsDict[selectionName] = sampleRateErr
                unscaledRatesDict[selectionName] = sampleUnscaledRate
                failRatesDict[selectionName] = sampleFailRate
                failRateErrsDict[selectionName] = sampleFailRateErr
                unscaledFailRatesDict[selectionName] = sampleFailUnscaledRate
                systematicsNominalDict[selectionName] = {"yield": sampleRate, "preselYield": None}
                # print("SICDEBUG: systematics nominalDict for sampleName={}, selectionName={} --> {} ".format(sampleName, selectionName, systematicsNominalDict[selectionName]["yield"]))
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
            systematics = GetSystematicsDict(scaledRootFile, sampleName, selectionNames, year)  # , sampleName == "LQToDEle_M-300_pair")
            if "zjet" in sampleName.lower():
                systematics.update({"DY_Norm": {sel: {"yield": (1+dyNormDeltaXOverX)*ratesDict[sel], "preselYield": None} for sel in selectionNames}})
            if "ttbar" in sampleName.lower() or "ttto2l2nu" in sampleName.lower():
                systematics.update({"TT_Norm": {sel: {"yield": (1+ttBarNormDeltaXOverX)*ratesDict[sel], "preselYield": None} for sel in selectionNames}})
        elif sampleName == "QCDFakes_DATA":
            # special handling for QCDFakes_DATA
            selsAndSysts = {sel: {"yield": (1+qcdNormDeltaXOverX)*ratesDict[sel], "preselYield": None} for sel in selectionNames}
            systematics = {"QCD_Norm": selsAndSysts}
            systematics.update({"nominal": systematicsNominalDict})
        currentSampleInfo = SampleInfo(sampleName, ratesDict, rateErrsDict, unscaledRatesDict, unscaledTotalEvts, failRatesDict, failRateErrsDict, unscaledFailRatesDict, systematics)
        sampleInfos[sampleName] = currentSampleInfo
        scaledRootFile.Close()
    if "mc" in bkgType.lower() and combineSingleTopAndDiboson and overrideSingleTopYields:
        for selectionName in selectionNames:
            # print("SICDEBUG: override SingleTop yield for sampleName={}, selectionName={}: was {} +/- {}, now {} +/- {}".format(otherBkgSampleName, selectionName, sampleInfos[otherBkgSampleName].rates[selectionName], sampleInfos[otherBkgSampleName].rateErrs[selectionName],
            #                                                                                                                     sampleInfos[dibosonSampleName].rates[selectionName] + sampleInfos[singleTopSampleName].rates[selectionName],
            #                                                                                                                     math.sqrt(pow(sampleInfos[dibosonSampleName].rateErrs[selectionName], 2) + pow(sampleInfos[singleTopSampleName].rateErrs[selectionName], 2))))
            # print("SICDEBUG: override SingleTop syst nominal yield for sampleName={}, selectionName={}: was {}, now {}".format(otherBkgSampleName, selectionName, sampleInfos[otherBkgSampleName].systematics["nominal"][selectionName]["yield"],
            #                                                                                                                    sampleInfos[dibosonSampleName].systematics["nominal"][selectionName]["yield"] + sampleInfos[singleTopSampleName].rates[selectionName]))
            sampleInfos[otherBkgSampleName].rates[selectionName] = sampleInfos[dibosonSampleName].rates[selectionName] + sampleInfos[singleTopSampleName].rates[selectionName]
            sampleInfos[otherBkgSampleName].rateErrs[selectionName] = math.sqrt(pow(sampleInfos[dibosonSampleName].rateErrs[selectionName], 2) + pow(sampleInfos[singleTopSampleName].rateErrs[selectionName], 2))
            # sampleInfos[otherBkgSampleName].systematics["nominal"][selectionName]["yield"] = sampleInfos[dibosonSampleName].systematics["nominal"][selectionName]["yield"] + sampleInfos[singleTopSampleName].rates[selectionName]
            # don't do anything about the others, like the unscaled rates (raw MC events).
            # the single top value will in any case be below the threshold specified where we take the fit value, so small
    return sampleInfos


def GetSystNames(year, forBackground):
    systematicsSet = set()
    for sampleName, systList in d_applicableSystematics[year].items():
        if not isinstance(systList, list):
            systList = [systList]
        if forBackground and sampleName not in background_names:
            continue
        if not forBackground and (IsBackground(sampleName) or IsComponentBackground(sampleName)):
            continue
        systematicsSet.update(systList)
    return list(systematicsSet)


def GetSignalSystNames(year):
    return GetSystNames(year, False)


def GetBackgroundSystNames(year):
    return GetSystNames(year, True)


def GetNTotalSysts(years):
    systNames = set()
    for year in years:
        systNames.update(GetBackgroundSystNames(year))
        systNames.update(GetSignalSystNames(year))
    return len(systNames), systNames


def ResetSystematic(sampleInfo, systName, resetWeights=False):
    for selection in selectionNames:
        for subStr in ["Up", "Down"]:
            if not systName+subStr in sampleInfo.systematics.keys():
                sampleInfo.systematics[systName+subStr] = {}
            if not selection in sampleInfo.systematics[systName+subStr].keys():
                sampleInfo.systematics[systName+subStr][selection] = {}
            sampleInfo.systematics[systName+subStr][selection]["yield"] = 0
            sampleInfo.systematics[systName+subStr][selection]["preselYield"] = 0
        if resetWeights and "weight" in systName.lower():
            systs = [syst for syst in sampleInfo.systematics.keys() if systName in syst]
            for syst in systs:
                sampleInfo.systematics[syst][selection]["yield"] = 0
                sampleInfo.systematics[syst][selection]["preselYield"] = 0


def AppendSystematics(systNames, selections, sampleInfo):
    for systName in systNames:
        sampleInfo.systematics[systName] = {}
        for selection in selections:
            sampleInfo.systematics[systName][selection] = {}


def AddSystematics(sampleInfo1, sampleInfo2, systNameStr):
    # print("SIC DEBUG3: AddSystematics(): about to add sample2={} to sample1={}; selectionNames={}".format(sampleInfo2.sampleName, sampleInfo1.sampleName, selectionNames))
    systs = [syst for syst in sampleInfo1.systematics.keys() if systNameStr in syst]
    # print("SIC DEBUG3: AddSystematics():systs={},  sampleInfo1.systematics.keys()={}".format(systs, sampleInfo1.systematics.keys()))
    for selection in selectionNames:
        # print("SIC DEBUG3: AddSystematics(): selection={}; do systs loop".format(selection))
        for syst in systs:
            # if "selection" in selection:
            #     print("SIC DEBUG3: add yield for syst={} selection={} from sample2={} of {} to sample1={} which has {}".format(syst, selection, sampleInfo2.sampleName, sampleInfo2.systematics[syst][selection]["yield"], sampleInfo1.sampleName, sampleInfo1.systematics[syst][selection]["yield"]))
            sampleInfo1.systematics[syst][selection]["yield"] += sampleInfo2.systematics[syst][selection]["yield"]
            hasPreselSystYield1, preselSystYield1 = sampleInfo1.HasPreselectionSystYield(syst, selection)
            hasPreselSystYield2, preselSystYield2 = sampleInfo2.HasPreselectionSystYield(syst, selection)
            if hasPreselSystYield2:
                if not hasPreselSystYield1:
                    sampleInfo1.systematics[syst][selection]["preselYield"] = 0
                sampleInfo1.systematics[syst][selection]["preselYield"] += sampleInfo2.systematics[syst][selection]["preselYield"]
        # print("SIC DEBUG3: AddSystematics(): selection={}; DONE with systs loop".format(selection))
    # print("SIC DEBUG3: AddSystematics(): DONE adding sample2={} to sample1={}".format(sampleInfo2.sampleName, sampleInfo1.sampleName))



def AddSystematic(sampleInfo1, sampleInfo2, systName):
    for selection in selectionNames:
        for subStr in ["Up", "Down"]:
            sampleInfo1.systematics[systName+subStr][selection]["yield"] += sampleInfo2.systematics[sysName+subStr][selection]["yield"]
            sampleInfo1.systematics[systName+subStr][selection]["preselYield"] += sampleInfo2.systematics[systName+subStr][selection]["preselYield"]


def AddSystematicQuad(sampleInfo1, sampleInfo2, systName):
    for selection in selectionNames:
        # here. we assume that it's a symmetric uncertainty
        delta = sampleInfo2.systematics[systName+"Up"][selection]["yield"] - sampleInfo2.systematics["nominal"][selection]["yield"]
        if not "deltaQuad" in sampleInfo1.systematics[systName+"Up"][selection].keys():
            sampleInfo1.systematics[systName+"Up"][selection]["deltaQuad"] = 0
        sampleInfo1.systematics[systName+"Up"][selection]["deltaQuad"] += delta**2
        for subStr in ["Up", "Down"]:
            sampleInfo1.systematics[systName+subStr][selection]["preselYield"] += sampleInfo2.systematics[systName+subStr][selection]["preselYield"]


def ComputeYieldFromDeltaQuad(sampleInfo, systName):
    for selection in selectionNames:
        # here. we assume that it's a symmetric uncertainty
        totalDelta = math.sqrt(sampleInfo.systematics[systName+"Up"][selection]["deltaQuad"])
        nominal = sampleInfo.systematics["nominal"][selection]["yield"]
        sampleInfo.systematics[systName+"Up"][selection]["yield"] = nominal + totalDelta
        sampleInfo.systematics[systName+"Down"][selection]["yield"] = nominal - totalDelta


def GetBranchTitle(sampleInfo, systName, doExcept=False):
    branchTitles = sampleInfo.systematics[systName]["branchTitles"]
    if doExcept and len(branchTitles) > 1:
        raise RuntimeError("Not sure how to handle multiple ({}) branch titles '{}' found for sample={}".format(len(branchTitles), branchTitles, sampleInfo.sampleName))
    return branchTitles


def GetPDFVariationTypeAndName(sampleInfo):
    branchTitles = GetBranchTitle(sampleInfo, "LHEPdfWeight_0")
    # print("SIC DEBUG3: GetPDFVariationTypeAndName for sample={} and branchTitles={} for LHEPdfWeight_0".format(sampleInfo.sampleName, branchTitles))
    if not isinstance(branchTitles, list):
        branchTitles = [branchTitles]
    pdfTypes = set()
    pdfNames = set()
    pdfNMembers = set()
    for bTitle in branchTitles:
        pdfType, pdfName, nMembers = cc.GetPDFVariationType(bTitle)
        pdfTypes.add(pdfType)
        pdfNames.add(pdfName)
        pdfNMembers.add(nMembers)
    if len(pdfTypes) > 1 or len(pdfNMembers) > 1:
        raise RuntimeError("Found multiple PDF types or different numbers of members from branches {} in sample={}: types={}, members={}".format(branchTitles, sampleInfo.sampleName, pdfTypes, pdfNMembers))
    return list(pdfTypes)[0], list(pdfNames)[0]


def ComputePDFVariationMC(sampleInfo, pdfKeys, verbose=False):
    pdfKeys = sorted(pdfKeys, key=lambda x: int(x[x.rfind("_")+1:]))
    # now, if we still have over 100, remove the last two
    if len(pdfKeys) > 100:
        pdfKeys = pdfKeys[:-2]
    elif len(pdfKeys) == 32:
        pdfKeys = pdfKeys[:-2]
    if verbose:
        print("INFO: CalculatePDFVariationMC() -- sampleName={}, we now have {} pdf variations to consider".format(sampleName, len(pdfKeys)), flush=True)
    for selection in selectionNames:
        nominal = sampleInfo.systematics["nominal"][selection]["yield"]
        pdfYieldsUnsorted = [sampleInfo.systematics[pdfKey][selection]["yield"] for pdfKey in pdfKeys]
        # Order the 100 yields and take the 84th and 16th.
        # See eq. 25 here: https://arxiv.org/pdf/1510.03865.pdf
        pdfYields = sorted(pdfYieldsUnsorted)
        if len(pdfKeys) == 100:
            pdfUp = pdfYields[83]
            pdfDown = pdfYields[15]
        elif len(pdfKeys) == 30:
            pdfUp = pdfYields[27]
            pdfDown = pdfYields[5]
        else:
            raise RuntimeError("Could not determine MC PDF variations as we have {} PDF keys".format(len(pdfKeys)))
        if verbose:
            print("INFO: CalculatePDFVariationMC() -- for sampleName={}, selection={}, nominal={}; pdfUp={}, pdfDown={}; pdfUp-nominal={}, pdfDown-nominal={}".format(
                    sampleInfo.sampleName, selection, nominal, pdfUp, pdfDown, pdfUp-nominal, pdfDown-nominal), flush=True)
            # if len(pdfKeys) == 100:
            #     print("INFO: CalculatePDFVariationMC() -- yield 15 = {}; yield 83 = {}".format(pdfYields[15], pdfYields[83]), flush=True)
            # elif len(pdfKeys) == 30:
            #     print("INFO: CalculatePDFVariationMC() -- yield 27 = {}; yield 5 = {}".format(pdfYields[27], pdfYields[5]), flush=True)
            # print("\tINFO: CalculatePDFVariationMC() -- pdfYields={}".format(pdfYields), flush=True)
        sampleInfo.systematics["LHEPdfWeightUp"][selection]["yield"] = pdfUp
        sampleInfo.systematics["LHEPdfWeightDown"][selection]["yield"] = pdfDown


def ComputePDFVariationHessian(sampleInfo, pdfKeys, verbose=False):
    # Sum in quadrature central - var, and use this as a symmetric uncertainty (both the up and down)
    pdfKeys = sorted(pdfKeys, key=lambda x: int(x[x.rfind("_")+1:]))
    # now, if we still have over 100, remove the last two
    if len(pdfKeys) > 100:
        pdfKeys = pdfKeys[:-2]
    elif len(pdfKeys) == 32:
        pdfKeys = pdfKeys[:-2]
    if verbose:
        print("INFO: sampleName={}, we now have {} pdf variations to consider".format(sampleInfo.sampleName, len(pdfKeys)))
    for selection in selectionNames:
        # if not sampleInfo.systematics["LHEPdfWeightHessianNominal"][selection]["yield"]:
        #     sampleInfo.systematics["LHEPdfWeightHessianNominal"][selection]["yield"] = sampleInfo.systematics["nominal"][selection]["yield"]
        # nominal = sampleInfo.systematics["LHEPdfWeightHessianNominal"][selection]["yield"]
        nominal = sampleInfo.systematics["nominal"][selection]["yield"]
        pdfYields = [sampleInfo.systematics[pdfKey][selection]["yield"] for pdfKey in pdfKeys]
        pdfVars = [pow(nominal-pdfYield, 2) for pdfYield in pdfYields]
        pdfDeltaX = math.sqrt(sum(pdfVars))
        if verbose:
            print("INFO: selection={}, nominal={}, pdfDeltaX={}; pdfYields={}; pdfVariations={}".format(selection, nominal, pdfDeltaX, pdfYields, pdfVars))
        if not "LHEPdfWeightHessianUp" in sampleInfo.systematics.keys():
            AppendSystematics(["LHEPdfWeightHessianUp", "LHEPdfWeightHessianDown"], selectionNames, sampleInfo)
        sampleInfo.systematics["LHEPdfWeightHessianUp"][selection]["yield"] = nominal + pdfDeltaX
        sampleInfo.systematics["LHEPdfWeightHessianDown"][selection]["yield"] = nominal - pdfDeltaX


def ComputePDFSystematic(sampleInfo):
    branchTitle = GetBranchTitle(sampleInfo, "LHEPdfWeight_0")
    pdfVariationType, pdfName = GetPDFVariationTypeAndName(sampleInfo)
    pdfKeys = [syst for syst in sampleInfo.systematics.keys() if syst[:syst.rfind("_")] == "LHEPdfWeight"]
    # print("SIC DEBUG3: ComputePDFSystematic for sample={}, branchTitle={}, pdfVariationType={}, pdfName={}; pdfKeys={}".format(sampleInfo.sampleName,
    #                                                                                                                            branchTitle, pdfVariationType,
    #                                                                                                                            pdfName, pdfKeys))
    if pdfVariationType != "mcNoCentral":
        pdfKeys.remove("LHEPdfWeight_0")  # don't consider index 0, central value
    verbose = False
    if "mc" in pdfVariationType:
        ComputePDFVariationMC(sampleInfo, pdfKeys, verbose)
    elif pdfVariationType == "hessian":
        ComputePDFVariationHessian(sampleInfo, pdfKeys, verbose)
    else:
        raise RuntimeError("Unknown PDF type '{}'. Can't calculate the PDF variations for this type (unimplemented).".format(pdfVariationType))


def CombineHessianAndMC(sampleInfo, sampleInfoHessian):
    for selection in selectionNames:
        systName = "LHEPdfWeight"
        totalMCDelta = 0
        if "deltaQuad" in sampleInfo.systematics[systName+"Up"][selection].keys():
            totalMCDelta = math.sqrt(sampleInfo.systematics[systName+"Up"][selection]["deltaQuad"])
        nominal = sampleInfo.systematics["nominal"][selection]["yield"]
        nominalHessian = sampleInfoHessian.systematics["nominal"][selection]["yield"]
        if systName+"HessianUp" in sampleInfoHessian.systematics.keys():
            hessianDelta = sampleInfoHessian.systematics[systName+"HessianUp"][selection]["yield"] - nominalHessian
            totalDelta = math.sqrt(totalMCDelta**2 + hessianDelta**2)
        else:
            totalDelta = totalMCDelta
        sampleInfo.systematics[systName+"Up"][selection]["yield"] = nominal + totalDelta
        sampleInfo.systematics[systName+"Down"][selection]["yield"] = nominal - totalDelta


def ComputeScaleSystematic(sampleInfo, verbose=False):
    branchTitle = GetBranchTitle(sampleInfo, "LHEScaleWeight_0", True)[0]
    validIndices, shapeTitles = cc.ParseShapeBranchTitle(branchTitle)
    shapeKeys = [syst for syst in sampleInfo.systematics.keys() if syst[:syst.rfind("_")] == "LHEScaleWeight"]
    validShapeKeys = [shapeKey for shapeKey in shapeKeys if shapeKey[shapeKey.rfind("_")+1:] in validIndices]
    for selection in selectionNames:
        nominal = sampleInfo.systematics["nominal"][selection]["yield"]
        shapeYields = [sampleInfo.systematics[shapeKey][selection]["yield"] for shapeKey in validShapeKeys]
        if verbose:
            print("\tINFO: For sampleName={}, systName=LHEScaleWeight, selection={}, found validShapeKeys={}, valid shapeTitles={}, shapeYields={}".format(
                    sampleInfo.sampleName, selection, validShapeKeys, shapeTitles, shapeYields))
        deltas = [shapeYield-nominal for shapeYield in shapeYields]
        deltasAbs = [math.fabs(delta) for delta in deltas]
        maxDeltaIdx = deltasAbs.index(max(deltasAbs))
        maxDelta = deltas[maxDeltaIdx]
        if verbose:
            print("\tINFO: For sampleName={}, systName=LHEScaleWeight, selection={}, found deltas={} with maxDelta={} and nominal={}".format(
                    sampleInfo.sampleName, selection, deltas, maxDelta, nominal))
        sampleInfo.systematics["LHEScaleWeightUp"][selection]["yield"] = nominal + math.fabs(maxDelta)
        sampleInfo.systematics["LHEScaleWeightDown"][selection]["yield"] = nominal - math.fabs(maxDelta)
        sampleInfo.systematics["LHEScaleWeightUp"][selection]["preselYield"] = sampleInfo.systematics[validShapeKeys[maxDeltaIdx]]["preselection"]["yield"]
        sampleInfo.systematics["LHEScaleWeightDown"][selection]["preselYield"] = sampleInfo.systematics[validShapeKeys[maxDeltaIdx]]["preselection"]["yield"]


def RecomputeBackgroundSystematic(systName, sampleInfos):
    for sample in sampleInfos.keys():
        if not IsBackground(sample):
            continue
        if sample not in dictSampleComponents.keys():
            continue  # this background only has one component
        if "LHEScaleWeight" == systName:
            ResetSystematic(sampleInfos[sample], systName, True)
            for component in dictSampleComponents[sample]:
                    ResetSystematic(sampleInfos[component], systName)
                    ComputeScaleSystematic(sampleInfos[component])
                    AddSystematicQuad(sampleInfos[sample], sampleInfos[component], systName)
            ComputeYieldFromDeltaQuad(sampleInfos[sample], systName) #FIXME: handle ZJets differently?
        if "LHEPdfWeight" == systName:
            ResetSystematic(sampleInfos[sample], systName, True)
            componentsWithMCVariations = []
            hessianComponent = None
            for component in dictSampleComponents[sample]:
                ResetSystematic(sampleInfos[component], systName)
                pdfType, pdfName = GetPDFVariationTypeAndName(sampleInfos[component])
                if "mc" in pdfType:
                    ComputePDFSystematic(sampleInfos[component])
                    componentsWithMCVariations.append(component)
                    AddSystematicQuad(sampleInfos[sample], sampleInfos[component], systName)
                else:
                    if not hessianComponent:
                        hessianComponent = copy.deepcopy(sampleInfos[component])
                        hessianComponent.sampleName = sample + " [hessian component]"
                    else:
                        hessianComponent += sampleInfos[component]
            if hessianComponent is not None:
                ComputePDFSystematic(hessianComponent)
            CombineHessianAndMC(sampleInfos[sample], hessianComponent)



def WriteDatacard(card_file_path, year):
    thresholdForAutoMCStats = 10
    yearStr = "# " + year + "\n"
    lumiStr = "# " + str(intLumi) + "\n\n"
    card_file = open(card_file_path, "w")
    card_file.write(yearStr)
    card_file.write(lumiStr)
    for i_signal_name, signal_name in enumerate(signal_names):
        doMassPointLoop = True
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
                txt_file_name = datacardHistoDir.rstrip("/") + "/datacard_" + year + "_" + fullSignalName + ".txt"
                indivMassPointCard = open(txt_file_name, "w")
                indivMassPointCard.write(yearStr)
                indivMassPointCard.write(lumiStr)
    
                datacardLines = []
                datacardLines.append("# " + txt_file_name + "\n\n")
                datacardLines.append("imax " + str(n_channels) + "\n")
                datacardLines.append("jmax " + str(n_background) + "\n")
                if doSystematics:
                    n_systs, systNames = GetNTotalSysts(years)
                    datacardLines.append("kmax " + str(n_systs) + "\n\n")
                else:
                    datacardLines.append("kmax 0\n\n")
                datacardLines.append("---------------\n")
                datacardLines.append("shapes * * "+shapeHistos_filePathBase.format(year)+" yieldFinalSelection_$MASS_$PROCESS\n")
                datacardLines.append("---------------\n")
                datacardLines.append("bin bin1\n\n")
    
                total_data, total_data_err = d_dataSampleInfos["DATA"].GetRateAndErr(selectionName)
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
                signalYield, signalYieldErr = d_signalSampleInfos[signalNameForFile].GetRateAndErr(selectionName)
                line = "rate {} ".format(signalYield)
                totalBackgroundYield = 0
                for background_name in background_names:
                    bkgYield, bkgYieldErr = d_backgroundSampleInfos[background_name].GetRateAndErr(selectionName)
                    line += "{} ".format(bkgYield)
                    totalBackgroundYield += bkgYield
                datacardLines.append(line.strip() + "\n")
                datacardLines.append("------------------------------\n")
                datacardLines.append("* autoMCStats "+str(thresholdForAutoMCStats)+"\n")
    
                # print signal_name, mass_point, total_signal, total_bkg, total_data
                # print signal_name+str(mass_point), total_signal, total_bkg
    
            # recall the form: systDict['PileupUp'/systematicFromHist]['ZJet_amcatnlo_ptBinned'/sampleName]['LQXXXX'/selection] = yield
            # for RPV, select proper signalSystDict based on ctau of signal
            # if doRPV:
            #     ctau = int(signal_name[signal_name.find("CTau") + 4:])
            #     signalSystDict = signalSystDictByCTau[ctau]
            if doSystematics:
                for syst in sorted(systematicsNamesBackground["all"]):
                    # if int(mass_point) > maxLQSelectionMass:
                    #     selectionNameSyst = "LQ"+str(maxLQSelectionMass)
                    # else:
                    #     selectionNameSyst = selectionName
                    selectionNameSyst = selectionName
                    systName = syst
                    # if systName == "Lumi":
                    #     systName += year.replace("preVFP", "").replace("postVFP", "")
                    line = systName + " lnN "
                    if syst in systematicsNamesSignal["all"] and selectionName != "preselection" and selectionName != "trainingSelection":
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
                        systEntry, deltaOverNominalUp, deltaOverNominalDown, systNomYield, systSelection, preselRatios = d_signalSampleInfos[signalNameForFile].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                        thisSigEvts, thisSigEvtsErr = d_signalSampleInfos[signalNameForFile].GetRateAndErr(selectionNameSyst)
                        thisSigSystUp = deltaOverNominalUp*systNomYield
                        thisSigSystDown = deltaOverNominalDown*systNomYield
                        d_systNominals[signalNameForFile][syst][selectionNameSyst] = systNomYield
                        _, d_systNominalErrs[signalNameForFile][syst][selectionNameSyst] = d_signalSampleInfos[signalNameForFile].GetRateAndErr(systSelection)
                        d_systUpDeltas[signalNameForFile][syst][selectionNameSyst] = thisSigSystUp
                        d_systDownDeltas[signalNameForFile][syst][selectionNameSyst] = thisSigSystDown
                        if EntryIsValid(systEntry):
                            d_signalSampleInfos[signalNameForFile].UpdateSystsApplied(syst)
                            line += str(systEntry) + " "
                            # print("INFO: Got valid syst entry For sample={} selection={} syst={}, systEntry={}, thisSigEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                            #         signalNameForFile, selectionNameSyst, syst, systEntry, thisSigEvts, thisSigSystUp, thisSigSystDown))
                        else:
                            print("INFO: Got invalid syst entry For sample={} selection={} syst={}, systEntry={}, thisSigEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                    signalNameForFile, selectionNameSyst, syst, systEntry, thisSigEvts, thisSigSystUp, thisSigSystDown))
                            line += "- "
                    else:
                        line += "- "
                    for ibkg, background_name in enumerate(background_names):
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
                        systEntry, deltaOverNominalUp, deltaOverNominalDown, systNomYield, systSelection, preselRatios = d_backgroundSampleInfos[background_name].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                        thisBkgEvts, thisBkgEvtsErr = d_backgroundSampleInfos[background_name].GetRateAndErr(selectionNameSyst)
                        thisBkgSystUp = deltaOverNominalUp*systNomYield
                        thisBkgSystDown = deltaOverNominalDown*systNomYield
                        d_systNominals[background_name][syst][selectionNameSyst] = systNomYield
                        _, d_systNominalErrs[background_name][syst][selectionNameSyst] = d_backgroundSampleInfos[background_name].GetRateAndErr(systSelection)
                        d_systUpDeltas[background_name][syst][selectionNameSyst] = thisBkgSystUp
                        d_systDownDeltas[background_name][syst][selectionNameSyst] = thisBkgSystDown
                        # bkg component info, for those backgrounds which have multiple components/subprocesses
                        if background_name in dictSampleComponents.keys():
                            for componentBkg in dictSampleComponents[background_name]:
                                if componentBkg not in list(d_systNominals.keys()):
                                    d_systNominals[componentBkg] = {}
                                    d_systNominalErrs[componentBkg] = {}
                                    d_systUpDeltas[componentBkg] = {}
                                    d_systDownDeltas[componentBkg] = {}
                                if syst not in list(d_systNominals[componentBkg].keys()):
                                    d_systNominals[componentBkg][syst] = {}
                                    d_systNominalErrs[componentBkg][syst] = {}
                                    d_systUpDeltas[componentBkg][syst] = {}
                                    d_systDownDeltas[componentBkg][syst] = {}
                                compSystEntry, compDeltaOverNominalUp, compDeltaOverNominalDown, compSystNomYield, compSystSelection, preselRatios = d_backgroundSampleInfos[componentBkg].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                                d_systNominals[componentBkg][syst][selectionNameSyst] = compSystNomYield
                                compBkgEvts, d_systNominalErrs[componentBkg][syst][selectionNameSyst] = d_backgroundSampleInfos[componentBkg].GetRateAndErr(compSystSelection)
                                d_systUpDeltas[componentBkg][syst][selectionNameSyst] = compDeltaOverNominalUp*compSystNomYield
                                d_systDownDeltas[componentBkg][syst][selectionNameSyst] = compDeltaOverNominalDown*compSystNomYield
                                # print("INFO: Got syst entry for sample={} selectionNameSyst={} systSelection={}, syst={}, systEntry={}, thisBkgEvts={}, systNomYield={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                #         componentBkg, selectionNameSyst, systSelection, syst, systEntry, compBkgEvts, systNomYield, deltaOverNominalUp*systNomYield, deltaOverNominalDown*systNomYield))
                        if EntryIsValid(systEntry):
                            d_backgroundSampleInfos[background_name].UpdateSystsApplied(syst)
                            line += str(systEntry) + " "
                        else:
                            # if "lumi" in syst.lower():
                            #     print("INFO2: Got invalid syst entry for sample={} selection={} syst={}, systEntry={}, thisBkgEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                            #             background_name, selectionNameSyst, syst, systEntry, thisBkgEvts, thisBkgSystUp, thisBkgSystDown))
                            line += "- "
                        # if "ZJet" in background_name and "800" in selectionNameSyst and "EES" in syst:
                        try:
                            if float(systEntry) < 0:
                                 #print("INFO: For sample={} selection={} syst={}, systEntry={}, thisBkgEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                 raise RuntimeError("For year={}, sample={} selection={} syst={}, systEntry={}, thisBkgEvts={}, d_systUpDeltas={}, d_systDownDeltas={}".format(
                                         "all", background_name, selectionNameSyst, syst, systEntry, thisBkgEvts, thisBkgSystUp, thisBkgSystDown))
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
combineSingleTopAndDiboson = True
overrideSingleTopYields = True
overrideThreshold = 2  # if < 2 raw MC events, then override the nominal yield with the fit
# doRPV = False
# forceGmNNormBkgStatUncert = False
# cc.finalSelectionName = "sT_eejj"  # "min_M_ej"
cc.finalSelectionName = "BDTOutput_"  # has LQ300 (for example) appended to it inside the code
# cc.finalSelectionName = "Meejj"
trainingSelectionCutName = "trainingSelection"
allYears = ["2016preVFP", "2016postVFP", "2017", "2018"]

sampleListForMerging = "$LQANA/config/sampleListForMerging_13TeV_eejj_{}.yaml"
sampleListsForMergingQCD = "$LQANA/config/sampleListForMerging_13TeV_QCD_dataDriven_{}.yaml"
eosBaseDir = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/"
qcdFile = eosBaseDir + "ultralegacy/analysis/{}/{}/output_qcdSubtractedYield/qcdSubtracted_plots.root"
eosDir = eosBaseDir + "ultralegacy/analysis/{}/{}/output_cutTable_lq_eejj_BDT/"

# LQ case
mass_points_int = [
    i for i in range(300, 3100, 100)
]  # go from 300-3000 in 100 GeV steps
mass_points = [str(i) for i in mass_points_int]
## mass_points.extend(["3500", "4000"])
# mass_points = [
#     str(i) for i in range(1500, 2100, 100)
# ]
# if doRPV:
#     mass_points = [
#         str(i) for i in range(200, 1250, 100)
#     ]  # go from 200-1200 in 100 GeV steps
#     signal_names = [
#         "Stop_M[masspoint]_CTau1000",
#         "Stop_M[masspoint]_CTau100",
#         "Stop_M[masspoint]_CTau10",
#         "Stop_M[masspoint]_CTau1",
#     ]
#     # signal_names = [ "Stop_M[masspoint]_CTau10","Stop_M[masspoint]_CTau1"]
#     # put in some more signals that don't fit the general pattern
#     # signal_names = ['Stop_M100_CTau100','Stop_M125_CTau100','Stop_M150_CTau100','Stop_M175_CTau100','Stop_M200_CTau50'] + signal_names
#     signalNameTemplate = "stop"

# maxLQSelectionMass = 1200  # max background selection point used
maxLQSelectionMass = 100000 # FIXME remove
# background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTBarFromDATA", "ZJet_amcatnlo_ptBinned", "WJet_amcatnlo_ptBinned", "DIBOSON","SingleTop"  ]
background_names = [
    "ZJet_amcatnlo_ptBinned_IncStitch",
    # "TTBarFromDATA",
    "TTTo2L2Nu"] + (["QCDFakes_DATA"] if doQCD else [])
singleTopSampleName = "SingleTop"
dibosonSampleName = "DIBOSON_nlo"
otherBkgSampleName = "OTHERBKG_dibosonNLO_singleTop"
if combineSingleTopAndDiboson:
    background_names.append(otherBkgSampleName)
else:
    background_names.extend([dibosonSampleName, singleTopSampleName])
background_fromMC_names = [bkg for bkg in background_names if "data" not in bkg.lower()]
background_QCDfromData = [bkg for bkg in background_names if "data" in bkg.lower() and "qcd" in bkg.lower()]
systTitleDict = OrderedDict([
        ("EleRecoSF", "Electron reconstruction"),
        ("EleIDSF", "Electron identification"),
        ("EES", "Electron energy scale"),
        ("EER", "Electron energy resolution"),
        ("JER", "Jet energy resolution"),
        ("JES", "Jet energy scale"),
        ("UnclusteredEne", "Unclustered energy (MET)"),
        ("EleTrigSF", "Trigger"),
        ("Pileup", "Pileup"),
        ("Prefire", "L1EcalPrefiring"),
        ("LHEPdfWeight", "PDF"),
        ("DY_Norm", "DYJ normalization"),
        ("DY_Shape", "DYJ shape"),
        ("TT_Norm", "TTbar normalization"),
        ("TT_Shape", "TTbar shape")
        ])
if combineSingleTopAndDiboson:
    systTitleDict["OtherBkg_Shape"] = "Other bkg. shape"
else:
     systTitleDict["Diboson_Shape"] = "Diboson shape"
     systTitleDict["ST_Shape"] = "SingleTop shape"
systTitleDict["QCD_Norm"] = "Multijet bkg. normalization"
systTitleDict["LumiCorrelated"] = "Lumi correlated"
otherBackgrounds = [otherBkgSampleName] if combineSingleTopAndDiboson else [dibosonSampleName, singleTopSampleName]

zjetsSampleName = GetSampleNameFromSubstring("ZJet", background_names)
ttbarSampleName = GetSampleNameFromSubstring("TTTo2L2Nu", background_names)

backgroundTitlesDict = {zjetsSampleName: "Z+jets", ttbarSampleName: "TTbar", "QCDFakes_DATA": "QCD(data)", "WJet_amcatnlo_ptBinned": "W+jets", "WJet_amcatnlo_jetBinned": "W+jets",
                        dibosonSampleName: "DIBOSON", "TRIBOSON": "TRIBOSON", "TTX": "TTX", singleTopSampleName: "SingleTop",
                        "PhotonJets_Madgraph": "Gamma+jets", otherBkgSampleName: "Other bkg. (VV+jets + SingleTop)"}
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
dyNormDeltaXOverX = 0.2
ttBarNormDeltaXOverX = 0.1
qcdNormDeltaXOverX = 0.4
lumiDeltaXOverX = {}
lumiCorrelatedDeltaXOverX = {}
lumi1718CorrelatedDeltaXOverX = {}
lumiDeltaXOverX["2016preVFP"] = 0.01
lumiCorrelatedDeltaXOverX["2016preVFP"] = 0.006
lumiDeltaXOverX["2016postVFP"] = 0.01
lumiCorrelatedDeltaXOverX["2016postVFP"] = 0.006
# from: https://arxiv.org/pdf/1506.04072.pdf tables 1-3: 5% is probably good enough for SingleTop PDF syst
# additionalBkgSystsDict["SingleTop"] = {"LHEPdfWeight": {sel: 0.10 for sel in selectionNames}}
lumiDeltaXOverX["2017"] = 0.02
lumiCorrelatedDeltaXOverX["2017"] = 0.009
lumi1718CorrelatedDeltaXOverX["2017"] = 0.006
lumiDeltaXOverX["2018"] = 0.015
lumiCorrelatedDeltaXOverX["2018"] = 0.02
lumi1718CorrelatedDeltaXOverX["2018"] = 0.002
d_flatSystematics = {}
for year in allYears:
    d_flatSystematics.update({year: {"Lumi"+year: lumiDeltaXOverX[year]} })
    d_flatSystematics[year].update({"LumiCorrelated": lumiCorrelatedDeltaXOverX[year]})
    if year == "2017" or year == "2018":
        d_flatSystematics[year].update({"LumiCorrelated1718": lumi1718CorrelatedDeltaXOverX[year]})

n_background = len(background_names)
# all bkg systematics, plus stat 'systs' for all bkg plus signal plus 3 backNormSysts
# update July/Sep 2020: no more data-driven TBar, so only 1 extra syst for eejj
n_channels = 1

backgroundsToRenormSystAtPresel = [zjetsSampleName, ttbarSampleName]
plots_filePath = "plots.root"
plotsDir = "makeDatacard_plots"
datacardHistoDir = "makeDatacard_separateMassPoints"
tablesDir = "makeDatacard_tables"
systematics_dictFilePath = "systematics_dict_{}.txt"
preselSystScaleFactorsJson = "preselectionSystSFs_{}.json"
systSelectionsUsedJson = "systSelectionsUsed_{}.json"
shapeHistos_filePathBase = "shapeHistos_{}.root"

###################################################################################################
# RUN
###################################################################################################
if len(sys.argv) < 4:
    print("ERROR: need to specify qcdAnalysisName, dataMCAnalysisName, year")
    print("Example: makeDatacard.py qcdAnalysisName dataMCAnalysisName 2016pre")
    raise RuntimeError("ERROR: need to specify qcdPath, dataMCPath, year")

qcdAnaName = sys.argv[1]
dataMCAnaName = sys.argv[2]
requestedYear = sys.argv[3]
if len(sys.argv) > 4:
    signalName = sys.argv[4]
    if "LQToDEle" in signalName:
        signalNameTemplate = "LQToDEle_M-{}_pair"
    elif "LQToBEle" in signalName:
        signalNameTemplate = "LQToBEle_M-{}_pair"
    else:
        raise RuntimeError("Couldn't understand provided signal name '{}'. Must contain LQToDEle or LQToBEle.".format(signalName))
else:
    print("INFO: Defaulting to LQToDEle_pair signal")
    signalNameTemplate = "LQToDEle_M-{}_pair"

if requestedYear != "all":
    years = requestedYear.split(",")
    yearsRequestedStr = ("_").join(years)
else:
    years = allYears
    yearsRequestedStr = "all"
for idx, year in enumerate(years):
    if year == "2016pre":
        year = "2016preVFP"
    elif year == "2016post":
        year = "2016postVFP"
    years[idx] = year


# if not signalNameTemplate.split("_")[0] in dataMCAnaName or not signalNameTemplate.split("_")[0] in qcdAnaName:
#     raise RuntimeError("signalNameTemplate specified is {}, while it does not appear to match the dataMC or qcd analysis names given.".format(signalNameTemplate))

signal_names = [signalNameTemplate]
signalNameList = [GetFullSignalName(signalNameTemplate, massPoint)[1] for massPoint in mass_points]

d_systTitles = OrderedDict()
intLumi = 0
validYear = False
for year in years:
    if year not in allYears:
        raise RuntimeError("Provided year '{}' is not one of 2016pre(VFP)/2016post(VFP)/2017/2018.".format(year))
    d_systTitles[year] = copy.deepcopy(systTitleDict)
    if "2016preVFP" == year:
        intLumi += 19501.601622000
        d_systTitles[year]["Lumi2016preVFP"] = "Lumi (2016preVFP)"
    elif "2016postVFP" == year:
        intLumi += 16812.151722000
        d_systTitles[year]["Lumi2016postVFP"] = "Lumi (2016postVFP)"
    elif "2017" == year:
        do2017 = True
        intLumi += 41477.877399
        d_systTitles[year]["Lumi2017"] = "Lumi (2017)"
        d_systTitles[year]["LumiCorrelated1718"] = "Lumi correlated 2017-2018"
    elif "2018" == year:
        do2018 = True
        intLumi += 59827.449483
        d_systTitles[year]["Lumi2018"] = "Lumi (2018)"
        d_systTitles[year]["LumiCorrelated1718"] = "Lumi correlated 2017-2018"
cc.intLumi = intLumi

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

bkgsToAddList = [singleTopSampleName] if len(singleTopSampleName) else []
bkgsToAddList.extend([dibosonSampleName] if len(dibosonSampleName) else [])
systematicsNamesBackground = OrderedDict((year, list(systTitlesDict.keys())) for year, systTitlesDict in d_systTitles.items())
allBkgSysts = {year: [syst for syst in systematicsNamesBg if "norm" not in syst.lower() and "shape" not in syst.lower()] for year, systematicsNamesBg in systematicsNamesBackground.items()}
d_applicableSystematics = {year: {bkg: list(allBkgSysts[year]) for bkg in background_fromMC_names + bkgsToAddList} for year in years}
for year in years:
    d_applicableSystematics[year].update({bkg: ["QCD_Norm"] for bkg in background_QCDfromData})
    d_applicableSystematics[year][zjetsSampleName].append("DY_Norm")
    d_applicableSystematics[year][zjetsSampleName].append("DY_Shape")
    d_applicableSystematics[year][ttbarSampleName].append("TT_Norm")
    d_applicableSystematics[year][ttbarSampleName].append("TT_Shape")
    if combineSingleTopAndDiboson:
        d_applicableSystematics[year][otherBkgSampleName].append("OtherBkg_Shape")
    else:
        d_applicableSystematics[year][dibosonSampleName].append("Diboson_Shape")
        d_applicableSystematics[year][singleTopSampleName].append("ST_Shape")
systematicsNamesSignal = {}
systematicsNamesSignal = {year: [syst for syst in systematicsNamesBg if "shape" not in syst.lower() and "norm" not in syst.lower() and "LHEPdfWeight" not in syst] for year, systematicsNamesBg in systematicsNamesBackground.items()}
# systematicsNamesSignal.remove("LHEPdfWeight")
for year in years:
    d_applicableSystematics[year].update({sig: systematicsNamesSignal[year] for sig in signalNameList})
n_systematics = {}
if doSystematics:
    n_systematics = {year:len(systematicsNames) for year, systematicsNames in systematicsNamesBackground.items()}
else:
    n_systematics = n_background + 1

# this has the TopPtReweight+updatedSF and the Z+jets St corrections at final selections
# filePath = os.environ["LQDATA"] + '/RunII/eejj_analysis_zJetsStCorrectionFinalSelections_21jul/output_cutTable_lq_eejj/'

# ttbar_data_filepath = ttbarFilePath + "analysisClass_lq_ttbarEst_plots.root"
# SIC 6 Jul 2020 remove
# ttbar_data_filepath = ""

sampleListsForMerging = [os.path.expandvars(sampleListForMerging.format(year)) for year in years]
# SIC 6 Jul 2020 remove
# sampleListForMergingTTBar = (
#     os.environ["LQANA"] + "/config/sampleListForMerging_13TeV_ttbarBkg_emujj.txt"
# )
qcdFilePaths = []
if doQCD:
    sampleListsForMergingQCD = [os.path.expandvars(sampleListsForMergingQCD.format(year)) for year in years]
    qcdFilePaths = [os.path.expandvars(qcdFile.format(year, qcdAnaName)) for year in years]
filePaths = [os.path.expandvars(eosDir.format(year, dataMCAnaName)) for year in years]
dataMC_filepaths = [filePath.rstrip("/") + "/" for filePath in filePaths]
# ---Check if sampleListForMerging files exist
for sampleListForMerging in sampleListsForMerging + sampleListsForMergingQCD:
    if not os.path.isfile(sampleListForMerging):
        print("ERROR: file " + sampleListForMerging + " not found")
        print("exiting...")
        sys.exit(1)

print("Launched like:")
for arg in sys.argv:
    print("\t" + arg)
print("Using input root paths/files:")
print("\t Data/MC:", dataMC_filepaths)
if doQCD:
    print("\t QCD(data):", qcdFilePaths)
else:
    print("\t No QCD used.")

dictSamplesByYear = {year: cc.GetSamplesToCombineDict(sampleListForMerging) for sampleListForMerging, year in zip(sampleListsForMerging, years)}
# TODO potentially rewrite to store components in SampleInfo class
dictSavedSamplesByYear = {}
dictSamplesContainingSampleByYear = {}
dictSampleComponentsByYear = {}
dictDesiredSamplesByYear = {}
dictDesiredBackgroundSamplesByYear = {}
for year, dictSamples in dictSamplesByYear.items():
    dictSavedSamples = {k: v for k, v in dictSamples.items() if dictSamples[k]["save"]}
    dictSamplesContainingSample = {}
    dictSampleComponents = {}
    dictDesiredSamples = {}
    savedSamples = list(dictSavedSamples.keys())
    for sample in background_fromMC_names + bkgsToAddList:
        dictDesiredSamples[sample] = dictSamples[sample]
        pieceList = dictSavedSamples[sample]["pieces"]
        if len(pieceList) < 2:
            continue  # in this case, we have only a single background component
        expandedPieces = cc.ExpandCompositePieces(pieceList, dictSamples) if len(pieceList) > 1 else [pieceList]
        dictSampleComponents[sample] = []
        for piece in expandedPieces:
            dictSamplesContainingSample[piece] = sample
            dictSampleComponents[sample].append(piece)
            dictDesiredSamples[piece] = dictSamples[piece]
            d_applicableSystematics[year][piece] = d_applicableSystematics[year][sample]
    dictDesiredBackgroundSamplesByYear[year] = dict(dictDesiredSamples)
    dictDesiredSamples.update(dictSavedSamples)
    dictSavedSamplesByYear[year] = dictSavedSamples
    dictSamplesContainingSampleByYear[year] = dictSamplesContainingSample
    dictSampleComponentsByYear[year] = dictSampleComponents
    dictDesiredSamplesByYear[year] = dictDesiredSamples

d_applicableSystematics["all"] = {}
allSamples = set()
for year in years:
    for sampleName in d_applicableSystematics[year].keys():
        if sampleName not in d_applicableSystematics["all"].keys():
            d_applicableSystematics["all"][sampleName] = []
        d_applicableSystematics["all"][sampleName].extend(d_applicableSystematics[year][sampleName])
        allSamples.add(sampleName)
for sampleName in allSamples:
    d_applicableSystematics["all"][sampleName] = set(d_applicableSystematics["all"][sampleName])

systematicsNamesSignal["all"] = []
for year in years:
    systematicsNamesSignal["all"].extend(systematicsNamesSignal[year])
systematicsNamesSignal["all"] = set(systematicsNamesSignal["all"])

d_systUpDeltas = {}
d_systDownDeltas = {}
d_systNominals = {}  # same as rates, but more conveniently indexed
d_systNominalErrs = {}
d_datacardStatErrs = {}

d_backgroundSampleInfos = {}
d_signalSampleInfos = {}
d_dataSampleInfos = {}
for idx, year in enumerate(years):
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

    # d_systUpDeltas[year] = {}
    # d_systDownDeltas[year] = {}
    # d_systNominals[year] = {}  # same as rates, but more conveniently indexed
    # d_systNominalErrs[year] = {}
    # d_datacardStatErrs[year] = {}

    # rates/etc.
    print("INFO: Filling background [MC] information for year {}...".format(year))
    dataMC_filepath = dataMC_filepaths[idx]
    # HACK!
    # if year == "2017" or year == "2018":
    #     dataMC_filepath = dataMC_filepath.replace("eejj_6mar2025_bdt_minBkg0p5", "eejj_6mar2025_bdt_minBkg0p5_orig")
    # END HACK!
    # dictDesiredSamples = dictDesiredSamplesByYear[year]
    dictDesiredBackgroundSamples = dictDesiredBackgroundSamplesByYear[year]
    dictSamples = dictSavedSamplesByYear[year]
    dictSamplesContainingSample = dictSamplesContainingSampleByYear[year]
    dictSampleComponents = dictSampleComponentsByYear[year]
    qcdFilePath = qcdFilePaths[idx]
    d_backgroundSampleInfos[year] = FillDicts(dataMC_filepath, list(dictDesiredBackgroundSamples.keys()), "MC", dictSavedSamples, dictSamplesContainingSample, year)
    if doQCD:
        print("INFO: Filling QCD[data] information for year {}...".format(year))
        bgFromDataSampleInfos = FillDicts(qcdFilePath, background_QCDfromData, "DATA", dictSavedSamples, dictSamplesContainingSample, year)
        d_backgroundSampleInfos[year].update(bgFromDataSampleInfos)
    # if doSystematics:
    #     for sampleName in additionalBkgSystsDict.keys():
    #         d_background_systs[year][sampleName].update(additionalBkgSystsDict[sampleName])
    # above would be similar for TTBarFromDATA
    print("INFO: Filling signal information for year {}...".format(year))
    d_signalSampleInfos[year] = FillDicts(dataMC_filepath, signalNameList, "signal", dictSavedSamples, dictSamplesContainingSample, year)
    print("INFO: Filling data information for year {}...".format(year))
    d_dataSampleInfos[year] = FillDicts(dataMC_filepath, ["DATA"], "DATA", dictSavedSamples, dictSamplesContainingSample, year)
    # print one of the systematics for checking
    # for syst in backgroundSystDict.keys():
    #    print 'Syst is:',syst
    #    print 'selection\t\tvalue'
    #    for selection in sorted(backgroundSystDict[syst].keys()):
    #        print selection+'\t\t'+str(backgroundSystDict[syst][selection])
    #    break
    # print signalSystDict
    # print backgroundSystDict

SumSampleInfoOverYears(d_backgroundSampleInfos, years)
SumSampleInfoOverYears(d_signalSampleInfos, years)
SumSampleInfoOverYears(d_dataSampleInfos, years)
if doSystematics:
    # now, need to recompute LHEScale/LHEPdfUncertainties
    RecomputeBackgroundSystematic("LHEScaleWeight", d_backgroundSampleInfos)
    RecomputeBackgroundSystematic("LHEPdfWeight", d_backgroundSampleInfos)

systematicsNamesBackground["all"] = []
d_systTitles["all"] = {}
for year in years:
    systematicsNamesBackground["all"].extend(systematicsNamesBackground[year])
    d_systTitles["all"].update(d_systTitles[year])
systematicsNamesBackground["all"] = list(dict.fromkeys((systematicsNamesBackground["all"])))

print("INFO: Preparing shape histograms...", end=' ')
CreateAndWriteHistograms(yearsRequestedStr)
print("Done")

# selectionNames = ["preselection"]
# selectionNames.extend(["LQ"+str(mass) for mass in mass_points])
print("INFO: Preparing datacard...", end=' ')
datacard_filePath = "tmp_card_file_{}_{}.txt".format(signalNameTemplate.split("_")[0], yearsRequestedStr)
WriteDatacard(datacard_filePath, yearsRequestedStr)
print("Done")

if not os.path.exists(tablesDir):
    os.makedirs(tablesDir)

if doSystematics:
    # tables
    columnNames = ["Systematic", "Signal (%)", "Background (%)"]
    for selectionName in selectionNames:
        print("INFO: total syst table for selection {}".format(selectionName))
        table = []
        for syst in systematicsNamesBackground["all"]:
            selectionNameSyst = selectionName
            if selectionName != "preselection" and selectionName != "trainingSelection":
                massPoint = selectionName.replace("LQ", "")
                if int(massPoint) > maxLQSelectionMass:
                    selectionNameSyst = "LQ"+str(maxLQSelectionMass)
            if syst in systematicsNamesSignal["all"] and selectionName != "preselection" and selectionName != "trainingSelection":
                for i_signal_name, signal_name in enumerate(signal_names):
                    fullSignalName, signalNameForFile = GetFullSignalName(signal_name, massPoint)
                    thisEntry, sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown, _, _, _ = d_signalSampleInfos[signalNameForFile].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
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
                thisEntry, bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown, _, _, _ = d_backgroundSampleInfos[background_name].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                thisBkgNominalEvts, _ = d_backgroundSampleInfos[background_name].GetRateAndErr(selectionName)
                # print("INFO: for background {}, systematic {}, selection {}, bkgYieldNominalAtSyst={}, bkgEvts={}".format(background_name, syst, selectionNameSyst, d_systNominals[background_name][syst][selectionNameSyst], thisBkgNominalEvts))
                # print("INFO: for background {}, systematic {}, selection {}, d_systNominals={}, bkgEvts={}".format(background_name, syst, selectionNameSyst, d_systNominals, thisBkgNominalEvts))
                totalBkgNominal += d_systNominals[background_name][syst][selectionNameSyst]
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
            table.append([d_systTitles["all"][syst], thisSigSystPercent, totalBkgSystPercent])
        print("Selection: {}".format(selectionName))
        print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".1f"))
        print(tabulate(table, headers=columnNames, tablefmt="latex", floatfmt=".1f"))
        print()
        with open(tablesDir+"/"+selectionName+".txt", "w") as table_file:
            table_file.write(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".1f"))
        with open(tablesDir+"/"+selectionName+".tex", "w") as table_file:
            table_file.write(tabulate(table, headers=columnNames, tablefmt="latex", floatfmt=".1f"))

    # print info on systematics used
    print(yearsRequestedStr)
    print("{0:40}\t{1}".format("sampleName", "systematics applied"))
    for sampleName, info in d_backgroundSampleInfos.items():
        if isinstance(info, dict):
            continue
        if sampleName in background_names:
            print("*", end="")
        print("{0:40}\t{1}".format(sampleName, sorted(list(info.systematicsApplied))))
    print("{0:40}\t{1}".format("sampleName", "systematics applied"))
    for sampleName, info in d_signalSampleInfos.items():
        if signalNameTemplate.format(300) not in sampleName:
            continue
        print("{0:40}\t{1}".format(sampleName, sorted(list(info.systematicsApplied))))
    print()

    # compute total syst for all backgrounds at each selection
    # also store info for systs renormalized at preselection
    # as well as the systematic selection used for all selections
    d_totalSystDeltaXOverX = {}
    d_preselectionSystScaleFactors = {}
    d_systematicSelectionsUsed = {}
    for selectionName in selectionNames:
        if selectionName not in list(d_totalSystDeltaXOverX.keys()):
            d_totalSystDeltaXOverX[selectionName] = {}
            d_preselectionSystScaleFactors[selectionName] = {}
            d_systematicSelectionsUsed[selectionName] = {}
        otherBkgTotalDeltaXOverXSqr = 0
        for background_name in background_names:
            if background_name not in list(d_preselectionSystScaleFactors.keys()):
                d_preselectionSystScaleFactors[selectionName][background_name] = {}
                d_systematicSelectionsUsed[selectionName][background_name] = {}
            thisBkgTotalSystOverNomSqr = 0
            thisBkgSystOverNom = 0
            for syst in systematicsNamesBackground["all"]:
                selectionNameSyst = selectionName
                if selectionName != "preselection" and selectionName != "trainingSelection":
                    massPoint = selectionName.replace("LQ", "")
                    if int(massPoint) > maxLQSelectionMass:
                        selectionNameSyst = "LQ"+str(maxLQSelectionMass)
                thisEntry, bkgSystDeltaOverNominalUp, bkgSystDeltaOverNominalDown, systNomYield, systSelection, preselRatios = d_backgroundSampleInfos[background_name].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                d_systematicSelectionsUsed[selectionName][background_name][syst] = systSelection
                if preselRatios is not None:
                    d_preselectionSystScaleFactors[selectionName][background_name][syst+"Up"] = preselRatios[0]
                    d_preselectionSystScaleFactors[selectionName][background_name][syst+"Down"] = preselRatios[1]
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
            fullSignalName, signalNameForFile = GetFullSignalName(signal_name, massPoint)
            if selectionName == "preselection" and selectionName != "trainingSelection":
                d_totalSystDeltaXOverX[selectionName][signalNameForFile] = 0.0
                continue
            thisSigTotalSystOverNomSqr = 0
            for syst in systematicsNamesSignal:
                thisEntry, sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown, _, _, _ = d_signalSampleInfos[signalNameForFile].GetSystematicEffectAbs("all", syst, selectionNameSyst, d_applicableSystematics)
                if thisEntry != "-":
                    thisSigSystOverNom = max(sigSystDeltaOverNominalUp, sigSystDeltaOverNominalDown)
                    print("INFO: selection={}, signal_name={}, syst={}, deltaX/X={}".format(selectionName, signalNameForFile, syst, thisSigSystOverNom))
                    thisSigTotalSystOverNomSqr += thisSigSystOverNom*thisSigSystOverNom
            print("INFO: selection={}, signal_name={}, total deltaX/X={}".format(selectionName, signalNameForFile, math.sqrt(thisSigTotalSystOverNomSqr)))
            d_totalSystDeltaXOverX[selectionName][signalNameForFile] = math.sqrt(thisSigTotalSystOverNomSqr)

    with open(systematics_dictFilePath.format(yearsRequestedStr), "w") as systematics_dictFile:
        systematics_dictFile.write(str(d_totalSystDeltaXOverX))
    print("systematics dict written to {}".format(systematics_dictFilePath.format(yearsRequestedStr)))

    with open(preselSystScaleFactorsJson.format(yearsRequestedStr), 'w') as fp:
        json.dump(d_preselectionSystScaleFactors, fp, indent=4)

    with open(systSelectionsUsedJson.format(yearsRequestedStr), 'w') as fp:
        json.dump(d_systematicSelectionsUsed, fp, indent=4)

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
            if not DoesSystematicApplyAnyYear(syst, sampleName, d_applicableSystematics):
                continue
            if not syst in d_totalBkgSysts.keys():
                d_totalBkgSysts[syst] = {}
            if not "totalYield" in d_totalBkgSysts.keys():
                d_totalBkgSysts["totalYield"] = {}
            tableRow = [d_systTitles["all"][syst]]
            systYieldRow = [d_systTitles["all"][syst]]
            systematicsByNameAndMass[d_systTitles["all"][syst]] = []
            systematicsSelectionsByNameAndMass[d_systTitles["all"][syst]] = []
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
                    # print("SIC DEBUG: for sample {}, syst {}, fill value (% deltaX/X) = {}".format(sampleName, syst, fillVal), flush=True)
                    # tableRow.append(fillVal)
                    if selection != "preselection" and selection != "trainingSelection":
                        mass = float(selection.replace("LQ", ""))
                        idxToFill = bisect(massRanges, mass) - 1
                        systHists[idxToFill].Fill(mass, fillVal)
                        valid = True
                if valid:
                    systematicsByNameAndMass[d_systTitles["all"][syst]].append(fillVal)
                    systematicsSelectionsByNameAndMass[d_systTitles["all"][syst]].append(mass)
                else:
                    systematicsByNameAndMass[d_systTitles["all"][syst]].append(0)
                    if selection != "preselection" and selection != "trainingSelection":
                        mass = float(selection.replace("LQ", ""))
                        systematicsSelectionsByNameAndMass[d_systTitles["all"][syst]].append(mass)
                    else:
                        systematicsSelectionsByNameAndMass[d_systTitles["all"][syst]].append(-1)
                tableRow.append(fillVal)
                # systYieldRow.append(systYieldVal+value)
                systYieldRow.append(systYieldVal)
                if sampleName in background_names or IsComponentBackground(sampleName):
                    toAdd = max(float(d_systUpDeltas[sampleName][syst][selection]), float(d_systDownDeltas[sampleName][syst][selection]))
                    rawEvts = d_backgroundSampleInfos[sampleName].GetUnscaledRate(selection)
                    val, err = d_backgroundSampleInfos[sampleName].GetRateAndErr(selection)
                if sampleName in background_names:
                    d_totalBkgSysts[syst][selection] += toAdd*toAdd
                    d_totalBkgSysts["totalYield"][selection] += value
                elif not IsComponentBackground(sampleName):
                    rawEvts = d_signalSampleInfos[sampleName].GetUnscaledRate(selection)
                    val, err = d_signalSampleInfos[sampleName].GetRateAndErr(selection)
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
                    # yieldRow.append("{:.4f} +/- {:.4f} [MC evts: {:.1f}]".format(value, err, thisBkgRawEvts))
                    yieldRow.append("{v:.2{c}} +/- {e:.2{c}}".format(v=val, e=err, c='e' if val < 1e-2 and math.fabs(val) > 0 else 'f'))
                    mcRow.append("{:.1f}".format(rawEvts))
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
                # print("DEBUG: compute max effect on sequence {} for sample={} and syst={}".format(effect, sampleName, systName))
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
    for syst in systematicsNamesBackground["all"]:
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
        # if "norm" in syst.lower():
        #     print("DEBUG: QCDFakes_DATA, syst={}, d_systNominals[{}][{}], total sum = {}".format(
        #         syst, "QCDFakes_DATA", syst, d_systNominals["QCDFakes_DATA"][syst], total))
    print(tabulate(totalBkgSystTable, headers=totalBkgSystsHeaders, tablefmt="fancy_grid", floatfmt=".2f"))

yieldTable = []
yieldTableHeaders = ["sample/sel."]
yieldTableHeaders.extend([selection for selection in selectionNames])
pctTable = []
pctTableHeaders = ["sample/sel. % of total"]
pctTableHeaders.extend([selection for selection in selectionNames])
sampleList = []
for sampleName in d_backgroundSampleInfos.keys():
    if sampleName in allYears or sampleName in sampleList:
        continue
    sampleList.append(sampleName)
    if sampleName in dictSampleComponents.keys():
        for componentBkg in dictSampleComponents[sampleName]:
            sampleList.append(componentBkg)
samplesAlwaysBelowOnePct = []
for sampleName in sampleList:
    yieldTableRow = [sampleName]
    pctTableRow = [sampleName]
    sampleBelowOnePct = True
    for selection in selectionNames:
        val, err = d_backgroundSampleInfos[sampleName].GetRateAndErr(selection)
        yieldTableRow.append("{v:.2{c}} +/- {e:.2{c}}".format(v=val, e=err, c='e' if val < 1e-2 and math.fabs(val) > 0 else 'f'))
        # total = sum([d_systNominals[sampleName][syst][selection] for sampleName in background_names])
        total = sum([d_backgroundSampleInfos[sampleName].GetRateAndErr(selection)[0] for sampleName in background_names])
        pctTableRow.append("{v:.2{c}} +/- {e:.2{c}}".format(v=100.*val/total, e=100.*err/total, c='f'))
        if val/total >= 0.01:
            sampleBelowOnePct = False
    if sampleBelowOnePct:
        samplesAlwaysBelowOnePct.append(sampleName)
    yieldTable.append(yieldTableRow)
    pctTable.append(pctTableRow)
print(tabulate(yieldTable, headers=yieldTableHeaders, tablefmt="fancy_grid", floatfmt=".2f"))
print(tabulate(pctTable, headers=pctTableHeaders, tablefmt="fancy_grid", floatfmt=".2f"))
print(tabulate(pctTable, headers=pctTableHeaders, tablefmt="latex", floatfmt=".2f"))
for sample in samplesAlwaysBelowOnePct:
    print("INFO: sample {} always below 1% of total background".format(sample))
print()

# make final selection tables
columnNames = [
    "MLQ",
    "signal"]
columnNames.extend(backgroundTitles)
columnNames.extend([
    "Total BG (stat) (syst)",
    "Data",
])
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
            thisSigEvts, thisSigEvtsErr = d_signalSampleInfos[signalNameForFile].GetRateAndErr(selectionName)
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
            thisBkgEvts, thisBkgEvtsErr = d_backgroundSampleInfos[background_name].GetRateAndErr(selectionName)
            # thisBkgEvtsErrUp, thisBkgEvtsErrDown = GetStatErrorsFromDatacard(d_datacardStatErrs[selectionName][background_name], thisBkgEvts)
            # print "GetStatErrorsFromDatacard for selection={}, background={}, thisBkgEvts={} + {} - {}".format(selectionName, background_name, thisBkgEvts, thisBkgEvtsErrUp, thisBkgEvtsErrDown)
            thisBkgEvtsErrUp = thisBkgEvtsErr
            thisBkgEvtsErrDown = thisBkgEvtsErr
            # thisBkgTotalEntries = d_background_unscaledRates[background_name][selectionName]
            thisBkgEffEntries = thisBkgEvts**2/thisBkgEvtsErr**2 if thisBkgEvtsErr != 0 else 0
            thisBkgSyst = 0
            if doSystematics:
                thisBkgSyst = d_backgroundSampleInfos[background_name].GetTotalSystDeltaOverNominal("all", selectionName, systematicsNamesBackground["all"], d_applicableSystematics)
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
                thisSigEvts, thisSigEvtsErrUp, thisSigEvtsErrDown, addDecimalsUntilNonzero=True
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
                thisSigEvts, thisSigEvtsErrUp, thisSigEvtsErrDown, addDecimalsUntilNonzero=True, latex=True
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
                latex=True,
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
    table_file.write(r"\hline")
    table_file.write("\n")
    ending = r"\end{tabular}}"  # extra } to end resizebox
    print(ending)
    table_file.write(ending+"\n")
    print()
    # table_file.write("\n")

with open(tablesDir + "/eventYieldsPaper.tex", "w") as table_file:
    print()
    print("Latex table: Paper")
    print()
    # latex table -- Paper
    for line in latexRowsPaper:
        print (line)
        table_file.write(line+"\n")
    print()
    # table_file.write("\n")

# acc * eff
plotSignalEfficiencyTimesAcceptance(dataMC_filepath, signalNameTemplate, [int(m) for m in mass_points], yearsRequestedStr)
print("datacard written to {}".format(datacard_filePath))
