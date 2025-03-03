#!/usr/bin/env python

import math
import copy
import os
from tabulate import tabulate
from ROOT import TFile, TH1


def GetHistForSample(fileName, sampleName):
    histName = histNameBase.format(sampleName)
    tfile = TFile.Open(fileName)
    hist = tfile.Get(histName)
    if not hist or hist is None:
        raise RuntimeError("Coudn't get hist {} from file {}".format(histName, filepath))
    hist = copy.deepcopy(hist)
    tfile.Close()
    return hist


def FormatFloat(val):
    return "{v:.2{c}}".format(v=val, c='e' if val < 1e-2 else 'f')


def GetReweightingFactor(hist, sampleName):
    verbose = False
    if verbose:
        print("INFO: GetReweightingFactor() - Going to attempt to find reweighting factor")
    reweight = 1.0
    foundReweight = False
    foundReweightBin = -1
    for xbin in range(0, hist.GetNbinsX()+1):
        binLabel = hist.GetXaxis().GetBinLabel(xbin)
        binContent = hist.GetBinContent(xbin)
        if "reweighting" in binLabel.lower():
            if binContent != 1.0:
                reweight = binContent
                if foundReweight:
                    raise RuntimeError("For hist {}, already found a reweighting factor from bin {} with label {} of value {}, but found another from bin {} with label {} of value {}. Cannot continue.".format(
                        hist.GetName(), foundReweightBin, hist.GetXaxis().GetBinLabel(foundReweightBin), reweight, xbin, binLabel, binContent
                        ))
                foundReweightBin = xbin
            foundReweight = True
    if not foundReweight:
        print("bin labels for hist=", list(hist.GetXaxis().GetLabels()))
        raise RuntimeError("WARN: Did not find reweighting factor from histogram")
    elif foundReweightBin == -1:
        print("WARN: For sample {}, only found a reweighting bin with content 1 (which might be OK)".format(sampleName))
    if verbose:
        print("INFO: GetReweightingFactor() - Found reweight factor of {} in bin {} with label {}".format(reweight, foundReweightBin, hist.GetXaxis().GetBinLabel(foundReweightBin)), flush=True)
    return reweight, foundReweightBin


def SumSamples(filepath, qcdfilepath, year, sampleNames):
    histTotal = None
    for sample in sampleNames:
        # print("INFO: SumSamples: consider sample={}".format(sample), flush=True)
        if "QCDFakes_DATA" in sample:
            histSample = GetHistForSample(qcdfilepath.format(year), sample)
        else:
            histSample = GetHistForSample(filepath.format(year, sample), sample)
        if histTotal is None:
            histTotal = histSample
        else:
            # print("INFO: SumSamples: Adding histo for sample {} to total".format(sample), flush=True)
            histTotal.Add(histSample)
    return histTotal


def MakeTables(sampleName, hist, txtFile, ltxFile):
    verbose = False
    hist2 = hist.Clone()
    hist2.Reset()
    for xbin in range(1, hist.GetNbinsX()):
        binContent = hist.GetBinContent(xbin+1)
        binLabel = hist.GetXaxis().GetBinLabel(xbin+1)
        hist2.SetBinContent(xbin, binContent)
        hist2.GetXaxis().SetBinLabel(xbin, binLabel)

    if verbose:
        for xbin in range(1, 10):
            print("INFO: pre division, numerator: bin={}, label={}, content={}, error={}".format(xbin, hist2.GetXaxis().GetBinLabel(xbin), hist2.GetBinContent(xbin), hist2.GetBinError(xbin)), flush=True)
    
    hist2.Divide(hist2, hist, 1, 1, "B")
    
    if verbose:
        for xbin in range(1, 10):
            print("INFO: pre-reweighting: division bin={}, label={}, content={}, error={}".format(xbin, hist2.GetXaxis().GetBinLabel(xbin), hist2.GetBinContent(xbin), hist2.GetBinError(xbin)), flush=True)
    
    totalEvents = hist.GetBinContent(1)
    if reweight:
        factor, reweightBin = GetReweightingFactor(hist2, sampleName)
        totalEvents *= factor
        for xbin in range(1, reweightBin+1):
            binContent = hist.GetBinContent(xbin) * factor
            hist.SetBinContent(xbin, binContent)
    
    if verbose:
        for xbin in range(1, 10):
            print("INFO: post-reweighting division bin={}, label={}, content={}, error={}".format(xbin, hist2.GetXaxis().GetBinLabel(xbin), hist2.GetBinContent(xbin), hist2.GetBinError(xbin)), flush=True)
    
    # make tables
    colNames = ["cutName", "Npass", "rel. eff. (%)", "abs. eff. (%)"]
    table = []
    cutNamesSeen = []
    for xbin in range(1, hist2.GetNbinsX()):
        cutName = hist.GetXaxis().GetBinLabel(xbin+1)
        nPass = hist.GetBinContent(xbin+1)
        relEff = hist2.GetBinContent(xbin) * 100.0
        absEff = nPass / totalEvents * 100.0
        # skip over cuts that were previously applied and have no effect here
        if relEff == 100.0 and cutName in cutNamesSeen:
            continue
        # skip reweighting cuts
        if reweight and "reweighting" in cutName.lower():
            continue
        # don't show the MejjLQXXX cuts
        if "meejj" in cutName.lower() and "lq" in cutName.lower():
            continue
        # for some reason, we have two empty bins at the end of these hists
        if cutName == "":
            continue
        nPass = FormatFloat(nPass)
        relEff = FormatFloat(relEff)
        absEff = FormatFloat(absEff)
        table.append([cutName, nPass, relEff, absEff])
        cutNamesSeen.append(cutName)
    print(tabulate(table, headers=colNames, tablefmt="fancy_grid", disable_numparse=True), file=txtFile)
    print(tabulate(table, headers=colNames, tablefmt="latex", disable_numparse=True), file=ltxFile)
    

###################################################################################################
# RUN
###################################################################################################
reweight = True
# includeQCD = True
allYears = ["2016preVFP", "2016postVFP", "2017", "2018"]
signalMasses = [500, 1000, 1500]
signalNameTemplate = "LQToDEle_M-{}_pair"
histNameBase = "histo1D__{}__EventsPassingCutsAllHist"
cutflowDir = "cutflows"
txtFileName = cutflowDir + "/cutflow_{}_{}.txt"
ltxFileName = cutflowDir + "/cutflow_{}_{}.tex"

filepathBase = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/{}/eejj_19feb2025_bdt_LQToDEle/output_cutTable_lq_eejj_BDT/analysisClass_lq_eejj_{}_plots.root"
qcdFilepathBase = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/{}/qcd_eejj_19feb2025_bdt_LQToDEle/output_qcdSubtractedYield/qcdSubtracted_plots.root"
bkgSamples = ["ZJet_amcatnlo_ptBinned_IncStitch", "TTTo2L2Nu", "OTHERBKG_dibosonNLO_singleTop"]
# FIXME: how to include QCD
# if includeQCD:
#     bkgSamples.append("QCDFakes_DATA")

if not os.path.exists(cutflowDir):
    os.makedirs(cutflowDir)

totalBkgHist = None
totalSigHists = [None]*len(signalMasses)
for year in allYears:
    bkgHist = SumSamples(filepathBase, qcdFilepathBase, year, bkgSamples)
    with open(txtFileName.format("allBkg", year), "w") as txtFile:
        with open(ltxFileName.format("allBkg", year), "w") as ltxFile:
            MakeTables("All background", bkgHist, txtFile, ltxFile)
    if totalBkgHist is None:
        totalBkgHist = bkgHist
    else:
        totalBkgHist.Add(bkgHist)
    for idx, mass in enumerate(signalMasses):
        sample = signalNameTemplate.format(mass)
        histSample = GetHistForSample(filepathBase.format(year, sample), sample)
        with open(txtFileName.format(sample, year), "w") as txtFile:
            with open(ltxFileName.format(sample, year), "w") as ltxFile:
                MakeTables(sample, histSample, txtFile, ltxFile)
        if totalSigHists[idx] is None:
            totalSigHists[idx] = histSample
        else:
            totalSigHists[idx].Add(histSample)
with open(txtFileName.format("allBkg", "RunII"), "w") as txtFile:
    with open(ltxFileName.format("allBkg", "RunII"), "w") as ltxFile:
        MakeTables("All background", totalBkgHist, txtFile, ltxFile)
for idx, mass in enumerate(signalMasses):
    sample = signalNameTemplate.format(mass)
    with open(txtFileName.format(sample, "RunII"), "w") as txtFile:
        with open(ltxFileName.format(sample, "RunII"), "w") as ltxFile:
            MakeTables(sample, totalSigHists[idx], txtFile, ltxFile)

print("Cutflows written into {}.".format(cutflowDir))
