from __future__ import print_function

from ROOT import (
    gROOT,
    gPad,
    TFile,
    TH1D,
    TGraphErrors,
    TCanvas,
    TColor,
    TLegend,
    kRed,
    kBlue,
    kCyan,
    kGreen,
    kAzure,
    kOrange,
    kGray,
    kViolet,
    kBlack,
    THStack,
    TGraph,
    TRatioPlot,
    TPad,
    TLine,
    TBox,
)
import ROOT
import numpy as np
import copy
import ctypes
import sys
import os
import math

minMLQ = int(sys.argv[1])
maxMLQ = int(sys.argv[2])
#year = str(sys.argv[3])
parameterized = False
#parameterized = True
folderName = ""#str(minMLQ)+"To"+str(maxMLQ)+"GeV"

base_folder = os.getenv("LQDATAEOS")+"/BDT_7FebSkim/LQToDEle"#{}".format(year)
optimizationFile = base_folder+"/optimizationPlots.root"
bdtPlotFile = base_folder+"/bdtPlots.root"

modelName = "LQToDEle"#_{}".format(year)

pdf_folder = base_folder+"/plots"
outFileName = base_folder+"/"+modelName+"Plots.root"
plotsForAN = pdf_folder+"/plotsForAN_{}.pdf".format(modelName)

gROOT.SetBatch(True)

if not os.path.isdir(pdf_folder):
    os.mkdir(pdf_folder)

outFile = TFile(outFileName,"recreate")
variables = ["sT_eejj", "PFMET_Type1_Pt", "M_e1e2", "M_e1j1", "M_e1j2", "M_e2j1", "M_e2j2", "Ele1_Pt", "Ele2_Pt", "MejMin", "MejMax", "Meejj","Jet1_Pt", "Jet2_Pt", "Jet3_Pt"]

bkgSamples = [
    "QCD",
    "SingleTop",
    "DIBOSON_nlo",
    "TTbar_powheg",
    "ZJet_amcatnlo_ptBinned",
]

bkgColors = [
    kCyan,
    kGreen+1,
    kBlue-7,
    kBlue,
    kRed,
]

LQmasses = list(range(minMLQ,maxMLQ+1,100))
#LQmasses = [1300, 2000, 2200, 2300]
print(LQmasses)
for mass in LQmasses:
    path = pdf_folder+"/"+str(mass)
    if not os.path.isdir(path):
        os.mkdir(path)

optTFile = TFile.Open(optimizationFile)
#get cut values
cutVsMLQPlot = optTFile.Get("optCutValVsLQMass")
cutValues = {}
for i in range(cutVsMLQPlot.GetN()):
    #This plot does not necessarily get filled in order of increasing mass
    mass = int(cutVsMLQPlot.GetPointX(i))
    cutValues[str(mass)] = cutVsMLQPlot.GetPointY(i)
c1 = TCanvas()
c1.SetGridy()
cutVsMLQPlot.GetYaxis().SetTitle("BDT cut")
cutVsMLQPlot.GetYaxis().SetRangeUser(-1,1)
cutVsMLQPlot.Draw("AP")
outFile.cd()
cutVsMLQPlot.Write()
c1.Print(pdf_folder+"/optCutValVsLQMass.pdf")
c1.Print(plotsForAN+"(","pdf")
#print(cutValues)

#plot max FOM vs. LQ mass
FOMVsMLQPlot = TGraph(len(LQmasses))
for i, mass in enumerate(LQmasses):
    #print("get fom plot ", "LQM"+str(mass)+"/fomValVsBDTCutGraphLQM"+str(mass))
    fomPlot = optTFile.Get("LQM"+str(mass)+"/fomValVsBDTCutGraphLQM"+str(mass))
    foms = fomPlot.GetY()
    cutVal = cutValues[str(mass)]
    index = int((cutVal - (-1))/0.02) #100 cuts between -1 and 1
    fomUsed = foms[index]
    FOMVsMLQPlot.SetPoint(i,mass,fomUsed)
    '''
    maxFOM = max(foms)
    if (mass<1000):
        FOMVsMLQPlot.SetPoint(i,mass,maxFOM)
    else:
        cutVal = cutValues[str(mass)]
        for j in range(fomPlot.GetN()):
            x = fomPlot.GetPointX(j)
            if (cutVal-0.01 < x and x < cutVal+0.01):
                fom = fomPlot.GetPointY(j)
                FOMVsMLQPlot.SetPoint(i,mass,fom)
    '''
    c = TCanvas()
    fomPlot.GetXaxis().SetTitle("BDT cut")
    fomPlot.Draw()
    #c.SetLogy()
    c.SetGridy()
    c.Print(pdf_folder+"/{}/FOMvsBDTCut.pdf".format(mass))
    if i==0:
        c.Print(pdf_folder+"/FOMPlots.pdf(","pdf")
    elif mass == LQmasses[-1]:
        c.Print(pdf_folder+"/FOMPlots.pdf)","pdf")
    else:
        c.Print(pdf_folder+"/FOMPlots.pdf","pdf")
    if mass == 2300:
        ymin = fomPlot.GetYaxis().GetBinLowEdge(fomPlot.GetYaxis().GetFirst())
        ymax = fomPlot.GetYaxis().GetBinUpEdge(fomPlot.GetYaxis().GetLast())
        lineMaxFOM = TLine(0.3804, ymin, 0.3804, ymax)
        lineMaxCut = TLine(0.6006, ymin, 0.6006, ymax)
        lineMaxFOM.SetLineColor(kBlue)
        lineMaxCut.SetLineColor(kRed)
        ltemp = TLegend(0.2,0.7,0.6,0.9)
        ltemp.AddEntry(lineMaxFOM, "BDT cut with max FOM and nB > 0.5. (FOM = 0.0283)","l")
        ltemp.AddEntry(lineMaxCut, "Max BDT cut with nB > 0.5 (FOM = 0.0277)", "l")
        ltemp.AddEntry(fomPlot, "Figure of merit", "p")
        lineMaxFOM.Draw('same')
        lineMaxCut.Draw('same')
        ltemp.Draw('same')
        fomPlot.Draw('psame')
        c.Print(pdf_folder+"/{}/FOM_withCuts_forAN.pdf".format(mass))


c2 = TCanvas()
c2.SetGridy()
c2.SetLogy()
FOMVsMLQPlot.SetName("maxFOM_vs_MLQ")
FOMVsMLQPlot.SetTitle("FOM of chosen BDT cut vs MLQ")
FOMVsMLQPlot.GetXaxis().SetTitle("MLQ (GeV)")
FOMVsMLQPlot.GetYaxis().SetTitle("FOM")
FOMVsMLQPlot.SetMarkerStyle(8)
FOMVsMLQPlot.Draw("ALP")
outFile.cd()
FOMVsMLQPlot.Write()
c2.Print(pdf_folder+"/maxFOMvsMLQ.pdf")
c2.Print(plotsForAN,"pdf")

#Plot Ns and Nb vs BDT cut for each mass
for mass in LQmasses:
    c = TCanvas()
    c.SetLogy()
    c.SetGridy()
    NsPlot = optTFile.Get("LQM"+str(mass)+"/nSigVsBDTCutGraphLQM"+str(mass))
    NsPlot.GetYaxis().SetTitle("")
    NsPlot.SetTitle("")
    NbPlot = optTFile.Get("LQM"+str(mass)+"/nBkgVsBDTCutGraphLQM"+str(mass))
    #print("use Nb plot ", NbPlot.GetName())
    NsList = np.array(NsPlot.GetY())
    NbList = np.array(NbPlot.GetY())
    NsMin = np.min(NsList[np.nonzero(NsList)])
    NbMin = np.min(NbList[np.nonzero(NbList)])
    plotMin = 0.1*min(NsMin,NbMin)
    NsMax = max(NsList)
    NbMax = max(NbList)
    plotMax = 2*max(NsMax,NbMax)
    NbPlot.SetMarkerColor(kRed)
    NbPlot.SetLineColor(kRed)
    NsPlot.SetLineColor(kBlue)
    NbPlot.GetYaxis().SetTitle("")
    NbPlot.SetTitle("Background Yield vs BDT Cut")
    #NbPlot.GetYaxis().SetRangeUser(plotMin,plotMax)
    NbPlot.GetXaxis().SetTitle("BDT cut")
    NsPlot.GetYaxis().SetRangeUser(plotMin,plotMax)
    NbPlot.SetMarkerSize(0.5)
    NsPlot.SetMarkerSize(0.5)
    NbPlot.SetMarkerColor(kBlack)
    NbPlot.SetLineColor(kBlack)
    NbPlot.Draw("AP")
    #NsPlot.Draw("Psame")
    l = TLegend(0.9,0.8,0.99,0.9)
    #l.AddEntry(NsPlot,"signal","p")
    l.AddEntry(NbPlot,"background yield", "p")
    #l.Draw("same")
    optCut = cutValues[str(mass)]
    line = TLine(optCut,NbPlot.GetYaxis().GetBinLowEdge(NbPlot.GetYaxis().GetFirst()),optCut,NbPlot.GetYaxis().GetBinUpEdge(NbPlot.GetYaxis().GetLast()))
    line.SetLineColor(kRed)
    line.Draw("same")
    if mass==2300:
        ltemp = TLegend(0.2,0.7,0.6,0.9)
        lineMaxFom = TLine(0.3804,NbPlot.GetYaxis().GetBinLowEdge(NbPlot.GetYaxis().GetFirst()),0.3804,NbPlot.GetYaxis().GetBinUpEdge(NbPlot.GetYaxis().GetLast()))
        lineMaxFom.SetLineColor(kBlue)
        lineMaxFom.Draw('same')
        ltemp.AddEntry(lineMaxFom, "BDT cut with max FOM and min nB > 0.5 (nB = 0.5087)","l")
        ltemp.AddEntry(line, "Max BDT cut with nB > 0.5 (nB = 0.5327)","l")
        ltemp.AddEntry(NbPlot, "background yield", "p")
        ltemp.Draw("same")
    NbPlot.Draw('psame')
    c.Print(pdf_folder+"/"+str(mass)+"/NBvsBDTCutMLQ"+str(mass)+".pdf")

    cEff = TCanvas()
    cEff.SetGridy()
    line = TLine(optCut, -0.1, optCut, 1.1)
    line.SetLineColor(kAzure+8)
    #signal efficiency and bkg rejection vs BDT cut
    sigEff = copy.deepcopy(NsPlot)
    for ipoint in range(sigEff.GetN()):
        Ns = sigEff.GetPointY(ipoint)
        BDTCut = sigEff.GetPointX(ipoint)
        sigEff.SetPoint(ipoint, BDTCut, Ns/NsMax)
    bkgRej = copy.deepcopy(NbPlot) #I'm going to overwrite all the points but this way I keep all the settings
    bkgRej.SetTitle("background rejection and signal efficiency for MLQ = "+str(mass))
    for ipoint in range(bkgRej.GetN()):
        Nb = bkgRej.GetPointY(ipoint)
        BDTCut = bkgRej.GetPointX(ipoint)
        bkgRej.SetPoint(ipoint, BDTCut, 1-Nb/NbMax)
    bkgRej.GetYaxis().SetRangeUser(-0.1,1.1)
    bkgRej.Draw("AXP")
    sigEff.Draw("XP")
    l.Draw("same")
    line.Draw("same")
    cEff.Print(pdf_folder+"/"+str(mass)+"/EffVsBDTCutMLQ"+str(mass)+".pdf")

#signal and background BDT output for each mass
c4 = TCanvas()
c4.SetLogy()
c4.SetGridy()
yieldAtCutSig = optTFile.Get("optNSVsLQMass")
yieldAtCutBkg = optTFile.Get("optNBVsLQMass")
totSigPlot = optTFile.Get("nSBeforeBDTCutVsLQMass")
totBkgPlot = optTFile.Get("nBBeforeBDTCutVsLQMass")
#The TGraphs don't get filled in order of increasing mass, so I can't just loop over all the masses
sigVsMLQ = {}
sigVsMLQ["yield"] = {}
sigVsMLQ["total"] = {}
bkgVsMLQ = {}
bkgVsMLQ["total"] = totBkgPlot.GetPointY(1) #all the same so it doesn't matter which one I get
for ipoint in range(yieldAtCutSig.GetN()):
    sigMass = int(yieldAtCutSig.GetPointX(ipoint))    
    sigVsMLQ["yield"][str(sigMass)] = yieldAtCutSig.GetPointY(ipoint)

    totSigMass = int(totSigPlot.GetPointX(ipoint))
    sigVsMLQ["total"][str(totSigMass)] = totSigPlot.GetPointY(ipoint)

    bkgMass = int(yieldAtCutBkg.GetPointX(ipoint))
    bkgVsMLQ[str(bkgMass)] = yieldAtCutBkg.GetPointY(ipoint)

sigEfficiency = TGraphErrors(len(LQmasses))
bkgRejection = TGraphErrors(len(LQmasses))
sigEffVsBkgRejMLQ = TGraph(len(LQmasses))
rocAreas = {}
for i,mass in enumerate(LQmasses):
    cBDTOutput = TCanvas()
    cBDTOutput.SetLogy()
    cBDTOutput.SetGridy()
    bkgPlot = optTFile.Get("LQM"+str(mass)+"/BDTOutputTotalBackgroundLQM"+str(mass))
    sigPlot = optTFile.Get("LQM"+str(mass)+"/hsig_BDTG_"+str(mass))
    #if mass==3000:
    #    bkgPlot.SetMarkerStyle(8)
     #   bkgPlot.SetMarkerSize(0.6)
     #   sigPlot.SetMarkerStyle(8)
     #   sigPlot.SetMarkerSize(0.6)
    #print("got hist ", sigPlot.GetName(), " for LQM ", mass)
    bkgPlotRebin = copy.deepcopy(bkgPlot)
    sigPlotRebin = copy.deepcopy(sigPlot)
    #bkgPlotRebin.Rebin(4)
    #sigPlotRebin.Rebin(4)
    plotMax = max(sigPlotRebin.GetMaximum(),bkgPlotRebin.GetMaximum())
    plotMin = min(sigPlotRebin.GetMinimum(0), bkgPlotRebin.GetMinimum(0))
    cut = cutValues[str(mass)]
    bkgPlotRebin.SetLineColor(kRed)
    bkgPlotRebin.SetMarkerColor(kRed)
    sigPlotRebin.SetLineColor(kBlue)
    sigPlotRebin.SetMarkerColor(kBlue)
    line = TLine(cut, 0, cut, plotMax)
    line.SetLineColor(kAzure+8)
    l = TLegend(0.4,0.7,0.6,0.9)
    l.AddEntry(sigPlotRebin,"BDT output on signal", "lp")
    l.AddEntry(bkgPlotRebin, "BDT output on background", "lp")
    l.AddEntry(line, "opt. cut value", "l")
    bkgPlotRebin.GetYaxis().SetRangeUser(0.2*plotMin,5*plotMax)
    #bkgPlotRebin.GetYaxis().SetRangeUser(-0.2,1)
    bkgPlotRebin.SetStats(0)
    bkgPlotRebin.GetXaxis().SetTitle("BDT output")
    bkgPlotRebin.GetXaxis().SetRangeUser(-1.1,1.1)
    bkgPlotRebin.SetTitle("BDT output MLQ = {} GeV".format(mass))
    sigPlotRebin.SetName("sigBDTOutput_MLQ"+str(mass))
    bkgPlotRebin.SetName("bkgBDTOutput_MLQ"+str(mass))
    outFile.cd()
    sigPlotRebin.Write()
    bkgPlotRebin.Write()
    bkgPlotRebin.Draw()
    sigPlotRebin.Draw("same")
    line.Draw("same")
    l.Draw("same")
    bkgPlotRebin.Draw("Same")
    sigPlotRebin.Draw("same")
    cBDTOutput.Print(pdf_folder+"/"+str(mass)+"/BDTOutputMLQ"+str(mass)+".pdf")
    if mass == LQmasses[0]:
        cBDTOutput.Print(pdf_folder+"/allBDTOutputs.pdf(","pdf")
    elif mass==LQmasses[-1]:
        cBDTOutput.Print(pdf_folder+"/allBDTOutputs.pdf)","pdf")
    else:
        cBDTOutput.Print(pdf_folder+"/allBDTOutputs.pdf","pdf")

    if mass==500 or mass==1500 or mass==2500:
        cBDTOutput.Print(plotsForAN,"pdf")
    #cBDTOutput.Print(pdf_folder+"/"+str(mass)+"/BDTOutputZoomMLQ"+str(mass)+".pdf")    

    #Make an ROC curve/auc by hand
    nBins = bkgPlot.GetNbinsX()
    ROCCurve = TGraph()
    ROCCurve.SetPointX(0,1.0)
    ROCCurve.SetPointY(0,0.0)
    for ibin in range(nBins): #zero negative bins out
        nBThisBin = bkgPlot.GetBinContent(ibin+1)
        if nBThisBin < 0:
            bkgPlot.SetBinContent(ibin+1,0)
    nSTot = sigPlot.Integral(1,bkgPlot.GetXaxis().GetLast())
    nBTot = bkgPlot.Integral(1,sigPlot.GetXaxis().GetLast())
    for ibin in range(nBins):
        nB = bkgPlot.Integral(ibin+1, bkgPlot.GetXaxis().GetLast())
        nS = sigPlot.Integral(ibin+1, sigPlot.GetXaxis().GetLast())
        if nB > 0:
            Y = 1.0 - nB/nBTot
        else:
            Y = 1.0
        X = nS/nSTot
        #print("MLQ{}: ({}, {})".format(mass,X,Y))
        ROCCurve.SetPointX(ibin,X)
        ROCCurve.SetPointY(ibin,Y)
        #print("X = {}/{} = {}".format(nS,nSTot,X))
        #print("Y = 1 - {}/{} = {}".format(nB,nBTot,Y))

    ROCCurve.SetPointX(nBins,0.0)
    ROCCurve.SetPointY(nBins,1.0)
    ROCToIntegrate = copy.deepcopy(ROCCurve)
    ROCToIntegrate.AddPoint(0,0)
    auc = ROCToIntegrate.Integral()
    rocAreas[str(mass)] = auc
    ROCCurve.SetTitle("ROC for MLQ = {} GeV, AUC = {}".format(mass,auc))
    ROCCurve.GetXaxis().SetTitle("signal efficiency")
    ROCCurve.GetYaxis().SetTitle("background rejection")
    c = TCanvas()
    c.SetGridx()
    c.SetGridy()
    ROCCurve.SetLineColor(kRed)
    ROCCurve.SetLineWidth(2)
    ROCCurve.Draw("AL")
    c.Print(pdf_folder+"/"+str(mass)+"/ROCFromBDTPlot"+str(mass)+".pdf")
    if mass == LQmasses[0]:
        c.Print(pdf_folder+"/allROCs.pdf(","pdf")
    elif mass == LQmasses[-1]:
        c.Print(pdf_folder+"/allROCs.pdf)","pdf")
    else:
        c.Print(pdf_folder+"/allROCs.pdf","pdf")

    if mass==500 or mass==1500 or mass==2500:
        c.Print(plotsForAN,"pdf")

    c = TCanvas()
    c.SetGridy()
    bkgPlotNormalized = bkgPlotRebin.DrawNormalized()
    sigPlotNormalized = sigPlotRebin.DrawNormalized("same")
    outFile.cd()
    sigPlotNormalized.SetName("sigNormalizedBDTOutput_MLQ"+str(mass))
    bkgPlotNormalized.SetName("bkgNormalizedBDTOutput_MLQ"+str(mass))
    sigPlotNormalized.Write()
    bkgPlotNormalized.Write()
    #plotMax = 1.2 * max(sigPlotNormalized.GetBinContent(sigPlotNormalized.GetMaximumBin()), bkgPlotNormalized.GetBinContent(bkgPlotNormalized.GetMaximumBin()))
    plotMax = 1
    bkgPlotNormalized.GetYaxis().SetRangeUser(0,plotMax)
    line.Draw("same")
    l.Draw("same")
    c.Print(pdf_folder+"/"+str(mass)+"/normalizedBDTOutputMLQ"+str(mass)+".pdf")
    '''    
    c.SetLogy()
    stack = THStack()
    l = TLegend(0.3,0.7,0.7,0.9)
    #combine QCD into a single hist
    qcdHist = optTFile.Get("LQM{}/BDTOutputQCDFakes_DYJLQM{}".format(mass,mass))
    qcd1FRHist = optTFile.Get("LQM{}/BDTOutputQCDFakes_DATALQM{}".format(mass,mass))
    qcd2FRHist = optTFile.Get("LQM{}/BDTOutputQCDFakes_DATA_2FRLQM{}".format(mass,mass))
    qcdHist.Add(qcd1FRHist)
    qcdHist.Add(qcd2FRHist)
    for i, sample in enumerate(bkgSamples):
        if "QCD" in sample:
            hist = qcdHist
        else:
            print("get hist LQM{}/BDTOutput{}LQM{}".format(mass, sample, mass))
            histName = "LQM{}/BDTOutput{}LQM{}".format(mass, sample, mass)
            hist = optTFile.Get(histName)
        if "QCD" in sample or "ZJet" in sample:
            nBins = hist.GetNbinsX()
            for ibin in range(1, nBins+1):
                if hist.GetBinContent(ibin) < 0:
                    hist.SetBinContent(ibin,0)
        hist.SetMinimum(1e6)
        hist.SetLineColor(bkgColors[i])
        hist.SetFillColor(bkgColors[i])
        l.AddEntry(hist, sample, "l")
        stack.Add(hist)
    sigPlot.SetLineColor(kBlack)
    sigPlot.SetMarkerColor(kBlack)
    sigPlot.SetStats(0)
    sigPlot.GetYaxis().SetRangeUser(1e-4, 1e4)
    sigPlot.SetTitle("BDT output for MLQ = {}".format(mass))
    line.SetLineColor(kGray+2)
    l.AddEntry(sigPlot, "signal", "l")
    l.AddEntry(line, "Cut value", "l")
    sigPlot.Draw()
    stack.Draw("histSame")
    sigPlot.Draw("SAME")
    line.Draw("same")
    #l.Draw("SAME")
    gPad.RedrawAxis()
    gPad.RedrawAxis("G")
    l.Draw("SAME")
    c.Print(pdf_folder+"/"+str(mass)+"/BDTOutputStack.pdf")
    '''
    #signal and bkg efficiencies at opt cut value vs. MLQ
    print(sigVsMLQ["yield"][str(mass)], sigVsMLQ["total"][str(mass)])
    sigEff = sigVsMLQ["yield"][str(mass)] / sigVsMLQ["total"][str(mass)]
    sigEfficiency.SetPoint(i,mass,sigEff)

    bkgRej = 1 - bkgVsMLQ[str(mass)] / bkgVsMLQ["total"]
    bkgRejection.SetPoint(i,mass,bkgRej)
    
    sigEffVsBkgRejMLQ.SetPoint(i,bkgRej,sigEff)
c = TCanvas()
c.SetGridy()
sigEffVsBkgRejMLQ.SetTitle("BDT cut efficiency vs bkg rejection")
sigEffVsBkgRejMLQ.SetName("effVsRej")
sigEffVsBkgRejMLQ.SetMarkerStyle(8)
sigEffVsBkgRejMLQ.SetMarkerSize(0.75)
sigEffVsBkgRejMLQ.GetYaxis().SetTitle("cut efficiency")
sigEffVsBkgRejMLQ.GetXaxis().SetTitle("bkg rejection")
sigEffVsBkgRejMLQ.Draw("ALP")
outFile.cd()
sigEffVsBkgRejMLQ.Write()
c.Print(pdf_folder+"/effVsRej.pdf")

cYield = TCanvas()
cYield.SetGridy()
cYield.SetLogy()
yieldAtCutSig.SetMarkerStyle(8)
yieldAtCutSig.SetMarkerSize(0.75)
yieldAtCutSig.GetXaxis().SetTitle("MLQ")
yieldAtCutSig.SetTitle("yield at opt. cut value vs LQ mass")
yieldAtCutSig.SetName("SigYieldAtOptCut")
yieldAtCutSig.SetMarkerColor(kBlue)
yieldAtCutBkg.SetMarkerStyle(8)
yieldAtCutBkg.SetMarkerSize(0.75)
yieldAtCutBkg.SetMarkerColor(kRed)
yieldAtCutBkg.SetName("BkgYieldAtOptCut")
yieldAtCutSig.GetYaxis().SetRangeUser(1e-4,50000)
yieldAtCutSig.GetXaxis().SetLimits(200,3100)
yieldAtCutBkg.GetXaxis().SetLimits(200,3100)
l = TLegend(0.7,0.8,0.9,0.9)
l.AddEntry(yieldAtCutSig,"signal","p")
l.AddEntry(yieldAtCutBkg,"background","p")
yieldAtCutSig.Draw("AP")
yieldAtCutBkg.Draw("P")
l.Draw("same")
outFile.cd()
yieldAtCutSig.Write()
yieldAtCutBkg.Write()
cYield.Print(pdf_folder+"/yieldAtOptCutvsMLQ.pdf")
cYield.Print(plotsForAN,"pdf")

cEff = TCanvas()
cEff.SetGridy()
sigEfficiency.SetMarkerStyle(8)
sigEfficiency.SetTitle("signal efficiency at opt cut vs MLQ")
sigEfficiency.GetYaxis().SetRangeUser(0.6,1.01)
sigEfficiency.GetXaxis().SetTitle("MLQ GeV")
sigEfficiency.GetXaxis().SetRangeUser(250,3050)
sigEfficiency.SetName("sigEffAtOptCutVsMLQ")
sigEfficiency.Draw("ALP")
outFile.cd()
sigEfficiency.Write()
cEff.Print(pdf_folder+"/sigEffVsMLQ.pdf")
bkgRejection.SetMarkerStyle(8)
bkgRejection.SetTitle("background rejection at opt cut vs MLQ")
bkgRejection.GetXaxis().SetTitle("MLQ GeV")
bkgRejection.GetXaxis().SetRangeUser(250,3050)
bkgRejection.SetName("bkgRejAtOptCutVsMLQ")
bkgRejection.Draw("ALP")
outFile.cd()
bkgRejection.Write()
cEff.Print(pdf_folder+"/bkgRejVsMLQ.pdf")

#bdtTFile = TFile.Open(bdtPlotFile)
#ROC curve and integral
rocaucVsMLQPlot = TGraph(len(LQmasses))
rocaucVsMLQPlot.SetName("rocaucVsMLQ")
oneMaucVsMLQPlot = TGraph(len(LQmasses))
oneMaucVsMLQPlot.SetName("1-rocaucVsMLQ")
for i, mass in enumerate(LQmasses):
    c3 = TCanvas()
#    rocPlot = bdtTFile.Get("LQM"+str(mass)+"/Graph")
    #plotForIntegral = copy.deepcopy(rocPlot)
    #plotForIntegral.AddPoint(0,0)
#    rocauc = plotForIntegral.Integral()
    rocauc = rocAreas[str(mass)]
    oneMauc = 1 - rocauc
    title = "ROC for MLQ = "+str(mass)+" (auc = "+str(rocauc)+")"
#    rocPlot.SetTitle(title)
#    rocPlot.SetName("ROC_MLQ"+str(mass))
#    rocPlot.Draw()
    outFile.cd()
#    rocPlot.Write()
#    c3.Print(pdf_folder+"/"+str(mass)+"/rocMLQ"+str(mass)+".pdf")
    rocaucVsMLQPlot.SetPoint(i,mass,rocauc)
    oneMaucVsMLQPlot.SetPoint(i,mass,oneMauc)
rocaucVsMLQPlot.SetMarkerStyle(8)
rocaucVsMLQPlot.GetXaxis().SetTitle("MLQ (GeV)")
rocaucVsMLQPlot.GetYaxis().SetTitle("ROC AUC")
rocaucVsMLQPlot.SetTitle("area under ROC vs MLQ")
c3.SetGridy()
rocaucVsMLQPlot.Draw("AP")
outFile.cd()
rocaucVsMLQPlot.Write()
c3.Print(pdf_folder+"/rocaucVsMLQ.pdf")
c3.Print(plotsForAN,"pdf")

oneMaucVsMLQPlot.SetMarkerStyle(8)
oneMaucVsMLQPlot.GetXaxis().SetTitle("MLQ (GeV)")
oneMaucVsMLQPlot.GetYaxis().SetTitle("1 - ROC AUC")
oneMaucVsMLQPlot.SetTitle("1 - area under ROC vs MLQ")
c3.SetLogy()
oneMaucVsMLQPlot.GetYaxis().SetRangeUser(1e-10,1)
oneMaucVsMLQPlot.Draw("AP")
outFile.cd()
oneMaucVsMLQPlot.Write()
c3.Print(pdf_folder+"/1-rocaucVsMLQ.pdf")
c3.Print(plotsForAN,"pdf")

#compare sig and bkg variable distributions for each mass
for i,mass in enumerate(LQmasses):
    if "2016preVFP" in base_folder:
        if "BEle" in pdf_folder:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToBEle_M-"+str(mass)+"_pair_TuneCP2_13TeV-madgraph-pythia8_APV.root")
        else:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToDEle_M-"+str(mass)+"_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8_APV.root")
    else:
        if "BEle" in pdf_folder:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToBEle_M-"+str(mass)+"_pair_TuneCP2_13TeV-madgraph-pythia8.root")
        else:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToDEle_M-"+str(mass)+"_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8.root")
    for var in variables:
       c5 = TCanvas()
       c5.SetLogy()
       #bkgPlot = bdtTFile.Get("LQM"+str(mass)+"/"+var+"_bkg")
       #sigPlot = bdtTFile.Get("LQM"+str(mass)+"/"+var+"_LQ"+str(mass))
       bkgPlot = TMVAFile.Get("dataset/InputVariables_Id/{}__Background_Id".format(var))
       sigPlot = TMVAFile.Get("dataset/InputVariables_Id/{}__Signal_Id".format(var))

    #   if not "Pt" in var and not "MET" in var and not "Mej" in var:
    #       bkgPlot.Rebin(25)
    #       sigPlot.Rebin(25)
    #   elif "Mej" in var:
    #       bkgPlot.Rebin(10)
    #       sigPlot.Rebin(10)
    #   else:
    #       bkgPlot.Rebin(5)
    #       sigPlot.Rebin(5)
       bkgPlot.SetLineColor(kRed)
       bkgPlot.SetFillColor(kRed)
       bkgPlot.SetFillStyle(3004)
       bkgPlot.SetMarkerColor(kRed)
       bkgPlot.SetStats(0)
       sigPlot.SetLineColor(kBlue)
       sigPlot.SetFillColor(kBlue)
       sigPlot.SetFillStyle(3005)
       #plotMax = max(bkgPlot.GetMaximum(), sigPlot.GetMaximum())
       #plotMin = min(bkgPlot.GetMinimum(0), sigPlot.GetMinimum(0))
       plotMax = 2*bkgPlot.GetMaximum()
       plotMin = min(bkgPlot.GetMinimum(0), sigPlot.GetMinimum(0))
       bkgPlot.GetYaxis().SetRangeUser(0.5*plotMin,2*plotMax)
       #bkgPlot.GetYaxis().SetRangeUser(0.1,1e5)
       #bkgPlot.GetXaxis().SetRangeUser(400,4700)
       bkgPlot.GetXaxis().SetTitle(var)
       bkgPlot.Draw()
       bkgPlot.GetYaxis().SetRangeUser(1e-3,2e5)
       plotMin = 0
       if "Pt" in var:
           plotMax = 700+100*i
           #bkgPlot.GetXaxis().SetLimits(0,2000)
           #sigPlot.GetXaxis().SetLimits(0,2000)
       elif "sT_eejj" in var:
           plotMax = 2400+150*i
           plotMin = 400
           #bkgPlot.GetXaxis().SetLimits(0,8000)
           #sigPlot.GetXaxis().SetLimits(0,8000)
       elif "Meejj" in var:
           plotMax = 3500+200*i
       elif "M_e1e2" in var:
           plotMax = 1400+200*i
           plotMin = 220
       else:
           plotMax = 1400+200*i
           #bkgPlot.GetXaxis().SetLimits(0,5000)
           #sigPlot.GetXaxis().SetLimits(0,5000)

       fakeHist = TH1D("forRange"+var+str(mass),"",1,plotMin,plotMax)
       fakeHist.SetTitle("signal and bkg for "+var+", MLQ="+str(mass)+" GeV")
       fakeHist.SetStats(0)
       fakeHist.GetYaxis().SetRangeUser(5e-3,2e5)
       bkgPlot.SetTitle("signal and bkg for "+var+", MLQ="+str(mass)+" GeV")
       fakeHist.Draw()
       bkgPlot.GetXaxis().SetRangeUser(0,bkgPlot.GetXaxis().GetXmax())
       bkgPlot.Draw("histSame")
       sigPlot.Draw("HISTsame") 
       l = TLegend(0.4,0.8,0.6,0.9)
       l.AddEntry(sigPlot,"signal","lp")
       l.AddEntry(bkgPlot,"background","lp")
       l.Draw("same")
       c5.Print(pdf_folder+"/"+str(mass)+"/"+var+"MLQ"+str(mass)+".pdf")
       if mass == LQmasses[0] and var == variables[0]:
           c5.Print(pdf_folder+"/allInputVars.pdf(","pdf")
       elif mass == LQmasses[-1] and var == variables[-1]:
           c5.Print(pdf_folder+"/allInputVars.pdf)","pdf")
       else:
           c5.Print(pdf_folder+"/allInputVars.pdf","pdf")

       if mass==500 or mass==1500 or mass==2500:
           c5.Print(plotsForAN,"pdf")

for i,mass in enumerate(LQmasses):
    #Also get overtraining plots 
    #see https://root-forum.cern.ch/t/tmva-signal-vs-background-bdt-response-plot-generated-with-macro-radically-different-from-manually-generated-plot/34002/4
    #cross-checked with the plots produced in the tmva gui and these are indeed the right plot names
    if "2016preVFP" in base_folder:
        if "BEle" in pdf_folder:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToBEle_M-"+str(mass)+"_pair_TuneCP2_13TeV-madgraph-pythia8_APV.root")
        else:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToDEle_M-"+str(mass)+"_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8_APV.root")
    else:
        if "BEle" in pdf_folder:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToBEle_M-"+str(mass)+"_pair_TuneCP2_13TeV-madgraph-pythia8.root")
        else:
            TMVAFile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToDEle_M-"+str(mass)+"_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8.root")
    bkgTrain = TMVAFile.Get("dataset/Method_BDT/BDTG/MVA_BDTG_Train_B")
    bkgTest = TMVAFile.Get("dataset/Method_BDT/BDTG/MVA_BDTG_B")
    sigTrain = TMVAFile.Get("dataset/Method_BDT/BDTG/MVA_BDTG_Train_S")
    sigTest = TMVAFile.Get("dataset/Method_BDT/BDTG/MVA_BDTG_S")

    bkgTrain.SetLineColor(kRed)
    bkgTrain.SetMarkerColor(kRed)
    sigTrain.SetLineColor(kBlue)
    sigTrain.SetMarkerColor(kBlue)
    bkgTest.SetLineColor(kRed)
    bkgTest.SetMarkerColor(kRed)
    sigTest.SetLineColor(kBlue)
    sigTest.SetMarkerColor(kBlue)

    sigTrain.SetMarkerStyle(8)
    sigTrain.SetMarkerSize(0.6)
    bkgTrain.SetMarkerStyle(8)
    bkgTrain.SetMarkerSize(0.6)

    sigTest.SetFillStyle(3354)
    sigTest.SetFillColor(kBlue)
    bkgTest.SetFillStyle(3345)
    bkgTest.SetFillColor(kRed)

    legend = TLegend(0.4,0.7,0.6,0.9)
    legend.AddEntry(sigTrain,"Training sample (sig)", "lp")
    legend.AddEntry(bkgTrain, "Training sample (bkg)", "lp")
    legend.AddEntry(sigTest, "Testing sample (sig)", "f")
    legend.AddEntry(bkgTest, "Testing sample (bkg)", "f")

    c = TCanvas()
    c.SetGridy()
    bkgTest.SetTitle("Overtraining check for MLQ = "+str(mass)+" GeV")
    bkgTest.SetStats(0)
    bkgTest.GetXaxis().SetTitle("BDT score")
    bkgTest.GetYaxis().SetTitle("(1/N) dN/dx")
    bkgTest.GetYaxis().SetRangeUser(0,21)
    bkgTest.Draw("HIST")
    sigTest.Draw("HISTSAME")
    bkgTrain.Draw("PSAME][")
    sigTrain.Draw("PSAME][")
    legend.Draw("SAME")
    c.Print(pdf_folder+"/"+str(mass)+"/BDTGOvertrainingPlot.pdf")
    if mass == LQmasses[0]:
        c.Print(pdf_folder+"/overtrainingPlots.pdf(","pdf")
    elif mass == LQmasses[-1]:
        c.Print(pdf_folder+"/overtrainingPlots.pdf)","pdf")
    else:
        c.Print(pdf_folder+"/overtrainingPlots.pdf","pdf")

    if mass==500 or mass==1500 or mass==2500:
        if mass==2500:
            plotsForAN+=")"
        c.Print(plotsForAN,"pdf")

#correlation matrices
if parameterized == True:
    exit()
exit()
cCorr = TCanvas()
#cCorr.SetLeftMargin(0.2)
#cCorr.SetBottomMargin(0.2)
if not os.path.isdir(pdf_folder+"/correlationMatrices"):
    os.mkdir(pdf_folder+"/correlationMatrices")
#Idk how to get the bin content based on the labels without hard-coding in the position like this.
jet1 = 10
jet2 = 11
jet3 = 12
ele1 = 8
ele2 = 9
mejMin = 13
mejMax = 14
jets = [jet1, jet2, jet3]
eles = [ele1, ele2]
mejs = [mejMin, mejMax]
jetLabels = ["jet1 pt","jet2 pt", "jet3 pt"]
eleLabels = ["ele1 pt", "ele2 pt"] 
mejLabels = ["Mej min", "Mej max"]
sigpoints = {}
bkgpoints = {}
for mass in LQmasses:
    sigpoints[str(mass)] = {}
    tfile = TFile.Open(base_folder+"/TMVA_ClassificationOutput_LQToDEle_M-"+str(mass)+"_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8_APV.root")
    sigMatrix = tfile.Get("dataset/CorrelationMatrixS")
    sigMatrix.GetXaxis().SetLabelSize(0.025)
    sigMatrix.GetYaxis().SetLabelSize(0.025)
    sigMatrix.SetTitle("Correlation Matrix (signal MLQ = "+str(mass)+"GeV)")
    sigMatrix.Draw("colzText")
    cCorr.Print(pdf_folder+"/correlationMatrices/MLQ"+str(mass)+".pdf")
    if mass == 300: #bkg ones are all the same so I only need it once
        bkgMatrix = tfile.Get("dataset/CorrelationMatrixB")
        bkgMatrix.GetXaxis().SetLabelSize(0.025)
        bkgMatrix.GetYaxis().SetLabelSize(0.025)
        bkgMatrix.Draw("colzText")
        for ijet in jets:
            bkgpoints[str(ijet)] = {}
            for iele in eles:
                bkgpoints[str(ijet)][str(iele)] = bkgMatrix.GetBinContent(bkgMatrix.GetBin(iele, ijet))
        for imej in mejs:
            bkgpoints[str(imej)] = {}
            for ijet in jets:
                bkgpoints[str(imej)][str(ijet)] = bkgMatrix.GetBinContent(bkgMatrix.GetBin(imej, ijet))
            for iele in eles:
                bkgpoints[str(imej)][str(iele)] = bkgMatrix.GetBinContent(bkgMatrix.GetBin(imej, iele))
        cCorr.Print(pdf_folder+"/correlationMatrices/background.pdf")
        

    for ijet in jets:
        sigpoints[str(mass)][str(ijet)] = {}
        for iele in eles:
                sigpoints[str(mass)][str(ijet)][str(iele)] = sigMatrix.GetBinContent(sigMatrix.GetBin(iele, ijet))
    for imej in mejs:
        sigpoints[str(mass)][str(imej)] = {}
        for ijet in jets:
            sigpoints[str(mass)][str(imej)][str(ijet)] = sigMatrix.GetBinContent(sigMatrix.GetBin(imej, ijet))
        for iele in eles:
            sigpoints[str(mass)][str(imej)][str(iele)] = sigMatrix.GetBinContent(sigMatrix.GetBin(imej, iele))

cCorr.SetGridy()
for ijet in jets:
    for iele in eles:
        plot = TGraph()
        plot.SetTitle("correlation between "+jetLabels[ijet-10]+" and "+eleLabels[iele-8]+" vs MLQ")
        for i,mass in enumerate(LQmasses):
            plot.SetPoint(i,mass,sigpoints[str(mass)][str(ijet)][str(iele)])
        bkgLine = TLine(250, bkgpoints[str(ijet)][str(iele)], 2050, bkgpoints[str(ijet)][str(iele)])
        bkgLine.SetLineColor(kRed)
        zeroLine = TLine(250,0,2050,0)
        zeroLine.SetLineColor(kGreen+1)
        plot.SetMarkerStyle(8)
        plot.SetMarkerSize(0.75)
        plot.GetXaxis().SetRangeUser(250,2050)
        plot.GetXaxis().SetTitle("MLQ (GeV)")
        plot.GetYaxis().SetRangeUser(-50,50)
        l = TLegend(0.7,0.7,0.9,0.9)
        l.AddEntry(plot, "signals", "p")
        l.AddEntry(bkgLine, "background value", "l")
        l.AddEntry(zeroLine, "zero", "l")
        plot.Draw("ALP")
        zeroLine.Draw("same")
        bkgLine.Draw("same")
        plot.Draw("LP")
        l.Draw("Same")
        cCorr.Print(pdf_folder+"/correlationMatrices/"+str(ijet)+"_"+str(iele)+".pdf")

eleAndJetLabels = eleLabels+jetLabels
print(bkgpoints["13"])
for imej in mejs:
    for iparticle in eles+jets:
        plot = TGraph()
        plot.SetTitle("correlation between "+mejLabels[imej-13]+" and "+eleAndJetLabels[iparticle-8]+" vs MLQ")
        for i, mass in enumerate(LQmasses):
            plot.SetPoint(i, mass, sigpoints[str(mass)][str(imej)][str(iparticle)])
        bkgLine = TLine(250, bkgpoints[str(imej)][str(iparticle)], 2050, bkgpoints[str(imej)][str(iparticle)])
        bkgLine.SetLineColor(kRed)
        zeroLine = TLine(250,0,2050,0)
        zeroLine.SetLineColor(kGreen+1)
        plot.SetMarkerStyle(8)
        plot.SetMarkerSize(0.75)
        plot.GetXaxis().SetRangeUser(250,2050)
        plot.GetXaxis().SetTitle("MLQ (GeV)")
        plot.GetYaxis().SetRangeUser(-50,50)
        l = TLegend(0.7,0.7,0.9,0.9)
        l.AddEntry(plot, "signals", "p")
        l.AddEntry(bkgLine, "background value", "l")
        l.AddEntry(zeroLine, "zero", "l")
        plot.Draw("ALP")
        zeroLine.Draw("same")
        bkgLine.Draw("same")
        plot.Draw("LP")
        l.Draw("Same")
        cCorr.Print(pdf_folder+"/correlationMatrices/"+str(imej)+"_"+str(iparticle)+".pdf")
