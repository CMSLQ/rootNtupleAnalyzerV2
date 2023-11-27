from __future__ import print_function

from ROOT import (
    gROOT,
    TFile,
    TH1D,
    TCanvas,
    TColor,
    TLegend,
    kRed,
    kGreen,
    kAzure,
    kOrange,
    kGray,
    kViolet,
    kBlack,
    THStack,
    TRatioPlot,
    TPad,
    TLine,
    TBox,
    TF1,
)
import ROOT
import numpy as np
import copy
import ctypes
import sys
import os
import math

def MakeStackAndRatioPlot(histDict,histoKeysForStack):
    bkgHist = 0
    stack = THStack()
    dataHist = histDict["data"]
    for key in histoKeysForStack:
        hist = histDict["bkg"][key]
        stack.Add(hist)
        if not bkgHist:
            bkgHist = copy.deepcopy(hist)
        else:
            bkgHist.Add(hist)
    stack.SetMaximum(1e4)
    stack.SetMinimum(1)
    ratioPlot = copy.deepcopy(dataHist)
    ratioPlot.Divide(bkgHist)
    ratioPlot.SetTitle("")
    ratioPlot.SetStats(0)
    ratioPlot.GetYaxis().SetRangeUser(0,2)
    ratioPlot.GetYaxis().SetLabelSize(0.08)
    ratioPlot.GetYaxis().SetTitle("1P1F data / bkg.")
    ratioPlot.GetYaxis().SetTitleSize(0.06)
    ratioPlot.GetXaxis().SetLabelSize(0.08)
    ratioPlot.GetXaxis().SetTitleSize(0.08)
    return stack, ratioPlot

#input_file = "$LQDATA/2016/qcdFRClosureTest/frClosureTest_2016pre_July23/FRCTCombined.root"
#input_file = "$LQDATA/2016/qcdFRClosureTest/frClosureTest_2016pre_Aug23/MCfix/FRCTCombined.root"
#input_file = "$LQDATA/2016/qcdFRClosureTest/frClosureTest_2016pre_sept2023/withSF/FRCTCombined.root"
#input_file = "$LQDATA/2017/qcdFRClosureTest/frClosureTest_2017_sept2023/gte1Jet/noSF/FRCTCombined.root"
input_file = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/FRCTCombined.root"
pdf_folder = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/plots"
fit_results = "/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/fitResults.txt"
if not os.path.isdir(pdf_folder) and pdf_folder != "":
    print("Making directory ", pdf_folder)
    os.mkdir(pdf_folder)
gROOT.SetBatch(True)


tfile = TFile.Open(input_file)
binsDict = {}

etaRegs = ["_BB","_BE","_EE",""]

#get histos
variable_names = [
    "Pt1stEle_PAS",
    "Me1j1_PAS",
    "Mee_PAS",
    "sT_PAS",
    "Me2j1_PAS",
    "MET_PAS",
    "Pt2ndEle_PAS",
    "Mt_MET_Ele1_PAS",
    "Mt_MET_Ele2_PAS",
    "Phi1stEle_PAS",
    "Phi2ndEle_PAS",
    "METPhi_PAS",
#    "HT",
    "Pt1stEle_tight",
    "Me1j1_tight",
    "Mee_tight",
    "sT_tight",
    "Me2j1_tight",
    "MET_tight",
    "Pt2ndEle_tight",
    "Mt_MET_Ele1_tight",
    "Mt_MET_Ele2_tight",
    "Phi1stEle_tight",
    "Phi2ndEle_tight",
    "METPhi_tight",
#    "HT_tight",
]
for name in variable_names:
    folderToMake = pdf_folder +"/"+ name
    if not os.path.isdir(folderToMake) and folderToMake != "":
        print("Making directory ", folderToMake)
        os.mkdir(folderToMake)
histDict = {}
for region in etaRegs:
    histDict[region] = {}
    for var in variable_names:
        histDict[region][var] = {}
        histDict[region][var]["bkg"] = {}
        histDict[region][var]["bkg"]["fakeRate"] = tfile.Get("histo1D__QCDFakes_DATA2F_WError__"+var+region)
        #histDict[region][var]["bkg"]["fakeRate"] = tfile.Get("histo1D__QCDFakes_DATA2F__"+var+region)
        histDict[region][var]["data"] = tfile.Get("histo1D__QCDFakes_DATA1P1F__"+var+region)
        histDict[region][var]["bkg"]["MCOnly"] = tfile.Get("histo1D__MCTotal-WHTBinnedIncStitch__"+var+region)

mcSamples = ["ZJet_NNLO_IncStitch", "WJet_HTBinned_IncStitch", "TTBar_powheg", "SingleTop", "GJets", "DIBOSON_nlo"]

mcShortNames = ["ZJets", "WJets", "TTBar", "ST", "GJets", "Diboson"]
allBkgNames = ["ZJets", "WJets", "TTBar", "ST", "GJets", "Diboson","fakeRate"]
allBkgNamesNoW = ["ZJets","TTBar", "ST", "GJets", "Diboson","fakeRate"]
mcShortNamesNoW = ["ZJets", "TTBar", "ST", "GJets", "Diboson"]
colors = [kRed, kViolet+6, 600, kGreen, kAzure+1, kOrange-3,kGray+1]

for index, sample in enumerate(mcSamples):
    for region in etaRegs:
        for variableName in variable_names:
            histo = tfile.Get("histo1D__"+sample+"__"+variableName+region)
            histDict[region][variableName]["bkg"][mcShortNames[index]] = histo


#comparison of WJets
WJetNames = ["WJet_amcatnlo_jetBinned", "WJet_HTBinned"]
c = TCanvas()
fPads1 = TPad("p1","",0.00,0.3,0.99,0.99)
fPads2 = TPad("p2","",0.00,0.00,0.99,0.301)
fPads1.SetFillColor(0)
fPads1.SetLineColor(0)
fPads2.SetFillColor(0)
fPads2.SetLineColor(0)
fPads1.SetBottomMargin(1e-2)
fPads2.SetTopMargin(3e-2)
fPads2.SetBottomMargin(3e-1)
fPads1.SetLogy()
fPads1.SetGridy()
fPads2.SetGridy()
fPads1.Draw()
fPads2.Draw()
for variable in variable_names:
    l = TLegend(0.6,0.8,0.9,0.9)
    for name in WJetNames:
        histo = tfile.Get("histo1D__"+name+"__"+variable)
        histo.SetLineWidth(2)
        #histo.Rebin(2)
        #histo.GetXaxis().SetRangeUser(0,2000)
        #print("found histo ", histo.GetName())
        if "HT" in name:
            histo.SetLineColor(kRed)
            histo.SetMarkerColor(kRed)
        #if "amcatnlo_jet" in name:
        #    histo.SetLineColor(kAzure+1)
        #    histo.SetMarkerColor(kAzure+1)
        histo.SetStats(0)
        histo.SetMinimum(0.01)
        var = variable.replace("_PAS","")
        var = var.replace("_tight","")
        histo.GetXaxis().SetTitle(var +" (GeV)")
        if "amcatnlo" in name:
            histoMax = 10*histo.GetBinContent(histo.GetMaximumBin())
            histo.GetYaxis().SetRangeUser(0.1,1000)#histoMax)
        if "MET_tight" in variable:
            histo.GetXaxis().SetRangeUser(0,100)
        if "HT" in variable:
            histo.GetXaxis().SetRange(0,500)
            histo.Rebin(2)
        if "Mt" in variable:
            histo.GetXaxis().SetRangeUser(0,1000)
        if "tight" in variable:
            histo.SetTitle("WJet samples "+var+" BDT selection")
        elif "PAS" in variable:
            histo.SetTitle("WJet samples "+var+" preselection")
        else:
            histo.SetTitle("WJet samples "+variable)
        if "HT" in name:
            l.AddEntry(histo, "WJet_LO_HTBinned","lp")
        else:
            l.AddEntry(histo, name, "lp")
        fPads1.cd()
        if "amcatnlo" in name:
            histo.Draw()
        else:
            histo.Draw("same")
    l.Draw("same")

    fPads2.cd()
    histoHT = tfile.Get("histo1D__WJet_HTBinned__"+variable)
    histoJet = tfile.Get("histo1D__WJet_amcatnlo_jetBinned__"+variable)
    #histoHT.Rebin(2)
    #histoJet.Rebin(2)
    ratio = copy.deepcopy(histoHT)
    ratio.Divide(histoJet)
    ratio.SetLineColor(kBlack)
    ratio.SetMarkerColor(kBlack)
    ratio.SetTitle("")
    ratio.SetStats(0)
    ratio.GetYaxis().SetRangeUser(0,2)
    ratio.GetYaxis().SetTitle("LO WJets / NLO WJets")
    ratio.Draw()
    xLow = ratio.GetXaxis().GetXmin()
    xHigh = ratio.GetXaxis().GetXmax()
    line = TLine(xLow,1,xHigh,1)
    line.SetLineStyle(7)
    line.SetLineColor(kGray+2)
    line.Draw("Same")
    ratio.Draw("same")
    c.Print("/eos/user/e/eipearso/LQ/lqData/2017/qcdFRClosureTest/frClosureTest_2017_nov2023-d/WJetComparisons/"+variable+".pdf")

#set colors and style
for region in etaRegs:
    for var in variable_names:
        histDict[region][var]["data"].SetLineWidth(2)
        for index, name in enumerate(allBkgNames):
            print("set color and style for "+name+" "+var+" "+region)
            if "fake" in name and "HT" in var:
                continue
            histo = histDict[region][var]["bkg"][name]
            histo.SetLineColor(colors[index])
            histo.SetFillColor(colors[index])
            histo.SetMarkerColor(colors[index])
            histo.SetLineWidth(2)
            histo.SetStats(0)
'''
for region in etaRegs:
    for var in variable_names:
        for name in allBkgNames:
            histo = copy.deepcopy(histDict[region][var]["bkg"][name])
            c = TCanvas()
            c.SetLogy()
            histo.GetXaxis().SetTitle("GeV")
            histo.GetYaxis().SetRangeUser(0.1, 3000)
            histo.SetTitle(name + " " +var+region)
            histo.Draw()
            c.Print(pdf_folder+"/"+var+"/"+name+"MC.pdf")
'''
#print(histDict)
#rebin histos
for variable in variable_names:
    binSize = 10
    lowestEdge = 0
    minWidth = 50

    if "Pt" in variable:
        plotRange = 1000
        lowestEdge = 50
    elif "MET_PAS" in variable or "MET_tight" in variable:
        plotRange = 100
        minWidth = 10
        binSize = 5
    elif "Mee_PAS" in variable:
        lowestEdge = 110
    elif "sT_PAS" in variable:
        lowestEdge = 200
    elif "Mee_tight" in variable:
        lowestEdge = 220
    elif "sT_tight" in variable:
        lowestEdge = 400
    elif "phi" in variable.lower():
        lowestEdge = -3.1416
        plotRange = 3.1416
        minWidth = 0.1
        binSize = (2*3.1416)/60 
    elif "Mt" in variable:
        plotRange = 1000
    else:
        plotRange = 2000
    nStartingBins = int(plotRange / binSize)
    if "Phi" in variable:
        nStartingBins = 60
    threshold = 100
    total = 0
    binEdges = [lowestEdge]
    currentBinWidth = 0
    hist = histDict[""][variable]["data"]
    if lowestEdge > 0:
        startBin = hist.GetXaxis().FindBin(lowestEdge+1)
    else: 
        startBin = 1 
    for i in range(startBin,nStartingBins+1):
        content = hist.GetBinContent(i)
        total += content
        currentBinWidth += binSize
        if currentBinWidth >= minWidth and total > threshold:
            lowEdge = hist.GetBinLowEdge(i)
            binEdges.append(lowEdge+binSize) #I need the high edge, which is low edge + bin width
            total = 0
            currentBinWidth = 0
            if lowEdge+binSize >= plotRange:
                break
    #If we get to the end of the for loop without reaching the end of the plot, there are < 100 events in the rest of the range.
    if lowEdge+binSize < plotRange:
        binEdges.pop()
        binEdges.append(plotRange)
        
    print("binning for ", variable)
    print(binEdges)
    binsDict[variable] = binEdges

histDictRebinned = {}
for region in etaRegs:
    histDictRebinned[region] = {}
    for var in variable_names:
        histDictRebinned[region][var] = {}
        bins = binsDict[var]
        #print("rebinning: ",var," Bins: ",bins)
        histDictRebinned[region][var]["data"] = histDict[region][var]["data"].Rebin(len(bins)-1,"data_rebinned",np.array(bins, dtype = float))
        #histDictRebinned[region][var]["data"] = copy.deepcopy(histDict[region][var]["data"])
        #histDictRebinned[region][var]["data"].Rebin(4)
        histDictRebinned[region][var]["bkg"] = {}
        histDictRebinned[region][var]["bkg"]["MCOnly"] = histDict[region][var]["bkg"]["MCOnly"].Rebin(len(bins)-1,"data_rebinned",np.array(bins, dtype = float))
        #histDictRebinned[region][var]["bkg"]["MCOnly"] = copy.deepcopy(histDict[region][var]["bkg"]["MCOnly"])
        #histDictRebinned[region][var]["bkg"]["MCOnly"].Rebin(4)
        for name in allBkgNames:
            histDictRebinned[region][var]["bkg"][name] = histDict[region][var]["bkg"][name].Rebin(len(bins)-1,name+"_rebinned",np.array(bins, dtype = float))
            #histDictRebinned[region][var]["bkg"][name] = copy.deepcopy(histDict[region][var]["bkg"][name])
            #histDictRebinned[region][var]["bkg"][name].Rebin(4)

#make background MC plots. One set as is and one rebinned, for all eta regions
c1 = TCanvas()
c1.SetGridy()
#c1.SetLogy()
c2 = TCanvas()
c2.SetGridy()
c2.SetLogy()
for reg in etaRegs:
    for var in variable_names:
        varPieces = var.split("_")
        axisTitle = varPieces[0] + " (GeV)"
        if "MET_tight" in var:
            l = TLegend(0.4,0.1,0.7,0.3)
        else:
            l = TLegend(0.5,0.7,0.8,0.9)
        c1.cd()
        for i,name in enumerate(mcShortNames):
            histo = histDict[reg][var]["bkg"][name]
            if i==0:
                histCopy = copy.deepcopy(histo)
                histCopy.GetXaxis().SetTitle(axisTitle)
                if "st" in var.lower():
                    histCopy.GetXaxis().SetRangeUser(150,2000)
                elif "mee" in var.lower():
                    histCopy.GetXaxis().SetRangeUser(50,2000)
                elif "me1j1" in var.lower():
                    histCopy.GetXaxis().SetRangeUser(0,2000)
                else:
                    histCopy.GetXaxis().SetRangeUser(0,1000)
                histCopy.Draw()
            else:
                histo.Draw("same")
            l.AddEntry(histo,name,"lp")
            l.Draw("same")
        c1.Print(pdf_folder+"/"+var+"/mcHists_"+var+reg+".pdf")
        c2.cd()
        for i,name in enumerate(mcShortNames):
            histo2 = histDictRebinned[reg][var]["bkg"][name]
            histo2.GetYaxis().SetRangeUser(0.1,1e5)
            #histo.Rebin(2)
            if i==0:
                hist2Copy = copy.deepcopy(histo2)
                hist2Copy.GetXaxis().SetTitle(axisTitle)
                if "st" in var.lower():
                    hist2Copy.GetXaxis().SetRangeUser(150,2000)
                elif "mee" in var.lower():
                    hist2Copy.GetXaxis().SetRangeUser(50,2000)
                elif "me1j1" in var.lower():
                    hist2Copy.GetXaxis().SetRangeUser(0,2000)
                elif "pt" in var.lower():
                    hist2Copy.GetXaxis().SetRangeUser(0,1000)
                hist2Copy.Draw()
            else:
                histo2.Draw("same")
        l.Draw("same")
        c2.Print(pdf_folder+"/"+var+"/mcHistsRebinned_"+var+reg+".pdf")

#make plot of un-rebinned data
for reg in etaRegs:
    for var in variable_names:
        histo = histDict[reg][var]["data"]
        c = TCanvas()
        if "MET_tight" in var:
            histo.GetXaxis().SetRangeUser(0,100)
        if "Mt" in var:
            histo.GetXaxis().SetRangeUser(0,1000)
        shortVarName = var.replace("_PAS", "")
        shortVarName.replace("_tight", "")
        histo.GetXaxis().SetTitle(shortVarName)
        if "tight" in var:
            histo.SetTitle(shortVarName+" BDT")
        elif "PAS" in var:
            histo.SetTitle(shortVarName+" preselection")
        else:
            histo.SetTitle(shortVarName)
        histo.Draw()
        c.Print(pdf_folder+"/"+var+"/dataHist_"+var+reg+".pdf")

#make stack and ratio plots
stackAllBkg = {}
ratioAllBkg = {}
for reg in etaRegs:
    stackAllBkg[reg] = {}
    ratioAllBkg[reg] = {}
    for var in variable_names:
        stackAllBkg[reg][var],ratioAllBkg[reg][var] = MakeStackAndRatioPlot(histDictRebinned[reg][var], allBkgNames)

#make canvas, pads and legend
c = TCanvas()
c.SetLogy()
fPads1 = TPad("p1","",0.00,0.3,0.99,0.99)
fPads2 = TPad("p2","",0.00,0.00,0.99,0.301)
fPads1.SetFillColor(0)
fPads1.SetLineColor(0)
fPads2.SetFillColor(0)
fPads2.SetLineColor(0)
fPads1.SetBottomMargin(1e-2)
fPads2.SetTopMargin(3e-2)
fPads2.SetBottomMargin(3e-1)
fPads1.SetLogy()
fPads2.SetGridy()
fPads1.Draw()
fPads2.Draw()
leg = TLegend(0.6,0.6,0.9,0.9)
legMETPlot = TLegend(0.3,0.1,0.6,0.3)
#Note: it doesn't matter which eta region and var name I use to make the legend bc it's the same samples for all of them
for name in allBkgNames:
    histo = histDictRebinned[""]["Pt1stEle_PAS"]["bkg"][name]
    leg.AddEntry(histo,name,"lp")
    legMETPlot.AddEntry(histo,name,"lp")
leg.AddEntry(histDictRebinned[""]["Pt1stEle_PAS"]["data"],"1P1F data","lp")
legMETPlot.AddEntry(histDictRebinned[""]["Pt1stEle_PAS"]["data"],"1P1F data","lp")
#draw stack and ratio plots
fitResults = {}
fitResults["stack"] = {}
for reg in etaRegs:
    fitResults["stack"][reg] = {}
    for var in variable_names:
        varPieces = var.split("_")
        axisTitle = varPieces[0] + " (GeV)"
        fPads1.cd()
        title = var+reg
        if "PAS" in var:
            title = title.replace("_PAS", "")
            title+=" preselection"
        if "tight" in var:
            title = title.replace("_tight", "")
            title+=" BDT selection"
        stackAllBkg[reg][var].SetTitle(title)
        stackAllBkg[reg][var].Draw("hist")
        histDictRebinned[reg][var]["data"].Draw("same")
        #histDict[reg]["data"].Draw("same")
        if "MET_PAS" in var or "MET_tight" in var:
            legMETPlot.Draw("same")
        else:
            leg.Draw("same")

        fPads2.cd()
        ratioAllBkg[reg][var].GetXaxis().SetTitle(axisTitle)
        ratioAllBkg[reg][var].Draw()
        line = TLine(0,1,2000,1)
        line.Draw("same")
        line.SetLineStyle(5)
        line.SetLineColorAlpha(13, 0.5)
        xAxisLow = ratioAllBkg[reg][var].GetXaxis().GetXmin()
        xAxisHigh = ratioAllBkg[reg][var].GetXaxis().GetXmax()
        errWindow = TBox(xAxisLow,0.75,xAxisHigh,1.25)
        errWindow.SetLineColor(kGray)
        errWindow.SetFillColorAlpha(kGray,0.50)
        errWindow.Draw("same")
        fitFunction = TF1("fit", "pol0", xAxisLow, xAxisHigh)
        fit = ratioAllBkg[reg][var].Fit(fitFunction, "S", "", xAxisLow, xAxisHigh)
        fitResults["stack"][reg][var] = {}
        fitResults["stack"][reg][var]["value"] = fit.Parameter(0)
        fitResults["stack"][reg][var]["error"] = fit.ParError(0)
        fitResults["stack"][reg][var]["chi2"] = fit.Chi2()
        fitResults["stack"][reg][var]["ndf"] = fit.Ndf()
        ratioAllBkg[reg][var].Draw("same")
        c.Print(pdf_folder+"/"+var+"/stack_plot_"+var+reg+".pdf")
    #c.Print(pdf_folder+"/stack_plot"+reg+"_not_rebinned.pdf")
#make the same stack/ratio plot without the 2F data added in
#fPads1.cd()
#fPads1.SetLogy()
#stackMCOnly.Draw("hist")
#histDict["rebinned"]["data"].Draw("same")
#fPads2.cd()
#fPads2.SetGridy()
#ratioMCOnly.Draw()
#line.Draw("same")
#c.Print(pdf_folder+"/stack_plot_without2F.pdf")

#make mc sub plot
c = TCanvas()
c.cd()
tpad1 = TPad("pad1","",0.00,0.3,0.99,0.99)
tpad2 = TPad("pad2","",0.00,0.00,0.99,0.301)
tpad1.SetLogy()
tpad1.SetGridy()
tpad2.SetGridy()
tpad1.Draw()
tpad2.Draw()
ratioHists = {}
fitResults["MCSub"] = {}
for reg in etaRegs:
    ratioHists[reg] = {}
    fitResults["MCSub"][reg] = {}
    for var in variable_names:
        varPieces = var.split("_")
        axisTitle = varPieces[0] +" (GeV)"
        histoMCSub = histDictRebinned[reg][var]["data"] - histDictRebinned[reg][var]["bkg"]["MCOnly"]
        #histoMCSub = histDict[reg]["data"] - histDict[reg]["bkg"]["MCOnly"]
        histoFR = histDictRebinned[reg][var]["bkg"]["fakeRate"]
        #histoFR = histDict[reg]["bkg"]["fakeRate"]
        histoMCSub.SetStats(0)
        title = var+reg
        if "PAS" in var:
            title = title.replace("_PAS", "")
            title+=" preselection"
        if "tight" in var:
            title = title.replace("_tight", "")
            title+=" BDT selection"
        histoMCSub.SetTitle(title)
        if "gte" in pdf_folder.lower(): #gte 1Jet has more stats so I need a bigger range here since I don't do this one log scale
            histoMCSub.GetYaxis().SetRangeUser(-4000,4000)
        else:
            histoMCSub.GetYaxis().SetRangeUser(0.1,1e4)
        l = TLegend(0.6,0.7,0.9,0.9)
        l.AddEntry(histoMCSub,"data - MC 1P1F region","lp")
        l.AddEntry(histoFR,"prediction by fake rate","lp")
        histoMCSub.GetXaxis().SetTitle(axisTitle)
        histoMCSub.Draw()
        histoFR.Draw("Esame")
        l.Draw("same")
        
        MC_FRratio = copy.deepcopy(histoFR)
        MC_FRratio.Divide(histoMCSub)
        tpad2.cd()
        MC_FRratio.GetYaxis().SetTitle("fake rate pred. / MCSub yield")
        MC_FRratio.SetTitle("")
        MC_FRratio.Draw()
        xAxisLow = MC_FRratio.GetXaxis().GetXmin()
        xAxisHigh = MC_FRratio.GetXaxis().GetXmax()
        errWindow = TBox(xAxisLow,0.75,xAxisHigh,1.25)
        errWindow.SetLineColor(kGray)
        errWindow.SetFillColorAlpha(kGray,0.50)
        errWindow.Draw("same")
        MC_FRratio.SetLineColor(kBlack)
        MC_FRratio.SetMarkerColor(kBlack)
        #ratioHists[reg][var] = ratio 
        fitFunction = TF1("fit", "pol0", xAxisLow, xAxisHigh)
        fitFunction.SetParameter(0,1)
        fit = MC_FRratio.Fit(fitFunction, "S", "", xAxisLow, xAxisHigh)
        fitResult = fit.Parameter(0)
        fitResultErr = fit.ParError(0)
        fitResults["MCSub"][reg][var] = {}
        fitResults["MCSub"][reg][var]["value"] = fitResult
        fitResults["MCSub"][reg][var]["error"] = fitResultErr
        fitResults["MCSub"][reg][var]["chi2"] = fit.Chi2()
        fitResults["MCSub"][reg][var]["ndf"] = fit.Ndf()
        MC_FRratio.GetYaxis().SetRangeUser(0,3)
        MC_FRratio.Draw("same")
        c.Print(pdf_folder+"/"+var+"/mcSub_fr_comparison_logScale_"+var+reg+".pdf")
    #c.Print(pdf_folder+"/mcSub_fr_comparison"+reg+"_not_rebinned.pdf")
with open(fit_results, "w") as resultsFile:
    resultsFile.write("MC subtraction ratio plot fit results:\n\n")
    for var in variable_names:
        if "PAS" in var:
            varToWrite = var.replace("PAS", "preselection")
        elif "tight" in var:
            varToWrite = var.replace("tight", "BDT")
        else:
            varToWrite = var
        line = str(fitResults["MCSub"][""][var]["value"]) + " +/- " + str(fitResults["MCSub"][""][var]["error"])
        resultsFile.write(varToWrite+"\n")
        resultsFile.write("    Chi2: "+str(fitResults["MCSub"][reg][var]["chi2"])+"\n")
        resultsFile.write("    NDF: "+str(fitResults["MCSub"][reg][var]["ndf"])+"\n")
        chi2OverNdf = fitResults["MCSub"][reg][var]["chi2"]/fitResults["MCSub"][reg][var]["ndf"]
        resultsFile.write("    Chi2 / NDF: "+str(chi2OverNdf)+"\n")
        resultsFile.write("    fit: "+line+"\n\n")
#put fit results in plot
fitResultsHistBDT = TH1D("fitResultsB", "fit results", 7,0,7)
fitResultsHistPresel = TH1D("fitResultsP", "fit results presel", 7,0,7)
binCounter = 1
varNamesPresel = ["Pt1stEle_PAS", "Mee_PAS", "Me1j1_PAS", "sT_PAS", "Pt2ndEle_PAS", "Me2j1_PAS", "MET_PAS"]
for var in varNamesPresel: 
    fitResultsHistPresel.SetBinContent(binCounter, fitResults["MCSub"][reg][var]["value"])
    fitResultsHistPresel.SetBinError(binCounter, fitResults["MCSub"][reg][var]["error"])
    fitResultsHistPresel.SetLineWidth(2)
    fitResultsHistPresel.GetXaxis().SetBinLabel(binCounter, var.replace("_PAS",""))
    fitResultsHistPresel.SetMarkerColor(kRed)
    fitResultsHistPresel.SetLineColor(kRed)
    varBDT = var.replace("PAS","tight")
    fitResultsHistBDT.SetBinContent(binCounter, fitResults["MCSub"][reg][varBDT]["value"])
    fitResultsHistBDT.SetBinError(binCounter, fitResults["MCSub"][reg][varBDT]["error"])
    fitResultsHistBDT.SetLineWidth(2)
    fitResultsHistBDT.GetXaxis().SetBinLabel(binCounter, varBDT.replace("_tight",""))
    binCounter+=1
print("fit of the fit results: ")
c = TCanvas()
c.SetGridy()
fitResultsHistBDT.SetStats(0)
fitResultsHistBDT.GetYaxis().SetRangeUser(0.5,1.0)
fitResultsHistBDT.SetNdivisions(25,"X")
f1 = TF1("fitResultsFit", "pol0", 0,7)
fitResultsHistBDT.Fit(f1,"","",0,7)
fitResultsHistBDT.Draw()
#fitResultsHistPresel.Draw("same")
l = TLegend(0.6,0.8,0.9,0.9)
l.AddEntry(fitResultsHistBDT, "BDT fit results", "lp")
l.AddEntry(fitResultsHistPresel, "preselection fit results", "lp")
#l.Draw("same")
c.Print(pdf_folder+"/fitResults.pdf")

#make plot of just FR prediction by itself for comparison w the AN
c = TCanvas()
c.cd()
c.SetLogy()
c.SetGridy()
oldPtBins = [50,100,150,220,300,400,500,600,800,1000]
histoFR = histDict[""]["Pt1stEle_PAS"]["bkg"]["fakeRate"].Rebin(len(oldPtBins)-1,"FR_with_old_ptBins",np.array(oldPtBins,dtype = float))
histoFR.SetStats(0)
histoFR.GetYaxis().SetRangeUser(0.1,1e5)
histoFR.GetXaxis().SetTitle("Pt1stEle (GeV)")
histoFR.Draw()
c.Print(pdf_folder+"/Pt1stEle_PAS/fakeRatePrediction.pdf")

#Get and print Pt1stEle, HEEPEle_Pt, LooseEle_Pt
'''
histoHEEPEle = tfile.Get("histo1D__QCDFakes_DATA1P1F__HEEPEle_Pt")
histoLooseEle = tfile.Get("histo1D__QCDFakes_DATA1P1F__LooseEle_Pt")
histo1stEle = tfile.Get("histo1D__QCDFakes_DATA1P1F__Pt1stEle_PAS")
histoHEEPEle.SetLineWidth(2)
histoHEEPEle.SetStats(0)
histoHEEPEle.SetTitle("pT, 1-pass-1-fail data")
histoHEEPEle.GetXaxis().SetTitle("pT (GeV)")
histoHEEPEle.GetYaxis().SetRangeUser(0.1,4500)
histoLooseEle.SetLineWidth(2)
histoLooseEle.SetLineColor(kGreen)
histoLooseEle.SetMarkerColor(kGreen)
histo1stEle.SetLineWidth(2)
histo1stEle.SetLineColor(kRed)
histo1stEle.SetMarkerColor(kRed)
leg = TLegend(0.6,0.7,0.9,0.9)
leg.AddEntry(histoHEEPEle, "pass HEEP electron pT", "lp")
leg.AddEntry(histoLooseEle, "fail HEEP electron pT", "lp")
leg.AddEntry(histo1stEle, "Pt1stEle_PAS", "lp")
c = TCanvas()
#c.SetLogy()
c.SetGridy()
histoHEEPEle.Draw()
histo1stEle.Draw("same")
histoLooseEle.Draw("same")
histoHEEPEle.Draw("same")
leg.Draw("same")
c.Print(pdf_folder+"/elePtComparison.pdf")
'''
#do MC subtraction and get yeild
frError = ctypes.c_double()
integralFR = histDict[""]["Mee_tight"]["bkg"]["fakeRate"].IntegralAndError(histDict[""]["Mee_tight"]["bkg"]["fakeRate"].GetXaxis().FindBin(110),histDict[""]["Mee_tight"]["bkg"]["fakeRate"].GetXaxis().GetLast(),frError)
#integralFR = histDict[""]["Mee_PAS"]["bkg"]["fakeRate"].Integral()
#histoFRErrSQ = tfile.Get("histo1D__QCDFakes_DATA2F__errFRsq_Mee")
#integralFRErrSQ = histoFRErrSQ.Integral()
#frError = math.sqrt(integralFRErrSQ)
data1P1FError = ctypes.c_double()
integral1P1F = histDict[""]["Mee_tight"]["data"].IntegralAndError(histDict[""]["Mee_tight"]["data"].GetXaxis().FindBin(110),histDict[""]["Mee_tight"]["data"].GetXaxis().GetLast(),data1P1FError)
MCError = ctypes.c_double()
integralMC = histDict[""]["Mee_tight"]["bkg"]["MCOnly"].IntegralAndError(histDict[""]["Mee_tight"]["bkg"]["MCOnly"].GetXaxis().FindBin(110),histDict[""]["Mee_tight"]["bkg"]["MCOnly"].GetXaxis().GetLast(),MCError)
MCSubError = math.sqrt(MCError.value**2 + data1P1FError.value**2)
integralMCSub = integral1P1F - integralMC
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("yield predicted by fake rate: ",integralFR," +/- ",frError.value)
print("data 1P1F region: ",integral1P1F," +/- ",data1P1FError.value)
print("MC 1P1F region: ",integralMC," +/- ",MCError.value)
print("actual yield (1P1F data - MC): ",integralMCSub," +/- ",MCSubError)
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

#print error histo
'''
histoFRErr = copy.deepcopy(histoFRErrSQ)
for i in range(histoFRErr.GetNbinsX()+1):
    errSQ = histoFRErr.GetBinContent(i)
    err = math.sqrt(errSQ)
    histoFRErr.SetBinContent(i,err)
c = TCanvas()
histoFRErr.Draw("hist")
c.Print(pdf_folder+"/test.pdf")
'''

print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("yield and uncertainty, first three bins, Mee BDT")
print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("sample    :    value +/- uncertainty")
for i in range(1,4):
    print("bin ", i)
    for name in allBkgNames:
        value = histDictRebinned[""]["Mee_tight"]["bkg"][name].GetBinContent(i)
        uncertainty = histDictRebinned[""]["Mee_tight"]["bkg"][name].GetBinError(i)
        print(name,"    :    ",value," +/- ",uncertainty)
    totBkg = histDictRebinned[""]["Mee_tight"]["bkg"]["MCOnly"].GetBinContent(i)
    totBkgErr = histDictRebinned[""]["Mee_tight"]["bkg"]["MCOnly"].GetBinError(i)
    print("total MC    :    ",totBkg," +/- ",totBkgErr)
    data = histDictRebinned[""]["Mee_tight"]["data"].GetBinContent(i)
    dataErr = histDictRebinned[""]["Mee_tight"]["data"].GetBinError(i)
    print("1P1F data    :    ",data," +/- ",dataErr)
    print("1P1F data - MC : ",data - totBkg, " +/- ", math.sqrt(totBkgErr**2 + dataErr**2))
    print("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

