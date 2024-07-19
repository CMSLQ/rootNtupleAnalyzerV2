from __future__ import print_function

from ROOT import (
    gROOT,
    gPad,
    TFile,
    TH1D,
    TH2D,
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
from tabulate import tabulate

def MakeStackAndRatioPlot(histDict1P1F, histDict2F, MCNames, year, variable):
    bkgHist = 0
    stack = THStack()
    dataHist = histDict1P1F[year]["data"][variable]
    for key in MCNames:
        hist = histDict1P1F[year][key][variable]
        stack.Add(hist)
        if not bkgHist:
            bkgHist = copy.deepcopy(hist)
        else:
            bkgHist.Add(hist)
    data2Fhist = histos2F[year][variable]
    stack.Add(data2Fhist)
    bkgHist.Add(data2Fhist)

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

def RebinTH2(hist2D):
    xAxisHigh = int(hist2D.GetXaxis().GetBinUpEdge(hist2D.GetXaxis().GetLast()))
    bins = [
        50,
        75,
        90,
        120,
        150,
        175,
        200,
        225,
        250,
        300,
        350,
        400,
        500,
        600,                                                            
        1000, 
    ]
    histRebinned = TH2D(hist2D.GetName()+"_rebinned", hist2D.GetTitle(), hist2D.GetNbinsX(), 0, xAxisHigh, len(bins)-1, np.array(bins,dtype=float))
    for xbin in range(1,hist2D.GetNbinsX()+1):
        for ybin in range(1,hist2D.GetNbinsY()+1):
            currentX = hist2D.GetXaxis().GetBinCenter(xbin)
            currentY = hist2D.GetYaxis().GetBinCenter(ybin)
            rebinnedX = histRebinned.GetXaxis().FindFixBin(currentX)
            rebinnedY = histRebinned.GetYaxis().FindFixBin(currentY)
            binToAdd = hist2D.GetBinContent(xbin, ybin)
            binErrToAdd = hist2D.GetBinErrorUp(xbin, ybin)
            rebinnedContent = histRebinned.GetBinContent(rebinnedX, rebinnedY)
            rebinnedContent+=binToAdd
            rebinnedErr = histRebinned.GetBinErrorUp(rebinnedX, rebinnedY)
            rebinnedErr = math.sqrt(rebinnedErr**2 + binErrToAdd**2)
            histRebinned.SetBinContent(rebinnedX, rebinnedY, rebinnedContent)
            histRebinned.SetBinError(rebinnedX, rebinnedY, rebinnedErr)
    return histRebinned

def AddFRError(dataHist, frHist, etaRegion):
    #Assumes dataHist is binned with the same pT binning as the frHist
    #fake rate is applied in the analysis class, here we just get the error and add that on
    for xbin in range(dataHist.GetNbinsX()+1):
        for ybin in range(dataHist.GetNbinsY()+1):
            yCoordinate = dataHist.GetYaxis().GetBinCenter(ybin)
            fakeRateHistBinY = frHist.GetYaxis().FindFixBin(yCoordinate)
            frXCoordinate = -1
            if etaRegion=="_Barrel":
                frXCoordinate = 0.75
            elif etaRegion=="_End1":
                frXCoordinate = 1.75
            elif etaRegion=="_End2":
                frXCoordinate = 2.25
            else:
                print("Cannot determine fake rate bin for eta region "+etaRegion)
            fakeRateHistBinX = frHist.GetXaxis().FindFixBin(frXCoordinate)
            fakeRate = frHist.GetBinContent(fakeRateHistBinX, fakeRateHistBinY)
            fakeRateErr = frHist.GetBinError(fakeRateHistBinX, fakeRateHistBinY)
            fakeRateEff = fakeRate / (1-fakeRate)
            fakeRateEffErr = fakeRateErr / ((1-fakeRate)**2)
            binContent = dataHist.GetBinContent(xbin,ybin)
            binErr = dataHist.GetBinErrorUp(xbin,ybin)
            if fakeRateEff>0:
                binErr = math.sqrt(binErr**2 + (fakeRateEffErr**2)*((binContent/fakeRateEff)**2))
            #err = sqrt( binErr**2 + (frEffErr**2)(N_unweighted)**2 )
            #Since the fake rate was already applied in the analysis class,
            #We have N_unweighted = binContent/fakeRateEff
            #Which only works if the dataHist has the same binning as the fake rate calc.
            dataHist.SetBinError(xbin, ybin, binErr)
    return dataHist

def Get2F1DHistFrom2DHists(dataHists, frHists, etaRegions, HEMRegions, var):
    for eta in etaRegions:
        for hem in HEMRegions:
            hist = dataHists[eta][hem]
            print(hist)
            frHist = frHists[hem]
            rebinnedHist = RebinTH2(hist)
            rebinnedHist = AddFRError(rebinnedHist, frHist, eta)
            dataHists[eta][hem] = rebinnedHist
    nbinsX = dataHists[etaRegions[0]][HEMRegions[0]].GetNbinsX()
    xAxisHigh = dataHists[etaRegions[0]][HEMRegions[0]].GetXaxis().GetBinUpEdge(dataHists[etaRegions[0]][HEMRegions[0]].GetXaxis().GetLast())
    hist1D = TH1D("histo1D__QCDFakes_DATA__{}".format(var),"histo1D__QCDFakes_DATA__{}".format(var), nbinsX, 0, xAxisHigh)
    for xbin in range(hist1D.GetNbinsX()+1):
        binContent = 0
        binErr = 0
        errList = []
        for eta in etaRegions:
            for hem in HEMRegions:
                hist2D = dataHists[eta][hem]
                for ybin in range(1, hist2D.GetNbinsY()+1):
                    binContent+=hist2D.GetBinContent(xbin, ybin)
                    binErrCurrent = hist2D.GetBinErrorUp(xbin, ybin)
                    errList.append(binErrCurrent)
        hist1D.SetBinContent(xbin, binContent)
        for ierr in errList:
            for jerr in errList:
                binErr += ierr * jerr
        binErr = math.sqrt(binErr)
        hist1D.SetBinError(xbin,binErr)
        #tfile.cd()
        #hist1D.Write()
    return hist1D


gROOT.SetBatch(True)

years = ["2016preVFP", "2016postVFP", "2017", "2018"]
#years = ["2017"]
filenames = {}
for year in years:
    filenames[year] = {}

for year in years:
    filenames[year]["2F"] = "$LQDATAEOS/qcdFRClosureTest_allYears_newErrorCalc/{}/2F/output_cutTable_lq_QCD_FakeRateClosureTest/qcdFRClosureTest_allYears_plots.root".format(year)
    filenames[year]["1P1F"] = "$LQDATAEOS/qcdFRClosureTest_allYears_newErrorCalc/{}/1P1F_SF/output_cutTable_lq_QCD_FakeRateClosureTest/qcdFRClosureTest_allYears_plots.root".format(year)
    filenames[year]["FR"] = os.getenv("LQINPUTS")+"/fakeRate/{}/fr2D{}.root".format(year, year.replace("2016",""))

pdf_folder = os.getenv("LQDATAEOS")+"/qcdFRClosureTest_allYears_newErrorCalc/plots_test"

if not os.path.isdir(pdf_folder):
    os.mkdir(pdf_folder)
    
for year in years + ["fullRunII"]:
    if not os.path.isdir(pdf_folder+"/"+year):
        os.mkdir(pdf_folder+"/"+year)

variableNameList = [
    "Pt1stEle_PAS",
    "Me1j1_PAS",
    "Mee_PAS",
    "sT_PAS",
    "MET_PAS",
    "Me2j1_PAS",
    "Pt2ndEle_PAS", 
#    "HT",
#    "Mt_MET_Ele1_PAS",
#    "Mt_MET_Ele2_PAS",
#    "Phi1stEle_PAS",
#    "Phi2ndEle_PAS",
#    "METPhi_PAS",
#    "nJet_PAS",
    "Pt1stEle_tight", 
    "Me1j1_tight", 
    "Mee_tight", 
    "sT_tight", 
    "MET_tight", 
    "Me2j1_tight", 
    "Pt2ndEle_tight",
#    "HT_tight",
#    "Mt_MET_Ele1_tight",
#    "Mt_MET_Ele2_tight",
#    "Phi1stEle_tight",
#    "Phi2ndEle_tight",
#    "METPhi_tight",
#    "MeeControlReg",
]

for var in variableNameList:
    for year in years+["fullRunII"]:
        f = pdf_folder+"/"+year+"/"+var
        if not os.path.isdir(f):
            os.mkdir(f)

histoNameData = "histo1D__QCDFakes_DATA__"
#mcSamples = {}
mcSamples = [
    "ZJet_amcatnlo_ptBinned_IncStitch",
    "WJet_HTBinned",
    "TTBar_powheg", 
    "SingleTop", 
    "GJets", 
    "DIBOSON_nlo",
]
'''
mcSamples["2018"] = [
    "ZJet_amcatnlo_ptBinnedIncStitch",
    "WJetHTBinned",
    "TTBar_powheg",
    "SingleTop",
    "GJets",
    "DIBOSON",
]
mcSamples["fullRunII"] = [
    "ZJet_amcatnlo_ptBinned_IncStitch",
    "WJet_HTBinned",
    "TTBar_powheg",
    "SingleTop",
    "GJets",
    "DIBOSON_nlo",
]
'''
mcShortNames = [
    "ZJets",
    "WJets",
    "TTBar",
    "ST",
    "GJets",
    "Diboson",
]

colors = [kRed, kViolet+6, 600, kGreen, kAzure+1, kOrange-3,kGray+1]

histos1P1F = {}
histos2F = {}

histos1P1F["fullRunII"] = {}
histos2F["fullRunII"] = {}
histos1P1F["fullRunII"]["data"] = {}
histos1P1F["fullRunII"]["MCTotal"] = {}

etaRegions = ["_Barrel", "_End1","_End2"]

histoNamesMC = []
for name in mcSamples:
    histoNamesMC.append("histo1D__"+name+"__")
    histos1P1F["fullRunII"]["histo1D__"+name+"__"] = {}
    histos2F["fullRunII"]["histo1D__"+name+"__"] = {}

data2DHists = {}
#Get hists 
for iyear,year in enumerate(years):
    histos1P1F[year] = {}
    for name in histoNamesMC:
        histos1P1F[year][name] = {}
    histos2F[year] = {}
    #1P1F histos:
    tfile = TFile.Open(filenames[year]["1P1F"])
    histos1P1F[year]["data"] = {}
    histos1P1F[year]["MCTotal"] = {}
    for var in variableNameList:
        histToGet = histoNameData+var #.replace("PAS","tight")
        print("Get 1P1F hist "+histToGet)
        histoData = tfile.Get(histToGet)
        histoData.SetLineWidth(2)
        histoData.SetStats(0)
        if "MET" in var:
            histoData.GetXaxis().SetRangeUser(0,100)
        histos1P1F[year]["data"][var] = copy.deepcopy(histoData)
        if iyear==0:
            histos1P1F["fullRunII"]["data"][var] = copy.deepcopy(histoData)
        else:
            histos1P1F["fullRunII"]["data"][var].Add(histoData)

        histoMCTotal = 0
        for i,name in enumerate(histoNamesMC):
            print("Get 1P1F hist "+name+var)
            histo = tfile.Get(name+var)#.replace("PAS","tight"))
            if "MET" in var:
                histo.GetXaxis().SetRangeUser(0,100)
     #       if year=="2017" and "Mee" in var:
     #           histo.Rebin(10)
            c = colors[i]
            histo.SetLineColor(c)
            histo.SetFillColor(c)
            histo.SetMarkerColor(c)
            histo.SetLineWidth(2)
            histo.SetStats(0)
            histos1P1F[year][name][var] = copy.deepcopy(histo)
            if iyear==0:
                histos1P1F["fullRunII"][name][var] = copy.deepcopy(histo)
            else:
                histos1P1F["fullRunII"][name][var].Add(histo)

            if i==0:
                histoMCTotal = copy.deepcopy(histo)
                histoMCTotal.SetName("histo1D__MCTotal__"+var)
            else:
                histoMCTotal.Add(histo)
        histos1P1F[year]["MCTotal"][var] = copy.deepcopy(histoMCTotal)
        if iyear==0:
            histos1P1F["fullRunII"]["MCTotal"][var] = copy.deepcopy(histoMCTotal)
        else:
            histos1P1F["fullRunII"]["MCTotal"][var].Add(histoMCTotal)

    #2F histos
    if year=="2018":
        HEMRegions = ["_pre319077","_post319077_HEMOnly","_post319077_noHEM"]
        print(HEMRegions)
    else:
        HEMRegions = [""]
    tfileFR = TFile.Open(filenames[year]["FR"])
    frHistBaseName = "fr2D_1Jet_TrkIsoHEEP7vsHLTPt_{}"
    frHists = {}
    if year == "2018":
        for hem in HEMRegions:
            if "pre" in hem:
                histName = frHistBaseName.format("pre319077")
            elif "HEMOnly" in hem:
                histName = frHistBaseName.format("HEMonly_post319077")
            elif "noHEM" in hem:
                histName = frHistBaseName.format("noHEM_post319077")
            else:
                print("cannot determine fake rate hist for HEM region "+HEMRegion)
            frHists[hem] = copy.deepcopy(tfileFR.Get(histName))
    else:
        histName = frHistBaseName.format("PAS")
        frHists[""] = copy.deepcopy(tfileFR.Get(histName))
    print(frHists)
    tfile2 = TFile.Open(filenames[year]["2F"])
    for var in variableNameList:
        for eta in etaRegions:
            data2DHists[eta] = {}
            for hem in HEMRegions:
                histoName2D = "histo2D__QCDFakes_DATA__{}{}{}".format(var, eta, hem)
                print("Get {} 2F hist ".format(year)+histoName2D)
                histoData2D = copy.deepcopy(tfile2.Get(histoName2D))
                data2DHists[eta][hem] = histoData2D
        histoData = Get2F1DHistFrom2DHists(data2DHists, frHists, etaRegions, HEMRegions, var)#, tfile) 
        if "MET" in var:
            histoData.GetXaxis().SetRangeUser(0,100)
      #  if year=="2017" and "Mee" in var:
      #      histoData.Rebin(10)
        c = colors[6]
        histoData.SetLineColor(c)
        histoData.SetMarkerColor(c)
        histoData.SetFillColor(c)
        histoData.SetLineWidth(2)
        histoData.SetStats(0)
        #print(histoData)
        histoNameErr = histoNameData+"errFRsq_"+var.replace("_PAS","")
        #print(histoNameErr)
        #histoErrSQ = tfile2.Get(histoNameErr+var)
        #for i in range(histoData.GetNbinsX()):
        #    errSQ = histoErrSQ.GetBinContent(i)
        #    err = math.sqrt(errSQ)
        #    histoData.SetBinError(i,err)
        histos2F[year][var] = copy.deepcopy(histoData)
        if iyear==0:
            histos2F["fullRunII"][var] = copy.deepcopy(histoData)
        else:
            histos2F["fullRunII"][var].Add(histoData)

#print(histos1P1F)
#print(histos2F)

binsDict = {}
#rebin histos
for year in years + ["fullRunII"]:
    binsDict[year] = {}
    for variable in variableNameList:
        #print("determine binning for {}, {}".format(year, variable))
        binSize = 10
        lowestEdge = 0
        minWidth = 50

        plotRange = 2000
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
            #binSize = 1
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
        threshold = 200
        total = 0
        binEdges = [lowestEdge]
        currentBinWidth = 0
        hist = histos1P1F[year]["data"][variable]
        #print(hist)
        if lowestEdge > 0:
            startBin = hist.GetXaxis().FindBin(lowestEdge+1)
            #if "Mee" in variable:
             #   startBin = 22
        else: 
            startBin = 1 
        lowEdge = lowestEdge
        #print("lowestEdge = ",lowestEdge)
        #print("startBin = {}, nStartingBins = {}".format(startBin, nStartingBins))
        for i in range(startBin,nStartingBins+1):
            content = hist.GetBinContent(i)
         #   print("content for bin {} = {}".format(i,content))
            total += content
            currentBinWidth += binSize
          #  print("currentBinWidth = ", currentBinWidth)
            if currentBinWidth >= minWidth and total > threshold:
                lowEdge = hist.GetBinLowEdge(i)
                binEdges.append(lowEdge+binSize) #I need the high edge, which is low edge + bin width
           #     print("append bin edge ", lowEdge+binSize)
                total = 0
                currentBinWidth = 0
                if lowEdge+binSize >= plotRange:
                    break
    #If we get to the end of the for loop without reaching the end of the plot, there are < 100 events in the rest of the range.
        if lowEdge+binSize < plotRange:
            binEdges.pop()
            binEdges.append(plotRange)
        
        #print("binning for ", variable)
       # print(binEdges)
        binsDict[year][variable] = binEdges
      
for year in years + ["fullRunII"]:
    for var in variableNameList:
        bins = binsDict["fullRunII"][var]
        histos1P1F[year]["data"][var] = histos1P1F[year]["data"][var].Rebin(len(bins)-1,histos1P1F[year]["data"][var].GetName(),np.array(bins, dtype = float))
        histos2F[year][var] = histos2F[year][var].Rebin(len(bins)-1,histos2F[year][var].GetName(),np.array(bins, dtype = float))
        if "Pt" in var:
            rebin = 5
        elif "MET" in var:
            rebin = 1
        else:
            rebin = 10
        #histos1P1F[year]["data"][var].Rebin(rebin)
        #histos2F[year][var].Rebin(rebin)
        for name in histoNamesMC:
           histos1P1F[year][name][var] = histos1P1F[year][name][var].Rebin(len(bins)-1,histos1P1F[year][name][var].GetName(),np.array(bins, dtype = float))
     #      histos1P1F[year][name][var].Rebin(rebin)
        histos1P1F[year]["MCTotal"][var] = histos1P1F[year]["MCTotal"][var].Rebin(len(bins)-1,histos1P1F[year]["MCTotal"][var].GetName(),np.array(bins, dtype=float))
     #   histos1P1F[year]["MCTotal"][var].Rebin(rebin)

#Do calc and make plots
stackAllBkg = {}
ratioAllBkg = {}
for year in years + ["fullRunII"]:
    stackAllBkg[year] = {}
    ratioAllBkg[year] = {}
    for var in variableNameList:
        stackAllBkg[year][var],ratioAllBkg[year][var] = MakeStackAndRatioPlot(histos1P1F, histos2F, histoNamesMC, year, var)

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
fPads1.SetGridy()
fPads2.SetGridy()
fPads1.Draw()
fPads2.Draw()
leg = TLegend(0.9,0.6,0.999,0.9)
#legMETPlot = TLegend(0.3,0.1,0.6,0.3)

leg.AddEntry(histos1P1F["2017"]["data"]["sT_tight"], "1-pass-1-fail data", "lp")
leg.AddEntry(histos2F["2017"]["sT_tight"],"predicted fakes", "lp")
for i,name in enumerate(histoNamesMC):
    leg.AddEntry(histos1P1F["2017"][name]["sT_tight"], mcShortNames[i], "lp")

legMETPlot = copy.deepcopy(leg)
#legMETPlot.SetX1NDC(0.3)
#legMETPlot.SetY1NDC(0.1)

for year in years + ["fullRunII"]:
    for var in variableNameList:
        title = var.replace("_tight"," BDT training region")
        title = var.replace("_PAS"," preselection")
        fPads1.cd()
        histos1P1F[year]["data"][var].SetTitle(title+" "+year)
        if "MET" in var:
            histos1P1F[year]["data"][var].GetXaxis().SetRangeUser(0,100)
        histos1P1F[year]["data"][var].GetYaxis().SetRangeUser(1,1e5)
        histos1P1F[year]["data"][var].Draw()
        stackAllBkg[year][var].Draw("histSame")
        histos1P1F[year]["data"][var].Draw("SAME")
        if "MET" in var:
            legMETPlot.Draw("SAME")
            gPad.Update()
            #legMETPlot.SetX1NDC(0.75)
            #legMETPlot.SetY1NDC(0.6)
            #legMETPlot.SetX2NDC(0.9)
            #legMETPlot.SetY2NDC(0.9)
            gPad.Modified()
        else:
            leg.Draw("SAME") 

        fPads2.cd()
        ratioAllBkg[year][var].GetXaxis().SetTitle(var.split("_")[0]+" GeV")
        if "MET" in var:
            ratioAllBkg[year][var].GetXaxis().SetRangeUser(0,100)
        ratioAllBkg[year][var].Draw()
        xAxisLow = ratioAllBkg[year][var].GetXaxis().GetXmin()
        xAxisHigh = ratioAllBkg[year][var].GetXaxis().GetXmax()
        errWindow = TBox(xAxisLow,0.75,xAxisHigh,1.25)
        errWindow.SetLineColor(kGray)
        errWindow.SetFillColorAlpha(kGray,0.50)
        errWindow.Draw("same")
        line = TLine(xAxisLow,1,xAxisHigh,1)
        line.SetLineStyle(5)
        line.SetLineColorAlpha(13, 0.5)
        line.Draw("SAME")
        ratioAllBkg[year][var].Draw("SAME")
        c.Print(pdf_folder+"/"+year+"/"+var+"/stack_plot_"+var+".pdf")
        if year == years[0] and var == variableNameList[0]:
            c.Print(pdf_folder+"/closureTestPlots.pdf(","pdf")
        else:
            c.Print(pdf_folder+"/closureTestPlots.pdf","pdf")

        histoMCSub = histos1P1F[year]["data"][var] - histos1P1F[year]["MCTotal"][var]
        histoFR = histos2F[year][var]
        histoMCSub.SetTitle(title+" "+year)
        histoMCSub.GetYaxis().SetRangeUser(0.1,1e4)
        fPads1.cd()
        l = TLegend(0.6,0.8,0.9,0.9)
        l.AddEntry(histoFR, "predicted fakes", "lp")
        l.AddEntry(histoMCSub,"observed fakes", "lp")
        histoMCSub.Draw()
        histoFR.Draw("ESAME")
        l.Draw("SAME")

        ratioPlot = copy.deepcopy(histoFR)
        ratioPlot.Divide(histoMCSub)
        fPads2.cd()
        if "MET" in var:
            ratioPlot.GetXaxis().SetRangeUser(0,100)
        ratioPlot.GetYaxis().SetTitle("predicted / observed")
        ratioPlot.SetTitle("")
        ratioPlot.SetLineColor(kBlack)
        ratioPlot.SetMarkerColor(kBlack)
        ratioPlot.GetYaxis().SetRangeUser(0,3)
        ratioPlot.GetYaxis().SetLabelSize(0.08)
        ratioPlot.GetYaxis().SetTitleSize(0.06)
        ratioPlot.GetXaxis().SetLabelSize(0.08)
        ratioPlot.GetXaxis().SetTitleSize(0.08)
        ratioPlot.GetXaxis().SetTitle(title+" GeV")
        ratioPlot.Draw()
        errWindow.Draw("SAME")
        ratioPlot.Draw("SAME")
        c.Print(pdf_folder+"/"+year+"/"+var+"/MCSub_"+var+".pdf")
        if year == "fullRunII" and var == variableNameList[-1]:
            c.Print(pdf_folder+"/closureTestPlots.pdf)","pdf")
        else:
            c.Print(pdf_folder+"/closureTestPlots.pdf","pdf")

#Do uncertainty calc
resultsFile = pdf_folder+"/results.txt"
with open(resultsFile, 'w') as f:
    f.write("fake rate closure test results \n\n")
variables = ["Mee_PAS", "Mee_tight"]
#variables = ["Mee_tight"]
for var in variables:
    for year in years + ["fullRunII"]:
        DataErr = ctypes.c_double()
        DataHist = histos1P1F[year]["data"][var]
        DataYield = DataHist.IntegralAndError(DataHist.GetXaxis().GetFirst(), DataHist.GetXaxis().GetLast(), DataErr)

        MCErr = ctypes.c_double()
        MCHist = histos1P1F[year]["MCTotal"][var]
        MCYield = MCHist.IntegralAndError(MCHist.GetXaxis().GetFirst(), MCHist.GetXaxis().GetLast(), MCErr)

        FRPredErr = ctypes.c_double()
        FRPredHist = histos2F[year][var]
        FRPred = FRPredHist.IntegralAndError(FRPredHist.GetXaxis().GetFirst(), FRPredHist.GetXaxis().GetLast(), FRPredErr)
        obsFakes = DataYield - MCYield
        obsFakesErr = math.sqrt(DataErr.value**2 + MCErr.value**2)

        ratio = FRPred / obsFakes
        ratioErr = math.sqrt(((1/obsFakes)*FRPredErr.value)**2 + ((FRPred/obsFakes**2)*obsFakesErr)**2)
        #Make a table
        headers = ["", "yield", "error"]
        table = [
            ["Predicted fakes",FRPred, FRPredErr.value],
            ["1-pass-1-fail data",DataYield, DataErr.value],
            ["1-pass-1-fail MC", MCYield, MCErr.value],
            ["Observed fakes", obsFakes, obsFakesErr],
            ["predicted / observed", ratio, ratioErr],
        ]
        if "PAS" in var:
            print(year+" preselection")
        elif "tight" in var:
            print(year+ " BDT training region")
        else:
            print(year+ " "+var)
        print("\npredicted and observed fakes: ")
        print(tabulate(table, headers=headers, stralign="left"))

        with open(resultsFile, 'a') as f:
            if "PAS" in var:
                f.write("Preselection")
            elif "tight" in var:
                f.write("BDT training region")
            else:
                f.write(var)
            f.write(year)
            f.write("\npredicted and observed fakes: \n")
            f.write(tabulate(table, headers=headers, stralign="left"))
            f.write("\n\nlatex table: \n")
            f.write(tabulate(table, headers=headers, tablefmt='latex',stralign="left"))
            f.write("\n\n")
        #Do breakdown of MC
        tableMC = []
        for iname, name in enumerate(histoNamesMC):
            err = ctypes.c_double()
            hist = histos1P1F[year][name][var]
            integral = hist.IntegralAndError(hist.GetXaxis().GetFirst(), hist.GetXaxis().GetLast(), err)
            percentage = round(integral / MCYield, 3)
            shortName = mcShortNames[iname]
            tableMC.append([shortName, integral, err.value, percentage])
        tableMC.append(["Total", MCYield, MCErr.value, ""])
        print("\nbreakdown of MC yield: ")
        print(tabulate(tableMC, headers=headers+["% of total MC"], stralign="left"))

        print("++++++++++++++++++++++++++++++++++++++++++++++++++\n")
        with open(resultsFile, 'a') as f:
            f.write("\nbreakdown of MC yield: \n")
            f.write(tabulate(tableMC, headers=headers+["% of total MC"], stralign="left"))
            f.write("\n\nlatex table:\n")
            f.write(tabulate(tableMC, headers=headers+["% of total MC"], tablefmt='latex', stralign="left"))
            f.write("\n++++++++++++++++++++++++++++++++++++++++++++++++++\n")
print("tables written to "+resultsFile)
