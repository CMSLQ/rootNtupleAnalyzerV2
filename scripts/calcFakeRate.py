from ROOT import (
    gROOT,
    TFile,
    TH2F,
    TCanvas,
    kRed,
    kBlue,
    kSpring,
    kAzure,
    TGraphAsymmErrors,
    TH1D,
    THStack,
    TLegend,
    kWhite,
    kGreen,
)
import numpy as np
import copy
import ctypes
import sys
import os
import math


def GetXBinsFromGraph(graph):
    bins = []
    nPoints = graph.GetN()
    for i in xrange(nPoints):
        bins.append(graph.GetPointX(i)-graph.GetErrorXlow(i))
    bins.append(graph.GetPointX(nPoints-1)+graph.GetErrorXhigh(nPoints-1))
    return bins


def GetCanvasTitle(varName, region, jetBin):
    titleStr = str(analysisYear)+" FR,"
    if "post" in varName.lower():
        titleStr += " Run >= "+varName[varName.lower().find("post")+4:]+","
    elif "pre" in varName.lower():
        titleStr += " Run < "+varName[varName.lower().find("pre")+3:]+","
    varNameSplit = varName.split("_")
    for token in varNameSplit:
        if "hem" in token.lower():
            titleStr += " "+token.replace("HEM", "HEM15and16")+","
    titleStr += " "+region.replace("End2", "Endcap2").replace("End1", "Endcap1").replace("Bar", "Barrel")
    if jetBin == "":
        titleStr += ", >= 0 jets"
    elif jetBin == "1Jet_":
        titleStr += ", >= 1 jets"
    elif jetBin == "2Jet_":
        titleStr += ", >= 2 jets"
    else:
        titleStr += ", "+jetBin
    # print "Created canvas title:", titleStr, "from var=", varName, "region=", reg, "jets=", jetBin
    return titleStr


def GetFakeRate(lowEnd, highEnd, reg, jets, histDict, verbose=False):
    verbose = True
    histo2D_Electrons = histDict[reg]["Electrons"][jets]
    histo2D_Jets = histDict[reg]["Jets"][jets]
    histo2D_Data = histDict[reg]["Total"][jets]
    if verbose:
        print "\t\tGetFakeRate: Considering region=", reg
        print "\t\tGetFakeRate: Considering jets=", jets
        print "\t\tGetFakeRate: Considering electron Pt:", str(lowEnd) + "-" + str(highEnd), "GeV"
        print "\t\tGetFakeRate: Considering histo for electrons=>", histo2D_Electrons.GetName(), "; histo for jets=>", histo2D_Jets.GetName(), "; histo for data=>", histo2D_Data.GetName()
        sys.stdout.flush()
    proj_Electrons = histo2D_Electrons.ProjectionY(
        "ElesTrkIso",
        histo2D_Electrons.GetXaxis().FindBin(lowEnd),
        histo2D_Electrons.GetXaxis().FindBin(highEnd) - 1,
    )
    proj_Jets = histo2D_Jets.ProjectionY(
        "JetsTrkIso",
        histo2D_Jets.GetXaxis().FindBin(lowEnd),
        histo2D_Jets.GetXaxis().FindBin(highEnd) - 1,
    )
    proj_Data = histo2D_Data.ProjectionY(
        "DataTrkIso",
        histo2D_Data.GetXaxis().FindBin(lowEnd),
        histo2D_Data.GetXaxis().FindBin(highEnd) - 1,
    )
    # proj_JetsEle_ = proj_Jets.Clone()
    # proj_JetsEle_.SetName('JetsEleTrkIso')
    eleErr = ctypes.c_double()
    ele = proj_Electrons.IntegralAndError(
        proj_Electrons.GetXaxis().FindBin(10), proj_Electrons.GetXaxis().FindBin(20) - 1, eleErr
    )
    # data = proj_Data_.Integral(proj_Data_.GetXaxis().FindBin(10),proj_Data_.GetXaxis().FindBin(20)-1)
    dataErr = ctypes.c_double()
    data = proj_Data.IntegralAndError(proj_Data.GetXaxis().GetFirst(), proj_Data.GetXaxis().GetLast(), dataErr)
    jets_SBErr = ctypes.c_double()
    jets_SB = proj_Jets.IntegralAndError(
        proj_Jets.GetXaxis().FindBin(10), proj_Jets.GetXaxis().FindBin(20) - 1, jets_SBErr
    )
    jets_SRErr = ctypes.c_double()
    jets_SR = proj_Jets.IntegralAndError(
        proj_Jets.GetXaxis().FindBin(0), proj_Jets.GetXaxis().FindBin(5) - 1, jets_SRErr
    )
    if jets_SB == 0:
        rJets = 0
        if jets_SR == 0:
            rJetsErr = math.sqrt(jets_SRErr.value**2+jets_SBErr.value**2)
        else:
            rJetsErr = jets_SRErr.value/jets_SR
    else:
        rJets = (jets_SR / jets_SB)
        rJetsErr = (jets_SR/jets_SB)*math.sqrt((jets_SRErr.value/jets_SR)**2+(jets_SBErr.value/jets_SB)**2)
    numerator = rJets * ele
    try:
        numeratorErr = numerator*math.sqrt((rJetsErr/rJets)**2+(eleErr.value/ele)**2)
    except ZeroDivisionError as e:
        if ele == 0:
            if rJets != 0:
                numeratorErr = numerator*math.sqrt((rJetsErr/rJets)**2)
                print "WARN: GetFakeRate:  Had a ZeroDivisionError in numeratorError; ele==0, so ignore it in error computation"
            else:
                numeratorErr = math.sqrt(rJetsErr**2+eleErr.value**2)
        elif rJets == 0:
            numeratorErr = numerator*math.sqrt((eleErr.value/ele)**2)
            print "WARN: GetFakeRate:  Had a ZeroDivisionError in numeratorError; rJets==0, so ignore it in error computation"
        else:
            print "ERROR in GetFakeRate:  Had a ZeroDivisionError in numeratorError; rJets=", rJets, "; ele=", ele
            sys.stdout.flush()
            raise e
    if verbose:
        print "\t\tGetFakeRate: Considering hists:", histo2D_Electrons.GetName(), ",", histo2D_Jets.GetName(), ",", histo2D_Data.GetName()
        # print 'endcap1 N_jet>=0: Nele=',ele,'data=',data,'contam=',ele/data
        # print 'barrel N_jet>=0: Nele=',ele,'data=',data,'contam=',ele/data
        print "\t\tGetFakeRate: nHEEPprime=", ele, " +/-", eleErr.value, "; jetsSR=", jets_SR, "+/-", jets_SRErr.value, "; jetsSB=", jets_SB, "+/-", jets_SBErr.value, ";data=", data, "+/-", dataErr.value
        sys.stdout.flush()
        print "\t\tGetFakeRate: FR=", numerator / data
        sys.stdout.flush()
    if writeOutput:
        suffix = "_"+reg+"_"+jets+"Pt"+str(lowEnd)+"To"+str(highEnd)
        if not outputFile.cd("TrkIso_Electrons"):
            outputFile.mkdir("TrkIso_Electrons").cd()
        proj_Electrons.SetName(proj_Electrons.GetName()+suffix)
        proj_Electrons.Write()
        if not outputFile.cd("TrkIso_Jets"):
            outputFile.mkdir("TrkIso_Jets").cd()
        proj_Jets.SetName(proj_Jets.GetName()+suffix)
        proj_Jets.Write()
        if not outputFile.cd("TrkIso_Data"):
            outputFile.mkdir("TrkIso_Data").cd()
        proj_Data.SetName(proj_Data.GetName()+suffix)
        proj_Data.Write()
    return numerator, numeratorErr, data, dataErr.value


def GetFakeRateMCSub(lowEnd, highEnd, reg, jets, histDict, verbose=False):
    if verbose:
        print "\t\tGetFakeRateMCSub"
        print "\t\t\tConsidering region=", reg
        print "\t\t\tConsidering jets=", jets
        print "\t\t\tConsidering electron Pt:", str(lowEnd) + "-" + str(highEnd), "GeV"
    histo2D_Electrons = histDict[reg]["Electrons"][jets]
    histo2D_MC = histDict[reg]["MC"][jets]
    mc2DHists = [histDict[reg][y][jets] for y in mcNames]
    histo2D_Data = histDict[reg]["Total"][jets]
    proj_Electrons = histo2D_Electrons.ProjectionY(
        "ElesTrkIso",
        histo2D_Electrons.GetXaxis().FindBin(lowEnd),
        histo2D_Electrons.GetXaxis().FindBin(highEnd) - 1,
    )
    proj_MC = histo2D_MC.ProjectionY(
        "TrkIsoMC",
        histo2D_MC.GetXaxis().FindBin(lowEnd),
        histo2D_MC.GetXaxis().FindBin(highEnd) - 1,
    )
    mcProjs = [
        histo2D.ProjectionY(
            "TrkIsoProj" + histo2D.GetName(),
            histo2D.GetXaxis().FindBin(lowEnd),
            histo2D.GetXaxis().FindBin(highEnd) - 1,
        )
        for histo2D in mc2DHists
    ]
    proj_Data = histo2D_Data.ProjectionY(
        "DataTrkIso",
        histo2D_Data.GetXaxis().FindBin(lowEnd),
        histo2D_Data.GetXaxis().FindBin(highEnd) - 1,
    )
    # number of HEEP electrons: pass trkIso < 5 GeV
    eleErr = ctypes.c_double()
    ele = proj_Electrons.IntegralAndError(
        proj_Electrons.GetXaxis().FindBin(0), proj_Electrons.GetXaxis().FindBin(5) - 1, eleErr
    )
    # denominator is all loose electrons in the region
    dataErr = ctypes.c_double()
    data = proj_Data.IntegralAndError(proj_Data.GetXaxis().GetFirst(), proj_Data.GetXaxis().GetLast(), dataErr)
    # MC for real electron subtraction: also trkIso < 5 GeV
    realEleErr = ctypes.c_double()
    realEle = proj_MC.IntegralAndError(
        proj_MC.GetXaxis().FindBin(0), proj_MC.GetXaxis().FindBin(5) - 1, realEleErr
    )
    realEleMCList = [
        proj.Integral(proj.GetXaxis().FindBin(0), proj.GetXaxis().FindBin(5) - 1)
        for proj in mcProjs
    ]
    numerator = ele - realEle
    numeratorErr = math.sqrt((eleErr.value)**2+(realEleErr.value)**2)
    if verbose:
        print "\t\t\tnumer:ele=", ele, "+/-", eleErr.value, ";numMC=", realEle, "+/-", realEleErr.value, "; numMCSubd=", numerator, "+/-", numeratorErr.value, "; denom:data=", data, "+/-", dataErr.value
        if realEle > ele:
            print "\t\t\tbreakdown of MC:"
            for i, sample in enumerate(mcNames):
                print "\t\t\t\t" + sample, realEleMCList[i]
            print "histo2D_Electrons has name:", histo2D_Electrons.GetName()
            print "histo2D_Electrons has entries:", histo2D_Electrons.GetEntries()
        print "\t\t\tFR=", numerator / data
    # make TrkIso for DY MC
    histo2D_DYElectrons = histDict[reg]["ZJets"][jets]
    proj_DYElectrons = histo2D_DYElectrons.ProjectionY(
        "DYElesTrkIso",
        histo2D_MC.GetXaxis().FindBin(lowEnd),
        histo2D_MC.GetXaxis().FindBin(highEnd) - 1,
    )
    if writeOutput:
        suffix = "_"+reg+"_"+jets+"Pt"+str(lowEnd)+"To"+str(highEnd)
        if not outputFile.cd("TrkIso_DYElectrons"):
            outputFile.mkdir("TrkIso_DYElectrons").cd()
        proj_DYElectrons.SetName(proj_DYElectrons.GetName()+suffix)
        proj_DYElectrons.Write()
    return numerator, numeratorErr, data, dataErr.value


def MakeFakeRatePlot(varName, reg, jets, bins, histDict, verbose=False, dataDriven=True):
    if verbose:
        print "MakeFakeRatePlot:varName=", varName, "reg=", reg, "jets=", jets, "dataDriven=", dataDriven
        sys.stdout.flush()
    # if reg == "Bar":
    #     bins = ptBinsBarrel
    # else:
    #     bins = ptBinsEndcap
    histNum = TH1D(
        varName+"_num_" + reg + "_" + jets + ("dataDriven" if dataDriven else "mcSub"),
        varName+"_num_" + reg + "_" + jets + ("dataDriven" if dataDriven else "mcSub"),
        len(bins) - 1,
        np.array(bins, dtype=float),
    )
    histDen = TH1D(
        varName+"_den_" + reg + "_" + jets + ("dataDriven" if dataDriven else "mcSub"),
        varName+"_den_" + reg + "_" + jets + ("dataDriven" if dataDriven else "mcSub"),
        len(bins) - 1,
        np.array(bins, dtype=float),
    )
    for index, binLow in enumerate(bins):
        if index >= (len(bins) - 1):
            break
        binHigh = bins[index + 1]
        if verbose:
            print "\tMakeFakeRatePlot:look at Pt:", str(binLow) + "-" + str(binHigh)
            sys.stdout.flush()
        if dataDriven:
            num, numErr, den, denErr = GetFakeRate(binLow, binHigh, reg, jets, histDict, verbose)
        else:
            num, numErr, den, denErr = GetFakeRateMCSub(binLow, binHigh, reg, jets, histDict, verbose)
        histNum.SetBinContent(histNum.GetXaxis().FindBin(binLow), num)
        histDen.SetBinContent(histDen.GetXaxis().FindBin(binLow), den)
        histNum.SetBinError(histNum.GetXaxis().FindBin(binLow), numErr)
        histDen.SetBinError(histDen.GetXaxis().FindBin(binLow), denErr)
        if verbose:
            print "\tMakeFakeRatePlot:num=", num, "den=", den
            print "\tMakeFakeRatePlot:set bin:", histNum.GetXaxis().FindBin(binLow), "to", num
            print "\tMakeFakeRatePlot:bin content is:", histNum.GetBinContent(
                histNum.GetXaxis().FindBin(binLow)
            )
            sys.stdout.flush()
    # c = TCanvas()
    # c.SetName('num'+reg+jets)
    # c.SetTitle('num'+reg+jets)
    # c.cd()
    # histNum.Draw()
    # c1 = TCanvas()
    # c1.SetName('den'+reg+jets)
    # c1.SetTitle('den'+reg+jets)
    # c1.cd()
    # histDen.Draw()

    graphFR = TGraphAsymmErrors()
    graphFR.BayesDivide(histNum, histDen)
    if analysisYear != 2018:
        graphFR.SetName(
            "fr{reg}_{jets}{dataDriven}".format(
                reg=reg, jets=jets, dataDriven="template" if dataDriven else "MCsub"
            )
        )
    else:
        # a bit nasty to hardcode the varName in the below
        graphFR.SetName(
            "fr{reg}_{var}_{jets}{dataDriven}".format(
                reg=reg, var=varName.replace("TrkIsoHEEP7vsHLTPt_", ""),
                jets=jets, dataDriven="template" if dataDriven else "MCsub"
            )
        )
    return graphFR, histNum, histDen


def MakeFRCanvas(plotList, titleList, canTitle):
    if len(titleList) != len(plotList):
        print "ERROR: titleList and plotList passed into MakeFRCanvas are different lengths!"
    can = TCanvas()
    can.cd()
    can.SetGridy()
    can.SetTitle(canTitle)
    can.SetName("y"+canTitle.lower().replace(',', '').replace(' ', '_').replace('>=', 'gte').replace('<', 'lt'))
    colorList = [1, kSpring-1, kAzure+1, kBlue, kGreen]
    markerStyleList = [8, 25, 23, 22]
    for i, plot in enumerate(plotList):
        if i == 0:
            plot.Draw("ap0")
            plot.GetXaxis().SetRangeUser(0, 1000)
            plot.GetXaxis().SetTitle("E_{T} [GeV]")
            plot.GetYaxis().SetRangeUser(0, 0.1)
            plot.GetYaxis().SetNdivisions(510)
            plot.GetYaxis().SetTitle("fake rate")
            plot.SetTitle(canTitle)
        if "graph" not in plot.ClassName().lower():
            plot.SetStats(0)
        plot.SetMarkerColor(colorList[i])
        plot.SetLineColor(colorList[i])
        plot.SetMarkerStyle(markerStyleList[i])
        plot.SetMarkerSize(0.9)
        if i == 0:
            plot.Draw("ap0")
        else:
            plot.Draw("psame0")
    # hsize=0.20
    # vsize=0.35
    # leg = TLegend(0.19,0.18,0.44,0.35)
    leg = TLegend(0.38, 0.71, 0.63, 0.88)
    if len(plotList) > 1:
        leg.SetFillColor(kWhite)
        leg.SetFillStyle(1001)
        leg.SetBorderSize(0)
        leg.SetShadowColor(10)
        leg.SetMargin(0.2)
        leg.SetTextFont(132)
        for i, title in enumerate(titleList):
            leg.AddEntry(plotList[i], title, "lp")
            # print "For canvas titled:", canTitle, ", added entry in legend for plot:", plotList[
            #     i
            # ].GetName(), "setting title to:", title
        leg.Draw()
    can.Modified()
    can.Update()
    return can, leg


def LoadHistosData(histos, varName, sample):
    for reg in detectorRegions:
        histos[reg] = {}
        for eType in electronTypes:
            histos[reg][eType] = {}
            for jet in jetBins:
                hist = tfile.Get(
                    histoBaseName.format(
                        sample=sample, type=eType, region=reg, jets=jet, varName=varName
                    )
                )
                histos[reg][eType][jet] = hist
                # print 'added histo:',hist.GetName(),'with entries:',hist.GetEntries(),'to','['+reg+']['+eType+']['+jet+']'


def LoadHistosMC(histos, varName):
    # we only care about electrons in the MC
    eType = "Electrons"
    for reg in detectorRegions:
        histos[reg]["MC"] = {}
        for name in mcNames:
            histos[reg][name] = {}
        for jet in jetBins:
            mcTotalHist = 0
            for i, name in enumerate(mcNames):
                histName = histoBaseName.format(
                    type=eType,
                    region=reg,
                    jets=jet,
                    sample=mcSamples[i],
                    varName=varName,
                )
                hist = tfile.Get(histName)
                if not hist:
                    print "ERROR: could not get hist '"+histName+"' from file: "+tfile.GetName()
                    exit(-1)
                histos[reg][name][jet] = hist
                if not mcTotalHist:
                    # mcTotalHist = hist.Clone()
                    mcTotalHist = copy.deepcopy(hist)
                    mcTotalHist.SetName(
                        histoBaseName.format(
                            type=eType,
                            region=reg,
                            jets=jet,
                            sample="MCTotal",
                            varName=varName,
                        )
                    )
                else:
                    mcTotalHist.Add(hist)
            histos[reg]["MC"][jet] = mcTotalHist


def GetJetBin(histName):
    if "Jet" in histName:
        return histName[histName.find("Jet") - 1: histName.find("Jet") + 3]
    else:
        return "0Jets"


# make FR 2D plot from FR TGraphAsymmErrors returned from MakeFakeRatePlot (frGraph)
# frGraph has x-axis: Pt and y-axis FR
# so we need to know the eta region here
def MakeFR2D(frGraph, reg, bins):
    if reg == "Bar":
        # bins = ptBinsBarrel
        etaToFill = 1
    elif reg == "End1":
        # bins = ptBinsEndcap
        etaToFill = 1.7
    elif reg == "End2":
        # bins = ptBinsEndcap
        etaToFill = 2
    else:
        print "ERROR: could not understand region given:", reg, "; return empty hist"
        return
    etaBins = [-2.5, -2.0, -1.566, -1.4442, 0, 1.4442, 1.566, 2.0, 2.5]
    frHist2d = TH2F(
        "test",
        "test",
        len(etaBins) - 1,
        np.array(etaBins, dtype=float),
        len(bins) - 1,
        np.array(bins, dtype=float),
    )
    # print 'consider frGraph with name {}'.format(frGraph.GetName())
    jets = GetJetBin(frGraph.GetName())
    name = "fr2D_{reg}_{jets}".format(reg=reg, jets=jets)
    frHist2d.SetNameTitle(name, name)
    frHist2d.GetXaxis().SetTitle("SuperCluster #eta")
    frHist2d.GetYaxis().SetTitle("p_{T} [GeV]")
    for iPoint in range(frGraph.GetN()):
        x = ctypes.c_double()
        y = ctypes.c_double()
        frGraph.GetPoint(iPoint, x, y)
        frHist2d.SetBinContent(
            frHist2d.GetXaxis().FindBin(etaToFill), frHist2d.GetYaxis().FindBin(x.value), y.value
        )
        # TODO set error
    return frHist2d


####################################################################################################
# RUN
####################################################################################################
# filename = "$LQDATA/nanoV6/2016/qcdFakeRateCalc/17jul2020/output_cutTable_lq_QCD_FakeRateCalculation/analysisClass_lq_QCD_FakeRateCalculation_plots.root"
filename = "$LQDATA/nanoV6/2017/qcdFakeRateCalc/20jul2020/output_cutTable_lq_QCD_FakeRateCalculation/analysisClass_lq_QCD_FakeRateCalculation_plots.root"
#filename = "$LQDATA/nanoV6/2018/qcdFakeRateCalc/20jul2020/output_cutTable_lq_QCD_FakeRateCalculation/analysisClass_lq_QCD_FakeRateCalculation_plots.root"

print "Opening file:", filename
tfile = TFile.Open(filename)
if not tfile or tfile.IsZombie():
    print "ERROR: TFile is zombie. Quitting here."
    exit(-1)

if '2016' in filename:
    analysisYear = 2016
elif '2017' in filename:
    analysisYear = 2017
elif '2018' in filename:
    analysisYear = 2018

outputFileName = "plots.root"
pdf_folder = "pdf"

gROOT.SetBatch(True)
writeOutput = True
doMuz = False  # do Muzamil's plots
doMCSubFR = True

histoBaseName = "histo2D__{sample}__{type}_{region}_{jets}{varName}"
dataSampleName = "QCDFakes_DATA"
electronTypes = ["Jets", "Electrons", "Total"]
# probably eventually expand to BarPlus, BarMinus, etc.
detectorRegions = ["Bar", "End1", "End2"]
regTitleDict = {}
# jetBins = ["", "1Jet_", "2Jet_", "3Jet_"]
jetBins = ["", "1Jet_", "2Jet_"]
# for MC
if '2016' in filename:
    mcSamples = [
         "ZJet_amcatnlo_ptBinned",
         # "ZJet_amcatnlo_Inc",
         # "WJet_amcatnlo_ptBinned",
         # "WJet_Madgraph_Inc",
         "WJet_amcatnlo_Inc",
         "TTbar_powheg",
         "SingleTop",
         "PhotonJets_Madgraph",
         "DIBOSON_nlo",
    ]
else:
    mcSamples = [
        "ZJet_jetAndPtBinned",
        "WJet_amcatnlo_jetBinned",
        "TTbar_powheg",
        "SingleTop",
        "PhotonJets_Madgraph",
        "DIBOSON_nlo",
    ]
mcNames = ["ZJets", "WJets", "TTBar", "ST", "GJets", "Diboson"]

# varNameList = ["TrkIsoHEEP7vsHLTPt_PAS", "TrkIsoHEEP7vsMTenu_PAS"]
if analysisYear != 2018:
    varNameList = ["TrkIsoHEEP7vsHLTPt_PAS"]
else:
    varNameList = [
            "TrkIsoHEEP7vsHLTPt_PAS", "TrkIsoHEEP7vsHLTPt_pre319077", "TrkIsoHEEP7vsHLTPt_post319077",
            "TrkIsoHEEP7vsHLTPt_noHEM_post319077", "TrkIsoHEEP7vsHLTPt_HEMonly_post319077"]

# ptBinsBarrel = [
# #    35,
# #    40,
#     45,
#     50,
#     60,
#     70,
#     80,
#     90,
#     100,
#     110,
#     130,
#     150,
#     170,
#     200,
#     250,
#     300,
#     400,
#     500,
#     600,
#     800,
#     1000,
# ]
# 2016: Photon22, 30, 36, 50, 75, 90, 120, 175
# 2017: Photon25, 33, 50, 75, 90, 120, 150, 175, 200
# 2018: Photon    33, 50, 75, 90, 120, 150, 175, 200
# ptBinsBarrel = [
#      45,
#      60,
#      75,
#      90,
#      105,
#      120,
#      135,
#      150,
#      170,
#      200,
#      250,
#      300,
#      400,
#      500,
#      600,
#      800,
#      1000,
#  ]
# ptBinsEndcap = [35, 50, 75, 100, 125, 150, 200, 250, 300, 350, 500, 1000]
# Z' binning -- 2016
# ptBinsEndcap = [36, 47, 50, 62, 75, 82, 90, 105, 120, 140, 175, 200, 225, 250, 300, 350, 400, 500, 600, 1000]
if analysisYear == 2016:
    ptBinsEndcap = [36, 50, 75, 90, 120, 140, 175, 200, 225, 250, 300, 350, 400, 500, 600, 1000]
else:
    ptBinsEndcap = [36, 50, 75, 90, 120, 150, 175, 200, 225, 250, 300, 350, 400, 500, 600, 1000]
# 2018
ptBinsEndcapHEM1516Only = [36, 50, 75, 90, 120, 150, 175, 200, 225, 250, 300, 350, 400, 500, 1000]
# ptBinsBarrel = ptBinsEndcap
ptBinsDict = {}
for varName in varNameList:
    ptBinsDict[varName] = {}
    for reg in detectorRegions:
        ptBinsDict[varName][reg] = ptBinsEndcap
if analysisYear == 2018:
    ptBinsDict["TrkIsoHEEP7vsHLTPt_HEMonly_post319077"]["End1"] = ptBinsEndcapHEM1516Only
    ptBinsDict["TrkIsoHEEP7vsHLTPt_HEMonly_post319077"]["End2"] = ptBinsEndcapHEM1516Only

allHistos = {}
for varName in varNameList:
    allHistos[varName] = {}
    LoadHistosData(allHistos[varName], varName, dataSampleName)
    LoadHistosMC(allHistos[varName], varName)

histList = []
myCanvases = []
if writeOutput:
    outputFile = TFile(outputFileName, "recreate")
    outputFile.cd()
# TEST end2 FR plots
# if doMuz:
#     # get Muzamil's hists
#     tfileMuzamilTwoJ = TFile(
#         "/afs/cern.ch/user/m/mbhat/work/public/Fakerate_files_2016/FR2D2JetScEt.root"
#     )
#     tfileMuzamilZeroJ = TFile(
#         "/afs/cern.ch/user/m/mbhat/work/public/Fakerate_files_2016/FR0JET_HLT.root"
#     )
#     muzHist2DZeroJ = tfileMuzamilZeroJ.Get("Endcap_Fake_Rate")
#     muzHist2DZeroJBar = tfileMuzamilZeroJ.Get("Barrel_Fake_Rate")
#     muzHistEnd2ZeroJ = muzHist2DZeroJ.ProjectionX("projMuzEnd2_0Jets", 2, 2)
#     muzHistEnd1ZeroJ = muzHist2DZeroJ.ProjectionX("projMuzEnd1_0Jets", 1, 1)
#     muzHistBarZeroJ = muzHist2DZeroJBar.ProjectionX("projMuzBar_0Jets")
# get Sam's hists
tfileZPrime = TFile(
    "/afs/cern.ch/user/s/scooper/work/public/Leptoquarks/QCDFakeRate/heepV7p0_2016_reminiAOD_noEleTrig_fakerate.root"
)
if tfileZPrime.IsZombie():
    print "ERROR: zprime tfile is zombie:", tfileZPrime.GetName()
    exit(-1)
zprimeHistEnd2ZeroJ = tfileZPrime.Get("frHistEEHigh")
zprimeHistEnd1ZeroJ = tfileZPrime.Get("frHistEELow")
zprimeHistBarZeroJ = tfileZPrime.Get("frHistEB")
histList.extend([zprimeHistEnd2ZeroJ, zprimeHistEnd1ZeroJ, zprimeHistBarZeroJ])
zprimeHistDict = {}
zprimeHistDict["End2"] = {}
zprimeHistDict["End2"][""] = zprimeHistEnd2ZeroJ
zprimeHistDict["End1"] = {}
zprimeHistDict["End1"][""] = zprimeHistEnd1ZeroJ
zprimeHistDict["Bar"] = {}
zprimeHistDict["Bar"][""] = zprimeHistBarZeroJ
print zprimeHistBarZeroJ.GetEntries()

# make list of histos to use for FR
# histos = [allHistos[varNameList[0]]]
#if '2018' in filename:
#    histos = allHistos

histDict = {}
numerHistDict = {}
denomHistDict = {}
for varName in allHistos:
    histos = allHistos[varName]
    histDict[varName] = {}
    numerHistDict[varName] = {}
    denomHistDict[varName] = {}
    for reg in detectorRegions:
        histDict[varName][reg] = {}
        numerHistDict[varName][reg] = {}
        denomHistDict[varName][reg] = {}
        bins = ptBinsDict[varName][reg]
        for jetBin in jetBins:
            histDict[varName][reg][jetBin] = {}
            numerHistDict[varName][reg][jetBin] = {}
            denomHistDict[varName][reg][jetBin] = {}
            histFR, histNum, histDen = MakeFakeRatePlot(
                    varName, reg, jetBin, bins, histos
            )
            histDict[varName][reg][jetBin]["data"] = histFR
            numerHistDict[varName][reg][jetBin]["data"] = histNum
            denomHistDict[varName][reg][jetBin]["data"] = histDen
            if doMCSubFR:
                histFRMC, histNumMC, histDenMC = MakeFakeRatePlot(
                        varName, reg, jetBin, bins, histos, dataDriven=False
                )
                histDict[varName][reg][jetBin]["mc"] = histFRMC
                numerHistDict[varName][reg][jetBin]["mc"] = histNumMC
                denomHistDict[varName][reg][jetBin]["mc"] = histDenMC

for varName in allHistos:
    for reg in detectorRegions:
        for jetBin in jetBins:
            varRegJetHistList = [histDict[varName][reg][jetBin]["data"]]
            titleList = ["Data-driven"]
            if doMCSubFR:
                varRegJetHistList.append(histDict[varName][reg][jetBin]["mc"])
                titleList.append("MCSub")
            if jetBin == "" and varName == "TrkIsoHEEP7vsHLTPt_PAS":
                varRegJetHistList.append(zprimeHistDict[reg][jetBin])
                titleList.append("2016 Zprime (E_{T}^{HLT})")
            myCanvases.append(
                MakeFRCanvas(
                    varRegJetHistList,
                    titleList,
                    GetCanvasTitle(varName, reg, jetBin)
                )
            )
            histList.extend(varRegJetHistList)

# # for drawing num/den hists
# numDenHistList = [histNumEnd2ZeroJ,histDenEnd2ZeroJ,histNumEnd1ZeroJ,histDenEnd1ZeroJ,histNumBarZeroJ,histDenBarZeroJ]
# for idx,hist in enumerate(numDenHistList):
#    if idx%2 != 0:
#        continue
#    if idx > len(numDenHistList)-2:
#        break
#    print 'consider hist:',hist.GetName(),'and',numDenHistList[idx+1].GetName()
#    c = TCanvas()
#    c.SetName(hist.GetName())
#    c.SetTitle(hist.GetTitle())
#    c.cd()
#    c.SetLogx()
#    c.SetLogy()
#    c.SetGridx()
#    c.SetGridy()
#    numDenHistList[idx+1].SetLineWidth(2)
#    numDenHistList[idx+1].SetLineColor(4)
#    numDenHistList[idx+1].SetMarkerColor(4)
#    numDenHistList[idx+1].Draw()
#    numDenHistList[idx+1].GetYaxis().SetRangeUser(2e-1,1e12)
#    hist.SetLineWidth(2)
#    hist.SetLineColor(2)
#    hist.SetMarkerColor(2)
#    hist.Draw('sames')

if writeOutput:
    outputFile.cd()
    if not os.path.isdir(pdf_folder) and pdf_folder != "":
        "Making directory", pdf_folder
        os.mkdir(pdf_folder)
for canLeg in myCanvases:
    canv = canLeg[0]
    canv.Draw()  # canvas
    # canLeg[-1][1].Draw() #legend
    if writeOutput:
        canv.Write()
        canv.Print(pdf_folder + "/" + canv.GetName()+".pdf")


##################################################
# make Et plot
##################################################
# histoNameZ = 'histo2D__ZJet_amcatnlo_ptBinned__Electrons_{region}_{jets}TrkIsoHEEP7vsPt_PAS'
# histoNameW = 'histo2D__WJet_amcatnlo_ptBinned__Electrons_{region}_{jets}TrkIsoHEEP7vsPt_PAS'
# histoNameData = 'histo2D__QCDFakes_DATA__Electrons_{region}_{jets}TrkIsoHEEP7vsPt_PAS'
histoNameDataLoose = "histo2D__QCDFakes_DATA__Total_{region}_{jets}TrkIsoHEEP7vsPt_PAS"
histosPt = {}
for reg in detectorRegions:
    histosPt[reg] = {}
    histosPt[reg]["ZElectrons"] = {}
    histosPt[reg]["WElectrons"] = {}
    histosPt[reg]["DataLooseElectrons"] = {}
    histosPt[reg]["DataElectrons"] = {}
    for jet in jetBins:
        # histo2D_MC = tfile.Get(histoNameZ.format(region=reg,jets=jet))
        # histo2D_MC = histos[reg]["ZJets"][jet]
        # proj_MC = histo2D_MC.ProjectionX(
        #     "EtZ",
        #     histo2D_MC.GetYaxis().FindBin(0),
        #     histo2D_MC.GetYaxis().FindBin(5) - 1,
        # )
        # histosPt[reg]["ZElectrons"][jet] = proj_MC
        # # histo2D_MC = tfile.Get(histoNameW.format(region=reg,jets=jet))
        # histo2D_MC = histos[reg]["WJets"][jet]
        # proj_MC = histo2D_MC.ProjectionX(
        #     "EtW",
        #     histo2D_MC.GetYaxis().FindBin(0),
        #     histo2D_MC.GetYaxis().FindBin(5) - 1,
        # )
        # histosPt[reg]["WElectrons"][jet] = proj_MC
        histo2D_data = tfile.Get(histoNameDataLoose.format(region=reg, jets=jet))
        proj_data = histo2D_data.ProjectionX(
            "EtDataLoose",
            histo2D_data.GetYaxis().FindBin(0),
            histo2D_data.GetYaxis().FindBin(5) - 1,
        )
        histosPt[reg]["DataLooseElectrons"][jet] = proj_data
        # histo2D_data = tfile.Get(histoNameData.format(region=reg,jets=jet))
        histo2D_data = histos[reg]["Electrons"][jet]
        proj_data = histo2D_data.ProjectionX(
            "EtData",
            histo2D_data.GetYaxis().FindBin(0),
            histo2D_data.GetYaxis().FindBin(5) - 1,
        )
        histosPt[reg]["DataElectrons"][jet] = proj_data
        # print 'added histo:',hist.GetName(),' to ','['+reg+']['+eType+']['+jet+']'

canvasEt = TCanvas()
canvasEt.cd()
canvasEt.SetLogy()
reg = "Bar"
jets = ""
rebinFactor = 20
histData = histosPt[reg]["DataElectrons"][jets]
histData.Rebin(rebinFactor)
histDataLoose = histosPt[reg]["DataLooseElectrons"][jets]
histDataLoose.Rebin(rebinFactor)
histDataLoose.SetMarkerColor(kRed)
histDataLoose.SetLineColor(kRed)
histDataLoose.SetLineStyle(2)
#
# stack = THStack()
# histZ = histosPt[reg]["ZElectrons"][jets]
# histZ.Rebin(rebinFactor)
# histZ.SetLineColor(7)
# histZ.SetLineWidth(2)
# histZ.SetFillColor(7)
# histZ.SetMarkerColor(7)
# stack.Add(histZ)
# histW = histosPt[reg]["WElectrons"][jets]
# # print 'W entries:',histW.GetEntries()
# histW.Rebin(rebinFactor)
# histW.SetLineColor(kBlue)
# histW.SetLineWidth(2)
# histW.SetFillColor(kBlue)
# histW.SetMarkerColor(kBlue)
# stack.Add(histW)
# stack.Draw("hist")
# stack.SetMaximum(1.2 * histData.GetMaximum())
# stack.SetMinimum(1e-1)
# stack.GetXaxis().SetTitle("Et [GeV]")
# stack.Draw("hist")
histData.Draw("psame")
histDataLoose.Draw("psame")
legEt = TLegend(0.38, 0.71, 0.63, 0.88)
legEt.SetFillColor(kWhite)
legEt.SetFillStyle(1001)
legEt.SetBorderSize(0)
legEt.SetShadowColor(10)
legEt.SetMargin(0.2)
legEt.SetTextFont(132)
# legEt.AddEntry(histZ, "ZJets", "lp")
# legEt.AddEntry(histW, "WJets", "lp")
legEt.AddEntry(histData, "Data", "lp")
legEt.AddEntry(histDataLoose, "Data (loose e)", "lp")
legEt.Draw()
canvasEt.Modified()
canvasEt.Update()
low = 800
high = 1000
print "integrals in 800-1000:"
print "dataLoose=", histDataLoose.Integral(
    histDataLoose.FindBin(low), histDataLoose.FindBin(high) - 1
)
print "dataEle=", histData.Integral(histData.FindBin(low), histData.FindBin(high) - 1)
# print "WEle=", histW.Integral(histW.FindBin(low), histW.FindBin(high) - 1)
# print "ZEle=", histZ.Integral(histZ.FindBin(low), histZ.FindBin(high) - 1)


if writeOutput:
    endcap2dHists = []
    for hist in histList:
        hist.Write()
        if "template" in hist.GetName():
            reg = hist.GetName()[2: hist.GetName().find("_")]
            bins = GetXBinsFromGraph(hist)
            hist2d = MakeFR2D(hist, reg, bins)
            if "End" in hist2d.GetName():
                endcap2dHists.append(hist2d)
            else:
                hist2d.Write()
    histNamesDone = []
    for i in range(len(endcap2dHists)):
        hist = endcap2dHists[i]
        name = hist.GetName()
        if name in histNamesDone:
            continue
        reg = name[2: name.find("_")]
        jets = GetJetBin(name)
        for j in range(len(endcap2dHists)):
            hist2 = endcap2dHists[j]
            name2 = hist2.GetName()
            if name == name2:
                continue
            reg2 = name2[2: name2.find("_")]
            jets2 = GetJetBin(name2)
            if reg == reg2 and jets == jets2:
                hist.Add(hist2)
                histNamesDone.extend([name, name2])
                name = name.replace("End1", "End")
                name = name.replace("End2", "End")
                hist.SetNameTitle(name, name)
                hist.Write()
                break

    outputFile.Close()

# histo2D_DY_zeroJ = tfile.Get('histo2D__ZJet_amcatnlo_ptBinned__Electrons_End1_TrkIsoHEEP7vsPt_PAS')
# histo2D_Jets_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Jets_End1_TrkIsoHEEP7vsPt_PAS')
# histo2D_Electrons_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Electrons_End1_TrkIsoHEEP7vsPt_PAS')
# histo2D_Data_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Total_End1_TrkIsoHEEP7vsPt_PAS')
##histo2D_DY_zeroJ = tfile.Get('histo2D__ZJet_amcatnlo_ptBinned__Electrons_bar_TrkIsoHEEP7vsPt_PAS')
##histo2D_Jets_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Jets_Bar_TrkIsoHEEP7vsPt_PAS')
##histo2D_Electrons_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Electrons_Bar_TrkIsoHEEP7vsPt_PAS')
##histo2D_Data_zeroJ = tfile.Get('histo2D__QCDFakes_DATA__Total_Bar_TrkIsoHEEP7vsPt_PAS')
#
# lowEnd=350 # GeV
# highEnd=500
##lowEnd=60 # GeV
##highEnd=70
# print
# print
# print 'Template plot section'
# print
# print 'Considering electron Pt:',str(lowEnd)+'-'+str(highEnd),'GeV'
# proj_DY_zeroJ = histo2D_DY_zeroJ.ProjectionY('DYZeroJTrkIso',histo2D_DY_zeroJ.GetXaxis().FindBin(lowEnd),histo2D_DY_zeroJ.GetXaxis().FindBin(highEnd))
# proj_Electrons_zeroJ = histo2D_Electrons_zeroJ.ProjectionY('ElesZeroJTrkIso',histo2D_Electrons_zeroJ.GetXaxis().FindBin(lowEnd),histo2D_Electrons_zeroJ.GetXaxis().FindBin(highEnd)-1)
# proj_Jets_zeroJ = histo2D_Jets_zeroJ.ProjectionY('JetsZeroJTrkIso',histo2D_Jets_zeroJ.GetXaxis().FindBin(lowEnd),histo2D_Jets_zeroJ.GetXaxis().FindBin(highEnd)-1)
# proj_Data_zeroJ = histo2D_Data_zeroJ.ProjectionY('DataZeroJTrkIso',histo2D_Data_zeroJ.GetXaxis().FindBin(lowEnd),histo2D_Data_zeroJ.GetXaxis().FindBin(highEnd)-1)
# print 'Using bin range:',histo2D_Data_zeroJ.GetXaxis().FindBin(lowEnd),'to',histo2D_Data_zeroJ.GetXaxis().FindBin(highEnd)-1
# can1 = TCanvas()
# can1.cd()
# proj_Data_zeroJ.Draw()
#
# proj_DY_zeroJ.SetLineColor(kBlue)
# proj_DY_zeroJ.SetMarkerColor(kBlue)
# proj_DY_zeroJ.SetLineStyle(2)
# proj_DY_zeroJ.SetLineWidth(2)
# proj_Jets_zeroJ.SetLineColor(kRed)
# proj_Jets_zeroJ.SetMarkerColor(kRed)
# proj_Jets_zeroJ.SetLineStyle(2)
# proj_Jets_zeroJ.SetLineWidth(2)
# proj_Data_zeroJ.SetLineColor(kBlue)
# proj_Data_zeroJ.SetMarkerColor(kBlue)
#
## ugly copy paste
# histo2D_DY_twoJ = tfile.Get('histo2D__ZJet_amcatnlo_ptBinned__Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS')
# histo2D_Jets_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Jets_End1_2Jet_TrkIsoHEEP7vsPt_PAS')
# histo2D_Electrons_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Electrons_End1_2Jet_TrkIsoHEEP7vsPt_PAS')
# histo2D_Data_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Total_End1_2Jet_TrkIsoHEEP7vsPt_PAS')
##histo2D_Jets_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Jets_Bar_2Jet_TrkIsoHEEP7vsPt_PAS')
##histo2D_Electrons_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Electrons_Bar_2Jet_TrkIsoHEEP7vsPt_PAS')
##histo2D_Data_twoJ = tfile.Get('histo2D__QCDFakes_DATA__Total_Bar_2Jet_TrkIsoHEEP7vsPt_PAS')
#
# proj_DY_twoJ = histo2D_DY_twoJ.ProjectionY('DYTwoJTrkIso',histo2D_DY_twoJ.GetXaxis().FindBin(lowEnd),histo2D_DY_twoJ.GetXaxis().FindBin(highEnd))
# proj_Electrons_twoJ = histo2D_Electrons_twoJ.ProjectionY('ElesTwoJTrkIso',histo2D_Electrons_twoJ.GetXaxis().FindBin(lowEnd),histo2D_Electrons_twoJ.GetXaxis().FindBin(highEnd)-1)
# proj_Jets_twoJ = histo2D_Jets_twoJ.ProjectionY('JetsTwoJTrkIso',histo2D_Jets_twoJ.GetXaxis().FindBin(lowEnd),histo2D_Jets_twoJ.GetXaxis().FindBin(highEnd)-1)
# proj_Data_twoJ = histo2D_Data_twoJ.ProjectionY('DataTwoJTrkIso',histo2D_Data_twoJ.GetXaxis().FindBin(lowEnd),histo2D_Data_twoJ.GetXaxis().FindBin(highEnd)-1)
#
# proj_DY_twoJ.SetLineColor(kBlue)
# proj_DY_twoJ.SetMarkerColor(kBlue)
# proj_DY_twoJ.SetLineStyle(2)
# proj_DY_twoJ.SetLineWidth(2)
# proj_Jets_twoJ.SetLineColor(kRed)
# proj_Jets_twoJ.SetMarkerColor(kRed)
# proj_Jets_twoJ.SetLineStyle(2)
# proj_Jets_twoJ.SetLineWidth(2)
# proj_Data_twoJ.SetLineColor(kBlue)
# proj_Data_twoJ.SetMarkerColor(kBlue)
#
#
# dyZeroJ = proj_DY_zeroJ.Integral(proj_DY_zeroJ.GetXaxis().FindBin(10),proj_DY_zeroJ.GetXaxis().FindBin(20)-1)
# eleZeroJ = proj_Electrons_zeroJ.Integral(proj_Electrons_zeroJ.GetXaxis().FindBin(10),proj_Electrons_zeroJ.GetXaxis().FindBin(20)-1)
# dataZeroJ_SB = proj_Data_zeroJ.Integral(proj_Data_zeroJ.GetXaxis().FindBin(10),proj_Data_zeroJ.GetXaxis().FindBin(20)-1)
# dataZeroJ = proj_Data_zeroJ.Integral()
# jetsZeroJ_SB = proj_Jets_zeroJ.Integral(proj_Jets_zeroJ.GetXaxis().FindBin(10),proj_Jets_zeroJ.GetXaxis().FindBin(20)-1)
# jetsZeroJ_SR = proj_Jets_zeroJ.Integral(proj_Jets_zeroJ.GetXaxis().FindBin(0),proj_Jets_zeroJ.GetXaxis().FindBin(5)-1)
##print 'endcap1 N_jet>=0: Nele=',eleZeroJ,'data=',dataZeroJ,'contam=',eleZeroJ/dataZeroJ
# print 'barrel N_jet>=0: Nele=',eleZeroJ,'data=',dataZeroJ,'contam=',eleZeroJ/dataZeroJ_SB
# print 'barrel N_jet>=0: nHEEPprime=',eleZeroJ,'jetsSR=',jetsZeroJ_SR,'jetsSB=',jetsZeroJ_SB,'data=',dataZeroJ
# print 'FR=',((jetsZeroJ_SR/jetsZeroJ_SB) * eleZeroJ)/dataZeroJ
#
# dyTwoJ = proj_DY_twoJ.Integral(proj_DY_twoJ.GetXaxis().FindBin(10),proj_DY_twoJ.GetXaxis().FindBin(20)-1)
# eleTwoJ = proj_Electrons_twoJ.Integral(proj_Electrons_twoJ.GetXaxis().FindBin(10),proj_Electrons_twoJ.GetXaxis().FindBin(20)-1)
# dataTwoJ_SB = proj_Data_twoJ.Integral(proj_Data_twoJ.GetXaxis().FindBin(10),proj_Data_twoJ.GetXaxis().FindBin(20)-1)
# dataTwoJ = proj_Data_twoJ.Integral()
# jetsTwoJ_SB = proj_Jets_twoJ.Integral(proj_Jets_twoJ.GetXaxis().FindBin(10),proj_Jets_twoJ.GetXaxis().FindBin(20)-1)
# jetsTwoJ_SR = proj_Jets_twoJ.Integral(proj_Jets_twoJ.GetXaxis().FindBin(0),proj_Jets_twoJ.GetXaxis().FindBin(5)-1)
##print 'endcap1 N_jet>=2: Nele=',eleTwoJ,'data=',dataTwoJ,'contam=',eleTwoJ/dataTwoJ
# print 'barrel N_jet>=2: Nele=',eleTwoJ,'data=',dataTwoJ,'contam=',eleTwoJ/dataTwoJ_SB
# print 'barrel N_jet>=2: nHEEPprime=',eleTwoJ,'jetsSR=',jetsTwoJ_SR,'jetsSB=',jetsTwoJ_SB,'data=',dataTwoJ
# print 'FR=',((jetsTwoJ_SR/jetsTwoJ_SB) * eleTwoJ)/dataTwoJ
#
#
# can = TCanvas()
# can.cd()
# can.SetLogy()
# proj_Jets_zeroJ.Scale(eleZeroJ/jetsZeroJ_SB)
# proj_JetsEle_zeroJ = proj_Jets_zeroJ.Clone()
# proj_JetsEle_zeroJ.SetName('JetsEleZeroJTrkIso')
# proj_JetsEle_zeroJ.Add(proj_Electrons_zeroJ)
# proj_JetsEle_zeroJ.SetLineStyle(2)
# proj_JetsEle_zeroJ.SetLineWidth(1)
# proj_JetsEle_zeroJ.SetLineColor(1)
# proj_JetsEle_zeroJ.SetMarkerColor(1)
# proj_Data_zeroJ.Draw()
##proj_Data_zeroJ.GetXaxis().SetRangeUser(0,50)
# proj_Data_zeroJ.GetYaxis().SetRangeUser(1e0,2e7)
# proj_DY_zeroJ.Draw('samehist')
# proj_Jets_zeroJ.Draw('samehist')
# proj_JetsEle_zeroJ.Draw('samehist')
# proj_Electrons_zeroJ.Draw('samehist')
#
# can2 = TCanvas()
# can2.cd()
# can2.SetLogy()
# proj_Jets_twoJ.Scale(eleTwoJ/jetsTwoJ_SB)
# proj_JetsEle_twoJ = proj_Jets_twoJ.Clone()
# proj_JetsEle_twoJ.SetName('JetsEleTwoJTrkIso')
# proj_JetsEle_twoJ.Add(proj_Electrons_twoJ)
# proj_JetsEle_twoJ.SetLineStyle(2)
# proj_JetsEle_twoJ.SetLineWidth(1)
# proj_JetsEle_twoJ.SetLineColor(1)
# proj_JetsEle_twoJ.SetMarkerColor(1)
# proj_Data_twoJ.Draw()
##proj_Data_twoJ.GetXaxis().SetRangeUser(0,50)
# proj_Data_twoJ.GetYaxis().SetRangeUser(1e0,2e7)
# proj_DY_twoJ.Draw('samehist')
# proj_Jets_twoJ.Draw('samehist')
# proj_JetsEle_twoJ.Draw('samehist')
# proj_Electrons_twoJ.Draw('samehist')
#
