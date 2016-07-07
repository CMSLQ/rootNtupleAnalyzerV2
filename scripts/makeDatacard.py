#!/usr/bin/env python

import os, sys, math
import subprocess as sp
from optparse import OptionParser
from prettytable import PrettyTable
from ROOT import *

from combineCommon import *


def GetSystDictFromFile(filename):
    # go custom text parsing :`(
    # format is like:
    # LQ300  :     0.0152215
    # selection point, 100*(deltaX/X) [rel. change in %]
    systDict = {}
    with open(filename,'r') as thisFile:
        for line in thisFile:
            line = line.strip()
            #print 'line=',line,'; with length=',len(line)
            if len(line)==0:
                continue
            #print 'line.strip()="'+line.strip()+'"'
            #print 'line.strip().split(":")=',line.strip().split(':')
            items = line.split(':')
            #print 'items[0].strip()='+items[0].strip()
            #print 'items[1].strip()='+items[1].strip()
            selectionPoint = items[0].strip()
            if '_' in selectionPoint:
                bkgName = selectionPoint.split('_')[1]
                if not bkgName in syst_background_names:
                    print 'ERROR: unknown background named:',bgkName,'not found in list of systematics background names:',syst_background_names
                selectionPoint = selectionPoint.split('_')[0]
                if not bkgName in systDict.keys():
                    systDict[bkgName] = {}
                systDict[bkgName][selectionPoint] = float(items[1].strip())/100.0
            # signal
            systDict[selectionPoint] = float(items[1].strip())/100.0
    return systDict

def FillSystDicts(systNames,isBackground=True):
    systDict = {}
    for syst in systNames:
        if isBackground:
          filePath = systematics_filepath+syst+'_sys.dat'
        else:
          filePath = systematics_filepath+'LQ'+syst+'_sys.dat'
        thisSystDict = GetSystDictFromFile(filePath)
        # this will give the form (for background):
        #   systDict['Trigger'][bkgname]['LQXXXX'] = value
        systDict[syst] = thisSystDict
    return systDict

def RoundToN(x, n):
    #if n < 1:
    #    raise ValueError("can't round to less than 1 sig digit!")
    ## number of digits given by n
    #return "%.*e" % (n-1, x)
    if isinstance(x,float):
        return round(x,n)
    else:
        return x

def GetTableEntryStr(evts,errStat,errSyst=0):
    if evts=='-':
      return evts
    # rounding
    evts = RoundToN(evts,2)
    errStat = RoundToN(errStat,2)
    if errSyst==0:
      return str(evts)+' +/- '+str(errStat)
    else:
      errSyst = RoundToN(errSyst,2)
      return str(evts)+' +/- '+str(errStat)+' +/- '+str(errSyst)

def GetXSecTimesIntLumi(sampleNameFromDataset):
    #print 'GetXSecTimesIntLumi(',sampleNameFromDataset+')'
    xsection = float(lookupXSection(sampleNameFromDataset,xsectionDict))
    intLumiF = float(intLumi)
    return xsection*intLumiF

def CalculateScaledRateError(sampleNameFromDataset, N_unscaled_tot, N_unscaled_pass_entries, N_unscaled_pass_integral, doScaling=True):
    #print 'CalculateScaledRateError(',sampleNameFromDataset, N_unscaled_tot, N_unscaled_pass_entries, N_unscaled_pass_integral,')'
    # binomial error
    p = N_unscaled_pass_entries/N_unscaled_tot
    q = 1-p
    w = N_unscaled_pass_integral/N_unscaled_pass_entries if N_unscaled_pass_entries != 0 else 0.0
    unscaledRateError = N_unscaled_tot*w*math.sqrt(p*q/N_unscaled_tot)
    if doScaling:
        xsecTimesIntLumi = GetXSecTimesIntLumi(sampleNameFromDataset)
        scaledRateError=unscaledRateError*(xsecTimesIntLumi/N_unscaled_tot)
    else:
        scaledRateError=unscaledRateError
    return scaledRateError


def FindUnscaledSampleRootFile(sampleName, isQCD=False):
  #print 'FindUnscaledSampleRootFile('+sampleName+')'
  #filePath
  analysisCode = 'analysisClass_lq_eejj' if not isQCD else 'analysisClass_lq_eejj_QCD'
  reducedSkimStrings = ['_reduced_skim','_pythia8_reduced_skim']
  for redSkimStr in reducedSkimStrings:
      rootFilename = filePath + "/" + analysisCode + "___" + sampleName + redSkimStr + ".root"
      if os.path.isfile(rootFilename):
         return rootFilename
  print "ERROR:  could not find unscaled root file for sample=",sampleName
  print 'Tried:',[filePath + "/" + analysisCode + "___" + sampleName + redSkimStr + ".root" for redSkimStr in reducedSkimStrings]
  print "Exiting..."
  exit(-1)
  #rootFilename = filePath + "/" + analysisCode + "___" + sampleName + reducedSkimStrings[0] + ".root"
  #while not os.path.isfile(rootFilename) and index<len(reducedSkimStrings):
  #    index+=1
  #    rootFilename = filePath + "/" + analysisCode + "___" + sampleName + reducedSkimStrings[index] + ".root"
  #    #if not os.path.isfile(rootFilename):
  #    #    print "ERROR: file " + rootFilename + " not found"
  #    #    print "exiting..."
  #    #    exit(-1)
  #return rootFilename


def GetRatesAndErrors(unscaledRootFile,combinedRootFile,unscaledTotalEvts,sampleName,selection,isDataOrQCD=False):
    #print 'GetRatesAndErrors(',unscaledRootFile,combinedRootFile,unscaledTotalEvts,sampleName,selection,')'
    if selection=='preselection':
        selection = 'PAS'
    #mejHist = combinedRootFile.Get('histo1D__'+sampleName+'__Mej_selected_min_'+selection)
    #if not mejHist:
    #  print 'ERROR: could not find hist','histo1D__'+sampleName+'__Mej_selected_min_'+selection,' in file:',combinedRootFile.GetName()
    #  print 'EXIT'
    #  exit(-1)
    #rate = mejHist.Integral()
    mejUnscaledHist = unscaledRootFile.Get('Mej_selected_min_'+selection)
    if not mejUnscaledHist:
      print 'ERROR: could not find hist','Mej_selected_min_'+selection,' in file:',unscaledRootFile.GetName()
      print 'EXIT'
      exit(-1)
    unscaledInt = mejUnscaledHist.Integral()
    unscaledRate = mejUnscaledHist.GetEntries()
    xsecTimesIntLumi = GetXSecTimesIntLumi(sampleName)
    if not isDataOrQCD:
        rate = unscaledInt*xsecTimesIntLumi/unscaledTotalEvts
        rateErr = CalculateScaledRateError(sampleName,unscaledTotalEvts,unscaledRate,unscaledInt)
    else:
        rate = unscaledInt
        rateErr = CalculateScaledRateError(sampleName,unscaledTotalEvts,unscaledRate,unscaledInt,False)
    #print 'INFO: hist','Mej_selected_min_'+selection,' in file:',unscaledRootFile.GetName()
    #print 'unscaledRate=',unscaledRate,'unscaled entries=',mejUnscaledHist.GetEntries()
    #print 'xsecTimesIntLumi=',xsecTimesIntLumi,'unscaledInt=',unscaledInt,'unscaledRate=',unscaledRate,'unscaledTotalEvts=',unscaledTotalEvts,'rate=unscaledInt*xsecTimesIntLumi/unscaledTotalEvts=',rate
    return rate,rateErr,unscaledRate

def GetUnscaledTotalEvents(unscaledRootFile):
    unscaledEvtsHist = unscaledRootFile.Get('EventsPassingCuts')
    unscaledTotalEvts = unscaledEvtsHist.GetBinContent(1)
    return unscaledTotalEvts

def FillDicts(rootFilename,qcdRootFilename):
    qcdTFile = TFile(qcdRootFilename)
    tfile = TFile(rootFilename)

    # backgrounds
    for i_bkg,bkg_name in enumerate(background_names):
        scaledRootFile = ''
        if not 'QCD' in bkg_name:
          scaledRootFile = tfile
          isQCD=False
        else:
          scaledRootFile = qcdTFile
          isQCD=True
        sampleList = dictSamples[bkg_name]
        sampleRate = 0
        sampleRateErr = 0
        sampleUnscaledRate = 0
        sampleUnscaledTotalEvts = 0
        #print 'PRESELECTION bkg_bame=',bkg_name
        for bkgSample in sampleList:
            bkgUnscaledRootFilename = FindUnscaledSampleRootFile(bkgSample,isQCD)
            bkgUnscaledRootFile = TFile.Open(bkgUnscaledRootFilename)
            if not bkgUnscaledRootFile:
              print 'ERROR: something happened when trying to open the file:',bkgUnscaledRootFilename
              exit(-1)
            unscaledTotalEvts = GetUnscaledTotalEvents(bkgUnscaledRootFile)
            sampleUnscaledTotalEvts+=unscaledTotalEvts
            # preselection
            #print '------>Call GetRatesAndErrors for sampleName=',bkgSample
            rate,rateErr,unscaledRate = GetRatesAndErrors(bkgUnscaledRootFile,scaledRootFile,unscaledTotalEvts,bkgSample,'preselection',isQCD)
            #print '------>rate=',rate,'rateErr=',rateErr,'unscaledRate=',unscaledRate
            sampleRate+=rate
            sampleUnscaledRate+=unscaledRate
            sampleRateErr+=(rateErr*rateErr)
            bkgUnscaledRootFile.Close()
        sampleRateErr = math.sqrt(sampleRateErr)
        #print 'sampleRate:',sampleRate,'sampleRateErr=',sampleRateErr,'sampleUnscaledRate=',sampleUnscaledRate
        bkgRatesDict = {}
        bkgRatesDict['preselection'] = sampleRate
        bkgRateErrsDict = {}
        bkgRateErrsDict['preselection'] = sampleRateErr
        bkgUnscaledRatesDict = {}
        bkgUnscaledRatesDict['preselection'] = sampleUnscaledRate
        bkgTotalEvts = sampleUnscaledTotalEvts
        # final selections
        for i_signal_name, signal_name in enumerate(signal_names):
            for i_mass_point, mass_point in enumerate(mass_points):
                selectionName = 'LQ'+mass_point
                sampleList = dictSamples[bkg_name]
                sampleRate = 0
                sampleRateErr = 0
                sampleUnscaledRate = 0
                #print selectionName,'bkg_bame=',bkg_name
                for bkgSample in sampleList:
                    bkgUnscaledRootFilename = FindUnscaledSampleRootFile(bkgSample,isQCD)
                    bkgUnscaledRootFile = TFile.Open(bkgUnscaledRootFilename)
                    unscaledTotalEvts = GetUnscaledTotalEvents(bkgUnscaledRootFile)
                    sampleUnscaledTotalEvts+=unscaledTotalEvts
                    # preselection
                    #print '------>Call GetRatesAndErrors for sampleName=',bkgSample
                    rate,rateErr,unscaledRate = GetRatesAndErrors(bkgUnscaledRootFile,scaledRootFile,unscaledTotalEvts,bkgSample,selectionName,isQCD)
                    #print '------>rate=',rate,'rateErr=',rateErr,'unscaledRate=',unscaledRate
                    #if isQCD:
                    #  print 'for sample:',bkgSample,'got unscaled entries=',unscaledRate
                    sampleRate+=rate
                    sampleUnscaledRate+=unscaledRate
                    sampleRateErr+=(rateErr*rateErr)
                    bkgUnscaledRootFile.Close()
                sampleRateErr = math.sqrt(sampleRateErr)
                #print 'sampleRate:',sampleRate,'sampleRateErr=',sampleRateErr,'sampleUnscaledRate=',sampleUnscaledRate
                bkgRatesDict[selectionName] = sampleRate
                bkgRateErrsDict[selectionName] = sampleRateErr
                bkgUnscaledRatesDict[selectionName] = sampleUnscaledRate
        # fill full dicts
        d_background_rates[bkg_name] = bkgRatesDict
        d_background_rateErrs[bkg_name] = bkgRateErrsDict
        d_background_unscaledRates[bkg_name] = bkgUnscaledRatesDict
        d_background_totalEvents[bkg_name] = bkgTotalEvts

    # signals
    for i_signal_name, signal_name in enumerate(signal_names):
        for i_mass_point, mass_point in enumerate(mass_points):
          if 'BetaHalf' in signal_name:
              signalNameForFile = 'LQToUE_' #FIXME
          else:
              signalNameForFile = 'LQToUE_M-'+mass_point+'_BetaOne'
          fullSignalName = signal_name+mass_point
          unscaledRootFilename = FindUnscaledSampleRootFile(signalNameForFile)
          unscaledRootFile = TFile.Open(unscaledRootFilename)
          unscaledTotalEvts = GetUnscaledTotalEvents(unscaledRootFile)
          # preselection
          rate,rateErr,unscaledRate = GetRatesAndErrors(unscaledRootFile,tfile,unscaledTotalEvts,signalNameForFile,'preselection')
          sigRatesDict = {}
          sigRatesDict['preselection'] = rate
          sigRateErrsDict = {}
          sigRateErrsDict['preselection'] = rateErr
          sigUnscaledRatesDict = {}
          sigUnscaledRatesDict['preselection'] = unscaledRate
          sigTotalEvts = unscaledTotalEvts
          # final selection
          for imp, mp in enumerate(mass_points):
              signalSelName = signal_name+mp
              selectionName = 'LQ'+mp
              rate,rateErr,unscaledRate = GetRatesAndErrors(unscaledRootFile,tfile,unscaledTotalEvts,signalNameForFile,selectionName)
              sigRatesDict[selectionName] = rate
              sigRateErrsDict[selectionName] = rateErr
              sigUnscaledRatesDict[selectionName] = unscaledRate
          unscaledRootFile.Close()

          # fill full dicts
          signalFullName = signal_name + mass_point
          d_signal_rates[signalFullName] = sigRatesDict
          d_signal_rateErrs[signalFullName] = sigRateErrsDict
          d_signal_unscaledRates[signalFullName] = sigUnscaledRatesDict
          d_signal_totalEvents[signalFullName] = sigTotalEvts

    qcdTFile.Close()
    tfile.Close()
  

###################################################################################################
# CONFIGURABLES
###################################################################################################

#signal_names = [ "LQ_BetaHalf_M", "LQ_M" ] 
signal_names = [ "LQ_M_" ] 
#FIXME add 200 GeV point
#mass_points = [str(i) for i in range(300,2050,50)] # go from 300-2000 in 50 GeV steps
#mass_points = [str(i) for i in range(300,1550,50)] # go from 300-1500 in 50 GeV steps
mass_points = [str(i) for i in range(200,1550,50)] # go from 300-1500 in 50 GeV steps
#systematics = [ "jes", "ees", "shape", "norm", "lumi", "eer", "jer", "pu", "ereco", "pdf" ]
systematicsNamesBackground = [ "Trigger", "Reco", "PU", "PDF", "Lumi", "JER", "JEC", "HEEP", "E_scale", "DYShape", "TTShape" ]
systematicsNamesSignal = [ "Trigger", "Reco", "PU", "PDF", "Lumi", "JER", "JEC", "HEEP", "E_scale" ]
#FIXME systematics
systematics = []
background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTbar_Madgraph", "WJet_Madgraph_HT", "ZJet_Madgraph_HT", "DIBOSON","SingleTop"  ]
# background names for systs
syst_background_names = ['GJets', 'QCDFakes_DATA', 'TTbar', 'WJets', 'DY', 'Diboson', 'Singletop']
# XXX FIXME TEST
maxLQselectionBkg = 'LQ1500' # max background selection point used

# add DYnorm and TTnorm by hand
# for z+jets, we have 0.026 stat error on the scale factor of 1
# we need to add 10% additional
# deltaX/X = (1+0.026+0.1*1)/1 - 1
extraZJetNormSyst = 0.1
zJetNormDeltaXOverX = (1+0.026+extraZJetNormSyst*1)/1 - 1
# for ttbar, we have 0.037 stat error on the scale factor of 0.815
# we need to add 10% additional
# deltaX/X = (0.815+0.037+0.1*0.815)/0.815 - 1
#extraTTbarNormSyst = 0.1
extraTTbarNormSyst = 0.1
ttBarNormDeltaXOverX = (0.815+0.037+extraTTbarNormSyst*0.815)/0.815 - 1
# QCDNorm is 0.40
qcdNormDeltaXOverX = 0.40

n_background = len ( background_names  )
#n_systematics = len ( systematics ) + n_background + 1
# all bkg systematics, plus stat 'systs' for all bkg plus signal plus 3 backNormSysts
n_systematics = len ( systematicsNamesBackground ) + n_background + 1 + 3
n_channels = 1

d_background_rates = {}
d_background_rateErrs = {}
d_background_unscaledRates = {}
d_background_totalEvents = {}
d_signal_rates = {}
d_signal_rateErrs = {}
d_signal_unscaledRates = {}
d_signal_totalEvents = {}

inputList = os.environ["LQANA"]+'/config/TestCombinationMay4/inputListAllCurrent.txt'
sampleListForMerging = os.environ["LQANA"]+'/config/sampleListForMerging_13TeV_eejj.txt'
sampleListForMergingQCD = os.environ["LQANA"]+'/config/sampleListForMerging_13TeV_eejj_QCD.txt'
xsection = os.environ["LQANA"]+'/versionsOfAnalysis_eejj/1jun_ttbarRescale/xsection_13TeV_2015_TTbarRescale.txt'
intLumi = 2570

#filePath = os.environ["LQDATA"] + '/RunII/eejj_analysis_finalSelsUnbugged_24may2016/output_cutTable_lq_eejj/'
# nominal
filePath = os.environ["LQDATA"] + '/RunII//eejj_analysis_ttbarRescaleFinalSels_2jun2016/output_cutTable_lq_eejj/'
#filePath = os.environ["LQDATA"] + '/RunII/eejj_analysis_ttbarRescaleFinalSels__asymptoticOpt_19jun2016/output_cutTable_lq_eejj_asymptoticOpt/'
#filePath = os.environ["LQDATA"] + '/RunII/eejj_analysis_ttbarRescaleFinalSels_zStBiasCorrDYJ_28jun2016/output_cutTable_lq_eejj/'
dataMC_filepath   = filePath+'analysisClass_lq_eejj_plots.root'
qcd_data_filepath = filePath+'analysisClass_lq_eejj_QCD_plots.root'
systematics_filepath = '/afs/cern.ch/user/m/mbhat/work/public/Systematics_txtfiles_07_06_2016/'


###################################################################################################
# RUN
###################################################################################################

#---Check if sampleListForMerging file exist
if(os.path.isfile(sampleListForMerging) == False):
    print "ERROR: file " + sampleListForMerging + " not found"
    print "exiting..."
    sys.exit()

#---Check if sampleListForMergingQCD file exist
if(os.path.isfile(sampleListForMergingQCD) == False):
    print "ERROR: file " + sampleListForMergingQCD + " not found"
    print "exiting..."
    sys.exit()

#---Check if xsection file exist
if(os.path.isfile(xsection) == False):
    print "ERROR: file " + xsection + " not found"
    print "exiting..."
    sys.exit()

print 'Launched like:'
for arg in sys.argv:
  print '\t'+arg
print 'Using tables:'
print '\t Data/MC:',dataMC_filepath
print '\t QCD(data):',qcd_data_filepath
print 'Using systematics from:',systematics_filepath

# get xsections
xsectionDict = ParseXSectionFile(xsection)
dictSamples = GetSamplesToCombineDict(sampleListForMerging)
dictSamplesQCD  = GetSamplesToCombineDict(sampleListForMergingQCD)
dictSamples.update(dictSamplesQCD)

# check to make sure we have xsections for all samples
for lin in open( inputList ):
    lin = string.strip(lin,"\n")
    if lin.startswith('#'):
      continue
    dataset_fromInputList = string.split( string.split(lin, "/" )[-1], ".")[0]
    xsection_val = lookupXSection(SanitizeDatasetNameFromInputList(dataset_fromInputList),xsectionDict)


# rates/etc.
FillDicts(dataMC_filepath,qcd_data_filepath)
# systematics
backgroundSystDict = FillSystDicts(systematicsNamesBackground)
signalSystDict = FillSystDicts(systematicsNamesSignal,False)
# print one of them for checking
#for syst in backgroundSystDict.keys():
#    print 'Syst is:',syst
#    print 'selection\t\tvalue'
#    for selection in sorted(backgroundSystDict[syst].keys()):
#        print selection+'\t\t'+str(backgroundSystDict[syst][selection])
#    break
#print signalSystDict
#print backgroundSystDict

card_file_path = "tmp_card_file.txt"
card_file = open ( card_file_path, "w" ) 

for i_signal_name, signal_name in enumerate(signal_names):
    for i_mass_point, mass_point in enumerate(mass_points):
        fullSignalName = signal_name + mass_point
        selectionName = 'LQ'+mass_point
        
        txt_file_name = fullSignalName + ".txt\n"

        card_file.write ( txt_file_name + "\n\n" )
        card_file.write ( "imax " + str ( n_channels    ) + "\n" ) 
        card_file.write ( "jmax " + str ( n_background  ) + "\n" ) 
        card_file.write ( "kmax " + str ( n_systematics ) + "\n\n" ) 
        
        card_file.write ( "bin 1\n\n" )

        if "BetaHalf" in signal_name: 
            #FIXME for unblinding
            #total_data = enujj_data["data"][i_mass_point]
            #card_file.write ( "observation " + str ( enujj_data["data"][i_mass_point] ) + "\n\n" )
            card_file.write ( "observation " + str ( -1 ) + "\n\n" )
        else : 
            #FIXME for unblinding
            #total_data = eejj_data["data"][i_mass_point]
            #card_file.write ( "observation " + str ( eejj_data["data"][i_mass_point] ) + "\n\n" )
            card_file.write ( "observation " + str ( -1 ) + "\n\n" )
        
        line = "bin " 
        for i_channel in range (0, n_background + 1) :
            line = line + "1 " 
        card_file.write (line + "\n") 

        line = "process " + signal_name + mass_point + " "
        for background_name in background_names:
            line = line + background_name + " "
        card_file.write (line + "\n") 

        line = "process 0 "
        for background_name in background_names:
            line = line + "1 "
        card_file.write (line + "\n\n") 

        # rate line
        line = "rate "
        total_bkg = 0.0
        total_signal = d_signal_rates[fullSignalName][selectionName]
        line = line + str(total_signal) + " "
        for ibkg,background_name in enumerate(background_names):
            thisBkgEvts = d_background_rates[background_name][selectionName]
            line += str(thisBkgEvts) + " "
            total_bkg += float(thisBkgEvts) 
        card_file.write ( line + "\n\n")

        #print signal_name, mass_point, total_signal, total_bkg, total_data
        print signal_name, mass_point, total_signal, total_bkg

        # recall the form: signal --> sysDict['Trigger']['LQXXXX'] = value
        #             backgrounds --> sysDict['Trigger'][bkgName]['LQXXXX'] = value
        for syst in signalSystDict.keys():
            line = syst + ' lnN '
            if selectionName not in signalSystDict[syst].keys():
                selectionNameSigSyst = maxLQselectionBkg
            else:
                selectionNameSigSyst = selectionName
            line += str(1+signalSystDict[syst][selectionNameSigSyst])
            line += ' '
            #else:
            #    print 'ERROR: could not find syst "',syst,'" in signalSystDict.keys():',signalSystDict.keys()
            for ibkg,background_name in enumerate(syst_background_names):
                #print 'try to lookup backgroundSystDict['+syst+']['+background_name+']['+selectionName+']'
                #print 'syst="'+syst+'"'
                if background_name=='' or 'QCD' in background_name:
                    #print 'empty background_name; use - and continue'
                    line += ' - '
                    continue
                if selectionName not in backgroundSystDict[syst][background_name].keys():
                    selectionNameBkgSyst = maxLQselectionBkg
                else:
                    selectionNameBkgSyst = selectionName
                try:
                  line += str(1+backgroundSystDict[syst][background_name][selectionNameBkgSyst])+' '
                except KeyError:
                    print 'Got a KeyError with: backgroundSystDict['+syst+']['+background_name+']['+selectionNameBkgSyst+']'
            card_file.write(line+'\n')

        # background-only special systs: "DYShape", "TTShape"
        for syst in ["DYShape","TTShape"]:
            line = syst + ' lnN - '
            for ibkg,background_name in enumerate(syst_background_names):
                if syst=='DYShape' and not 'DY' in background_name or syst=='TTShape' and not 'TT' in background_name:
                    #print 'empty background_name; use - and continue'
                    line += ' - '
                    continue
                if selectionName not in backgroundSystDict[syst][background_name].keys():
                    selectionNameBkgSyst = maxLQselectionBkg
                else:
                    selectionNameBkgSyst = selectionName
                try:
                  line += str(1+backgroundSystDict[syst][background_name][selectionNameBkgSyst])+' '
                except KeyError:
                    print 'Got a KeyError with: backgroundSystDict['+syst+']['+background_name+']['+selectionNameBkgSyst+']'
            card_file.write(line+'\n')
        
        # background norm systs
        foundTTBar = False
        foundZJet = False
        foundQCD = False
        for ibkg,background_name in enumerate(syst_background_names):
            # XXX WARNING: hardcoded background name (ick); some checking is done at least
            if 'TTbar' in background_name and not foundTTBar:
                line = 'norm_ttbar lnN - '
                line += ' - '*(ibkg)
                line += str(1+ttBarNormDeltaXOverX)+' '
                line += ' - '*(len(syst_background_names)-ibkg-1)+'\n'
                card_file.write(line)
                foundTTBar = True
            elif 'DY' in background_name and not foundZJet:
                line = 'norm_zjet lnN - '
                line += ' - '*(ibkg)
                line += str(1+zJetNormDeltaXOverX)+' '
                line += ' - '*(len(syst_background_names)-ibkg-1)+'\n'
                card_file.write(line)
                foundZJet = True
            elif 'QCD' in background_name and not foundQCD:
                line = 'norm_QCD lnN - '
                line += ' - '*(ibkg)
                line += str(1+qcdNormDeltaXOverX)+' '
                line += ' - '*(len(syst_background_names)-ibkg-1)+'\n'
                card_file.write(line)
                foundQCD = True
        if not foundTTBar or not foundZJet or not foundQCD:
            print 'ERROR: could not find one or more of [ttbar,zjet,QCD] background names for normalization syst; check background names'
            exit(-1)

        card_file.write("\n")

        # background stat error part
        for i_background_name ,background_name in enumerate(background_names):
            thisBkgEvts = d_background_rates[background_name][selectionName]
            thisBkgEvtsErr = d_background_rateErrs[background_name][selectionName]
            thisBkgTotalEntries = d_background_unscaledRates[background_name][selectionName]

            if thisBkgEvts != 0.0: 
                lnN_f = 1.0 + 1.0/math.sqrt(thisBkgTotalEntries+1) # Poisson becomes Gaussian, approx by logN with this kappa
                gmN_weight = thisBkgEvts / thisBkgTotalEntries # for small uncertainties, use gamma distribution with alpha=(factor to go to signal region from control/MC)
            else: 
                # for small uncertainties, use gamma distribution with alpha=(factor to go to signal region from control/MC)
                # since we can't compute evts/entries, we use it from the preselection (following LQ2)
                gmN_weight = d_background_rates[background_name]['preselection'] / d_background_unscaledRates[background_name]['preselection']
            
            line_ln = "stat_" + background_name + " lnN -"
            line_gm = "stat_" + background_name + " gmN " + str(int(thisBkgTotalEntries)) + " -"
            for i_tmp in range ( 0, i_background_name ):
                line_ln = line_ln + " -"
                line_gm = line_gm + " -"
            line_ln = line_ln + " " + str(lnN_f)
            line_gm = line_gm + " " + str(gmN_weight)
            for i_tmp in range ( i_background_name, len(background_names) -1 ):
                line_ln = line_ln + " -"
                line_gm = line_gm + " -"

            if thisBkgTotalEntries > 10:
                card_file.write (line_ln + "\n")
            else:
                card_file.write (line_gm + "\n")
            #if background_name=='TTbar_Madgraph':
            #    print 'selectionName=',selectionName
            #    print 'thisBkgEvts=',thisBkgEvts
            #    print 'thisBkgEvtsErr=',thisBkgEvtsErr
            #    print 'thisBkgTotalEntries=',thisBkgTotalEntries
            #    print 'line_gm=',line_gm
            
        # signal stat error part
        # always use lnN error
        thisSigEvts = d_signal_rates[fullSignalName][selectionName]
        thisSigEvtsErr = d_signal_rateErrs[fullSignalName][selectionName]
        thisSigTotalEntries =d_signal_unscaledRates[fullSignalName][selectionName]
        #if thisSigEvts == 0.0: 
        #  print 'ERROR: signal events for this signal (',fullSignalName,'came out to be zero...stat error not supported. Quitting!'
        #  exit(-1)
        if thisSigEvts != 0.0:
          lnN_f = 1.0 + 1.0/math.sqrt(thisSigTotalEntries+1) # Poisson becomes Gaussian, approx by logN with this kappa
          gmN_weight = thisSigEvts / thisSigTotalEntries
        else:
          gmN_weight = d_signal_rates[background_name]['preselection'] / d_signal_unscaledRates[background_name]['preselection']
        line_ln = "stat_Signal lnN " + str(lnN_f)
        line_gm = "stat_Signal gmN " + str(int(thisBkgTotalEntries)) + " " + str(gmN_weight)
        for i_background_name ,background_name in enumerate(background_names):
            line_ln = line_ln + " -"
            line_gm = line_ln + " -"
        if thisSigTotalEntries > 10:
            card_file.write (line_ln + "\n")
        else:
            card_file.write (line_gm + "\n")

        # DONE!
        card_file.write("\n\n\n")

print 'datacard written to:',card_file_path

# make final selection tables
columnNames = ['MLQ','signal','Z+jets','ttbar','QCD(data)','Other','Data','Total BG']
otherBackgrounds = ['PhotonJets_Madgraph','WJet_Madgraph_HT','DIBOSON','SingleTop']
#background_names =  [ "PhotonJets_Madgraph", "QCDFakes_DATA", "TTbar_Madgraph", "WJet_Madgraph_HT", "ZJet_Madgraph_HT", "DIBOSON","SingleTop"  ]
latexRows = []
t = PrettyTable(columnNames)
t.float_format = "4.3"
selectionNames = ['LQ'+mass_point for mass_point in mass_points]
selectionNames.insert(0,'preselection')
for i_signal_name, signal_name in enumerate(signal_names):
    for selectionName in selectionNames:
        massPoint = selectionName.replace('LQ','')
        fullSignalName = signal_name + massPoint
        # signal events
        thisSigEvts = '-'
        thisSigEvtsErr = '-'
        #print 'selectionName=',selectionName
        if selectionName!='preselection':
            thisSigEvts = d_signal_rates[fullSignalName][selectionName]
            thisSigEvtsErr = d_signal_rateErrs[fullSignalName][selectionName]
        backgroundEvts = {}
        backgroundEvtsErr = {}
        totalBackground = 0.0
        totalBackgroundErrStat = 0.0
        totalBackgroundErrSyst = 0.0 #FIXME TODO
        otherBackground = 0.0
        otherBackgroundErrStat = 0.0
        for i_background_name ,background_name in enumerate(background_names):
            thisBkgEvts = d_background_rates[background_name][selectionName]
            thisBkgEvtsErr = d_background_rateErrs[background_name][selectionName]
            thisBkgTotalEntries = d_background_unscaledRates[background_name][selectionName]
            totalBackground+=thisBkgEvts
            totalBackgroundErrStat+=(thisBkgEvtsErr*thisBkgEvtsErr)
            if background_name in otherBackgrounds:
              otherBackground+=thisBkgEvts
              otherBackgroundErrStat+=(thisBkgEvtsErr*thisBkgEvtsErr)
            backgroundEvts[background_name] = thisBkgEvts
            backgroundEvtsErr[background_name] = thisBkgEvtsErr
        totalBackgroundErrStat = math.sqrt(totalBackgroundErrStat)
        otherBackgroundErrStat = math.sqrt(otherBackgroundErrStat)
        row = [selectionName,GetTableEntryStr(thisSigEvts,thisSigEvtsErr),
            GetTableEntryStr(backgroundEvts['ZJet_Madgraph_HT'],backgroundEvtsErr['ZJet_Madgraph_HT']),
            GetTableEntryStr(backgroundEvts['TTbar_Madgraph'],backgroundEvtsErr['TTbar_Madgraph']),
            GetTableEntryStr(backgroundEvts['QCDFakes_DATA'],backgroundEvtsErr['QCDFakes_DATA']),
            GetTableEntryStr(otherBackground,otherBackgroundErrStat),
            GetTableEntryStr(-1,-1),
            GetTableEntryStr(totalBackground,totalBackgroundErrStat,totalBackgroundErrSyst),
            ]
        latexRow = ''
        for i,entry in enumerate(row):
            entry = entry.replace('+/-','$\pm$')
            entry = entry.replace('LQ','')
            latexRow+=entry
            if i<len(row)-1:
                latexRow+=' & '
            else:
                latexRow+=' \\\\ '
        latexRows.append(latexRow)
        if selectionName=='preselection':
          latexRows.append('\\hline')
        t.add_row(row)
print t

print
# latex table
for line in latexRows:
  print line
print

exit(0)

