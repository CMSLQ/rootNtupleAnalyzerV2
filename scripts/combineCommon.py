#!/usr/bin/env python

#---Import
import sys
import string
import math


def SanitizeDatasetNameFromInputList(dataset_fromInputList):
  # "hack" for data-driven QCD samples: name is created by the createInputList script
  # do this first, since it's at the very end of the filename
  # XXX FIXME special hacks for datasets
  #if dataset_fromInputList.contains('_reduced_skim'):
  #  #dataset_fromInputList = dataset_fromInputList[0:dataset_fromInputList.find('_reduced_skim')]
  #  dataset_fromInputList.replace('_reduced_skim','')
  dataset_fromInputList = dataset_fromInputList.replace('_reduced_skim','')
  #XXX FIXME
  ## special hack for handling repated madgraphMLM samples
  #if dataset_fromInputList.endswith('_madgraphMLM'):
  #  dataset_fromInputList = dataset_fromInputList[0:dataset_fromInputList.find('_madgraphMLM')]
  ##XXX FIXME
  ## special hack for handling repated amcatnloFXFX samples
  #elif dataset_fromInputList.endswith('_amcatnloFXFX'):
  #  dataset_fromInputList = dataset_fromInputList[0:dataset_fromInputList.find('_amcatnloFXFX')]
  if dataset_fromInputList.endswith('_pythia8'):
    dataset_fromInputList = dataset_fromInputList[0:dataset_fromInputList.find('_pythia8')]
  #if '__' in dataset_fromInputList:
  #  dataset_fromInputList = dataset_fromInputList[0:dataset_fromInputList.find('__')]
  return dataset_fromInputList


def SanitizeDatasetNameFromFullDataset(dataset):
  #print 'SanitizeDatasetNameFromFullDataset: dataset looks like:'+dataset
  # this logic is somewhat copied from the submission script for the ntuples:
  #    https://github.com/CMSLQ/submitJobsWithCrabV2/blob/master/createAndSubmitJobsWithCrab3.py
  if not 'Run2015' in dataset:
    outputFileNames = []
    outputFileNames.append(dataset[1:dataset.find('_Tune')])
    outputFileNames.append(dataset[1:dataset.find('_13TeV')])
    try:
      outputFileNames.append(dataset.split('/')[1])
    except IndexError:
      print "Had an IndexError trying to split('/') dataset:",dataset,'; this can happen if this is a piece (not a full dataset) containing multiple samples that has not been defined earlier in the sampleListToCombineFile'
      exit(-1)

    # use the one with the shortest filename
    outputFile = sorted(outputFileNames, key=len)[0]
    if 'ext' in dataset:
       extN = dataset[dataset.find('_ext')+4]
       outputFile = outputFile+'_ext'+extN
    if 'madgraphMLM' in dataset:
      outputFile+='_madgraphMLM'
    elif 'amcatnloFXFX' in dataset:
      outputFile+='_amcatnloFXFX'
  else:
    outputFile = dataset[1:].replace('/','__')
    outputFile = outputFile.split('__')[0]+'__'+outputFile.split('__')[1]
  return outputFile


def GetSamplesToCombineDict(sampleListForMerging):
  dictSamples = {}
  for l,line in enumerate( open( sampleListForMerging ) ):
    # ignore comments
    if line.startswith('#'):
      continue
    line = string.strip(line,"\n")
    # ignore empty lines
    if len(line) <= 0:
      continue

    #print 'line from samplesToCombineFile looks like:"'+line+'"'
    # line looks like: "ZJet_Madgraph_Inc    DYJetsToLL_M-5to50 DYJetsToLL_M-50"
    # or "LQ_M200   /LQToUE_M-200_BetaOne_TuneCUETP8M1_13TeV-pythia8/RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2/MINIAODSIM"

    # the rule is that the name of each 'piece' here must match the inputList name and filename 
    for i,piece in enumerate(line.split()):
      #print "i=", i , "  piece= " , piece
      if (i == 0):
        key = piece
        dictSamples[key] = []
      elif piece in dictSamples:
        dictSamples[key].extend(dictSamples[piece])
      else:
        #print 'GetSamplesToCombineDict: SanitizeDatasetNameFromFullDataset(',piece,')'
        piece = SanitizeDatasetNameFromFullDataset(piece)
        dictSamples[key].append(piece)
  return dictSamples


def ParseXSectionFile(xsectionFile):
  xsectionDict = {}

  for line in open( xsectionFile ):

    # ignore comments
    if line.startswith('#'):
      continue
    line = string.strip(line,"\n")
    # ignore empty lines
    if len(line) <= 0:
      continue

    try:
      dataset, xsection_val = string.split(line)
    except ValueError:
      print 'ERROR: could not split line "',line,'"'
      exit(-1)

    outputFile = SanitizeDatasetNameFromFullDataset(dataset)
    xsectionDict[outputFile] = xsection_val
    #print outputFile + " " + xsection_val
  
  return xsectionDict


def lookupXSection(datasetNameFromInputList,xsectionDict):
  if len(xsectionDict) <= 0:
    print
    print 'ERROR: xsectionDict is empty. Cannot lookupXSection for',datasetNameFromInputList
    exit(-1)
  for dataset in xsectionDict.keys():
    if dataset.startswith(datasetNameFromInputList):
      return xsectionDict[dataset]
  print 'ERROR: lookupXSection(): xsectionDict does not have an entry for',datasetNameFromInputList
  exit(-1)
  #try:
  #  #return xsectionDict[datasetNameFromInputList]
  #except KeyError:
  #  print
  #  #for key,val in xsectionDict.iteritems():
  #  #  print 'sample=',key,'xsection=',val
  #  print 'ERROR: lookupXSection(): xsectionDict does not have an entry for',datasetNameFromInputList
  #  exit(-1)


#def lookupXSection(datasetFromAnalysis,xsectionFile):
#  dataset_forXsec = datasetFromAnalysis
#  if datasetFromAnalysis.endswith('_reduced_skim'):
#    dataset_forXsec = datasetFromAnalysis[0:datasetFromAnalysis.find('_reduced_skim')]
#
#  xsectionDataset = ''
#  xsectionVal = -999
#  for lin1 in open( xsectionFile ):
#
#    lin1 = string.strip(lin1,"\n")
#
#    (dataset , xsection_val) = string.split(lin1)
#    #print dataset + " " + xsection_val
#
#    dataset_mod_1 = dataset[1:].replace('/','__')
#    #print dataset_mod_1 + " " + xsection_val
#
#    # TODO: fix this hack!
#    if ( dataset_mod_1.startswith(dataset_forXsec+'_Tune') or
#         dataset_mod_1.startswith(dataset_forXsec+'_13TeV') or
#         #('Run20' in dataset_mod_1 and dataset_mod_1.startswith(dataset_forXsec)) or
#         ('Run20' in dataset_mod_1 and dataset[1:].replace('/','_').startswith(dataset_forXsec)) or
#         (dataset_forXsec=='DYJetsToLL_Mbin_M-50' and dataset_mod_1.startswith('DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8')) or
#         (dataset_forXsec=='TTJets_SingleLeptFromTbar_ext1' and dataset_mod_1.startswith('TTJets_SingleLeptFromTbar_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9-v2')) or
#         (dataset_forXsec=='TTJets_DiLept_ext1' and dataset_mod_1.startswith('TTJets_DiLept_TuneCUETP8M1_13TeV-madgraphMLM-pythia8__RunIISpring15DR74-Asympt25ns_MCRUN2_74_V9_ext1-v1'))):
#      # TODO: fix this hack
#      # special handling of DYJetsToLL_M-50
#      # for madgraph; Mbin has _Mbin in it
#      if dataset_forXsec=='DYJetsToLL_M-50':
#        if not 'madgraph' in dataset:
#          continue
#      #elif dataset_forXsec=='DYJetsToLL_Mbin_M-50':
#      #  if not 'amcatnloFXFX' in dataset:
#      #    continue
#      #
#      if len(xsectionDataset) <= 0:
#        xsectionDataset = dataset
#        xsectionVal = xsection_val
#      elif xsectionVal != xsection_val:
#        print 'ERROR: Two datasets in xsection file start with',dataset_forXsec
#        print '1)',xsectionDataset,xsectionVal
#        print '2)',dataset,xsection_val
#        print 'Cannot figure out which is correct; exiting'
#        sys.exit()
#
#  xsection_val = xsectionVal
#  dataset = xsectionDataset
#  if xsection_val < -1: # -1 being the value for data
#    print "ERROR: xsection for dataset " + dataset_forXsec + " not found in " + xsectionFile
#    print "Expected a line in xsection file to start with "+dataset_forXsec+"_Tune but couldn't find one."
#    print "exiting..."
#    sys.exit()
#  return dataset,xsection_val
#
#

#def AddHisto(inputHistoName, outputHisto, inputRootFileName, currentWeight,
#             rebin=int(1), currentColor=int(1), currentFillStyle=int(1001), currentMarker=int(1)):

def UpdateTable(inputTable, outputTable):
    if not outputTable:
        for j,line in enumerate( inputTable ):
            outputTable[int(j)]={'variableName': inputTable[j]['variableName'],
                                 'min1': inputTable[j]['min1'],
                                 'max1': inputTable[j]['max1'],
                                 'min2': inputTable[j]['min2'],
                                 'max2': inputTable[j]['max2'],
                                 'N':       float(inputTable[j]['N']),
                                 'errN':    pow( float(inputTable[j]['errN']), 2 ),
                                 'Npass':       float(inputTable[j]['Npass']),
                                 'errNpass':    pow( float(inputTable[j]['errNpass']), 2 ),
                                 'EffRel':      float(0),
                                 'errEffRel':   float(0),
                                 'EffAbs':      float(0),
                                 'errEffAbs':   float(0),
                                 }
    else:
        for j,line in enumerate( inputTable ):
            #print 'outputTable[int(',j,')][N]=',outputTable[int(j)]['N'],'inputTable[',j,']','[N]=',inputTable[j]['N']
            outputTable[int(j)]={'variableName': inputTable[j]['variableName'],
                                 'min1': inputTable[j]['min1'],
                                 'max1': inputTable[j]['max1'],
                                 'min2': inputTable[j]['min2'],
                                 'max2': inputTable[j]['max2'],
                                 'N':       float(outputTable[int(j)]['N']) + float(inputTable[j]['N']),
                                 'errN':    float(outputTable[int(j)]['errN']) + pow( float(inputTable[j]['errN']), 2 ),
                                 'Npass':       float(outputTable[int(j)]['Npass']) + float(inputTable[j]['Npass']),
                                 'errNpass':    float(outputTable[int(j)]['errNpass']) + pow( float(inputTable[j]['errNpass']), 2 ),
                                 'EffRel':      float(0),
                                 'errEffRel':   float(0),
                                 'EffAbs':      float(0),
                                 'errEffAbs':   float(0),
                                 }
    return
            

def SubtractTables(inputTable, outputTable, zeroNegatives = False):
    # subtract the inputTable from the outputTable
    if not outputTable:
        print 'ERROR: no outputTable found! cannot subtract input from nothing; FATAL'
        exit(-1)
    else:
        for j,line in enumerate( inputTable ):
            #print 'outputTable[int(',j,')][N]=',outputTable[int(j)]['N'],'inputTable[',j,']','[N]=',inputTable[j]['N']
            newN = float(outputTable[int(j)]['N']) - float(inputTable[j]['N'])
            newNpass = float(outputTable[int(j)]['Npass']) - float(inputTable[j]['Npass'])
            if newN < 0.0 and zeroNegatives:
                newN = 0.0
            if newNpass < 0.0 and zeroNegatives:
                newNpass = 0.0
            outputTable[int(j)]={'variableName': inputTable[j]['variableName'],
                                 'min1': inputTable[j]['min1'],
                                 'max1': inputTable[j]['max1'],
                                 'min2': inputTable[j]['min2'],
                                 'max2': inputTable[j]['max2'],
                                 'N':       newN,
                                 'errN':    math.sqrt(pow(float(outputTable[int(j)]['errN']),2) + pow( float(inputTable[j]['errN']), 2 )),
                                 'Npass':       newNpass,
                                 'errNpass':    math.sqrt(pow(float(outputTable[int(j)]['errNpass']),2) + pow( float(inputTable[j]['errNpass']), 2 )),
                                 'EffRel':      float(0),
                                 'errEffRel':   float(0),
                                 'EffAbs':      float(0),
                                 'errEffAbs':   float(0),
                                 }
    return
            

def ScaleTable(inputTable, scaleFactor, errScaleFactor):
    if not inputTable:
        print 'ERROR: no inputTable found! cannot scale nothing; FATAL'
        exit(-1)
    else:
        for j,line in enumerate( inputTable ):
            nOrig = float(inputTable[int(j)]['N'])
            errNorig = float(inputTable[int(j)]['errN'])
            nNew =  nOrig * scaleFactor
            if nOrig > 0.0:
              errNnew = nNew*math.sqrt(pow(errNorig/nOrig,2)+pow(errScaleFactor/scaleFactor,2))
            else:
              errNnew = nNew*math.sqrt(pow(errNorig/nOrig,2)+pow(errScaleFactor/scaleFactor,2))
            nPassOrig = float(inputTable[int(j)]['Npass'])
            errNPassOrig = float(inputTable[j]['errNpass'])
            nPassNew =  nPassOrig * scaleFactor
            errNpassNew = nPassNew*math.sqrt(pow(errNPassOrig/nPassOrig,2)+pow(errScaleFactor/scaleFactor,2))
            
            inputTable[int(j)]={'variableName': inputTable[j]['variableName'],
                                 'min1': inputTable[j]['min1'],
                                 'max1': inputTable[j]['max1'],
                                 'min2': inputTable[j]['min2'],
                                 'max2': inputTable[j]['max2'],
                                 'N':           nNew,
                                 'errN':        errNnew,
                                 'Npass':       nPassNew,
                                 'errNpass':    errNpassNew,
                                 'EffRel':      float(0),
                                 'errEffRel':   float(0),
                                 'EffAbs':      float(0),
                                 'errEffAbs':   float(0),
                                 }
    return
            

def SquareTableErrorsForEfficiencyCalc(table):
    if not table:
        print 'ERROR: no inputTable found! cannot convert nothing; FATAL'
        exit(-1)
    else:
        for j,line in enumerate( table ):
            table[int(j)]={'variableName': table[j]['variableName'],
                                 'min1': table[j]['min1'],
                                 'max1': table[j]['max1'],
                                 'min2': table[j]['min2'],
                                 'max2': table[j]['max2'],
                                 'N':          float(table[j]['N']) ,
                                 'errN':       pow(float(table[j]['errN']),2) , 
                                 'Npass':      float(table[j]['Npass']) ,
                                 'errNpass':   pow(float(table[j]['errNpass']),2) , 
                                 'EffRel':      float(0),
                                 'errEffRel':   float(0),
                                 'EffAbs':      float(0),
                                 'errEffAbs':   float(0),
                                 }
    return
            

def CalculateEfficiency(table):
    # this also (sneakily) converts 'errors' in the tables (which are really errSqr) to sqrt(errors)
    for j,line in enumerate( table ):
        if( j == 0):
            table[int(j)] = {'variableName':       table[int(j)]['variableName'],
                             'min1':        table[int(j)]['min1'],
                             'max1':        table[int(j)]['max1'],
                             'min2':        table[int(j)]['min2'],
                             'max2':        table[int(j)]['max2'],
                             'N':          float(table[j]['N']) ,
                             'errN':       int(0), 
                             'Npass':      float(table[j]['Npass']) ,
                             'errNpass':   int(0), 
                             'EffRel':     int(1),
                             'errEffRel':  int(0),
                             'EffAbs':     int(1),
                             'errEffAbs':  int(0),
                             }
        else:
            N = float(table[j]['N']) 
            errN = math.sqrt(float(table[j]["errN"]))
            if( float(N) > 0 ):
                errRelN = errN / N 
            else:
                errRelN = float(0)

            Npass = float(table[j]['Npass']) 
            errNpass = math.sqrt(float(table[j]["errNpass"]))
            if( float(Npass) > 0 ):
                errRelNpass = errNpass / Npass
            else:
                errRelNpass = float(0)

            if(Npass > 0  and N >0 ):
                EffRel = Npass / N
                errRelEffRel = math.sqrt( errRelNpass*errRelNpass + errRelN*errRelN )
                errEffRel = errRelEffRel * EffRel
            
                EffAbs = Npass / float(table[0]['N'])
                errEffAbs = errNpass / float(table[0]['N'])
            else:
                EffRel = 0
                errEffRel = 0
                EffAbs = 0
                errEffAbs = 0 
            
            table[int(j)]={'variableName': table[int(j)]['variableName'],
                           'min1': table[int(j)]['min1'],
                           'max1': table[int(j)]['max1'],
                           'min2': table[int(j)]['min2'],
                           'max2': table[int(j)]['max2'],
                           'N':       N,
                           'errN':    errN, 
                           'Npass':       Npass,
                           'errNpass':    errNpass, 
                           'EffRel':      EffRel,
                           'errEffRel':   errEffRel,
                           'EffAbs':      EffAbs,
                           'errEffAbs':   errEffAbs,
                           }
            #print table[j]
    return


#--- TODO: FIX TABLE FORMAT (NUMBER OF DECIMAL PLATES AFTER THE 0)

def WriteTable(table, name, file):
    print >>file, name
    print >>file, "variableName".rjust(25),
    print >>file, "min1".rjust(15),
    print >>file, "max1".rjust(15),
    print >>file, "min2".rjust(15),
    print >>file, "max2".rjust(15),
    print >>file, "Npass".rjust(17),
    print >>file, "errNpass".rjust(17),
    print >>file, "EffRel".rjust(15),
    print >>file, "errEffRel".rjust(15),
    print >>file, "EffAbs".rjust(15),
    print >>file, "errEffAbs".rjust(15)

    for j, line in enumerate(table):
        print >>file, table[j]['variableName'].rjust(25),
        print >>file, table[j]['min1'].rjust(15),
        print >>file, table[j]['max1'].rjust(15),
        print >>file, table[j]['min2'].rjust(15),
        print >>file, table[j]['max2'].rjust(15),
        ###
        if( table[j]['Npass'] >= 0.1 ):
            print >>file, ("%.04f" % table[j]['Npass']).rjust(17),
        else:
            print >>file, ("%.04e" % table[j]['Npass']).rjust(17),
        ### 
        if( table[j]['errNpass'] >= 0.1):    
            print >>file, ("%.04f" % table[j]['errNpass']).rjust(17),
        else:
            print >>file, ("%.04e" % table[j]['errNpass']).rjust(17),
        ### 
        if( table[j]['EffRel'] >= 0.1 ):
            print >>file, ("%.04f" % table[j]['EffRel']).rjust(15),
        else:
            print >>file, ("%.04e" % table[j]['EffRel']).rjust(15),
        ### 
        if( table[j]['errEffRel'] >= 0.1 ):
            print >>file, ("%.04f" % table[j]['errEffRel']).rjust(15),    
        else:
            print >>file, ("%.04e" % table[j]['errEffRel']).rjust(15),
        ### 
        if( table[j]['EffAbs'] >= 0.1 ):
            print >>file, ("%.04f" % table[j]['EffAbs']).rjust(15),
        else:
            print >>file, ("%.04e" % table[j]['EffAbs']).rjust(15),
        ### 
        if( table[j]['errEffAbs'] >= 0.1 ):
            print >>file, ("%.04f" % table[j]['errEffAbs']).rjust(15)
        else:
            print >>file, ("%.04e" % table[j]['errEffAbs']).rjust(15)         
        ###
            
    print >>file, "\n"

    #--- print to screen
    
    print "\n"
    print name
    print "variableName".rjust(25),
    print "min1".rjust(15),
    print "max1".rjust(15),
    print "min2".rjust(15),
    print "max2".rjust(15),
    print "Npass".rjust(17),
    print "errNpass".rjust(17),
    print "EffRel".rjust(15),
    print "errEffRel".rjust(15),
    print "EffAbs".rjust(15),
    print "errEffAbs".rjust(15)

    for j, line in enumerate(table):
        print table[j]['variableName'].rjust(25),
        print table[j]['min1'].rjust(15),
        print table[j]['max1'].rjust(15),
        print table[j]['min2'].rjust(15),
        print table[j]['max2'].rjust(15),
        ###
        if( table[j]['Npass'] >= 0.1 ):
            print ("%.04f" % table[j]['Npass']).rjust(17),
        else:
            print ("%.04e" % table[j]['Npass']).rjust(17),
        ### 
        if( table[j]['errNpass'] >= 0.1):    
            print ("%.04f" % table[j]['errNpass']).rjust(17),
        else:
            print ("%.04e" % table[j]['errNpass']).rjust(17),
        ### 
        if( table[j]['EffRel'] >= 0.1 ):
            print ("%.04f" % table[j]['EffRel']).rjust(15),
        else:
            print ("%.04e" % table[j]['EffRel']).rjust(15),
        ### 
        if( table[j]['errEffRel'] >= 0.1 ):
            print ("%.04f" % table[j]['errEffRel']).rjust(15),    
        else:
            print ("%.04e" % table[j]['errEffRel']).rjust(15),
        ### 
        if( table[j]['EffAbs'] >= 0.1 ):
            print ("%.04f" % table[j]['EffAbs']).rjust(15),
        else:
            print ("%.04e" % table[j]['EffAbs']).rjust(15),
        ### 
        if( table[j]['errEffAbs'] >= 0.1 ):
            print ("%.04f" % table[j]['errEffAbs']).rjust(15)
        else:
            print ("%.04e" % table[j]['errEffAbs']).rjust(15)         
        ###

    return


