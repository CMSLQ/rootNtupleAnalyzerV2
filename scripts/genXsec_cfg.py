# Instructions:
# 1) Need CMSSW_12 or higher release area
# 2) cmsenv inside release area
# 3) if needed (first time only), install DBS client with 'pip3 install dbs3-client'
# 4) run this script: 'cmsRun genXsec_cfg.py dataset=/Full/NanoAOD/dataset'
import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
try:
    from dbs.apis.dbsClient import DbsApi
    from dbs.exceptions.dbsClientException import dbsClientException
except ImportError:
    print()
    print("ERROR: Could not load dbs APIs")
    print("You need to run this script using cmsRun after doing 'cmsenv' inside a CMSSW release area with major version >= 12, and then (if necessary-first time only) installing the DBS client with 'pip3 install dbs3-client'.")
    exit(-1)

options = VarParsing ('analysis')
options.register ('dataset',
                  'default',
                  VarParsing.multiplicity.singleton,
                  VarParsing.varType.string,
                  "Full name/path of dataset")

options.parseArguments()

datasetpath = options.dataset
if datasetpath=='default':
    print('ERROR: must specify full dataset name of NANOAODSIM dataset: dataset=/Primary/Secondary/NANOAODSIM or dataset=test to skip DBS querying')
    print( 'Exiting.')
    print()
    exit(-1)

doQuery = True
if datasetpath=='test':
    doQuery = False

if doQuery:
    if 'MINIAOD' not in datasetpath:
        print('DBS: query for parent of dataset='+datasetpath+'...', end='')
        
        try:
            api = DbsApi(url = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
            parents = api.listDatasetParents(dataset=datasetpath)
        except dbsClientException as ex:
            print()
            print("Caught API Exception %s: %s "  % (ex.name,ex))
            exit(-1)
        
        #for parent in parents:
        #    print parent
        if len(parents) > 1:
            print()
            print('ERROR: got multiple parents for dataset:')
            for parent in parents:
                print(parent['parent_dataset'])
            print('ERROR: not sure how to continue.')
            exit(-2)
        print('Done')
        parentDataset = parents[0]['parent_dataset']
        print('DBS: Found parent dataset='+parentDataset)
    else:
        parentDataset = datasetpath
    print('DBS: query for files in dataset...', end ='')
    
    try:
        api = DbsApi(url = 'https://cmsweb.cern.ch/dbs/prod/global/DBSReader')
        #dbsFiles = api.listFiles(dataset=datasetpath)
        dbsFiles = api.listFiles(dataset=parentDataset)
    except dbsClientException as ex:
        print()
        print("Caught API Exception %s: %s "  % (ex.name,ex))
        exit(-3)
    
    print('Done')
    #for dbsFile in dbsFiles:
    #    print dbsFile['logical_file_name']
    # get first N LFNs
    filesToUse = 100
    #filesToUse = 25
    firstNFiles = [dbsFile['logical_file_name'] for dbsFile in dbsFiles]
    # if needed
    filesToSkip = ["6B651A93-8182-564A-A5C5-80D0D08AB7F1"]
    # filesToSkip = []
    for skipFile in filesToSkip:
        firstNFiles[:] = [x for x in firstNFiles if skipFile not in x]
    if len(firstNFiles) > filesToUse:
        firstNFiles = firstNFiles[-1*filesToUse:]
    # print firstNFiles
    print('INFO: using last',len(firstNFiles),'out of',len(dbsFiles),'files in parent dataset as input to GenXSecAnalyzer')
else:
    print("INFO: not doing DBS query.")

print('Run.')
print()

process = cms.Process('XSec')

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(options.maxEvents)
        )

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1) )
secFiles = cms.untracked.vstring() 
if doQuery:
    process.source = cms.Source ("PoolSource",
            fileNames = cms.untracked.vstring(firstNFiles),
            secondaryFileNames = secFiles)
else:
    process.source = cms.Source ("PoolSource",
            fileNames = cms.untracked.vstring(
                # can put filenames here manually
                '/store/mc/RunIISummer16MiniAODv2/GJets_Pt-20To100_13TeV-sherpa/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/70000/F87154C7-C4CC-E611-BB3D-5065F3819221.root',
              ),
    )

process.xsec = cms.EDAnalyzer("GenXSecAnalyzer")

process.ana = cms.Path(process.xsec)
process.schedule = cms.Schedule(process.ana)
