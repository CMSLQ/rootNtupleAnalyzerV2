#!/usr/bin/env python3

# ---Import
import sys
import string
from optparse import OptionParser
import os
from ROOT import TFile, gROOT, SetOwnership, TObject
import re
import traceback
import subprocess
import shlex
import shutil
import copy
import time
from graphlib import TopologicalSorter
from collections import OrderedDict
import pprint
import tempfile
from termcolor import colored
from pathlib import Path

import combineCommon
from combinePlotsBatch import SeparateArgs, FillDictFromOptionByYear


gROOT.SetBatch(True)

ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')


def RunCommand(cmd, workingDir=None, env=None, suppressOutput=False):
    origCmd = cmd
    if not suppressOutput:
        print(colored("\t{}".format(origCmd), "green"), flush=True)
    try:
        useShell = False
        # if env is None:
        #     useShell = True
        #     cmd = shlex.split(cmd)
        cmd = shlex.split(cmd)
        process = subprocess.run(cmd, shell=useShell, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env=env)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        if suppressOutput:
            print(colored("\t{}".format(origCmd), "green"), flush=True)
        print(colored("RunCommand had an error; output: {}".format(result), "red"), flush=True)
        raise e
    if not suppressOutput:
        print(process.stdout.decode(), flush=True)


def CheckForFile(filename):
    filename = filename.replace("root://eoscms/", "/eos/cms/").replace("root://eosuser/", "/eos/user/")
    if not filename.startswith("/eos"):
        if os.path.isfile(filename):
            return True
    else:
        if not filename.startswith("/eos/user"):
            serverName = filename.split('/store')[0].strip('/').replace('/', '')
        else:
            serverName = "eosuser"
        my_env = os.environ.copy()
        my_env["EOS_MGM_URL"] = "root://{}.cern.ch/".format(serverName)
        result = subprocess.run(["eos", "ls", filename], env=my_env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if result == 0:
            return True
    return False


def RemoveFile(filename):
    filename = filename.replace("root://eoscms/", "/eos/cms/").replace("root://eosuser/", "/eos/user/")
    if not os.path.isfile(filename):
        return
    if not filename.startswith("/eos"):
        os.remove(filename)
    else:
        if not filename.startswith("/eos/user"):
            serverName = filename.split('/store')[0].strip('/').replace('/', '')
        else:
            serverName = "eosuser"
        my_env = os.environ.copy()
        my_env["EOS_MGM_URL"] = "root://{}.cern.ch/".format(serverName)
        result = subprocess.check_call(["eos", "rm", filename], env=my_env, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def CreateAndSubmitJobs(sampleDict):
    condorSubfileDir = "combinePlotsCondor"
    Path(condorSubfileDir).mkdir(parents=True, exist_ok=True)
    shFilename = condorSubfileDir + "/condor.sh"
    samplesSet = set()
    for year, sampleList in sampleDict.items():
        samplesToCombine = [sample for sample, keys in sampleList.items() if keys["save"]]
        samplesSet.update(samplesToCombine)
    WriteCondorShFile(samplesSet, shFilename)
    print("INFO: wrote sh file to {}".format(shFilename))
    subFilename = shFilename.replace(".sh", ".sub")
    WriteCondorSubFile(len(samplesSet), subFilename)
    print("INFO: wrote sub file to {}".format(subFilename))
    SubmitCondorJob(str(Path(subFilename).name), condorSubfileDir)


def WriteCondorShFile(sampleList, filename):
    execName = "combinePlotsBatch.py"
    with open(filename, "w") as outputfile:
        outputfile.write("#!/bin/bash\n")
        outputfile.write("echo \"Running at site: $GLIDEIN_CMSSite\"\n")
        # hardcoded root is a bit nasty FIXME
        outputfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh\n')
        outputfile.write("echo \"Sourced LCG env.\"\n")
        # ROOT likes HOME set
        outputfile.write('[ -z "$HOME" ] && export HOME='+os.getenv('HOME')+'\n')
        outputfile.write("echo \"Set HOME=$HOME.\"\n")
        outputfile.write('export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH\n')
        outputfile.write("echo \"exported LD_LIBRARY_PATH.\"\n")
        outputfile.write("echo \"Running executable.\"\n")
        execArgs = CreateExecArgs()
        for idx, sample in enumerate(sampleList):
            outputfile.write("if [ $1 -eq {} ]; then\n".format(idx))
            outputfile.write('    ./' + execName + ' --sample {} '.format(sample) + ' '.join(execArgs)+"\n")
            outputfile.write('    retVal=$?\n')
            outputfile.write('    if [ $retVal -ne 0 ]; then\n')
            outputfile.write('      echo "./'+execName+' return error code=$retVal; quitting here."\n')
            outputfile.write('      exit $retVal\n')
            outputfile.write('    fi\n')
            outputfile.write('fi\n')


def GetNameRelativeToSecondParent(path):
    return str(Path(path).parent.parent.name) + "/" + GetNameRelativeToParent(path)


def GetNameRelativeToParent(path):
    return str(Path(path).parent.name) + "/" + str(Path(path).name)


def CreateExecArgs():
   args = ['-i ' + ",".join([GetNameRelativeToParent(inputList) for idx, inputList in enumerate(SeparateArgs(options.inputList))])]
   args += ['-c ' + options.analysisCode]
   args += ['-d ' + ",".join([str(Path(sampleList).resolve()) for idx, sampleList in enumerate(SeparateArgs(options.inputDir))])]
   args += ['-l ' + options.intLumi]
   args += ['-y ' + options.years]
   args += ['-x ' + ",".join([GetNameRelativeToSecondParent(xsection) for idx, xsection in enumerate(SeparateArgs(options.xsection))])]
   args += ['-o ' + options.outputDir]
   args += ['-s ' + ",".join([str(Path(sampleList).name) for idx, sampleList in enumerate(SeparateArgs(options.sampleListForMerging))])]
   if options.fitDiagFilepath is not None or options.postFitJSON is not None:
       if options.fitDiagFilepath:
           args += ['--fitDiagFilepath ' + str(Path(options.fitDiagFilepath).name)]
       elif options.postFitJSON:
           args += ['--postFitJSON ' + str(Path(options.postFitJSON).name)]
       if options.preFit:
           args += ['--preFit']
       elif options.postFit:
           args += ['--postFit']
           args += ['--fitType ' + options.fitType]
   return args


def WriteCondorSubFile(nJobs, condorFilename):
    shFilepath = condorFilename.replace(".sub", ".sh")
    outputmain = str(Path(shFilepath).parent.resolve())
    shFilename = str(Path(shFilepath).name)
    Path(outputmain + "/output").mkdir(exist_ok=True)
    Path(outputmain + "/error").mkdir(exist_ok=True)
    Path(outputmain + "/log").mkdir(exist_ok=True)
    with open(condorFilename, 'w') as condorFile:
        condorFile.write('executable  = '+ shFilename + '\n')
        condorFile.write("arguments = $(ProcId)\n")
        condorFile.write('output      = '+outputmain+'/output/$(Process).out\n')
        condorFile.write('error       = '+outputmain+'/error/$(Process).err\n')
        condorFile.write('log         = '+outputmain+'/log/$(Process).log\n')
        # http://batchdocs.web.cern.ch/batchdocs/local/submit.html
        #  - cms connect shouldn't use JobFlavor or the requirements
        #  - assume this is lxbatch if queue option specified
        if options.queue is None:
            condorFile.write('+REQUIRED_OS = rhel9\n')
            condorFile.write('use_x509userproxy = true\n')
            # outputRootFile = outputPrefix+"_$(Process).root"
            # outputDatFile = outputPrefix+"_$(Process).dat"
            # outputPath = outputmain+"/output/"
            # condorFile.write('transfer_output_files = '+outputRootFile+','+outputDatFile+'\n')
            # condorFile.write('transfer_output_remaps = "'+outputRootFile+' = '+outputPath+outputRootFile+'; '+outputDatFile+' = '+outputPath+outputDatFile+'"\n')
            condorFile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.\n")
            condorFile.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
            condorFile.write('transfer_output_files = ""\n')
        else:
            condorFile.write('+JobFlavour = "'+options.queue+'"\n')
            # require EL9
            condorFile.write('MY.WantOS = "el9"\n')
            condorFile.write('transfer_output_files = ""\n')
            condorFile.write('RequestCpus = 4\n')
        # make sure the job finishes with exit code 0
        condorFile.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        condorFile.write('max_retries = 3\n')
        condorFile.write('should_transfer_files = YES\n')
        lqAnaBase = os.getenv("LQANA")
        inputFilesToTransfer = [lqAnaBase + "/scripts/combinePlotsBatch.py", lqAnaBase + "/scripts/combineCommon.py"]
        inputFilesToTransfer.extend([str(Path(xsection).resolve().parent.parent) for idx, xsection in enumerate(SeparateArgs(options.xsection))])
        inputFilesToTransfer.extend([str(Path(inputList).resolve().parent) for idx, inputList in enumerate(SeparateArgs(options.inputList))])
        inputFilesToTransfer.extend([str(Path(sampleList).resolve()) for idx, sampleList in enumerate(SeparateArgs(options.sampleListForMerging))])
        if options.fitDiagFilepath is not None:
            inputFilesToTransfer.extend([str(Path(options.fitDiagFilepath).resolve())])
        elif options.postFitJSON is not None:
            inputFilesToTransfer.extend([str(Path(options.postFitJSON).resolve())])
        filesToTransfer = ",".join(inputFilesToTransfer)
        # filesToTransfer += ","
        # filesToTransfer += options.filesToTransfer
        condorFile.write('transfer_input_files = '+filesToTransfer+'\n')
        condorFile.write("queue {}\n".format(nJobs))


def SubmitCondorJob(subFilename, condorSubfileDir):
    if not options.dryRun:
        cmd = "condor_submit {}".format(subFilename)
        runCommandArgs = [cmd, condorSubfileDir, None, False]
        trial = 0
        maxTries = 5
        listFailedCommands = []
        while trial < maxTries:
            try:
                RunCommand(*runCommandArgs)
            except subprocess.CalledProcessError as e:
                if trial < maxTries-1:
                    print(colored("Caught exception running condor_submit; retry", "yellow"), flush=True)
                trial += 1
                continue
            break
        if trial >= maxTries:
            listFailedCommands.append(" ".join(runCommandArgs[0:1]))
            return " ".join(runCommandArgs[0:1])
    else:
        print("INFO: dry run enabled; not submitting {} to condor".format(str(Path(subFilename))))


####################################################################################################
# RUN
####################################################################################################
# ---Option Parser
usage = "usage: %prog [options] \nExample: \n./combinePlotsBatchDriver.py -i /home/santanas/Workspace/Leptoquarks/rootNtupleAnalyzer/config/inputListAllCurrent.txt -c analysisClass_genStudies -d /home/santanas/Workspace/Leptoquarks/rootNtupleAnalyzer/data/output -l 100 -x /home/santanas/Data/Leptoquarks/RootNtuples/V00-00-06_2008121_163513/xsection_pb_default.txt -o /home/santanas/Workspace/Leptoquarks/rootNtupleAnalyzer/data/output -s /home/santanas/Workspace/Leptoquarks/rootNtupleAnalyzer/config/sampleListForMerging.txt"

parser = OptionParser(usage=usage)

parser.add_option(
    "-i",
    "--inputList",
    dest="inputList",
    help="list of all datasets to be used (full path required)",
    metavar="LIST",
)

parser.add_option(
    "-c",
    "--code",
    dest="analysisCode",
    help="name of the CODE.C code used to generate the rootfiles (which is the beginning of the root file names before ___)",
    metavar="CODE",
)

parser.add_option(
    "-d",
    "--inputDir",
    dest="inputDir",
    help="the directory INDIR contains the rootfiles with the histograms to be combined (full path required)",
    metavar="INDIR",
)

parser.add_option(
    "-l",
    "--intLumi",
    dest="intLumi",
    help="results are rescaled to the integrated luminosity INTLUMI (in pb-1)",
    metavar="INTLUMI",
)

parser.add_option(
    "-x",
    "--xsection",
    dest="xsection",
    help="the file XSEC contains the cross sections (in pb) for all the datasets (full path required). Use -1 as cross section value for no-rescaling",
    metavar="XSEC",
)

parser.add_option(
    "-o",
    "--outputDir",
    dest="outputDir",
    help="the directory OUTDIR contains the output of the program (full path required)",
    metavar="OUTDIR",
)

parser.add_option(
    "-s",
    "--sampleListForMerging",
    dest="sampleListForMerging",
    help="put in the file SAMPLELIST the name of the sample with the associated strings which should  match with the dataset name (full path required)",
    metavar="SAMPLELIST",
)

parser.add_option(
    "-t",
    "--tablesOnly",
    action="store_true",
    dest="tablesOnly",
    default=False,
    help="only combine tables, do not do plots",
    metavar="TABLESONLY",
)

parser.add_option(
    "-b",
    "--ttbarBkg",
    action="store_true",
    dest="ttbarBkg",
    default=False,
    help="do the ttbar background prediction from data; don't write out any other plots",
    metavar="TTBARBKG",
)

parser.add_option(
    "-q",
    "--qcdClosure",
    action="store_true",
    dest="qcdClosure",
    default=False,
    help="do the QCD observed 1 pass HEEP 1 fail HEEP yield (cej), subtracting non-QCD processes from data using MC; don't write out any other plots",
    metavar="QCDCLOSURE",
)

parser.add_option(
    "-k",
    "--keepInputFiles",
    action="store_true",
    dest="keepInputFiles",
    default=False,
    help="Don't delete individual sample root files (intermediate outputs); defaults to False",
    metavar="KEEPINPUTFILES",
)

parser.add_option(
    "-z", #all the letters that I thought were sensible to use for this were already taken. -Emma
    "--histInclusionList",
    dest="histInclusionList",
    default="",
    help="Use an inclusion list to specify which histograms make it into the output file. Defaults to False",
    metavar="HISTINCLUSIONLIST",
)

parser.add_option(
    "--queue",
    dest="queue",
    default="tomorrow",
    help="name of the lxbatch queue",
    metavar="QUEUE"
)

parser.add_option(
    "--dryRun",
    dest="dryRun",
    help="don't submit jobs to condor, just write the submission files (dry run)",
    metavar="dryRun",
    action="store_true",
    default=False,
)

parser.add_option(
    "-y",
    "--years",
    dest="years",
    help="analysis years",
    metavar="YEARS",
)

parser.add_option(
    "--fitDiagFilepath",
    dest="fitDiagFilepath",
    default=None,
    help="FitDiagnostics root files path (from combine output)",
    metavar="FITDIAGFILEPATH",
)

parser.add_option(
    "--postFitJSON",
    dest="postFitJSON",
    default=None,
    help="Post-fit JSON file containing separated stat/syst uncertainties (fitdiagnostics after processing)",
    metavar="POSTFITJSON",
)

parser.add_option(
    "--postFit",
    dest="postFit",
    default=False,
    action="store_true",
    help="do postfit final selection plots",
    metavar="POSTFIT",
)

parser.add_option(
    "--preFit",
    dest="preFit",
    default=False,
    action="store_true",
    help="do prefit final selection plots",
    metavar="PREFIT",
)

parser.add_option(
    "--fitType",
    dest="fitType",
    default=None,
    help="fit type for postfit final selection plots",
    metavar="FITTYPE",
)

(options, args) = parser.parse_args()

requiredOpts = [options.inputList, options.analysisCode, options.inputDir, options.intLumi, options.xsection, options.outputDir, options.sampleListForMerging, options.years]
requiredOptNames = ["inputList", "analysisCode", "inputDir", "intLumi", "xsection", "outputDir", "sampleListForMerging", "years"]
missingOpts = []
for idx, opt in enumerate(requiredOpts):
    if opt is None:
        missingOpts.append(requiredOptNames[idx])
if len(missingOpts):
    print("ERROR: one or more required options not given:", missingOpts)
    raise RuntimeError(usage)

mergingFiles = SeparateArgs(options.sampleListForMerging)
# ---Check if sampleListForMerging files exist
for mergingFile in mergingFiles:
    if os.path.isfile(mergingFile) is False:
        raise RuntimeError("File " + options.sampleListForMerging + " not found")

xsectionFiles = SeparateArgs(options.xsection)
# ---Check if xsection files exist
for xsectionFile in xsectionFiles:
    if os.path.isfile(xsectionFile) is False:
        raise RuntimeError("File " + options.xsection + " not found")
xsectionFileDict = FillDictFromOptionByYear(options.xsection, options.years)
xsectionDict = {}
for year, xsFile in xsectionFileDict.items():
    xsectionDict[year] = combineCommon.ParseXSectionFile(xsFile)

useInclusionList = False
histoNamesToUse = []
if os.path.isfile(options.histInclusionList) is True:
    print("using inclusion list: ",options.histInclusionList)
    useInclusionList = True
    with open(options.histInclusionList) as f:
        histoNamesToUse = [line.rstrip() for line in f]
#print("including histos: ",histoNamesToUse)

if options.postFit and options.preFit:
    raise RuntimeError("Can't specify both preFit and postFit options.")
if options.postFit or options.preFit:
    if options.fitDiagFilepath is None and options.postFitJSON is None:
        raise RuntimeError("With options preFit or postFit, must provide either (1) path to FitDiagnostics root files with prefit or postfit results using --fitDiagFilepath or (2) post-fit JSON file using --postFitJSON")
    if options.tablesOnly:
        raise RuntimeError("With options preFit or postFit, doesn't make sense to do tables only as tables are not affected (at least at the present time).")
if options.postFit and options.fitType is None:
    raise RuntimeError("With option postFit, must provide fit type (b or sb).")

print("Launched like:")
print("python ", end=' ')
for arg in sys.argv:
    print(" " + arg, end=' ')
print()

doPDFReweight2016LQSignals = False
if doPDFReweight2016LQSignals:
    print("Doing PDF reweighting for 2016 LQ B/D signal samples")

if not os.path.exists(options.outputDir):
    os.makedirs(options.outputDir)

sampleLists = FillDictFromOptionByYear(options.sampleListForMerging, options.years)
dictSamples = {}
for year in SeparateArgs(options.years):
    sampleList = sampleLists[year]
    dictSamples[year] = combineCommon.GetSamplesToCombineDict(sampleList)
    
# --- Declare efficiency tables
dictFinalTables = {}
# --- Declare histograms
dictFinalHisto = {}
# --- Samples to save in final histo dict
samplesToSave = {}
sampleFiles = {}

if options.ttbarBkg or options.qcdClosure:
    raise RuntimeError("TTBar Bkg and QCD closure options not yet implemented here")

# if options.ttbarBkg:
#     ttbarDataRawSampleName = "TTBarUnscaledRawFromDATA"
#     nonTTbarAMCBkgSampleName = "NONTTBARBKG_amcatnloPt_amcAtNLODiboson_emujj"
#     samplesToSave.extend([ttbarDataRawSampleName, nonTTbarAMCBkgSampleName])
# if options.qcdClosure:
#     qcdDataSampleName = "SinglePhoton_all"
#     nonQCDBkgSampleName = "ALLBKG_powhegTTBar_ZJetWJetPt_amcAtNLODiboson"  # Z, W, TTBar, SingleTop, Diboson, gamma+jets
#     samplesToSave.extend([qcdDataSampleName, nonQCDBkgSampleName])

# check to make sure we have xsections for all samples
inputLists = FillDictFromOptionByYear(options.inputList, options.years)
inputDirs = FillDictFromOptionByYear(options.inputDir, options.years)
dictDatasetsFileNames = {}
for year in SeparateArgs(options.years):
    inputList = inputLists[year]
    inputDir = inputDirs[year]
    for lin in open(inputList):
        lin = lin.strip("\n")
        if lin.startswith("#"):
            continue
        dataset_fromInputList = lin.split("/")[-1].split(".")[0]
        dataset_fromInputList = dataset_fromInputList.replace("_tree", "")
        xsection_val = combineCommon.lookupXSection(
            combineCommon.SanitizeDatasetNameFromInputList(
                dataset_fromInputList.replace("_tree", "")
            ),
            xsectionDict[year]
        )
    foundAllFiles, dictDatasetsFileNames[year] = combineCommon.FindInputFiles(inputList, options.analysisCode, inputDir)
    if not foundAllFiles:
        raise RuntimeError("Some files not found for year {}.".format(year))
    else:
        print("\bDone.  All root/dat files are present for year {}.".format(year))
print()


CreateAndSubmitJobs(dictSamples)

# now handle special backgrounds
# FIXME: will need special handling of these
# if options.ttbarBkg:
#     # special actions for TTBarFromData
#     # subtract nonTTbarBkgMC from TTbarRaw
#     # FIXME: we hardcode the sample names for now
#     ttbarDataPredictionTable = dictFinalTables[ttbarDataRawSampleName]
#     # nonTTbarAMCBkgSampleName = 'NONTTBARBKG_amcatnloPt_emujj'
#     # move to amcAtNLO diboson
#     nonTTbarAMCBkgTable = dictFinalTables[nonTTbarAMCBkgSampleName]
#     ttBarPredName = "TTBarFromDATA"
#     # Mar17 fixing muon pt and eta-->2.4
#     Rfactor = 0.418559  # Ree,emu = Nee/Nemu[TTbarMC]
#     errRfactor = 0.002474
#     print("TTBar data-driven: Using Rfactor =", Rfactor, "+/-", errRfactor)
#     print("TTBar data-driven: Using non-ttbar background sample:", nonTTbarAMCBkgSampleName)
#     # print '0) WHAT DOES THE RAW DATA TABLE LOOK LIKE?'
#     # WriteTable(ttbarDataPredictionTable, ttbarDataRawSampleName, outputTableFile)
#     # remove the x1000 from the nonTTbarBkgMC
#     combineCommon.ScaleTable(nonTTbarAMCBkgTable, 1.0 / 1000.0, 0.0)
#     # print '1) WHAT DOES THE SCALED MC TABLE LOOK LIKE?'
#     # WriteTable(nonTTbarMCBkgTable, nonTTbarMCBkgSampleName, outputTableFile)
#     # subtract the nonTTBarBkgMC from the ttbarRawData, NOT zeroing entries where we run out of data
#     combineCommon.SubtractTables(nonTTbarAMCBkgTable, ttbarDataPredictionTable)
#     # print '2) WHAT DOES THE SUBTRACTEDTABLE LOOK LIKE?'
#     # WriteTable(ttbarDataPredictionTable, ttBarPredName, outputTableFile)
#     # scale by Ree,emu
#     combineCommon.ScaleTable(ttbarDataPredictionTable, Rfactor, errRfactor)
#     # print '3) WHAT DOES THE RfactroCorrectedTABLE LOOK LIKE?'
#     # WriteTable(ttbarDataPredictionTable, ttBarPredName, outputTableFile)
#     combineCommon.SquareTableErrorsForEfficiencyCalc(ttbarDataPredictionTable)
#     combineCommon.CalculateEfficiency(ttbarDataPredictionTable)
#     # print '4) WHAT DOES THE SCALEDTABLE AFTER EFF CALCULATION LOOK LIKE?'
#     combineCommon.WriteTable(ttbarDataPredictionTable, ttBarPredName, outputTableFile)
# 
# if options.qcdClosure:
#     # special actions for the QCD closure test
#     # subtract nonQCD from QCDData yield
#     qcdDataTable = dictFinalTables[qcdDataSampleName]
#     nonQCDBkgTable = dictFinalTables[nonQCDBkgSampleName]
#     qcdClosureSampleName = "QCDClosureObserved"
#     # print "TTBar data-driven: Using non-ttbar background sample:", nonTTbarAMCBkgSampleName
#     # print '0) WHAT DOES THE RAW DATA TABLE LOOK LIKE?'
#     # WriteTable(ttbarDataPredictionTable, ttbarDataRawSampleName, outputTableFile)
#     # remove the x1000 from the nonQCDBkgMC
#     combineCommon.ScaleTable(nonQCDBkgTable, 1.0 / 1000.0, 0.0)
#     # print '1) WHAT DOES THE SCALED MC TABLE LOOK LIKE?'
#     # WriteTable(nonTTbarMCBkgTable, nonTTbarMCBkgSampleName, outputTableFile)
#     # subtract the nonTTBarBkgMC from the ttbarRawData, NOT zeroing entries where we run out of data
#     combineCommon.SubtractTables(nonQCDBkgTable, qcdDataTable)
#     # print '2) WHAT DOES THE SUBTRACTEDTABLE LOOK LIKE?'
#     # WriteTable(ttbarDataPredictionTable, ttBarPredName, outputTableFile)
#     combineCommon.SquareTableErrorsForEfficiencyCalc(qcdDataTable)
#     combineCommon.CalculateEfficiency(qcdDataTable)
#     # print '4) WHAT DOES THE SCALEDTABLE AFTER EFF CALCULATION LOOK LIKE?'
#     combineCommon.WriteTable(qcdDataTable, qcdClosureSampleName, outputTableFile)
# 


# FIXME for ttbar/QCD
# if options.ttbarBkg:
#     # special actions for TTBarFromData
#     # subtract nonTTbarBkgMC from TTbarRaw
#     ttbarDataPredictionHistos = dictFinalHisto[ttbarDataRawSampleName]
#     # print 'ttbarDataPredictionHistos:',ttbarDataPredictionHistos
#     for n, histo in ttbarDataPredictionHistos.items():
#         # subtract the nonTTBarBkgMC from the ttbarRawData
#         # find nonTTbarMCBkg histo; I assume they are in the same order here
#         histoToSub = dictFinalHisto[nonTTbarAMCBkgSampleName][n]
#         ## also write histos that are subtracted
#         # histToSub.Write()
#         # print 'n=',n,'histo=',histo
#         outputTfile.cd()
#         histoTTbarPred = histo.Clone()
#         histoTTbarPred.Add(histoToSub, -1)
#         # scale by Rfactor
#         histoTTbarPred.Scale(Rfactor)
#         histoTTbarPred.SetName(
#             re.sub(
#                 "__.*?__",
#                 "__" + ttBarPredName + "__",
#                 histoTTbarPred.GetName(),
#                 flags=re.DOTALL,
#             )
#         )
#         histoTTbarPred.Write()
# 
# if options.qcdClosure:
#     # special actions for QCDClosure observed
#     # subtract nonQCDBkgMC from data
#     qcdClosureHistos = dictFinalHisto[qcdDataSampleName]
#     # print 'qcdClosureHistos:',qcdClosureHistos
#     for n, histo in qcdClosureHistos.items():
#         # find nonTTbarMCBkg histo; assume they are in the same order here
#         histoToSub = dictFinalHisto[nonQCDBkgSampleName][n]
#         ## also write histos that are subtracted
#         # histToSub.Write()
#         # print 'n=',n,'histo=',histo
#         histoQCDClosure = histo.Clone()
#         histoQCDClosure.Add(histoToSub, -1)
#         histoQCDClosure.SetName(
#             re.sub(
#                 "__.*?__",
#                 "__" + qcdClosureSampleName + "__",
#                 histoQCDClosure.GetName(),
#                 flags=re.DOTALL,
#             )
#         )
#         histoQCDClosure.Write()
# 
