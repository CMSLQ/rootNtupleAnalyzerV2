#!/usr/bin/env python3

# ---Import
import sys
import string
from optparse import OptionParser
import os
from ROOT import TFile, gROOT, SetOwnership, TObject
import re
import glob
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


def GetURLs(fileList):
    fileList = [filename.replace("/eos/cms/", "root://eoscms//").replace("/eos/user/", "root://eosuser//") for filename in fileList]
    return fileList


def CreateAndSubmitJobs(sampleDict):
    condorSubfileDir = "combineRunIIPlotsCondor"
    Path(condorSubfileDir).mkdir(parents=True, exist_ok=True)
    samplesToCombine = sampleDict.keys()
    shFilename = condorSubfileDir + "/condor.sh"
    WriteCondorShFile(sampleDict, shFilename)
    print("INFO: wrote sh file to {}".format(shFilename))
    subFilename = shFilename.replace(".sh", ".sub")
    WriteCondorSubFile(len(samplesToCombine), subFilename)
    print("INFO: wrote sub file to {}".format(subFilename))
    SubmitCondorJob(str(Path(subFilename).name), condorSubfileDir)


def WriteCondorShFile(sampleDict, filename):
    execName = "hadd"
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
        for idx, sample in enumerate(sampleDict.keys()):
            outputFile = "analysisClass_{}_combinedRunII_plots.root".format(sample)
            execArgs = " -fk207 {} {}".format(outputFile, ' '.join(GetURLs(sampleDict[sample])))
            outputfile.write("if [ $1 -eq {} ]; then\n".format(idx))
            outputfile.write('    ' + execName + execArgs + "\n")
            outputfile.write('    retVal=$?\n')
            outputfile.write('    if [ $retVal -ne 0 ]; then\n')
            outputfile.write('      echo "'+execName+' return error code=$retVal; quitting here."\n')
            outputfile.write('      exit $retVal\n')
            outputfile.write('    fi\n')
            outputfile.write('    xrdcp -f ' + outputFile + ' {}/{}'.format(options.outputDir, outputFile) + '\n')
            outputfile.write('fi\n')


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
        # make sure the job finishes with exit code 0
        # condorFile.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        condorFile.write('max_retries = 3\n')
        # condorFile.write('should_transfer_files = YES\n')
        # lqAnaBase = os.getenv("LQANA")
        # inputFilesToTransfer = [lqAnaBase + "/scripts/combinePlotsBatch.py", lqAnaBase + "/scripts/combineCommon.py"]
        # inputFilesToTransfer.extend([str(Path(options.inputList).resolve()), str(Path(options.sampleListForMerging).resolve()), str(Path(options.xsection).resolve())])
        # filesToTransfer = ",".join(inputFilesToTransfer)
        # filesToTransfer += ","
        # filesToTransfer += options.filesToTransfer
        # condorFile.write('transfer_input_files = '+filesToTransfer+'\n')
        # condorFile.write('transfer_input_files = '+filesToTransfer+'\n')
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


def GetSamplesAndFiles(sampleList, analysisName, yearList, basePath, plotFileTemplate="analysisClass_lq_eejj_{}_plots.root"):
    dictSamplesAndFiles = {}
    for sample in sampleList:
        fileList = []
        for year in yearList:
            path = basePath.format(year, analysisName)
            globStr = path + plotFileTemplate.format(sample)
            files = glob.glob(globStr)
            if len(files) != 1:
                raise RuntimeError("Found {} files ({}) while globbing for {}, while 1 and only 1 was expected".format(len(files), files, globStr))
            fileList.extend(files)
        dictSamplesAndFiles[sample] = fileList
    return dictSamplesAndFiles


####################################################################################################
# RUN
####################################################################################################
eosBasePath = "/eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/{}/{}/"
includeQCD = True
qcdFilePath = "output_qcdSubtractedYield/"  # after the first analysisName dir, and assuming the first analysisNameDir needs 'qcd_' prepended
qcdFileName = "qcdSubtracted_plots.root"
samplesToCombine = ["ZJet_amcatnlo_ptBinned_IncStitch", "TTTo2L2Nu", "SingleTop", "DIBOSON_nlo", "DATA"]
samplesToCombine.extend(["LQToBEle_M-{}_pair".format(mass) for mass in range(300, 3100, 100)])
samplesToCombine.extend(["LQToDEle_M-{}_pair".format(mass) for mass in range(300, 3100, 100)])
# ---Option Parser
usage = "usage: %prog [options] \nExample: \ncombineRunIIPlots.py -a eejj_13nov2024_bdt_LQToDEle/output_cutTable_lq_eejj_BDT -o /eos/cms/store/group/phys_exotica/leptonsPlusJets/LQ/scooper/ultralegacy/analysis/2016preVFP/eejj_13nov2024_bdt_LQToDEle/combinedRunIIPlots -y 2016preVFP,2016postVFP"

parser = OptionParser(usage=usage)

parser.add_option(
    "-a",
    "--analysisName",
    dest="analysisName",
    help="the name of the analysis directory (e.g., eejj_13nov2024_bdt_LQToDEle)",
    default = None,
    metavar="ANALYSISNAME",
)
parser.add_option(
    "-o",
    "--outputDir",
    dest="outputDir",
    help="the directory OUTDIR contains the output of the program (full path required)",
    default = None,
    metavar="OUTDIR",
)
parser.add_option(
    "-y",
    "--years",
    dest="years",
    help="comma-separated list of analysis years to combine, or 'all'",
    default="all",
    metavar="YEARS",
)
parser.add_option(
    "-t",
    "--dryRun",
    dest="dryRun",
    help="don't submit jobs to condor, just write the submission files (dry run)",
    metavar="dryRun",
    action="store_true",
    default=False,
)
parser.add_option(
    "--queue",
    dest="queue",
    default="microcentury",
    help="name of the lxbatch queue",
    metavar="QUEUE"
)


(options, args) = parser.parse_args()

if options.analysisName is None:
    raise RuntimeError("Must specify analysis dir name with '-a' option")
if options.outputDir is None:
    raise RuntiemError("Must specify output dir name with '-o' option")

if not os.path.exists(options.outputDir):
    os.makedirs(options.outputDir)

yearList = ["2016preVFP", "2016postVFP","2017", "2018"] if options.years == "all" else options.years.split(",")
analysisName = options.analysisName
analysisNameBase = analysisName
if analysisName.find("/") >= 0:
    analysisNameBase = analysisName[0:analysisName.find("/")+1]

dictSamplesAndFiles = GetSamplesAndFiles(samplesToCombine, analysisName, yearList, eosBasePath)
if includeQCD:
    dictSamplesAndFiles.update(GetSamplesAndFiles(["qcd"], "qcd_"+analysisNameBase+qcdFilePath, yearList, eosBasePath, qcdFileName))
CreateAndSubmitJobs(dictSamplesAndFiles)
