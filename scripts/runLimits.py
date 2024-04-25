#!/usr/bin/env python3
import os
from optparse import OptionParser
import subprocess
import shlex
import glob
from termcolor import colored
import time
import datetime
import numpy as np
import re
import math
import multiprocessing
import traceback
from pathlib import Path
from bisect import bisect
from tabulate import tabulate
from ROOT import TFile, TGraph, TSpline3
from combineCommon import SeparateDatacards
from BR_Sigma_EE_vsMass import BR_Sigma_EE_vsMass
from ComboPlotLQ1 import ComboPlot

ansi_escape = re.compile(r'\x1B(?:[@-Z\\-_]|\[[0-?]*[ -/]*[@-~])')

def GetExecTimeStr(startTime, stopTime):
    delta = datetime.timedelta(seconds=stopTime - startTime)
    hours, rem = divmod(delta.seconds, 3600)
    minutes, seconds = divmod(rem, 60)
    execTimeStr = "{}s".format(seconds)
    if minutes > 0:
        execTimeStr = "{}m".format(minutes) + execTimeStr
    if hours > 0:
        execTimeStr = "{}h".format(hours) + execTimeStr
    return execTimeStr


def ExtractCMSSWEnv(cmsswPreambleCmds):
    cmd = cmsswPreambleCmds + 'echo ~~~~START_ENVIRONMENT_HERE~~~~ && set'
    try:
        process = subprocess.run(cmd, check=True, shell=True, env={}, stdout=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        print(colored("ExtractCMSSWEnv had an error; output: {}".format(result), "red"), flush=True)
        raise e
    env = process.stdout.decode().splitlines()
    envDict = {}
    record = False
    idx = 0
    while idx < len(env):
        e = env[idx]
        if record:
            # print("e={}".format(e))
            e = e.strip().split('=')
            # print("e.strip().split('=')={}".format(e))
            if e[1] == "'":
                nextLine = env[idx+1].strip()
                if "=" not in nextLine and nextLine == "'":
                    e[1] += nextLine
                    idx+=2
                else:
                    raise RuntimeError("Don't know how to process {} and {}".format(e, env[idx+1]))
            # print("{}={}".format(e[0], e[1]))
            envDict[e[0]] = e[1]
        elif e.strip() == '~~~~START_ENVIRONMENT_HERE~~~~': 
            record = True
        idx += 1
    if "SHELLOPTS" in envDict.keys():
        del envDict["SHELLOPTS"]
    return envDict


# def RunCommand(cmd, workingDir=None, cleanEnv=False, shell=False):
def RunCommand(cmd, workingDir=None, env=None, suppressOutput=False):
    print(colored("\t{}".format(cmd), "green"), flush=True)
    try:
        # if not shell:
        #     cmd = shlex.split(cmd)
        # if cleanEnv:
        #     process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env={}, shell=shell)
        # else:
        #     process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, shell=shell)

        useShell = False
        if env is None:
            useShell = True
        #     cmd = shlex.split(cmd)
        cmd = shlex.split(cmd)
        process = subprocess.run(cmd, shell=useShell, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env=env)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        print(colored("RunCommand had an error; output: {}".format(result), "red"), flush=True)
        raise e
    if not suppressOutput:
        print(process.stdout.decode(), flush=True)
        #print(process.stderr.decode())


def ConvertDatacardToWorkspace(datacard, mass):
    workspaceFileName = datacard+".m{}.root".format(mass)
    if doBatch:
        cmd = 'text2workspace.py {} -o {} -m {}'.format(datacard, workspaceFileName, mass)
        # print(cmsswEnv, flush=True)
        RunCommand(cmd, None, cmsswEnv)
    else:
        cmd = 'text2workspace.py {} -o {} -m {}'.format(datacard, workspaceFileName, mass)
        RunCommand(cmd)
    return Path(workspaceFileName).resolve()


def FindCardWorkspace(cardFile, mass, exceptOnAbsence=False):
    workspaceFileName = cardFile+".m{}.root".format(mass)
    listOfWorkspaceFiles = sorted(glob.glob(workspaceFileName), key=os.path.getmtime)
    if exceptOnAbsence and len(listOfWorkspaceFiles) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(workspaceFileName, listOfWorkspaceFiles))
    if len(listOfWorkspaceFiles) > 0:
        return Path(listOfWorkspaceFiles[-1]).resolve()


def GenerateAsimovToyData(args):
    workspace, mass, toysDir, dictAsimovToysByScaleFactor, signalScaleFactor = args
    if doBatch:
        cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.{} -m {} {}'.format(workspace, signalScaleFactor, mass, commonCombineArgs.format(signalScaleFactor))
        RunCommand(cmd, toysDir, cmsswEnv, True)
    else:
        cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov.{} -m {} {}'.format(workspace, signalScaleFactor, mass, commonCombineArgs.format(signalScaleFactor))
        RunCommand(cmd, toysDir)
    toyFile = Path(sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)[-1]).resolve()
    dictAsimovToysByScaleFactor[signalScaleFactor] = toyFile
    return toyFile


def FindAsimovToyData(mass, signalScaleFactor, toysDir):
    listOfToyFiles = sorted(glob.glob(toysDir+'/higgsCombine.asimov.{}.GenerateOnly.mH{}.*.root'.format(signalScaleFactor, mass)), key=os.path.getmtime)
    if len(listOfToyFiles) > 0:
        return Path(listOfToyFiles[-1]).resolve()


def FindFile(globString):
    fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    if len(fileList) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(globString, fileList))
    return Path(fileList[-1]).resolve()


def GetFileList(globString):
    fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    if len(fileList) < 1:
        raise RuntimeError("Globbing for {} did not result in any files.".format(fileList))
    return fileList


def MakeCondorDirs(dirName):
    condorDir = dirName.strip("/")+"/condor"
    condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    condorOutDir = dirName.strip("/")+"/condor/out"
    condorErrDir = dirName.strip("/")+"/condor/error"
    condorLogDir = dirName.strip("/")+"/condor/log"
    dirsToHave = [condorDir, condorSubfileDir, condorOutDir, condorErrDir, condorLogDir]
    for idir in dirsToHave:
        if not os.path.isdir(idir):
            print("INFO: Making directory", idir, flush=True)
            Path(idir).mkdir(exist_ok=True)


def WriteCondorSubFile(filename, shFilename, combinedOutputFile=""):
    basename = shFilename.replace("condor_", "").replace(".sh", "")
    with open(filename, "w") as subfile:
        subfile.write("executable = subFiles/" + shFilename + "\n")
        subfile.write("arguments = $(ProcId)\n")
        subfile.write("output                = out/"   + basename + ".$(ClusterId).$(ProcId).out\n")
        subfile.write("error                 = error/" + basename + ".$(ClusterId).$(ProcId).err\n")
        subfile.write("log                   = log/"   + basename + ".$(ClusterId).log\n")
        subfile.write("\n")
        if combinedOutputFile != "":
            subfile.write("transfer_output_files = {}\n".format(combinedOutputFile))
            subfile.write("\n")
        subfile.write("# Send the job to Held state on failure.\n")
        subfile.write("on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)\n")
        subfile.write("\n")
        subfile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.\n")
        subfile.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
        subfile.write("\n")
        subfile.write("+JobFlavour=\"workday\"\n")
        subfile.write("MY.WantOS = \"el7\"\n")
        subfile.write("queue 1\n")


def WriteCondorShFile(filename, combineCmd, mass, signalScaleFactor=1.0, quantile=-1, rMin=-1, rMax=-1):
    with open(filename, "w") as shfile:
        shfile.write("#!/bin/sh\n")
        shfile.write("ulimit -s unlimited\n")
        shfile.write("set -e\n")
        shfile.write("cd " + cmsswDir + "\n")
        # shfile.write("export SCRAM_ARCH=slc7_amd64_gcc900\n")  # try without this
        shfile.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
        shfile.write("eval `scramv1 runtime -sh`\n")
        # shfile.write("cd {}\n".format(condorDir))
        shfile.write("cd -\n")
        shfile.write("\n")
        shfile.write("if [ $1 -eq 0 ]; then\n")
        if rMin > 0 and rMax > 0:
            sigSFRound = round(signalScaleFactor, 6)
            scanPoints = 50
            stepSize = (rMax-rMin)/scanPoints
            print("INFO: Writing commands to compute grid of {} limits in the range r=[{}, {}] by steps of {}".format(scanPoints, rMin, rMax, stepSize))
            for rVal in np.arange(rMin, rMax+stepSize, stepSize):
                thisStepCmd = combineCmd
                thisStepCmd += ' -n .signalScaleFactor{}.POINT.{}'.format(sigSFRound, rVal)
                thisStepCmd += ' --singlePoint {}'.format(rVal)
                shfile.write("  {}\n".format(thisStepCmd))
            quant = "quant{}.".format(quantile) if quantile > 0 else ""
            combinedOutputFile = "higgsCombine.signalScaleFactor{}.HybridNew.mH{}.{}grid.root".format(sigSFRound, mass, quant)
            shfile.write("  hadd -fk207 {} higgsCombine.*.root\n".format(combinedOutputFile))
        else:
            shfile.write("  {}\n".format(combineCmd))
            combinedOutputFile = ""
        shfile.write("fi\n")
    return combinedOutputFile


def GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor, batch, rMin, rMax):
    rAbsAcc = 0.0001
    clsAcc = 0 if batch else 0.001
    toys = 10000
    # toys = 500  # reduced for shape-based limits as test, but still took hours
    cmd = '{}'.format(workspace)
    cmd += ' -v1'
    cmd += ' -M HybridNew'
    cmd += ' --LHCmode LHC-limits'
    cmd += ' --saveHybridResult'
    cmd += ' --saveToys'
    cmd += ' --seed -1'
    if not batch:
        cmd += ' --iterations 2'
        cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += ' --clsAcc {}'.format(clsAcc)
    if not batch:
        cmd += ' --rMin {}'.format(rMin)
        cmd += ' --rMax {}'.format(rMax)
    cmd += ' -T {}'.format(toys)
    cmd += ' -m {}'.format(mass)
    if not batch:
        cmd += ' -n .signalScaleFactor{}'.format(round(signalScaleFactor, 6))
    # cmd += ' -H AsymptoticLimits'
    # cmd += ' --fork 4'
    cmd += commonCombineArgs.format(signalScaleFactor)
    if quantile > 0:
        cmd += ' --expectedFromGrid {}'.format(quantile)
        if genAsimovToyFile != "":
            cmd += ' -D {}:toys/toy_asimov'.format(genAsimovToyFile)
    return cmd


def SubmitHybridNewBatch(args):
    workspace, mass, dirName, listFailedCommands, quantile, genAsimovToyFile, signalScaleFactor, rMin, rMax = args
    condorDir = dirName.strip("/")+"/condor"
    condorSubfileDir = dirName.strip("/")+"/condor/subFiles"
    # condorOutDir = dirName.strip("/")+"/condor/out"
    # condorErrDir = dirName.strip("/")+"/condor/error"
    # condorLogDir = dirName.strip("/")+"/condor/log"
    # dirsToHave = [condorDir, condorSubfileDir, condorOutDir, condorErrDir, condorLogDir]
    # for idir in dirsToHave:
    #     if not os.path.isdir(idir):
    #         print("INFO: Making directory", idir, flush=True)
    #         Path(idir).mkdir(exist_ok=True)
    if cmsswDir is None:
        raise RuntimeError("Need to specify valid CMSSW area with combine checked out.")
    #cwd = os.getcwd()
    #os.chdir(condorDir)
    cmdArgs = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor, True, rMin, rMax)
    cmd = "combine " + cmdArgs
    quantileStr = str(quantile).replace(".", "p")
    taskName = 'hybridNewLimits.M{}.{}'.format(mass, quantileStr)
    if signalScaleFactor != 1.0:
        taskName += '.signalScaleFactor{}'.format(round(signalScaleFactor, 6))
    shFilename = condorSubfileDir + "/condor_" + taskName + ".sh"
    condorFile = shFilename.replace(".sh", ".sub")
    combinedOutputFile = WriteCondorShFile(shFilename, cmd, mass, signalScaleFactor, quantile, rMin, rMax)
    WriteCondorSubFile(condorFile, Path(shFilename).name, combinedOutputFile)
    # # cmd = cmsswPreambleCmd
    # # cmd += "cd {} && combineTool.py ".format(cwd + "/" + condorDir)
    # cmd = "combineTool.py "
    # cmd += cmdArgs
    # cmd += ' --job-mode condor'
    # cmd += ' --task-name {}'.format(taskName)
    # cmd += ' --sub-opts=\'+JobFlavour="workday"\\nMY.WantOS = "el7"\''
    # # cmd += ' --dry-run'
    # # RunCommand(cmd, dirName, True, True)
    # runCommandArgs = [cmd, condorDir, cmsswEnv]
    # # now submit
    # # condorFile = "condor_hybridNewLimits.M{}.{}.sub".format(mass, quantileStr)
    # # cmd = "source /etc/profile && k5start -f {}".format(kerberosKeytab) + " && condor_submit {}".format(condorFile)
    # # runCommandArgs = [cmd, cwd + "/" + condorDir, True, True]
    # try:
    #     RunCommand(*runCommandArgs)
    # except subprocess.CalledProcessError as e:
    #     try:
    #         print(colored("Caught exception running combineTool batch submission command; retry", "yellow"), flush=True)
    #         RunCommand(*runCommandArgs)
    #     except subprocess.CalledProcessError as e:
    #         listFailedCommands.append( runCommandArgs[0:1])
    #         return runCommandArgs[0:1]
    # # cmd = "condor_submit {}".format(condorFile)
    # # RunCommand(cmd, cwd + "/" + condorDir)
    # # os.chdir(cwd)
    cmd = "condor_submit subFiles/{}".format(Path(condorFile).name)
    runCommandArgs = [cmd, condorDir, cmsswEnv, True]
    try:
        RunCommand(*runCommandArgs)
    except subprocess.CalledProcessError as e:
        try:
            print(colored("Caught exception running condor_submit; retry", "yellow"), flush=True)
            RunCommand(*runCommandArgs)
        except subprocess.CalledProcessError as e:
            listFailedCommands.append( runCommandArgs[0:1])
            return runCommandArgs[0:1]


# http://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/part3/nonstandard/#nuisance-parameter-impacts
def MakeImpacts(workspace, mass, dirName, signalScaleFactor):
    impactsDir = dirName.strip("/")+"/impacts"
    if not os.path.isdir(impactsDir):
        print("INFO: Making directory", impactsDir, flush=True)
        Path(impactsDir).mkdir(exist_ok=True)
    if cmsswDir is None:
        raise RuntimeError("Need to specify valid CMSSW area with combine checked out.")
    cwd = os.getcwd()
    cmd = cmsswPreambleCmd
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} --doInitialFit --robustFit 1".format(workspace, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, cmsswEnv, True)
    cmd = cmsswPreambleCmd
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} --robustFit 1 --doFits".format(workspace, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, cmsswEnv, True)
    cmd = cmsswPreambleCmd
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + impactsDir)
    cmd += " -M Impacts -d {} -m {} -o impacts.m{}.json".format(workspace, mass, mass)
    cmd += commonCombineArgs.format(signalScaleFactor)
    RunCommand(cmd, dirName, cmsswEnv, True)
    cmd = cmsswPreambleCmd
    cmd += "cd {} && plotImpacts.py ".format(cwd + "/" + impactsDir)
    cmd += " -i impacts.m{}.json -o impacts.m{}".format(mass, mass)
    RunCommand(cmd, dirName, cmsswEnv, True)

    
def RunHybridNewInteractive(workspace, mass, dirName, quantile=-1, genAsimovToyFile="", signalScaleFactor=1.0):
    cmd = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile, signalScaleFactor)
    cmd = 'combine ' + cmd
    RunCommand(cmd, dirName)
    quantileString = "quant{:.3f}".format(quantile) if quantile > 0 else "."
    globString = dirName+'/higgsCombineTest.HybridNew.mH{}.*.{}.root'.format(mass, quantileString)
    return FindFile(globString)


def RunAsymptotic(workspace, mass, dirName, batch, blinded=True, signalScaleFactor=1.0):
    rAbsAcc = 0.0001
    cmd = 'combine'
    cmd += ' {}'.format(workspace)
    cmd += ' -v1'
    cmd += ' -M AsymptoticLimits'
    cmd += ' --seed -1'
    cmd += ' -m {}'.format(mass)
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += commonCombineArgs.format(signalScaleFactor)
    if blinded:
        cmd += ' --run blind --expectSignal 0'
    RunCommand(cmd, dirName, cmsswEnv)
    return Path(sorted(glob.glob(dirName+'/higgsCombineTest.AsymptoticLimits.mH{}.*.root'.format(mass)), key=os.path.getmtime)[-1]).resolve()


def ComputeLimitsFromGrid(workspace, mass, dirName, filename, quantile):
    rAbsAcc = 0.00001
    rRelAcc = 0.005
    plotFilename = "limit_scan_m{}_quant{}.pdf".format(mass, quantile)
    cmd = 'combine'
    cmd += ' {}'.format(workspace)
    # cmd += ' -v1'
    cmd += ' -M HybridNew'
    cmd += ' --LHCmode LHC-limits'
    cmd += ' --readHybridResults'
    cmd += ' --grid={}'.format(filename)
    cmd += ' --expectedFromGrid {}'.format(quantile)
    cmd += ' -m {}'.format(mass)
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += ' --rRelAcc {}'.format(rRelAcc)
    cmd += ' --plot={}'.format(plotFilename)
    RunCommand(cmd, dirName, cmsswEnv)
    return Path(sorted(glob.glob(dirName+'/higgsCombineTest.HybridNew.mH{}.{}.root'.format(mass, "quant{:.3f}".format(quantile))), key=os.path.getmtime)[-1]).resolve()


def ExtractLimitResult(rootFile):
    # higgsCombineTest.MethodName.mH$MASS.[word$WORD].root
    # higgsCombineTest.HybridNew.mH120.quant0.500.root
    if not os.path.isfile(rootFile):
        raise RuntimeError("ERROR: Did not find the root file {}. Exiting.".format(rootFile))
    tfile = TFile.Open(str(rootFile))
    limitTree = tfile.Get("limit")
    bytesRead = limitTree.GetEntry(0)
    if bytesRead <= 0:
        raise RuntimeError("ERROR: Something went wrong: read {} bytes from 'limit' tree in file {}. Exiting.".format(bytesRead, tfile))
    limit = limitTree.limit
    limitErr = limitTree.limitErr
    quantile = limitTree.quantileExpected
    signalScaleParam = 0
    if "trackedParam_signalScaleParam" in limitTree.GetListOfBranches():
        # print("Scale original limit r={}+/-{} by signalScaleParam={}; final limit={}".format(limit, limitErr, limitTree.trackedParam_signalScaleParam, limit/limitTree.trackedParam_signalScaleParam))
        limit /= limitTree.trackedParam_signalScaleParam
        signalScaleParam = limitTree.trackedParam_signalScaleParam
    tfile.Close()
    return limit, limitErr, quantile, signalScaleParam


def ExtractAsymptoticLimitResult(rootFile):
    if not os.path.isfile(rootFile):
        raise RuntimeError("ERROR: Did not find the root file {}. Exiting.".format(rootFile))
    tfile = TFile.Open(str(rootFile))
    limitTree = tfile.Get("limit")
    bytesRead = limitTree.GetEntry(0)
    if bytesRead <= 0:
        raise RuntimeError("ERROR: Something went wrong: read {} bytes from 'limit' tree in file {}. Exiting.".format(bytesRead, tfile))
    limits = []
    limitErrs = []
    quantiles = []
    for iEntry in range(0, limitTree.GetEntries()):
        limitTree.GetEntry(iEntry)
        limits.append(limitTree.limit)
        limitErrs.append(limitTree.limitErr)
        quantiles.append(round(limitTree.quantileExpected, 3))
    tfile.Close()
    return limits, limitErrs, quantiles


def ReadXSecFile(filename):
    xsThByMass = {}
    yPDFupByMass = {}
    yPDFdownByMass = {}
    with open(os.path.expandvars(filename), "r") as xsecFile:
        for line in xsecFile:
            line = line.strip()
            if line.startswith("#"):
                continue
            split = line.split()
            if len(split) != 7:
                raise RuntimeError("length of this line is not 7; don't know how to handle it. Quitting.  Line looks like '"+line+"'")
            mass = float(split[0])
            xs = float(split[1])
            xsThByMass[mass] =  xs
            yPDFupByMass[mass] = xs*(1+float(split[5])/100.)
            yPDFdownByMass[mass] = xs*(1-float(split[6])/100.)
    return xsThByMass, yPDFupByMass, yPDFdownByMass


def CreateArraysForPlotting(xsecLimitsByMassAndQuantile):
    massList = sorted(list(xsecLimitsByMassAndQuantile.keys()))
    shadeMassList = []
    xs_medExpList = []
    xs_oneSigmaExpList = []
    xs_twoSigmaExpList = []
    xs_obsList = []
    if str(-1) in xsecLimitsByMassAndQuantile[list(xsecLimitsByMassAndQuantile.keys())[0]].keys():
        hasObserved = True
    else:
        hasObserved = False
    for mass in massList:
        try:
            xs_medExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.5)])
        except KeyError:
            xs_medExpList.append(xsecLimitsByMassAndQuantile[mass]["0.500"])
        if hasObserved:
            xs_obsList.append(xsecLimitsByMassAndQuantile[mass][str(-1)])
        try:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.16)])
        except KeyError:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass]["0.160"])
        xs_twoSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.025)])
    for mass in reversed(massList):
        try:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.84)])
        except KeyError:
            xs_oneSigmaExpList.append(xsecLimitsByMassAndQuantile[mass]["0.840"])
        xs_twoSigmaExpList.append(xsecLimitsByMassAndQuantile[mass][str(0.975)])
    masses = np.array(massList, dtype="f")
    shadeMasses = np.concatenate([masses, np.flip(masses)])
    xsMedExp = np.array(xs_medExpList, dtype="f")
    xsObs = np.array(xs_obsList, dtype="f")
    xsOneSigmaExp = np.array(xs_oneSigmaExpList, dtype="f")
    xsTwoSigmaExp = np.array(xs_twoSigmaExpList, dtype="f")
    return masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs


def GetBetaRangeToAttempt(mass):
    # get a restricted beta range to test for each mass point
    minBetasPerMass = {300: 0, 400: 0, 500: 0.05, 600: 0.05, 700: 0.1, 800: 0.12, 900: 0.15, 1000: 0.2, 1100: 0.25, 1200: 0.3, 1300: 0.45, 1400: 0.55, 1500: 0.65}
    maxBetasPerMass = {300: 0.1, 400: 0.1, 500: 0.15, 600: 0.2, 700: 0.2, 800: 0.25, 900: 0.3, 1000: 0.35, 1100: 0.45, 1200: 0.55, 1300: 0.75, 1400: 1.0, 1500: 1.0}
    if mass > 1500:
        return 0.9, 1.0
    else:
        return minBetasPerMass[mass], maxBetasPerMass[mass]


def CheckErrorFile(errFileName, throwException=True):
    # print("\tChecking error file {}...".format(errFileName), end = "", flush=True)
    with open(errFileName, "r") as errFile:
        for line in errFile:
            if len(line):
                print(colored("Found unexpected content '{}' in {}".format(line, errFileName), "red"))
                if throwException:
                    raise RuntimeError("Found unexpected content '{}' in {}".format(line, errFileName))
    # print("OK")


def ReadBatchResults(massList, condorDir):
    for mass in massList:
        rLimitsByMassAndQuantile[mass] = {}
        xsecLimitsByMassAndQuantile[mass] = {}
        signalScaleFactorsByMassAndQuantile[mass] = {}
        cardFile = combinedDatacard if doShapeBasedLimits else dirName + "/datacards/tmpDatacard_m{}_card0_combCardFile.txt".format(mass)
        cardWorkspace = FindCardWorkspace(cardFile, mass)
        for quantileExp in quantilesExpected:
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantileStr = str(quantileExp).replace(".", "p")
            globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            errFileName = FindFile(globString)
            CheckErrorFile(errFileName, False)  # don't throw exception
            # now run HybridNew again to get limits from grid
            globString = condorDir+'/higgsCombine.signalScaleFactor*.HybridNew.mH{}.quant{}.grid.root'.format(mass, quantileExp)
            rootFileName = FindFile(globString)
            resultFileName = ComputeLimitsFromGrid(cardWorkspace, mass, condorDir, rootFileName, quantileExp)
            # now extract limits from the root file produced above
            # globString = condorDir+'/higgsCombineTest.HybridNew.mH{}.quant{}.root'.format(mass, quantileExp)
            # rootFileName = FindFile(globString)
            # limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
            limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(resultFileName)
            rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
            xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
            signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = signalScaleFactor
    return xsecLimitsByMassAndQuantile, signalScaleFactorsByMassAndQuantile


def ReadBatchResultsBetaScan(massList, condorDir):
    rLimitsByMassAndQuantile = {}
    xsecLimitsByMassAndQuantile = {}
    signalScaleFactorsByMassAndQuantile = {}
    for quantileExp in quantilesExpected:
        rLimitsByMassAndQuantile[str(quantileExp)] = {}
        xsecLimitsByMassAndQuantile[str(quantileExp)] = {}
        signalScaleFactorsByMassAndQuantile[str(quantileExp)] = {}
        for mass in massList:
            quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
            quantileStr = str(quantileExp).replace(".", "p")
            globString = condorDir+'/error/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
            errorFiles = GetFileList(globString)
            for errFile in errorFiles:
                CheckErrorFile(errFile)
                lastPart = errFile.split("{}.".format(quantileStr))[-1]
                lastPartCut = lastPart.rstrip(".0.err")
                sigScaleFact = lastPartCut[0:lastPartCut.rfind(".")].strip("signalScaleFactor")
                sigScaleFactRound = round(float(sigScaleFact), 6)
                rootGlobString = condorDir+'/higgsCombine.signalScaleFactor{}.HybridNew.mH{}.*.{}.root'.format(sigScaleFactRound, mass, quantile)
                rootFileName = FindFile(rootGlobString)
                limit, limitErr, quantileFromFile, signalScaleFactor = ExtractLimitResult(rootFileName)
                betaVal = math.sqrt(signalScaleFactor)
                betaValRound = round(betaVal, 6)
                if not betaValRound in rLimitsByMassAndQuantile[str(quantileExp)].keys():
                    rLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                    signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound] = {}
                rLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit
                xsecLimitsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = limit * xsThByMass[float(mass)]
                signalScaleFactorsByMassAndQuantile[str(quantileExp)][betaValRound][mass] = signalScaleFactor
                # print("\tINFO: ReadBatchResultsBetaScan() - fill xsecLimitsByMassAndQuantile[{}][{}][{}]".format(quantileExp, betaValRound, mass))
    return xsecLimitsByMassAndQuantile


def sigmaval(spline, mval):
    return spline.Eval(mval)
    
 
def mval(spline, sigma):
    testm = 150
    oldtestm = 1500
    inc = 50
    dif = 55
    olddif = 0
    while abs(oldtestm - testm) > 0.01:
        testsigma = sigmaval(spline, testm)
        olddif = dif
        dif = testsigma - sigma
        if testm > 1500:
            break
        if dif*olddif <= 0.0:
            inc = -inc/2.3
        oldtestm = testm
        #print '**' + str(testm) + '  ' + str(testsigma) +'  ' +str(dif) + '   ' + str(dif*olddif)
        testm = testm + inc
    return testm


def loggraph(inputarrayX, inputarrayY):
    logarray = []
    for j in inputarrayY:
        logarray.append(math.log10(j))
        # logarray.append(math.log(j))
    x = np.array(inputarrayX, dtype="f")
    y = np.array(logarray, dtype="f")
    # print("\tINFO: loggraph() - x={}, y={}".format(x, y))
    g = TGraph(len(x), x, y)
    return g


def logspline(inputarrayX,inputarrayY):
    logarray = []
    for j in inputarrayY:
        logarray.append(math.log(j))
    x = array("d",inputarrayX)
    y = array("d",logarray)
    g = TGraph(len(x),x,y)
    outspline = TSpline3("",g)
    return outspline


def get_simple_intersection(graph1, graph2, xmin, xmax):
    num = (xmax-xmin)*10
    inc = (xmax - xmin)/(1.0*num)
    dif = []
    sdif = []
    x = xmin +0.1
    xvals = []
    xx = []
    yy = []
    xvals = []
    while x < (xmax-.1):
        # print("\tINFO: get_simple_intersection() - x=", str(x) + '   xmax-.1='+ str(xmax-.1))
        # print("\tINFO: get_simple_intersection() - graph1.Eval(x) =", graph1.Eval(x))
        # print("\tINFO: get_simple_intersection() - math.exp(graph1.Eval(x)) =", math.exp(graph1.Eval(x)))
        # print("\tINFO: get_simple_intersection() - graph2.Eval(x) =", graph2.Eval(x))
        # print("\tINFO: get_simple_intersection() - math.exp(graph2.Eval(x)) =", math.exp(graph2.Eval(x)))
        # print("\tINFO: get_simple_intersection() - math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)) =", math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)))
        thisdif = (math.exp(graph1.Eval(x)) - math.exp(graph2.Eval(x)))
        xx.append(math.exp(graph1.Eval(x)))
        yy.append(math.exp(graph2.Eval(x)))
        sdif.append(thisdif)
        dif.append(abs(thisdif))
        xvals.append(x)
        # print("\tINFO: get_simple_intersection() [1] -", str(x) + '   ' +str(math.exp(graph1.Eval(x))) + '    '+str(math.exp(graph2.Eval(x))) + '    ' + str(thisdif))
        x = x+inc
	#print 'Done Looping for Difs'
    mindif = min(dif)
    bestmass = 0	
    
    for x in range(len(dif)-2):
       a = sdif[x]
       b = sdif[x+1]	
       # print("\tINFO: get_simple_intersection() [2] -", str(xvals[x+1]) +'    '+str(a)  + '     ' +str(b))
       if ((a/abs(a))*(b/abs(b))) < 0.0 and a > 0.0 :
           # print('\tINFO: get_simple_intersection() - Limit found at: '+ (str(xvals[x])))
           bestmass = xvals[x]
           break
    return [bestmass, mindif]


def ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile):
    mTh = np.array(list(xsThByMass.keys()), dtype="f")
    xsTh = np.array(list(xsThByMass.items()), dtype="f")
    g = TGraph(len(mTh), mTh, xsTh);
    # spline = TSpline3("xsection", g)
    print("INFO: ComputeBetaLimits() - make theory graph")
    logtheory = loggraph(list(xsThByMass.keys()), list(xsThByMass.values()))  # make graph with log10 of y values
    massLimitsByQuantileAndBetaVal = {}
    for quantile in xsecLimitsByMassAndQuantile.keys():
        print("INFO: ComputeBetaLimits() - examine quantile={}".format(quantile))
        massLimitsByQuantileAndBetaVal[quantile] = {}
        for betaVal in xsecLimitsByMassAndQuantile[quantile].keys():
            print("\tINFO: examine betaVal={}".format(betaVal))
            limit_set = list(xsecLimitsByMassAndQuantile[quantile][betaVal].values())
            massList = list(xsecLimitsByMassAndQuantile[quantile][betaVal].keys())
            if len(massList) < 2:
                print("\tWARN: Skipping beta value={} as we only have one mass point tested here! Need to adjust beta scan range.".format(betaVal))
                continue
            print("\tINFO: examine massList={}, limit_set={}".format(massList, limit_set))
            # print("\tINFO: ComputeBetaLimits() - make loggraph for fitted_limits")
            fitted_limits = loggraph(massList, limit_set)
            #FIXME: get_simple_intersection seems to assume both graphs are in log base e rather than log10; the spline (really, graph?) is log10, though.
            # print("\tINFO: ComputeBetaLimits() - get_simple_intersection")
            goodm = get_simple_intersection(logtheory, fitted_limits, min(mTh), max(mTh))
            massLimitsByQuantileAndBetaVal[quantile][betaVal] = round(goodm[0], 3)
            print("\tINFO: for betaVal={}: bestmass, mindif={}".format(betaVal, goodm))
            if goodm[0] < min(massList) or goodm[0] > max(massList):
                print("\tWARN: For beta value={}, intersection mass is outside of massList range {}! Need to adjust beta scan range.".format(betaVal, massList))
    return massLimitsByQuantileAndBetaVal



def CreateComboArrays(xsecLimitsByMassAndQuantile):
    retVal = {}
    for quantile in xsecLimitsByMassAndQuantile.keys():
        retVal[quantile] = {}
        retVal[quantile]["betas"] = []
        retVal[quantile]["massLimits"] = []
        sortedBetas = list(sorted(xsecLimitsByMassAndQuantile[quantile].keys()))
        # for betaVal, mass in xsecLimitsByMassAndQuantile[quantile].items():
        for betaVal in sortedBetas:
            retVal[quantile]["betas"].append(betaVal)
            retVal[quantile]["massLimits"].append(xsecLimitsByMassAndQuantile[quantile][betaVal])
    return retVal


def MakeResultTable(masses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp):
    doObs = False
    if len(xsObs):
        doObs = True
    table = []
    for idx, mass in enumerate(list(masses)):
        tableLine = [int(mass), xsMedExp[idx], xsOneSigmaExp[idx], xsTwoSigmaExp[idx]]
        if doObs:
            tableLine.append(xsObs[idx])
        table.append(tableLine)
    columnNames = ["MLQ", "xs exp [pb]", "xs 1sigma exp [pb]", "xs 2sigma exp [pb]"]
    if doObs:
        columnNames.append("xs obs [pb]")
    return table, columnNames


####################################################################################################
# Run
####################################################################################################
if __name__ == "__main__":
    cmsswDir = "/afs/cern.ch/user/s/scooper/work/private/cmssw/11_3_4/LQLimitSetting/src"
    kerberosKeytab = "/afs/cern.ch/user/s/scooper/scooper.keytab"
    doBatch = True
    doShapeBasedLimits = False
    doAsymptoticLimits = False
    blinded = True
    # massList = list(range(300, 3100, 100))
    massList = list(range(300, 1900, 100))
    # massList = list(range(1000, 2000, 100))
    betasToScan = list(np.linspace(0.002, 1, 500))[:-1] + [0.9995]
    
    quantilesExpected = [0.025, 0.16, 0.5, 0.84, 0.975]
    xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    if doBatch:
        cmsswPreambleCmd = "source /etc/profile && k5start -f {}".format(kerberosKeytab) + " && cd {}".format(cmsswDir) + " && eval $(/cvmfs/cms.cern.ch/common/scram ru -sh) && cd - && "
        cmsswEnv = ExtractCMSSWEnv(cmsswPreambleCmd)
    commonCombineArgs = " --setParameters signalScaleParam={} --freezeParameters signalScaleParam --trackParameters signalScaleParam"
    
    parser = OptionParser(
        usage="%prog datacard",
    )
    parser.add_option(
        "-d",
        "--datacard",
        dest="datacard",
        help="combined datacard",
        metavar="datacard",
    )
    parser.add_option(
        "-r",
        "--readResults",
        action="store_true",
        help="read limit results from batch",
        metavar="readResults",
        default=False,
    )
    parser.add_option(
        "-n",
        "--name",
        dest="name",
        help="name of limit calculation (for bookkeeping)",
        metavar="name",
    )
    parser.add_option(
        "-i",
        "--doImpacts",
        dest="doImpacts",
        help="produce impact plots",
        metavar="doImpacts",
        action="store_true",
        default=False,
    )
    parser.add_option(
        "-b",
        "--doBetaScan",
        dest="doBetaScan",
        help="do scan of beta values",
        metavar="doBetaScan",
        action="store_true",
        default=False,
    )
    (options, args) = parser.parse_args()
    if options.datacard is None and options.readResults is None:
        raise RuntimeError("Need either option -d to specify datacard, or option r to specify reading limit results from batch")
    if options.name is None:
        raise RuntimeError("Option -n to specify name of limit results dir is required")
    if options.doBetaScan and not doBatch:
        raise RuntimeError("Won't do beta scan without batch submission enabled (see doBatch parameter inside script).")
    
    manager = multiprocessing.Manager()
    dictAsimovToysByScaleFactor = manager.dict()
    listFailedCommands = manager.list()
    dirName = options.name
    xsThByMass, yPDFUpByMass, yPDFDownByMass = ReadXSecFile(xsThFilename)
    rLimitsByMassAndQuantile = {}
    xsecLimitsByMassAndQuantile = {}
    signalScaleFactorsByMassAndQuantile = {}
    failedBatchCommands = []
    combinedDatacard = options.datacard
    if not options.readResults:
        startTime = time.time()
        #TODO: check for previous datacards?
        separateDatacardsDir = dirName+"/datacards"
        asimovToysDir = dirName+"/asimovData"
        if not os.path.isdir(dirName):
            print("INFO: Making directory", dirName, flush=True)
            Path(dirName).mkdir(exist_ok=True)
        if not os.path.isdir(separateDatacardsDir):
            print("INFO: Making directory", separateDatacardsDir, flush=True)
            Path(separateDatacardsDir).mkdir(exist_ok=True)
        if not os.path.isdir(asimovToysDir):
            print("INFO: Making directory", asimovToysDir, flush=True)
            Path(asimovToysDir).mkdir(exist_ok=True)
        
        if not doShapeBasedLimits:
            massListFromCards, cardFilesByMass = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
        for mass in massList:
            cardFile = combinedDatacard if doShapeBasedLimits else cardFilesByMass[mass]
            print("INFO: Computing limits for mass {}".format(mass), flush=True)
            rLimitsByMassAndQuantile[mass] = {}
            xsecLimitsByMassAndQuantile[mass] = {}
            signalScaleFactorsByMassAndQuantile[mass] = {}
            cardWorkspace = FindCardWorkspace(cardFile, mass)
            if cardWorkspace is not None:
                print("INFO: Using previously-generated card workspace: {}".format(cardWorkspace), flush=True)
            else:
                cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)
            if doAsymptoticLimits:
                print("INFO: Doing interactive AsymptoticLimits for mass {}".format(mass), flush=True)
                asStartTime = time.time()
                rootFileName = RunAsymptotic(cardWorkspace, mass, dirName, blinded, 1.0)
                #TODO: should we also run again, scaling the signal yields, as per below? it might not matter.
                asStopTime = time.time()
                execTimeStr = GetExecTimeStr(asStartTime, asStopTime)
                print("Asymptotic calculation execution time:", execTimeStr, flush=True)
                limits, limitErrs, quantiles = ExtractAsymptoticLimitResult(rootFileName)
                for index, quantileExp in enumerate(quantiles):
                    limit = limits[index]
                    rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                    xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
            else:
                print("INFO: Doing interactive AsymptoticLimits for mass {} to estimate scan range".format(mass), flush=True)
                asymptoticRootFileName = RunAsymptotic(cardWorkspace, mass, dirName, blinded, 1.0)
                limits, limitErrs, quantiles = ExtractAsymptoticLimitResult(asymptoticRootFileName)
                rValuesByQuantile = dict(zip(quantiles, limits))
                for quantileExp in quantilesExpected:
                    signalScaleFactor = 1.0
                    asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysDir)
                    if asimovToyFile is not None:
                        print("INFO: Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                    else:
                        asimovToyFile = GenerateAsimovToyData((cardWorkspace, mass, asimovToysDir, dictAsimovToysByScaleFactor, signalScaleFactor))
                    if doBatch:
                        MakeCondorDirs(dirName)
                        if quantileExp == 0.025 and mass > 800:  # adjust scan range upwards for lowest quantile and higher masses
                            rMax = rValuesByQuantile[quantileExp]*1.8
                            rMin = rValuesByQuantile[quantileExp]*0.8
                        else:
                            rMax = rValuesByQuantile[quantileExp]*1.3
                            rMin = rValuesByQuantile[quantileExp]*0.75
                        failedCmd = SubmitHybridNewBatch((cardWorkspace, mass, dirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor, rMin, rMax))
                        if failedCmd is not None:
                            failedBatchCommands.append(failedCmd)
                        # beta scan jobs
                        if options.doBetaScan:
                            betaDirName = dirName+"/betaScan"
                            if not os.path.isdir(betaDirName):
                                print("INFO: Making directory", betaDirName, flush=True)
                                Path(betaDirName).mkdir(exist_ok=True)
                            MakeCondorDirs(betaDirName)
                            minBeta, maxBeta = GetBetaRangeToAttempt(int(mass))
                            minBetaPosition = bisect(betasToScan, minBeta)
                            maxBetaPosition = bisect(betasToScan, maxBeta)
                            print("INFO: for mass {}, scan beta range {} to {} or indices {} to {}".format(mass, minBeta, maxBeta, minBetaPosition, maxBetaPosition))
                            betasToSubmit = betasToScan[minBetaPosition : maxBetaPosition]
                            ncores = 8
                            with multiprocessing.Pool(ncores) as pool:
                                asimovJobs = 0
                                for beta in betasToSubmit:
                                    signalScaleFactor = beta*beta
                                    asimovToysBetaDir = betaDirName+"/asimovToys"
                                    if not os.path.isdir(asimovToysBetaDir):
                                        print("INFO: beta scan - Making directory", asimovToysBetaDir, flush=True)
                                        Path(asimovToysBetaDir).mkdir(exist_ok=True)
                                    asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                                    if asimovToyFile is not None:
                                        # print("INFO: beta scan - Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                                        pass
                                    else:
                                        try:
                                            pool.apply_async(GenerateAsimovToyData, [[cardWorkspace, mass, asimovToysBetaDir, dictAsimovToysByScaleFactor, signalScaleFactor]])
                                            asimovJobs += 1
                                        except KeyboardInterrupt:
                                            print("\n\nCtrl-C detected: Bailing.")
                                            pool.terminate()
                                            sys.exit(1)
                                        except Exception as e:
                                            print("ERROR: caught exception in asimov toy job for scale factor: {}; exiting".format(signalScaleFactor))
                                            traceback.print_exc()
                                            pool.terminate()
                                            exit(-2)
                                # now close the pool and wait for jobs to finish
                                pool.close()
                                if asimovJobs > 0:
                                    print("Waiting for {} Asimov toy generation jobs to finish for all beta values for mass={}, quantile={}...".format(asimovJobs, mass, quantileExp))
                                pool.join()
                            # HybridNew
                            with multiprocessing.Pool(ncores) as pool:
                                hybridNewJobs = 0
                                for beta in betasToSubmit:
                                    signalScaleFactor = beta*beta
                                    try:
                                        asimovToyFile = FindAsimovToyData(mass, signalScaleFactor, asimovToysBetaDir)
                                        pool.apply_async(SubmitHybridNewBatch, [[cardWorkspace, mass, betaDirName, listFailedCommands, quantileExp, asimovToyFile, signalScaleFactor]])
                                        hybridNewJobs += 1
                                    except KeyboardInterrupt:
                                        print("\n\nCtrl-C detected: Bailing.")
                                        pool.terminate()
                                        sys.exit(1)
                                    except Exception as e:
                                        print("ERROR: caught exception in hybridnew job for scale factor: {}; exiting".format(signalScaleFactor))
                                        traceback.print_exc()
                                        pool.terminate()
                                        exit(-2)
                                # now close the pool
                                pool.close()
                                print("Waiting for submission of {} HybridNew jobs for mass={}, quantile={}...".format(hybridNewJobs, mass, quantileExp))
                                pool.join()
                                # if failedCmd is not None:
                                #     failedBatchCommands.append(failedCmd)
                    else:
                        hnStartTime = time.time()
                        rootFileName = RunHybridNewInteractive(cardWorkspace, mass, dirName, quantile=quantileExp, genAsimovToyFile=asimovToyFile, signalScaleFactor=signalScaleFactor)
                        hnStopTime = time.time()
                        execTimeStr = GetExecTimeStr(hnStartTime, hnStopTime)
                        print("HybridNew calculation execution time:", execTimeStr, flush=True)
                        limit, limitErr, quantile, signalScaleFactor = ExtractLimitResult(rootFileName)
                        rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                        xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
                        signalScaleFactorsByMassAndQuantile[mass][str(quantileExp)] = signalScaleFactor
        
    if len(failedBatchCommands):
        print(colored("batch commands failed:", "red"))
        for item in failedBatchCommands:
            print(item)
    if len(listFailedCommands):
        print(colored("batch commands failed:", "red"))
        for item in listFailedCommands:
            print(item)

    if not doBatch or options.readResults:
        if options.readResults:
            condorDir = dirName.rstrip("/")+"/condor"
            print("INFO: Reading results from batch from {}...".format(condorDir), flush=True)
            xsecLimitsByMassAndQuantile, signalScaleFactorsByMassAndQuantile = ReadBatchResults(massList, condorDir)
            print("DONE", flush=True)
        masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs = CreateArraysForPlotting(xsecLimitsByMassAndQuantile)
        # print("mData =", list(masses))
        # print("x_shademasses =", list(shadeMasses))
        # print("xsUp_expected =", list(xsMedExp))
        # print("xsUp_observed =", list(xsObs))
        # print("y_1sigma =", list(xsOneSigmaExp))
        # print("y_2sigma =", list(xsTwoSigmaExp))
        table, columnNames = MakeResultTable(masses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
        print(tabulate(table, headers=columnNames, tablefmt="fancy_grid", floatfmt=".6f"))
        stopTime = time.time()
        if not doBatch:
            execTimeStr = GetExecTimeStr(startTime, stopTime)
            print("Total limit calculation execution time:", execTimeStr)
        print("Make plot and calculate mass limit")
        BR_Sigma_EE_vsMass(dirName, masses, shadeMasses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
        if options.doBetaScan:
            mainDirName = dirName
            betaDirName = dirName+"/betaScan"
            condorDir = betaDirName.strip("/")+"/condor"
            print("INFO: Reading results from batch from {}...".format(condorDir), flush=True, end="")
            xsecLimitsByMassAndQuantile = ReadBatchResultsBetaScan(massList, condorDir)
            print("DONE", flush=True)
            # print("xsecLimitsByMassAndQuantile={}".format(xsecLimitsByMassAndQuantile))
            print("INFO: Compute beta limits...", flush=True, end="")
            massLimitsByQuantileAndBetaVal = ComputeBetaLimits(xsThByMass, xsecLimitsByMassAndQuantile)
            print("massLimitsByQuantileAndBetaVal={}".format(massLimitsByQuantileAndBetaVal))
            print("DONE", flush=True)
            comboArraysByQuantile = CreateComboArrays(massLimitsByQuantileAndBetaVal)
            print("comboArraysByQuantile['0.5']['betas']={}".format(comboArraysByQuantile["0.5"]["betas"]))
            print("comboArraysByQuantile['0.5']['massLimits']={}".format(comboArraysByQuantile["0.5"]["massLimits"]))
            ComboPlot(mainDirName, comboArraysByQuantile["0.5"]["betas"], comboArraysByQuantile["0.5"]["massLimits"]) # , m_observed_lljj)
        if options.doImpacts:
            print("INFO: Making nuisance parameter impact plots...", flush=True, end="")
            for mass in massList:
                datacardDir = dirName.strip("/")+"/datacards"
                cardWorkspace = FindCardWorkspace(datacardDir + "/*", mass, True)
                MakeImpacts(cardWorkspace, mass, dirName, signalScaleFactorsByMassAndQuantile[mass][str(quantilesExpected[0])])
            print("DONE", flush=True)
