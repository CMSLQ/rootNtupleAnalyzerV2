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
from pathlib import Path
from ROOT import TFile
from combineCommon import SeparateDatacards
from BR_Sigma_EE_vsMass import BR_Sigma_EE_vsMass

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


def RunCommand(cmd, workingDir=None, cleanEnv=False, shell=False):
    print(colored("\t{}".format(cmd), "green"), flush=True)
    try:
        if not shell:
            cmd = shlex.split(cmd)
        if cleanEnv:
            process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, env={}, shell=shell)
        else:
            process = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=workingDir, shell=shell)
    except subprocess.CalledProcessError as e:
        result = ansi_escape.sub('', e.output.decode())
        print(colored("RunCommand had an error; output: {}".format(result), "red"), flush=True)
        raise e
    print(process.stdout.decode(), flush=True)
    #print(process.stderr.decode())


def ConvertDatacardToWorkspace(datacard, mass):
    workspaceFileName = datacard+".m{}.root".format(mass)
    if doBatch:
        cmd = cmsswPreambleCmd
        cmd += 'text2workspace.py {} -o {} -m {}'.format(datacard, workspaceFileName, mass)
        RunCommand(cmd, None, True, True)
    else:
        cmd = 'text2workspace.py {} -o {} -m {}'.format(datacard, workspaceFileName, mass)
        RunCommand(cmd)
    return Path(workspaceFileName).resolve()


def FindCardWorkspace(cardFile, mass):
    workspaceFileName = cardFile+".m{}.root".format(mass)
    listOfWorkspaceFiles = sorted(glob.glob(workspaceFileName), key=os.path.getmtime)
    if len(listOfWorkspaceFiles) > 0:
        return Path(listOfWorkspaceFiles[-1]).resolve()


def GenerateAsimovToyData(workspace, mass, toysDir):
    if doBatch:
        cmd = cmsswPreambleCmd
        cmd += 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov -m {}'.format(workspace, mass)
        RunCommand(cmd, toysDir, True, True)
    else:
        cmd = 'combine -M GenerateOnly {} -t -1 --seed -1 --saveToys -n .asimov -m {}'.format(workspace, mass)
        RunCommand(cmd, toysDir)
    return Path(sorted(glob.glob(toysDir+'/higgsCombine.asimov.GenerateOnly.mH{}.*.root'.format(mass)), key=os.path.getmtime)[-1]).resolve()


def FindAsimovToyData(mass, toysDir):
    listOfToyFiles = sorted(glob.glob(toysDir+'/higgsCombine.asimov.GenerateOnly.mH{}.*.root'.format(mass)), key=os.path.getmtime)
    if len(listOfToyFiles) > 0:
        return Path(listOfToyFiles[-1]).resolve()


def FindFile(globString):
    fileList = sorted(glob.glob(globString), key=os.path.getmtime)
    if len(fileList) != 1:
        raise RuntimeError("Globbing for {} did not result in one file as expected, but {}".format(globString, fileList))
    return Path(fileList[-1]).resolve()


def GetHybridNewCommandArgs(workspace, mass, dirName, quantile=-1, genAsimovToyFile=""):
    rAbsAcc = 0.0001
    toys = 4000
    # toys = 500  # reduced for shape-based limits as test, but still took hours
    cmd = '{}'.format(workspace)
    cmd += ' -M HybridNew'
    cmd += ' --LHCmode LHC-limits'
    cmd += ' --saveHybridResult'
    cmd += ' --seed -1'
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    cmd += ' -T {}'.format(toys)
    cmd += ' -m {}'.format(mass)
    cmd += ' -H AsymptoticLimits'
    cmd += ' --fork 4'
    if quantile > 0:
        cmd += ' --expectedFromGrid {}'.format(quantile)
        if genAsimovToyFile != "":
            cmd += ' -D {}:toys/toy_asimov'.format(genAsimovToyFile)
    return cmd


def SubmitHybridNewBatch(workspace, mass, dirName, quantile, genAsimovToyFile):
    condorDir = dirName.strip("/")+"/condor"
    if not os.path.isdir(condorDir):
        print("INFO: Making directory", condorDir, flush=True)
        Path(condorDir).mkdir(exist_ok=True)
    if cmsswDir is None:
        raise RuntimeError("Need to specify valid CMSSW area with combine checked out.")
    cwd = os.getcwd()
    cmdArgs = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile)
    quantileStr = str(quantile).replace(".", "p")
    cmd = cmsswPreambleCmd
    cmd += "cd {} && combineTool.py ".format(cwd + "/" + condorDir)
    cmd += cmdArgs
    cmd += ' --job-mode condor'
    cmd += ' --task-name hybridNewLimits.M{}.{}'.format(mass, quantileStr)
    cmd += ' --sub-opts=\'+JobFlavour="workday"\''
    cmd += ' --dry-run'
    RunCommand(cmd, dirName, True, True)
    # now submit
    condorFile = "condor_hybridNewLimits.M{}.{}.sub".format(mass, quantileStr)
    # cmd = "source /etc/profile && k5start -f {} && cd {}".format(kerberosKeytab, cmsswDir) + " && eval $(/cvmfs/cms.cern.ch/common/scram ru -sh) && cd - && condor_submit {}".format(condorFile)
    cmd = "source /etc/profile && k5start -f {}".format(kerberosKeytab) + " && condor_submit {}".format(condorFile)
    RunCommand(cmd, cwd + "/" + condorDir, True, True)
    # cmd = "condor_submit {}".format(condorFile)
    # RunCommand(cmd, cwd + "/" + condorDir)


def RunHybridNewInteractive(workspace, mass, dirName, quantile=-1, genAsimovToyFile=""):
    cmd = GetHybridNewCommandArgs(workspace, mass, dirName, quantile, genAsimovToyFile)
    cmd = 'combine -v1 ' + cmd
    RunCommand(cmd, dirName)
    quantileString = "quant{:.3f}".format(quantile) if quantile > 0 else "."
    globString = dirName+'/higgsCombineTest.HybridNew.mH{}.*.{}.root'.format(mass, quantileString)
    return FindFile(globString)


def RunAsymptotic(workspace, mass, dirName, blinded=True):
    rAbsAcc = 0.0001
    cmd = 'combine'
    # cmd += ' -v1'
    cmd += ' -M AsymptoticLimits {}'.format(workspace)
    cmd += ' --seed -1'
    cmd += ' -m {}'.format(mass)
    cmd += ' --rAbsAcc {}'.format(rAbsAcc)
    if blinded:
        cmd += ' --run blind --expectSignal 0'
    RunCommand(cmd, dirName)
    return Path(sorted(glob.glob(dirName+'/higgsCombineTest.AsymptoticLimits.mH{}.*.root'.format(mass)), key=os.path.getmtime)[-1]).resolve()


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
    tfile.Close()
    return limit, limitErr, quantile


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
    massList = list(range(1000, 2100, 100))
    
    quantilesExpected = [0.025, 0.16, 0.5, 0.84, 0.975]
    xsThFilename = "$LQANA/config/xsection_theory_13TeV_scalarPairLQ.txt"
    cmsswPreambleCmd = "source /etc/profile && cd {}".format(cmsswDir) + " && eval $(/cvmfs/cms.cern.ch/common/scram ru -sh) && cd - && "
    
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
        default=False,
        help="read limit results from batch",
        metavar="readResults",
    )
    parser.add_option(
        "-n",
        "--name",
        dest="name",
        help="name of limit calculation (for bookkeeping)",
        metavar="name",
    )
    (options, args) = parser.parse_args()
    if options.datacard is None and options.readResults is None:
        raise RuntimeError("Need either option -d to specify datacard, or option r to specify reading limit results from batch")
    if options.name is None:
        raise RuntimeError("Option -n to specify name of limit results dir is required")
    
    dirName = options.name
    xsThByMass, yPDFUpByMass, yPDFDownByMass = ReadXSecFile(xsThFilename)
    rLimitsByMassAndQuantile = {}
    xsecLimitsByMassAndQuantile = {}
    if not options.readResults:
        startTime = time.time()
        combinedDatacard = options.datacard
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
            massList, cardFilesByMass = SeparateDatacards(combinedDatacard, 0, separateDatacardsDir)
        for mass in massList:
            cardFile = combinedDatacard if doShapeBasedLimits else cardFilesByMass[mass]
            print("INFO: Computing limits for mass {}".format(mass), flush=True)
            rLimitsByMassAndQuantile[mass] = {}
            xsecLimitsByMassAndQuantile[mass] = {}
            cardWorkspace = FindCardWorkspace(cardFile, mass)
            if cardWorkspace is not None:
                print("INFO: Using previously-generated card workspace: {}".format(cardWorkspace), flush=True)
            else:
                cardWorkspace = ConvertDatacardToWorkspace(cardFile, mass)
            if doAsymptoticLimits:
                print("INFO: Doing interactive AsymptoticLimits for mass {}".format(mass), flush=True)
                asStartTime = time.time()
                rootFileName = RunAsymptotic(cardWorkspace, mass, dirName, blinded)
                asStopTime = time.time()
                execTimeStr = GetExecTimeStr(asStartTime, asStopTime)
                print("Asymptotic calculation execution time:", execTimeStr, flush=True)
                limits, limitErrs, quantiles = ExtractAsymptoticLimitResult(rootFileName)
                for index, quantileExp in enumerate(quantiles):
                    limit = limits[index]
                    rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                    xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
            else:
                asimovToyFile = FindAsimovToyData(mass, asimovToysDir)
                if asimovToyFile is not None:
                    print("INFO: Using previously-generated Asimov toy file: {}".format(asimovToyFile), flush=True)
                else:
                    asimovToyFile = GenerateAsimovToyData(cardWorkspace, mass, asimovToysDir)
                for quantileExp in quantilesExpected:
                    hnStartTime = time.time()
                    if doBatch:
                        SubmitHybridNewBatch(cardWorkspace, mass, dirName, quantile=quantileExp, genAsimovToyFile=asimovToyFile)
                    else:
                        rootFileName = RunHybridNewInteractive(cardWorkspace, mass, dirName, quantile=quantileExp, genAsimovToyFile=asimovToyFile)
                        hnStopTime = time.time()
                        execTimeStr = GetExecTimeStr(hnStartTime, hnStopTime)
                        print("HybridNew calculation execution time:", execTimeStr, flush=True)
                        limit, limitErr, quantile = ExtractLimitResult(rootFileName)
                        rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                        xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
        
    if not doBatch or options.readResults:
        if options.readResults:
            condorDir = dirName.strip("/")+"/condor"
            print("INFO: Reading results from batch from {}...".format(condorDir), flush=True)
            for mass in massList:
                rLimitsByMassAndQuantile[mass] = {}
                xsecLimitsByMassAndQuantile[mass] = {}
                for quantileExp in quantilesExpected:
                    quantile = "quant{:.3f}".format(quantileExp) if quantileExp > 0 else "."
                    # quantile = quantileString.replace("quant0.", "0p")
                    quantileStr = str(quantileExp).replace(".", "p")
                    globString = condorDir+'/hybridNewLimits.M{}.{}.*.0.err'.format(mass, quantileStr)
                    errFileName = FindFile(globString)
                    # print("\tChecking error file {}...".format(errFileName), end = "", flush=True)
                    with open(errFileName, "r") as errFile:
                        for line in errFile:
                            if len(line):
                                raise RuntimeError("Found unexpected content '{}' in {}".format(line, errFileName))
                    # print("OK")
                    globString = condorDir+'/higgsCombine.Test.HybridNew.mH{}.*.{}.root'.format(mass, quantile)
                    rootFileName = FindFile(globString)
                    limit, limitErr, quantile = ExtractLimitResult(rootFileName)
                    rLimitsByMassAndQuantile[mass][str(quantileExp)] = limit
                    xsecLimitsByMassAndQuantile[mass][str(quantileExp)] = limit * xsThByMass[float(mass)]
            print("DONE", flush=True)
        masses, shadeMasses, xsMedExp, xsOneSigmaExp, xsTwoSigmaExp, xsObs = CreateArraysForPlotting(xsecLimitsByMassAndQuantile)
        print("mData =", list(masses))
        print("x_shademasses =", list(shadeMasses))
        print("xsUp_expected =", list(xsMedExp))
        print("xsUp_observed =", list(xsObs))
        print("y_1sigma =", list(xsOneSigmaExp))
        print("y_2sigma =", list(xsTwoSigmaExp))
        stopTime = time.time()
        if not doBatch:
            execTimeStr = GetExecTimeStr(startTime, stopTime)
            print("Total limit calculation execution time:", execTimeStr)
        print("Make plot and calculate mass limit")
        BR_Sigma_EE_vsMass(dirName, masses, shadeMasses, xsMedExp, xsObs, xsOneSigmaExp, xsTwoSigmaExp)
