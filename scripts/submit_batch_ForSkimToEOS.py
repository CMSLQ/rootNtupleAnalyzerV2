#! /usr/bin/env python3

import os
import sys
import string
import re
from optparse import OptionParser
import subprocess
import shlex
from pathlib import Path


def GetSiteListFromDAS(fileNameList, dataset):
    # the best we can easily do is to just submit to all sites which host at least part of the dataset
    firstName = fileNameList[0]
    query = "dataset file=/"+firstName.split("//")[-1].strip()
    cmd = 'dasgoclient --query="{}"'.format(query)
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    output = proc.communicate()[0].decode().strip()
    query = "site dataset="+output
    cmd = 'dasgoclient --query="{}"'.format(query)
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    output = proc.communicate()[0].decode().strip()
    siteList = list(output.replace("_Disk", "").split())
    return siteList


def GetSitesHostingAllFiles(dataset):
    query = "site dataset={} | grep site.name, site.block_completion".format(dataset)
    cmd = 'dasgoclient --query="{}"'.format(query)
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    output = proc.communicate()[0].decode().strip()
    siteList = list(output.replace("_Disk", "").split("\n"))
    sitesHostingAllFiles = [line.split()[0].strip() for line in siteList if "100.00%" in line]
    return sitesHostingAllFiles


def GetDefaultCondorSites():
    # cmd = 'get_condor_sites --default'
    # proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE)
    # output = proc.communicate()[0].decode().strip()
    # siteList = list(output.split())
    siteList = ""
    with open("/etc/ciconnect/config.ini", "r") as siteConfFile:
        for line in siteConfFile:
            if "DefaultSites" in line:
                siteList = line.split("=")[-1].strip().split(",")
    return siteList
    

def PrepareJobScript(outputname):
    with open(outputname, "w") as outputfile:
        outputfile.write("#!/bin/bash\n")
        outputfile.write("echo \"Running at site: $GLIDEIN_CMSSite\"\n")
        # hardcoded root is a bit nasty FIXME
        outputfile.write("echo \"Now source LCG env.\"\n")
        outputfile.write('source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc13-opt/setup.sh\n')
        outputfile.write("echo \"-->Sourced LCG env.\"\n")
        # ROOT likes HOME set
        outputfile.write('[ -z "$HOME" ] && export HOME='+os.getenv('HOME')+'\n')
        outputfile.write("echo \"-->Set HOME=$HOME.\"\n")
        outputfile.write('export LD_LIBRARY_PATH=.:$LD_LIBRARY_PATH\n')
        outputfile.write("echo \"-->exported LD_LIBRARY_PATH.\"\n")
        inputList = inputfilename.split('/')[-1]
        # merge files if nano skim or reduced skim requested
        if options.reducedSkim or options.nanoSkim:
            with open(outputmain+"/input/"+inputList) as f:
                outFileNames = f.read().splitlines()
            outputfile.write("echo \"Running haddnano.py.\"\n")
            outputfile.write("./haddnano.py inputTree.root %s\n" % (" ".join(outFileNames)))
            # overwrite original inputList to just have inputTree.root
            with open(inputfilename, "w") as newInputList:
                newInputList.write("inputTree.root\n")
            outputfile.write('retVal=$?\n')
            outputfile.write('if [ $retVal -ne 0 ]; then\n')
            outputfile.write('  echo "./haddnano.py return error code=$retVal; quitting here."\n')
            outputfile.write('  exit $retVal\n')
            outputfile.write('fi\n')
        outputfile.write("echo \"Running executable.\"\n")
        outputfile.write('./'+execName+' '+inputList+" "+cutfile.split('/')[-1]+" "+options.treeName+" "+outputPrefix+"_"+str(ijob)+" "+outputPrefix+"_"+str(ijob)+"\n")
        outputfile.write('retVal=$?\n')
        outputfile.write('if [ $retVal -ne 0 ]; then\n')
        outputfile.write('  echo "./'+execName+' return error code=$retVal; quitting here."\n')
        outputfile.write('  exit $retVal\n')
        outputfile.write('fi\n')
        outputfile.write("xrdfs "+options.eosHost+" mkdir -p \""+outputcondoreosdir+"\"\n")
        outputfile.write("xrdcp -fs "+"\""+outputPrefix+"_"+str(ijob)+".root\" \""+options.eosHost+outputcondoreosdir+"/"+dataset+"_"+str(ijob)+".root\"\n")
        outputfile.write("xrdcp -fs "+"\""+outputPrefix+"_"+str(ijob)+".dat\" \""+options.eosHost+outputcondoreosdir+"/"+dataset+"_"+str(ijob)+".dat\"\n")
        if options.reducedSkim:
            outputfile.write("if [ -f "+outputPrefix+"_"+str(ijob)+"_reduced_skim.root ]; then\n")
            outputfile.write("    xrdfs "+options.eosHost+" mkdir -p \""+outputskimeosdir+"\"\n")
            outputfile.write("    xrdcp -fs "+"\""+outputPrefix+"_"+str(ijob)+"_reduced_skim.root\" \""+options.eosHost+outputskimeosdir+"/"+dataset+"_"+str(ijob)+"_rsk.root\"\n")
        else:
            # flat skim
            outputfile.write("if [ -f "+outputPrefix+"_"+str(ijob)+"_skim.root ]; then\n")
            outputfile.write("    xrdfs "+options.eosHost+" mkdir -p \""+outputskimeosdir+"\"\n")
            outputfile.write("    xrdcp -fs "+"\""+outputPrefix+"_"+str(ijob)+"_skim.root\" \""+options.eosHost+outputskimeosdir+"/"+dataset+"_"+str(ijob)+"_sk.root\"\n")
        outputfile.write("fi\n")


def WriteSubmitFile(condorFileName):
    with open(condorFileName, 'w') as condorFile:
        condorFile.write('executable  = '+outputmain+'/src/submit_$(Process).sh\n')
        condorFile.write('N = '+str(ijobmax)+'\n')
        condorFile.write('output      = '+outputmain+'/output/$(Process).out\n')
        condorFile.write('error       = '+outputmain+'/error/$(Process).err\n')
        condorFile.write('log         = '+outputmain+'/log/$(Process).log\n')
        # http://batchdocs.web.cern.ch/batchdocs/local/submit.html
        #  - cms connect shouldn't use JobFlavor or the requirements
        #  - assume this is lxbatch if queue option specified
        if options.queue is None:
            condorFile.write('+REQUIRED_OS = rhel9\n')
            condorFile.write('use_x509userproxy = true\n')
            outputRootFile = outputPrefix+"_$(Process).root"
            outputDatFile = outputPrefix+"_$(Process).dat"
            outputPath = outputmain+"/output/"
            condorFile.write('transfer_output_files = '+outputRootFile+','+outputDatFile+'\n')
            condorFile.write('transfer_output_remaps = "'+outputRootFile+' = '+outputPath+outputRootFile+'; '+outputDatFile+' = '+outputPath+outputDatFile+'"\n')
            condorFile.write("# Periodically retry the jobs every 10 minutes, up to a maximum of 5 retries.\n")
            condorFile.write("periodic_release =  (NumJobStarts < 3) && ((CurrentTime - EnteredCurrentStatus) > 600)\n")
        else:
            condorFile.write('+JobFlavour = "'+queue+'"\n')
            # require CentOS7
            condorFile.write('MY.WantOS = "el9"\n')
            condorFile.write('transfer_output_files = ""\n')
            if os.getenv("LQINPUTS") is not None:
                condorFile.write('environment = "LQINPUTS={}"\n'.format(os.getenv("LQINPUTS")))
        # make sure the job finishes with exit code 0
        # condorFile.write('on_exit_remove = (ExitBySignal == False) && (ExitCode == 0)\n')
        condorFile.write('max_retries = 3\n')
        condorFile.write('should_transfer_files = YES\n')
        # condorFile.write('stream_output = True\n')
        # condorFile.write('stream_error = True\n')
        # exePath = os.path.dirname(os.path.abspath(options.executable))
        inputFilesToTransfer = [cutfile,options.executable,outputmain+'/input/input_$(Process).list']
        filesToTransfer = ",".join(inputFilesToTransfer)
        filesToTransfer += ","
        filesToTransfer += options.filesToTransfer
        condorFile.write('transfer_input_files = '+filesToTransfer+'\n')
        condorFile.write('queue $(N)\n')


usage = "usage: %prog [options] \nExample: ./scripts/submit_batch.py -i HeepStudies_v1/MinimumBias__Commissioning10-SD_EG-v9__RECO_short.txt -c HeepStudies_v1/cutFile_HeepElectronStudiesV1.txt -o TestFrancesco/Mydataset -t rootTupleTree/tree -n 2 -q 1nh -d /eos/cms/store/user/eberry/"

parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputlist", dest="inputlist",
                  help="list of all datasets to be used",
                  metavar="LIST")

parser.add_option("-c", "--cutfile", dest="cutfile",
                  help="name of the cut file",
                  metavar="CUTFILE")

parser.add_option("-o", "--output", dest="output",
                  help="the directory OUTDIR contains the output of the program",
                  metavar="OUTDIR")

parser.add_option("-t", "--treeName", dest="treeName",
                  help="name of the root tree",
                  metavar="TREENAME")

parser.add_option("-n", "--ijobmax", dest="ijobmax",
                  help="max number of jobs, limited automatically to the number of files in inputlist",
                  metavar="IJOBMAX")

# http://batchdocs.web.cern.ch/batchdocs/local/submit.html
parser.add_option("-q", "--queue", dest="queue",
                  help="name of the queue",
                  metavar="QUEUE")

parser.add_option("-d", "--eosDir", dest="eosDir",
                  help="full path of the eos directory for the skim output",
                  metavar="EOSDIR")

parser.add_option("-m", "--eosHost", dest="eosHost",
                  help="root:// MGM URL for the eos host for the skim output",
                  metavar="EOSHOST", default="root://eoscms.cern.ch/")

parser.add_option("-j", "--json", dest="jsonFileName",
                  help="json filename",
                  metavar="JSONFILE", default="")

parser.add_option("-b", "--branchSel", dest="branchSelFileName",
                  help="branch selection filename",
                  metavar="BRANCHSELFILE", default="")

parser.add_option("-e", "--exe", dest="executable",
                  help="executable",
                  metavar="EXECUTABLE", default="")

parser.add_option("-r", "--reducedSkim", dest="reducedSkim",
                  help="is this a reduced skim?",
                  metavar="REDUCEDSKIM", default=False, action="store_true")

parser.add_option("-s", "--nanoSkim", dest="nanoSkim",
                  help="is this a nanoAOD skim?",
                  metavar="NANOSKIM", default=False, action="store_true")

parser.add_option("-f", "--filesToTransfer", dest="filesToTransfer",
                  help="comma-separated list of files for condor to transfer",
                  metavar="FILESTOTRANSFER", default="")


(options, args) = parser.parse_args()

if len(sys.argv) < 14:
    print(usage)
    sys.exit()

################################################
pwd = os.getcwd()

outputmain = options.output.rstrip("/")
if not re.search("^/", outputmain):
    outputmain = pwd + "/" + outputmain

inputlist = options.inputlist
if not re.search("^/", inputlist):
    inputlist = pwd + "/" + inputlist

cutfile = options.cutfile
if not re.search("^/", cutfile):
    cutfile = pwd + "/" + cutfile

execName = options.executable.split('/')[-1]
################################################
# write on local disk
################################################
os.system("mkdir -p "+outputmain)
os.system("mkdir -p "+outputmain+"/log/")
os.system("mkdir -p "+outputmain+"/error/")
os.system("mkdir -p "+outputmain+"/input/")
os.system("mkdir -p "+outputmain+"/src/")
os.system("mkdir -p "+outputmain+"/output/")
# os.system("mkdir -p "+outputmain+"/skim/")
#################################################
# output prefix
outputPrefix = outputmain.split("/")[-1]
#################################################
# dataset
dataset = outputPrefix.split("___")[-1]
#################################################
# queue
queue = options.queue
# special override
if queue is not None and "DYJetsToEE_M-50_massWgtFix" in dataset:
    queue = "testmatch"
################################################
# create eos dir
################################################
outputeosdir = options.eosDir.rstrip('/') + '/' + dataset
outputcondoreosdir = options.eosDir.rstrip('/') + '/condor/' + dataset
outputskimeosdir = options.eosDir.rstrip('/') + '/skim/' + dataset
################################################
with open(inputlist, "r") as inputFile:
    allInputFiles = inputFile.readlines()
numfiles = len(allInputFiles)
ijobmax = int(options.ijobmax)
if ijobmax < 0:
    ijobmax = numfiles
    filesperjob = 1
else:
    if ijobmax > numfiles:
        ijobmax = numfiles
    filesperjob = int(numfiles/ijobmax)
    if numfiles % ijobmax != 0:
        filesperjob = filesperjob+1
        ijobmax = int(numfiles/filesperjob)
        if numfiles % filesperjob != 0:
            ijobmax = ijobmax+1
input = open(inputlist)
for ijob in range(ijobmax):
    # prepare the list file
    inputfilename = outputmain+"/input/input_"+str(ijob)+".list"
    inputfile = open(inputfilename, "w")
    filesThisJob = 0
    for i in range(filesperjob):
        line = input.readline()
        if line.startswith("#"):
            continue
        if line != "":
            inputfile.write(line)
            filesThisJob += 1
        continue
    inputfile.close()

    if filesThisJob < 1:
        os.remove(inputfilename)
        ijobmax -= 1
        continue
    # prepare the exec script
    outputname = outputmain+"/src/submit_"+str(ijob)+".sh"
    PrepareJobScript(outputname)

input.close()

# write condor submit file
condorFileName = outputmain+'/condorSubmit.sub'
WriteSubmitFile(condorFileName)

if options.reducedSkim or options.nanoSkim:
    fullDatasetName = Path(inputlist.replace(".txt", "_dataset.txt")).read_text()
    # siteList = GetSiteListFromDAS(allInputFiles, dataset)
    siteList = GetSitesHostingAllFiles(fullDatasetName)
    if len(siteList) < 1:
        siteList = GetDefaultCondorSites()
    # apply exclusions
    sitesToExclude = ["T2_TR_METU", "T2_UA_KIPT", "T1_ES_PIC", "T2_TW_NCHC", "T2_US_MIT", "T2_UA_KIPT", "T2_ES_IFCA", "T2_EE_Estonia", "T2_ES_CIEMAT"]
    for site in sitesToExclude:
        if site in siteList:
            siteList.remove(site)
else:
    siteList = []

failedToSub = False
# FIXME don't cd and use absolute paths in the condor submission instead
oldDir = os.getcwd()
os.chdir(outputmain)
print('submit jobs for', options.output.rstrip("/"), end="")
condorCmd = 'condor_submit '+condorFileName
# check site list for T1s or T2s
foundT1orT2 = False
for site in siteList:
    if site.startswith("T1") or site.startswith("T2"):
        foundT1orT2 = True
if foundT1orT2:
    condorSites = 'CONDOR_DEFAULT_DESIRED_SITES="'+",".join(siteList)
    condorSites = condorSites.strip(",")+'"'
    condorCmd = condorSites + " " + condorCmd
    print(" with", condorSites)
else:
    print()
exitCode = os.WEXITSTATUS(os.system(condorCmd))
# print 'from condor_submit '+condorFileName+',got exit code='+str(exitCode)
if exitCode != 0:
    print('\texited with '+str(exitCode)+'; try to resubmit')
    exitCode = os.WEXITSTATUS(os.system('condor_submit '+condorFileName))
    if exitCode != 0:
        failedToSub = True
os.chdir(oldDir)
if failedToSub:
    exit(-1)
