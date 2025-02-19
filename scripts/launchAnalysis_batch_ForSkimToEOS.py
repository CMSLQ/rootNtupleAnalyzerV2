#!/usr/bin/python3

from optparse import OptionParser
import os
import sys
import time
import subprocess as sp
import shlex
import re
import glob
import math

# --------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------


def GetTotalsFromEventsFile(eventsFile):
    totalEvents = 0
    totalFiles = 0
    with open(eventsFile, "r") as theFile:
        for line in theFile:
            totalEvents += int(line.strip())
            totalFiles += 1
    return totalEvents, totalFiles


def eos_isdir(path):
    command = "/usr/bin/eos ls " + path
    rfdir_stdout, rfdir_stderr = sp.Popen(
        command, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    ).communicate()
    if "No such file or directory" in rfdir_stderr:
        return False
    elif rfdir_stderr == "":
        return True
    else:
        print("Error: ")
        print("  command = " + command)
        print("  stdout  = " + stdout)
        print("  stderr  = " + stderr)
        sys.exit()


def eos_mkdir(path):

    if not eos_isdir(path):
        os.system("/usr/bin/eos mkdir -p " + path)
    else:
        print("  EOS folder already exists: " + path)

    if not eos_isdir(path):
        print("Error: Could not make this EOS folder: " + path)
        sys.exit()

    return path


def eos_getsize(path):
    nsls_stdout, nsls_stderr = sp.Popen(
        "/usr/bin/eos ls -l " + path, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, universal_newlines=True
    ).communicate()

    if "No such file or directory" in nsls_stderr:
        print("Error: No such file in EOS: " + path)
        sys.exit()

    file_size_in_bytes = int(nsls_stdout.split()[-5])

    return int(file_size_in_bytes)


def get_min_positive_integer_in_list(list):
    min_positive_integer = 1e20
    for item in list:
        if item < min_positive_integer:
            min_positive_integer = int(item)
    return int(min_positive_integer)


def get_n_largest_file_sizes_in_bytes_in_inputlist(n, inputlist):

    n_largest_file_sizes_in_bytes = []

    sublist_file = open(inputlist, "r")
    for sublist_line in sublist_file:

        file_name = sublist_line.strip()
        file_size_in_bytes = 0

        if file_name[:15] == "root://eoscms//":
            file_name = file_name[15:]
            file_size_in_bytes = int(eos_getsize(file_name))
        else:
            file_size_in_bytes = int(os.path.getsize(file_name))

        if len(n_largest_file_sizes_in_bytes) < n:
            n_largest_file_sizes_in_bytes.append(file_size_in_bytes)
        elif (
            get_min_positive_integer_in_list(n_largest_file_sizes_in_bytes)
            < file_size_in_bytes
        ):
            n_largest_file_sizes_in_bytes.append(file_size_in_bytes)
            n_largest_file_sizes_in_bytes.remove(
                int(get_min_positive_integer_in_list(n_largest_file_sizes_in_bytes))
            )

    sublist_file.close()

    return n_largest_file_sizes_in_bytes


def get_mean(list):
    mean = 0
    for item in list:
        mean = mean + item
    float_mean = float(mean) / float(len(list))
    return int(float_mean)


def get_n_files(inputlist):
    sublist_file = open(inputlist, "r")
    files = []
    for sublist_line in sublist_file:
        if ".root" in sublist_line:
            files.append(sublist_line)
    sublist_file.close()
    return len(files)


def FindInputFileAndModifyCutFile(localCutFile, keyword):
    found_file = False
    localFile = ""
    cutfile_lines = []
    with open(localCutFile, "r") as cutfile:
        for line in cutfile:
            if line.strip() == "":
                cutfile_lines.append(line)
                continue
            if line.split()[0] == "#":
                cutfile_lines.append(line)
                continue
            if line.strip().split()[0] == keyword:
                if found_file is True:
                    print("Error: You are only allowed to have one {} file per cut file.".format(keyword))
                    print("cut file = " + options.cutfile)
                    sys.exit()
                # if len(line.split()) != 2:
                #     print("Error: this line in your cut file does not make sense:")
                #     print(line)
                #     print("cut file = " + options.cutfile)
                #     sys.exit()
                inputFileInCutFile = line.split()[1]
                inputFile = os.path.expandvars(inputFileInCutFile)
                # print("INFO: found input file {} for keyword {}".format(inputFile, keyword))
                if os.path.isfile(inputFile):
                    print("Moving the {} file to the local output directory...".format(keyword), end=' ')
                    os.system("cp " + inputFile + " " + options.outputDir)
                    print("... done ")
                    found_file = True
                    localFile = options.outputDir + "/" + inputFile.split("/")[-1]
                    localFileBasePath = os.path.basename(localFile)
                    cutfile_lines.append(line.replace(inputFileInCutFile, localFileBasePath))
                elif os.path.isdir(inputFile):
                    print("Moving the {} dir to the local output directory...".format(keyword), end=' ')
                    os.system("cp -R " + inputFile + " " + options.outputDir)
                    print("... done ")
                    found_file = True
                    localFileBasePath = inputFile.strip("/").split("/")[-1]
                    localFile = options.outputDir + "/" + localFileBasePath
                    if localFileBasePath[-1] != "/":
                        localFileBasePath += "/"
                    print("INFO: replace inputFile={} in line {} with localFileBasePath={}".format(inputFileInCutFile, line, localFileBasePath))
                    cutfile_lines.append(line.replace(inputFileInCutFile, localFileBasePath))
                else:
                    print("Error: No {} file or dir here: {}".format(keyword, inputFile))
                    sys.exit()
            else:
                cutfile_lines.append(line)
    if found_file:
        with open(localCutFile, "w") as cutfile:
            for line in cutfile_lines:
                cutfile.write(line)
    return found_file, localFile

# --------------------------------------------------------------------------------
# Parse options
# --------------------------------------------------------------------------------

usage = "usage: %prog [options] \nExample: python ./scripts/launchAnalysis_batch_ForSkimToEOS.py <options>"

parser = OptionParser(usage=usage)

parser.add_option(
    "-i",
    "--inputlist",
    dest="inputlist",
    help="list of all datasets to be used",
    metavar="LIST",
)

parser.add_option(
    "-o",
    "--output",
    dest="outputDir",
    help="the directory OUTDIR contains the output of the program",
    metavar="OUTDIR",
)

parser.add_option(
    "-n",
    "--treeName",
    dest="treeName",
    help="name of the root tree",
    metavar="TREENAME",
)

parser.add_option(
    "-c", "--cutfile", dest="cutfile", help="name of the cut file", metavar="CUTFILE"
)

parser.add_option(
    "-j",
    "--ijobmax",
    dest="ijobmax",
    help="max number of jobs, limited automatically to the number of files in inputlist",
    metavar="IJOBMAX",
)

# http://batchdocs.web.cern.ch/batchdocs/local/submit.html
parser.add_option(
    "-q",
    "--queue",
    dest="queue",
    help="name of the queue (choose among espresso (20 min), microcentury (1 hr), longlunch (2 hrs), workday (8 hrs), etc.; see http://batchdocs.web.cern.ch/batchdocs/local/submit.html)",
    metavar="QUEUE",
)

parser.add_option(
    "-w",
    "--wait",
    dest="wait",
    help="number of seconds to wait between dataset submissions",
    metavar="WAIT",
)

parser.add_option(
    "-d",
    "--eosDir",
    dest="eosDir",
    help="the EOS directory where skims are stored",
    metavar="EOSDIR",
)

parser.add_option(
    "-m",
    "--eosHost",
    dest="eosHost",
    help="root:// MGM URL for the eos host for the skim output",
    metavar="EOSHOST",
    default="root://eoscms.cern.ch/",
)

parser.add_option(
    "-e",
    "--exe",
    dest="executable",
    help="executable",
    metavar="EXECUTABLE",
    default="build/main",
)

parser.add_option(
    "-r",
    "--reducedSkim",
    dest="reducedSkim",
    help="is this a reduced skim?",
    metavar="REDUCEDSKIM",
    default=False,
    action="store_true",
)

parser.add_option(
    "-s",
    "--nanoSkim",
    dest="nanoSkim",
    help="is this a nanoAOD skim?",
    metavar="NANOSKIM",
    default=False,
    action="store_true",
)

(options, args) = parser.parse_args()

if (
    not options.inputlist
    or not options.outputDir
    or not options.treeName
    or not options.cutfile
    or not options.ijobmax
    # or not options.queue
    # or not options.wait
    or not options.eosDir
):
    print("One of [outputDir,treeName,cutfile,ijobmax,eosDir] not specified")
    parser.print_help()
    sys.exit()
if options.reducedSkim and options.nanoSkim:
    print("options reducedSkim and nanoSkim are mutually exclusive")
    parser.print_help()
    sys.exit()

if "eos/user" in options.eosDir:
    options.eosHost = "root://eosuser.cern.ch/"
# set eos mgm url
print("INFO: Using", options.eosHost, "as eosHost")
os.environ["EOS_MGM_URL"] = options.eosHost

# --------------------------------------------------------------------------------
# Make the EOS dir
# --------------------------------------------------------------------------------

print("Making the EOS output directory...", end=' ')
sys.stdout.flush()

#eosPath = eos_mkdir(options.eosDir)
os.system("xrdfs "+options.eosHost+" mkdir -p \""+options.eosDir+"\"\n")
eosPath = options.eosDir

print("... done")

# --------------------------------------------------------------------------------
# Make the output directory
# --------------------------------------------------------------------------------

print("Making the local output directory...", end=' ')
sys.stdout.flush()

if not os.path.isdir(options.outputDir):
    os.system("mkdir -p " + options.outputDir)

if not os.path.isdir(options.outputDir):
    print("Error: I can't make this folder: " + options.outputDir)
    sys.exit(-2)

print("... done ")

# --------------------------------------------------------------------------------
# Look for the exe file.  If it exists, move it to the output directory
# --------------------------------------------------------------------------------
if os.path.isfile(options.executable):
    print("Moving the exe to the local output directory...", end=' ')
    os.system("cp " + options.executable + " " + options.outputDir + "/")
    print("... done ")
else:
    print("Warning: No file here: '" + options.executable + "'" + "; proceeding anyway")
    # sys.exit()

# --------------------------------------------------------------------------------
# Look for the cut file.  If it exists, move it to the output directory
# --------------------------------------------------------------------------------
print("Moving the cutfile to the local output directory...", end=' ')
if not os.path.isfile(options.cutfile):
    print("Error: No cut file here: " + options.cutfile)
    sys.exit(-1)
else:
    os.system("cp " + options.cutfile + " " + options.outputDir + "/")
print("... done ")
localCutFile = options.outputDir+"/"+os.path.basename(options.cutfile)

# --------------------------------------------------------------------------------
# Look for the inputList file
# --------------------------------------------------------------------------------

print("Moving the inputlist to the local output directory...", end=' ')

if not os.path.isfile(options.inputlist):
    print("Error: No input list here: " + options.inputlist)
    sys.exit()
else:
    os.system("cp " + options.inputlist + " " + options.outputDir + "/")

print("... done ")

# --------------------------------------------------------------------------------
# Get JSON and other files specified in the cut file and copy them to the output directory,
# also modifying the path in the local cut file
# --------------------------------------------------------------------------------
inputFileNames = ["JSON", "BranchSelection", "TriggerSFFileName", "EGMScaleJSONFileName", "JECTextFilePath", "JERTextFilePath"]
inputFilesToTransfer = []
for inputFile in inputFileNames:
    foundFile, filename = FindInputFileAndModifyCutFile(localCutFile, inputFile)
    if foundFile:
        inputFilesToTransfer.append(filename)

# ----------------------------------------------------------------------------------------
# For reduced skims: Check for haddnano.py.  If it exists, move it to the output directory
# ----------------------------------------------------------------------------------------
if options.reducedSkim or options.nanoSkim:
    haddnanoPath = os.path.expandvars("$LQANA/scripts/haddnano.py")
    print("Moving haddnano.py to the local output directory...", end=' ')
    if not os.path.isfile(haddnanoPath):
        print("Error: No haddnano.py here: " + haddnanoPath)
        sys.exit(-1)
    else:
        os.system("cp " + haddnanoPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/haddnano.py")

    print("... done ")
    # Also, copy CMSJME libs
    cmsjmelibPath = os.path.expandvars("$LQANA/build/_deps/cmsjmecalculators-build/libCMSJMECalculators.so")
    print("Moving CMSJME libs to the local output directory...", end=' ')
    if not os.path.isfile(cmsjmelibPath):
        print("Error: No libCMSJMECalculators.so here: " + cmsjmelibPath)
        sys.exit(-1)
    else:
        os.system("cp " + cmsjmelibPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/libCMSJMECalculators.so")
    cmsjmeDictlibPath = os.path.expandvars("$LQANA/build/_deps/cmsjmecalculators-build/libCMSJMECalculatorsDict.so")
    if not os.path.isfile(cmsjmeDictlibPath):
        print("Error: No libCMSJMECalculatorsDict.so here: " + cmsjmeDictlibPath)
        sys.exit(-1)
    else:
        os.system("cp " + cmsjmeDictlibPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/libCMSJMECalculatorsDict.so")
    print("... done ")
    # copy custom root libs
    dictLibPath = os.path.expandvars("$LQANA/build/libmyapp_dict_lib.so")
    print("Moving other libs to the local output directory...", end=' ')
    if not os.path.isfile(dictLibPath):
        print("Error: No libmyapp_dict_lib.so here: " + dictLibPath)
        sys.exit(-1)
    else:
        os.system("cp " + dictLibPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/libmyapp_dict_lib.so")
    dictLibPath = os.path.expandvars("$LQANA/build/libmyapp_dict_rdict.pcm")
    if not os.path.isfile(dictLibPath):
        print("Error: No libmyapp_dict_rdict.pcm here: " + dictLibPath)
        sys.exit(-1)
    else:
        os.system("cp " + dictLibPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/libmyapp_dict_rdict.pcm")
    dictLibPath = os.path.expandvars("$LQANA/build/libmyapp_dict.rootmap")
    if not os.path.isfile(dictLibPath):
        print("Error: No libmyapp_dict.rootmap here: " + dictLibPath)
        sys.exit(-1)
    else:
        os.system("cp " + dictLibPath + " " + options.outputDir + "/")
        inputFilesToTransfer.append(options.outputDir + "/libmyapp_dict.rootmap")
    print("... done ")

# --------------------------------------------------------------------------------
# Check if path is a link
# --------------------------------------------------------------------------------
nameToCheck = "analysisClass"
if "qcd" in options.executable.lower():
    nameToCheck = "analysisClassQCD"
print("Checking the link to {}...".format(nameToCheck), end=' ')

if not os.path.islink("src/{}.C".format(nameToCheck)):
    print()
    print("Error: src/{}.C is not a symbolic link".format(nameToCheck))
    sys.exit()
code_name = os.readlink("./src/{}.C".format(nameToCheck)).split("/")[-1].split(".C")[0]

print("... done")

# --------------------------------------------------------------------------------
# Launch
# --------------------------------------------------------------------------------

print("Launching jobs!")

submitted_jobs = 0
total_jobs = 0

failedCommands = list()

cmdsToRun = []
with open(options.inputlist, "r") as inputlist_file:
    for line in inputlist_file:
        if not len(line.strip()) > 0:
            continue
        if line.strip().startswith("#"):
            continue
    
        dataset = line.strip().split("/")[-1].split(".txt")[0]
        eventsFilename = line.strip().replace(".txt", "_nevents.txt")
        eventsFileList = glob.glob(eventsFilename)
        jobs_to_submit = int(options.ijobmax)
        if len(eventsFileList) == 1:
            # approx splitting by max events per job
            totalNumEvents, totalNumFiles = GetTotalsFromEventsFile(eventsFileList[0])
            # maxEventsDesiredPerJob = 2.5e6
            maxEventsDesiredPerJob = 2e6
            suggestedJobs = math.ceil(totalNumEvents/maxEventsDesiredPerJob)
            jobs_to_submit = min ( suggestedJobs, int(options.ijobmax))
    
        if jobs_to_submit > 0:
            total_jobs = total_jobs + jobs_to_submit
        else:
            with open(line.strip()) as inputFile:
                numLines = sum(1 for line in inputFile)
            total_jobs = total_jobs + numLines
    
        command = "./scripts/submit_batch_ForSkimToEOS.py"
        command = command + " -i " + line.strip()
        command = (
            command + " -c " + options.outputDir + "/" + options.cutfile.split("/")[-1]
        )
        command = command + " -t " + options.treeName
        command = command + " -o " + options.outputDir + "/" + code_name + "___" + dataset
        command = command + " -n " + str(jobs_to_submit)
        if options.queue is not None:
          command = command + " -q " + options.queue
        command = command + " -d " + eosPath
        command = (
            command + " -e " + options.outputDir + "/" + options.executable.split("/")[-1]
        )
        command = command + " -m " + options.eosHost
        if len(inputFilesToTransfer):
            command += " -f " + ",".join(inputFilesToTransfer)
        if options.reducedSkim:
            command = command + " -r "
        elif options.nanoSkim:
            command += " -s "
    
        cmdsToRun.append(command)

# use 6 processes to submit in parallel
# procsToUse = 6
# for i in range(0, len(cmdsToRun), procsToUse):
#     cmds = cmdsToRun[i:i+procsToUse]
#     for cmd in cmds:
#         print(cmd)
#     procs = [sp.Popen(shlex.split(i), stdout=sp.PIPE, stderr=sp.PIPE) for i in cmds]
#     for idx, p in enumerate(procs):
#         res = p.communicate()
#         for line in res[0].decode(encoding='utf-8').split('\n'):
#               print(line)
#         if p.returncode != 0:
#              print("stderr =", res[1])
#              print("Failed submitting jobs for this dataset; add to failedCommands list")
#              failedCommands.append(cmds[idx])
# one at a time
for cmd in cmdsToRun:
    try:
        proc = sp.run(shlex.split(cmd), check=True, universal_newlines=True, stdout=sp.PIPE, stderr=sp.PIPE)
    except sp.CalledProcessError as ex:
        print("stdout = ", ex.stdout)
        print("stderr = ", ex.stderr)
        print("Failed submitting jobs for this dataset; add to failedCommands list")
        failedCommands.append(cmd)
        # raise ex
        continue
    output = proc.stdout.strip()
    print(output)
    jobsThisDataset = int(re.search("(\d+)\sjob\(s\)\ssubmitted", output).group(1))  # parse from condor_submit output
    submitted_jobs += jobsThisDataset
    if proc.returncode != 0:
        print("stderr =", proc.stderr)
        print("Failed submitting jobs for this dataset; add to failedCommands list")
        failedCommands.append(cmd)

print("submitted {} jobs in total".format(submitted_jobs))

if len(failedCommands) > 0:
    print("list of failed commands:")
    for cmd in failedCommands:
        print(cmd)
    exit(1)
