#!/usr/bin/env python3

from optparse import OptionParser
import os,sys, errno, time
import subprocess as sp
import re
try:
  from CRABClient.UserUtilities import config, getUsernameFromSiteDB
except ImportError:
  print()
  print('ERROR: Could not load CRABClient.UserUtilities.  Please source the crab3 setup:')
  print('source /cvmfs/cms.cern.ch/crab3/crab.sh')
  exit(-1)

import deleteCrabSandboxes

#-------------------------------------------------------------------------------------------------------------
# Notes
# Use this script first, which then calls the submit_crab3 script to actually submit the jobs for each dataset
#-------------------------------------------------------------------------------------------------------------

#-----------------------------------------------------------------------------------------------------------------------
# Function to copy any files on afs specified in the cutfile to the local directory, and feed them into the crab sandbox
#-----------------------------------------------------------------------------------------------------------------------
# NB: If you put a '#' comment in the middle of the line, this may not handle it properly
def localizeCutFile(cutfileLines):
  global options
  extraInputFiles=[]
  replacementLine=''
  lineIdx=0
  for index,line in enumerate(cutfileLines):
      if line.strip() == "" : continue
      if line.startswith('#') : continue
      #if line.strip()[:4] == "JSON":
      lineSp = re.split(r'(\s+)', line)
      for itemIndex,item in enumerate(lineSp):
        if '/afs' in item:
          if not os.path.isfile( item ) : 
              print("Error: No file here: " + item)
              sys.exit(-1)
          else :
              os.system ( "cp " + item + " " + options.outputDir )
              lineIdx = index
              itemIdx = itemIndex
              lineSp[itemIndex] = item.split('/')[-1]
              extraInputFiles.append(lineSp[itemIndex])
              replacementLine = ''.join(lineSp)
              print('INFO: changing cutfile line with afs file specified from "'+cutfileLines[lineIdx].rstrip('\n')+'"')
              cutfileLines[lineIdx] = replacementLine
              print('\tto "'+cutfileLines[lineIdx].rstrip('\n')+'"')
              with open(options.outputDir+'/'+options.cutfile.split('/')[-1],'w') as myfile:
                myfile.writelines(cutfileLines)
  return extraInputFiles


#--------------------------------------------------------------------------------
# Parse options
#--------------------------------------------------------------------------------

usage = "usage: %prog [options] \nExample: python ./scripts/launchAnalysis_batch_ForSkimToEOS.py <options>"

parser = OptionParser(usage=usage)

parser.add_option("-i", "--inputlist", dest="inputlist",
                  help="list of all datasets to be used",
                  metavar="LIST")

parser.add_option("-o", "--output", dest="outputDir",
                  help="the directory OUTDIR contains the output of the program",
                  metavar="OUTDIR")

parser.add_option("-n", "--treeName", dest="treeName",
                  help="name of the root tree; defaults to rootTupleTree/tree",
                  metavar="TREENAME")

parser.add_option("-c", "--cutfile", dest="cutfile",
                  help="name of the cut file",
                  metavar="CUTFILE")

#parser.add_option("-j", "--ijobmax", dest="ijobmax",
#                  help="max number of jobs, limited automatically to the number of files in inputlist",
#                  metavar="IJOBMAX")

# TODO: eventually support this by using the crab scheduler?
#parser.add_option("-q", "--queue", dest="queue",
#                  help="name of the queue (choose among cmst3 8nm 1nh 8nh 1nd 1nw)",
#                  metavar="QUEUE")

# FIXME SIC: this doesn't work when using the API?
parser.add_option("-z", "--dryRun", dest="dryRun",
                  metavar="DRYRUN",default=False,action="store_true")

parser.add_option("-r", "--reducedSkim", dest="isReducedSkimTask",
                  metavar="REDUCEDSKIMTASK",default=False,action="store_true")

parser.add_option("-s", "--skim", dest="isSkimTask",
                  help="use this option to create skims",
                  metavar="SKIMTASK",default=False,action="store_true")

parser.add_option("-d", "--eosDir", dest="eosDir",
                  help="the EOS directory where output root/dat files are stored",
                  metavar="EOSDIR")

parser.add_option("-l", "--overrideOutputLength", dest="overrideOutputLength",
                  metavar="OVERRIDEOUTPUTLENGTH",default=False,action="store_true")

parser.add_option("-f", "--cernT2Only", dest="submitCERNT2only",
                  metavar="submitCERNT2only",default=False,action="store_true")

parser.add_option("-e", "--executable", dest="executable",
                  metavar="executable",default='main')

(options, args) = parser.parse_args()

if ( not options.inputlist 
     or not options.outputDir
     or not options.cutfile 
     #or not options.ijobmax 
     #or not options.queue 
     or not options.eosDir ) : 
    parser.print_help()
    sys.exit()

if not options.treeName:
  options.treeName='rootTupleTree/tree'

#--------------------------------------------------------------------------------
# Delete crab sandboxes
#--------------------------------------------------------------------------------
print('First, delete existing crab sandboxes...', end=' ')
try:
  deleteCrabSandboxes.main()
except deleteCrabSandboxes.Crab3ToolsException:
  print('WARNING: Something went wrong deleting the existing crab sandboxes; proceeding anyway but we might run out of quota')
print('... done ')

#--------------------------------------------------------------------------------
# Make the output directory
#--------------------------------------------------------------------------------

print("Making the local output directory...", end=' ')

if not os.path.isdir ( options.outputDir ) : 
    os.system ( "mkdir -p " + options.outputDir )

if not os.path.isdir ( options.outputDir ) : 
    print("Error: I can't make this folder: " + options.outputDir) 
    sys.exit() 

print("... done ")

#--------------------------------------------------------------------------------
# Look for the cut file.  If it exists, move it to the output directory
#--------------------------------------------------------------------------------

print("Moving the cutfile to the local output directory...", end=' ')

if not os.path.isfile ( options.cutfile ) : 
    print("Error: No cut file here: " + options.cutfile) 
    sys.exit() 
else : 
    os.system ( "cp " + options.cutfile + " " + options.outputDir + "/" )

print("... done ")

#--------------------------------------------------------------------------------
# Look for the inputList file
#--------------------------------------------------------------------------------

print("Moving the inputlist to the local output directory...", end=' ')

if not os.path.isfile ( options.inputlist ) : 
    print("Error: No input list here: " + options.inputlist) 
    sys.exit() 
else : 
    os.system ( "cp " + options.inputlist + " " + options.outputDir + "/" )

print("... done ")

#--------------------------------------------------------------------------------
# Get any /afs/... files from cut file and copy them to the output directory
#--------------------------------------------------------------------------------

print("Moving any /afs/... files specified in the cutfile to the local output directory...")

with open(options.outputDir+'/'+options.cutfile.split('/')[-1],'r') as cutfile:
  cutfileLines = cutfile.readlines()

additionalInputFiles = localizeCutFile(cutfileLines)
print('INFO: Found additional input files:',additionalInputFiles)

print("... done ")
                   
#--------------------------------------------------------------------------------
# Check if path is a link
#--------------------------------------------------------------------------------

print("Checking the link to analysisClass.C...", end=' ')

if not os.path.islink ( "src/analysisClass.C" ) :
    print("Error: src/analysisClass.C is not a symbolic link")
    sys.exit()
code_name = os.readlink ( "./src/analysisClass.C" ).split("/")[-1].split(".C")[0]

print("... done")

#--------------------------------------------------------------------------------
# Launch
#--------------------------------------------------------------------------------

print("Launching jobs...")
if options.isSkimTask or options.isReducedSkimTask:
  print('INFO: This is a SKIM task')
else:
  print('INFO: This is an ANA task')

inputlist_file = file ( options.inputlist,"r" )

total_jobs = 0
failedCommands = ''

for i,line in enumerate(inputlist_file): 
    if line[0] == "#" : continue
    dataset = line.strip().split("/")[-1].split(".txt")[0]
    
    sublist = line.strip()
    if len(sublist) <= 0:
      continue

    #jobs_to_submit = int(options.ijobmax)
    
    #total_jobs = total_jobs + jobs_to_submit

    command = "./scripts/submit_crab3.py"
    command = command + " -i " + sublist
    command = command + " -c " + options.outputDir+'/'+options.cutfile.split('/')[-1]
    command = command + " -t " + options.treeName 
    command = command + " -o " + options.outputDir + "/" + code_name + "___" + dataset
    #command = command + " -n " + str(jobs_to_submit)
    #command = command + " -q " + options.queue
    command = command + " -d " + options.eosDir
    #FIXME not 100% implemented because of using the link name as the code name
    command = command + " -e " + options.executable
    if options.isSkimTask:
      command+=" -s"
    elif options.isReducedSkimTask:
      command+=" -r"
    if options.dryRun:
      command+=" -z"
    if options.overrideOutputLength:
      command+=" -l"
    if options.submitCERNT2only:
      command+=" -f"
    if len(additionalInputFiles) > 0:
      command+=' -j '
      for item in additionalInputFiles:
        command+=options.outputDir+'/'+item+','
      command = command.rstrip(',')
    
    print(command)
    ret = os.system  ( command ) 
    if ret != 0:
      print('ERROR: something went wrong when running the command:')
      print('\t'+command)
      print('add to list')
      failedCommands+=command
      failedCommands+='\n'
    #time.sleep (  float( options.wait ) ) 
    ## Delete crab sandboxes every 20 submissions
    # Don't do this for now--it might screw up submissions...
    #if i % 20 == 0 and not i==0:
    #  print '----> Do periodic delete of crab sandboxes'
    #  try:
    #    deleteCrabSandboxes.main()
    #  except deleteCrabSandboxes.Crab3ToolsException:
    #    print 'WARNING: Something went wrong deleting the existing crab sandboxes; proceeding anyway but we might run out of quota'
    #  print '----> Done'


inputlist_file.close() 

print('list of failed commands:')
print(failedCommands)
#print "total jobs =", total_jobs

