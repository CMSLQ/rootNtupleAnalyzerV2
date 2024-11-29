#!/usr/bin/env python3
import os
from optparse import OptionParser
import subprocess
import re
import shutil
from pathlib import Path
from combineCommon import SeparateDatacards, DeleteTmpFiles, GetYearAndIntLumiFromDatacard


def SeparateCombineCardsOutput(output):
    perCardOutput = []
    split = re.split('(Combination)', output)[1:]
    for i in range(0, len(split)-1, 2):
        block = "".join(split[i:i+2])
        perCardOutput.append(block)
    return perCardOutput


def GetShapeRootFile(datacard):
    shapeLine = ""
    with open(datacard, "r") as datacardFile:
        for line in datacardFile:
            if "shapes" in line and "root" in line:
                shapeLine = line
                break  # take first occurance
    return str(Path(datacard).parent.resolve() / shapeLine.split()[3])


combinedOutputCardName = "combCardFile.txt"
# cmsswDir = os.path.expandvars("$LQLIMITS")
#FIXME: need to update this area with EL9 LCG, but seems that this works for now, since it's just python
combineDir = "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/limitSetting/HiggsAnalysis/CombinedLimit/"

cardsPerMassDir = "combinedCardsPerMass"

parser = OptionParser(
    usage="%prog [datacard1] [datacard2] ... [datacardN] [label1] [label2] ... [labelN]",
    )

(options, args) = parser.parse_args()
if not args:
    parser.error("No input datacards specified.")
if len(args)%2 != 0:
    parser.error("Odd number of arguments specified; must specify number of labels equal to number of datacards")

#FIXME: change this so that it uses the combineCards format of label1=card1 label2=card2
datacards = args[:len(args)//2]
labels = args[len(args)//2:]
for label in labels:
    if os.path.isfile(label):
        raise RuntimeError("Got a file {} where a label was expected instead. Must specify number of labels equal to number of datacards.".format(label))

massListByCombinedDatacard = {}
allTmpFilesByMass = {}
totalIntLumi = 0
years = []
shapesRootFiles = []
for index, combinedDatacard in enumerate(datacards):
    year, intLumi = GetYearAndIntLumiFromDatacard(combinedDatacard)
    years.append(year)
    totalIntLumi += intLumi
    shapesRootFiles.append(GetShapeRootFile(combinedDatacard))
    #FIXME add dirpath option to SeparateDatacards
    massList, tmpFilesByMass, cardContentsByMass = SeparateDatacards(combinedDatacard, index, dirPath=os.getcwd())
    massListByCombinedDatacard[combinedDatacard] = massList
    for mass, tmpFile in tmpFilesByMass.items():
        if mass not in allTmpFilesByMass:
            allTmpFilesByMass[mass] = []
        allTmpFilesByMass[mass].append(tmpFile)

print("INFO: Found total int lumi = {}/pb for years = {}".format(totalIntLumi, years))
#print(massListByCombinedDatacard)
referenceMassList = []
referenceDatacard = ""
for combinedCard, massList in massListByCombinedDatacard.items():
    if len(referenceMassList) < 1:
        referenceMassList = sorted(massList)
        referenceDatacard = combinedCard
        continue
    if sorted(massList) != referenceMassList:
        DeleteTmpFiles(allTmpFilesByMass)
        raise RuntimeError("mass list {} from parsing {} is not the same as mass list {} from parsing {}. Can't combine the datacards.".format(massList, combinedCard, referenceMassList, referenceDatacard))

Path(cardsPerMassDir).mkdir(exist_ok=True)

combineCardsCommands = []
with open(combinedOutputCardName, "w") as combCardFile:
    combCardFile.write("# {}\n".format(years))
    combCardFile.write("# {}\n".format(totalIntLumi))
    for mass in referenceMassList:
        cardsForMass = allTmpFilesByMass[mass]
        combineCardsArgs = [label+"="+card for label, card in zip(labels, cardsForMass)]
        #TODO: support args to combineCards
        combineCardsCommands.append("combineCards.py {}".format(" ".join(combineCardsArgs)))
    print("INFO: Creating combined cards for masses {}".format(referenceMassList), flush=True)
    
    # cmd = 'cd {} && eval `scram runtime -sh` && combineCards.py {}'.format(cmsswDir, " ".join(combineCardsArgs))
    # cmd = 'cd {} && eval `scram runtime -sh` && {}'.format(cmsswDir, "&& ".join(combineCardsCommands))
    cmd = 'cd {} && . env_lcg.sh && {}'.format(combineDir, "&& ".join(combineCardsCommands))
    #output = subprocess.check_output(cmd, env={})
    # print("INFO: combineCardsCommands looks like '{}'".format(combineCardsCommands))
    # print("INFO: Run command '{}'".format(cmd), flush=True)
    process = subprocess.Popen(['/bin/bash', '-l', '-c', cmd], env={}, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    errcode = process.returncode
    if errcode != 0:
        raise RuntimeError("Command {} failed with return code {}.\nStderr: {}\n Exiting here.".format("/bin/bash -l -c "+cmd, errcode, err.decode().strip()))
    combinedCardsPerMass = SeparateCombineCardsOutput(out.decode())
    for i, mass in enumerate(referenceMassList):
        with open(cardsPerMassDir.rstrip("/")+"/combCardFile_m{}.txt".format(mass), "w") as combCardFileMass:
            combCardFileMass.write(re.sub("\S*/(\w*\.root)", "\g<1>", combinedCardsPerMass[i]))
        print("INFO: Wrote combined file for mass {} to {}".format(mass, cardsPerMassDir.rstrip("/")+"/combCardFile_m{}.txt".format(mass)))
        # copy root files
        for rf in shapesRootFiles:
            shutil.copy2(rf, cardsPerMassDir)
        combCardFile.write("# LQ_M{}.txt\n".format(mass))
        # combCardFile.write(combinedCardsPerMass[i])
        # rewrite root file path to include the path to the copy we made
        combCardFile.write(re.sub("(\w*\.root)", cardsPerMassDir+"/\g<0>", combinedCardsPerMass[i]))
        combCardFile.write("\n\n")
    # combCardFile.write(out.decode())
    # print("Wrote combination for mass {}".format(mass))
    
DeleteTmpFiles(allTmpFilesByMass)
print("INFO: Wrote combined file for all masses to {}".format(combinedOutputCardName))
