import os
import shutil
import subprocess

#Script copies weight files to the EOS directory where we're keeping them, then updates cutfiles with the location of the weight files.
#Then, it gets the lines with the actual cuts from the log of the optimization step and updates the cutfiles with that too.
#If running for 2017 or 2018, use era = ""

years = ["2016preVFP", "2016postVFP","2017","2018"]
#era = "postVFP"
date = "7Feb2025"

#weight files
weightFiles = os.getenv("LQDATAEOS")+"/BDT_7FebSkim/LQToDEle/dataset/weights"#.format(year,era)
weightFilesDest = "/eos/user/e/eipearso/LQ_BDTWeightFiles/LQToDEle/HTLO-amcatnlo_{}/weights".format(date)

print("copying weight files to "+weightFilesDest)

#optimization results
optResults = os.getenv("LQDATAEOS")+"/BDT_7FebSkim/LQToDEle/minNB0p5/bdtOptimization.log"#.format(year,era)

#cut files
#filenameBase = os.getenv("LQMACRO")+"/config{}/Analysis/{}LQToBEle/".format(year,era+"/")
#filenameBase = os.getenv("LQMACRO")+"/config{}/Analysis/{}HTLO-amcatnlo/".format(year,era+"/")
filenames = [
    "cutTable_lq_eejj_BDT.txt",
    "cutTable_lq_eejj_QCD_doubleFR.txt",
    "cutTable_lq_eejj_QCD_singleFR.txt",
]
filesToUpdate = []
for year in years:
    for f in filenames:
        if year=="2016preVFP":
            directory = os.getenv("LQMACRO")+"/config2016/Analysis/preVFP/HTLO-amcatnlo/"
            #directory = os.getenv("LQMACRO")+"/config2016/Analysis/preVFP/LQToBEle/"
        elif year=="2016postVFP":
            directory = os.getenv("LQMACRO")+"/config2016/Analysis/postVFP/HTLO-amcatnlo/"
            #directory = os.getenv("LQMACRO")+"/config2016/Analysis/postVFP/LQToBEle/"
        else:
            directory = os.getenv("LQMACRO")+"/config{}/Analysis/HTLO-amcatnlo/".format(year)
            #directory = os.getenv("LQMACRO")+"/config{}/Analysis/LQToBEle/".format(year)
        fullPath = directory+f
        print("add file {}".format(fullPath))
        filesToUpdate.append(fullPath)

#Copy weight files
if not os.path.isdir(weightFilesDest):
    os.makedirs(weightFilesDest)
weightFileLocations = {}
for wf in os.listdir(weightFiles):
    if "xml" in wf:
       # print("copy file "+weightFiles+"/"+wf)
       # shutil.copy(weightFiles+"/"+wf, weightFilesDest)
        mass = wf.split("_")[2] #gets M-[mass]
        mass = mass.split("-")[1] #gets just the mass
        weightFileLocations[mass] = weightFilesDest+"/"+wf
    else:
        continue 
print(weightFileLocations)
#get weight file lines for cutfile
weightFilesList = []
for mass in list(range(300,3100,100)):
    line = "BDTWeightFileLQM{}   ".format(mass)+weightFileLocations[str(mass)]+" - - - -1"
    print(line)
    weightFilesList.append(line)

#get cut values for cutfile
cutValsTxt = optResults.replace(".log",".txt")
with open(cutValsTxt, 'w') as outFile:
    proc = subprocess.Popen(['tail', '-n', '160', optResults], stdout=outFile)
    proc.wait()
#Note: the text file here includes the LaTex formatted table, which is intentional.
cutVals = open(cutValsTxt,'r')
optLines = cutVals.readlines()
#I accidentally added a couple of useless lines at the end of the log file
optLines.pop()
optLines.pop()
nLinesToKeep = 112
nLines = len(optLines)
cutValLines = []
#for i in range(nLines-nLinesToKeep, nLines):
#    cutValLines.append(optLines[i])

cutValsDict = {}
for l in optLines:
    if "BDTOutput" in l:
        mass = l.split()[0]
        mass = mass.split("LQ")[1]
        l = l.replace("trainingSelection","MeejjLQ{}".format(mass))
        cutValsDict[mass] = l
    else:
        continue

for mass in range(300, 3100, 100):
    cutValLines.append("#------------------------------------------------------------------------------------------------------------------\n")
    cutValLines.append("# LQ M {} optimization\n".format(mass))
    cutValLines.append("#------------------------------------------------------------------------------------------------------------------\n")
    cutValLines.append("MeejjLQ{}                   {}                 +inf                    -                   -   2       500 0 5000 minselection:trainingSelection\n".format(mass, mass))
    cutValLines.append(cutValsDict[str(mass)])

for l in cutValLines:
    print(l)

#update cutfiles
for f in filesToUpdate:
    #path = filenameBase+f
    print(f)
    if "LQToBEle" in f:
        fToRead = f.replace("LQToBEle","HTLO-amcatnlo")
    else:
        fToRead = f
    with open(fToRead, 'r') as cutfile:
        cutfile.seek(0)
        lines = cutfile.readlines()

    iWeightFiles=0
    iCutVals = 0

    for i,l in enumerate(lines):
        if not "BDT" in l:
            continue #find weight file part
        iWeightFiles = i+2 #The first weight file is 2 lines after this
        break

    for i,l in enumerate(lines):
        if not "CUTS" in l:
            continue
        iCutVals = i+2
        break

    #replace weight files 
    j = 0
    for i in range(iWeightFiles, iWeightFiles+28): #28 masses
        lines[i] = weightFilesList[j]+"\n"
        j+=1

    #replace cut values
    '''
    j = 0
    for i in range(len(lines)-len(cutValLines), len(lines)):
        lines[i] = cutValLines[j]
        j+=1
    '''
    lines = lines[:iCutVals+1]
    for l in cutValLines:
        lines.append(l)

    with open(f, 'w') as cutfile:
        cutfile.writelines(lines)

