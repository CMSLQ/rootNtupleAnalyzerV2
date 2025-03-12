#!/usr/bin/env python3
import sys
import pathlib
import combineCommon

# XXX FIXME Feb. 2025
# WARNING: THIS SCRIPT SEEMS BROKEN!
# For some reason, the output of this does not agree with the output of the per-sample dat files we make in the combinePlots step!

def GetSampleName(datFile):
    with open(datFile, "r") as theFile:
        for line in theFile:
            return line.split()[-1]


tableDir = sys.argv[1].rstrip("/") + "/"
samplesToCombineFile = sys.argv[2]

files = [f for f in pathlib.Path(tableDir).glob("*.dat")]

sampleTables = {}
for datFile in files:
    datFile = str(datFile)
    data = combineCommon.ParseDatFile(datFile)
    sampleName = GetSampleName(datFile)
    # print("For file", datFile, "got sampleName=", sampleName)
    rootFilename = datFile.replace("tables.dat", "plots.root")
    data = combineCommon.FillTableErrors(data, rootFilename, sampleName)
    # sampleTable = combineCommon.UpdateTable(data, sampleTable)
    sampleTables[sampleName] = data


outputDatFile = "combinedTable.dat"

# combinedTable = combineCommon.CalculateEfficiency(sampleTable)

dictSamples = combineCommon.GetSamplesToCombineDict(samplesToCombineFile)
with open(outputDatFile, "w") as theFile:
    for sample in dictSamples.keys():
        if sample in sampleTables.keys():
            # combineCommon.WriteTable(table, sample, theFile, writeErrNPass=False)
            combineCommon.WriteTable(sampleTables[sample], sample, theFile)
