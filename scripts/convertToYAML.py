#!/usr/bin/env python3
import sys
from collections import OrderedDict
import ruamel.yaml as yaml

#sampleDict = OrderedDict()
sampleDict = {}
for l, line in enumerate(open(sys.argv[1])):
    line = line.strip()
    if line.startswith("#"):
        continue
    if len(line) <= 0:
        continue
    #print(line)
    pieces = line.split()[1:]
    pieceList = [piece for piece in pieces]
    name = line.split()[0]
    #sampleDict[name] = OrderedDict()
    sampleDict[name] = {}
    sampleDict[name]["pieces"] = pieceList

with open("sample.yaml", "w") as outfile:
    yaml.dump(sampleDict, outfile, Dumper=yaml.RoundTripDumper)
