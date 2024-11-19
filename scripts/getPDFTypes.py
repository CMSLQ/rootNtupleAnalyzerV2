#!/usr/bin/env python3

import sys
import glob
from pathlib import Path
from combineCommon import ExtractBranchTitles, GetPDFVariationType, GetBranchTitle
from ROOT import TFile


def get_files(path=Path('.')):
    files = []
    for entry in path.iterdir():
        if entry.is_file():
            files.append(entry)
        elif entry.is_dir():
            files.extend(get_files(entry))
    return files


filePath = sys.argv[1]
directory_path = Path(filePath)
rootFiles = [str(f) for f in get_files(directory_path) if ".root" in str(f)]

for rFile in rootFiles:
    tfile = TFile.Open(rFile)
    sampleTMap = tfile.Get("systematicNameToBranchesMap")
    sampleSystHist = tfile.Get("systematics")
    systNameToBranchTitleDict = ExtractBranchTitles(sampleSystHist, sampleTMap)
    sample = ""
    try:
        pdfVariationType, pdfName = GetPDFVariationType(GetBranchTitle("LHEPdfWeight", sample, systNameToBranchTitleDict)[0])
    except Exception as e:
        print("ERROR: could not get PDF info for file:", tfile, "because of", e)
        continue
    print("Found pdfVariationType={}, pdfName={}, in file={}".format(pdfVariationType, pdfName, tfile.GetName()))
    tfile.Close()
