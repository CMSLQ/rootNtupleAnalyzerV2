#!/bin/env python
import ROOT
import numpy
import sys
import time
import os
import shlex
import subprocess


if len(sys.argv) < 3 :
        print("Syntax: haddnano.py out.root input1.root input2.root ...")
ofname=sys.argv[1]
files=sys.argv[2:]


def zeroFill(tree,brName,brObj,allowNonBool=False) :
        #print("--> zeroFill(brName={}, brObj={}".format(brName, brObj))
        #print("\tbrObjName={}".format(brObj.GetName()))
        # typename: (numpy type code, root type code)
        branch_type_dict = {'Bool_t':('?','O'), 'Float_t':('f4','F'), 'UInt_t':('u4','i'), 'Long64_t':('i8','L'), 'Double_t':('f8','D')}
        brType = brObj.GetLeaf(brName).GetTypeName()
        if (not allowNonBool) and (brType != "Bool_t") :
                print("Did not expect to back fill non-boolean branches",tree,brName,brObj.GetLeaf(br).GetTypeName())
        else :
                if brType not in branch_type_dict: raise RuntimeError('Impossible to backfill branch of type %s'%brType)
                buff=numpy.zeros(1,dtype=numpy.dtype(branch_type_dict[brType][0]))
                b=tree.Branch(brName,buff,brName+"/"+branch_type_dict[brType][1])
                b.SetBasketSize(tree.GetEntries()*2) #be sure we do not trigger flushing
                for x in range(0,tree.GetEntries()):    
                        b.Fill()
                b.ResetAddress()


def RunCommand(cmd):
    print(cmd)
    proc = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = proc.communicate()
    output = out.decode()
    error = err.decode()
    returnCode = proc.returncode
    if returnCode != 0:
        print()
        print("cmd output:", output)
        print("cmd error:", error)
    print(flush=True)
    return returnCode


def GetTFile(lfn):
    print("Checking to see if this file {} is available at {}...".format(lfn, cmsSite))
    # # try this trick sent on mattermost to get the prefix
    # import json
    # try:
    #     data = json.load(open("/cvmfs/cms.cern.ch/SITECONF/{}/storage.json".format(cmsSite)))
    #     sites = [cmsSite for cmsSite in data if cmsSite["type"]=="DISK"]
    #     prefix = None
    #     for site in sites:
    #         theprotoc = next(protoc for protoc in site["protocols"] if protoc["protocol"]=="XRootD")
    #         if "prefix" not in theprotoc.keys():
    #             continue
    #         prefix = theprotoc["prefix"]
    #         host = prefix[prefix.find("//")+2:prefix.rfind("//")]
    #         if not len(host):
    #             # go to next slash instead
    #             host = prefix[prefix.find("//")+2:prefix.find("/", prefix.find("//")+2)]
    #         prefixEnd = prefix[prefix.find(host)+len(host):]
    #         # if not prefixEnd.endswith("/"):
    #         #     prefixEnd += "/"
    #         print("Found prefix from storage.json on cvmfs:", prefix, "with host:", host, "prefixEnd:", prefixEnd)
    #         break
    #     if prefix is None:
    #         raise RuntimeError("Could not find a site specification in the storage json with 'prefix' key in XRootD protocol")
    #     # cmd = 'xrdfs {} ls {}'.format(host, lfnTest)
    #     # returnCode = RunCommand(cmd)
    #     # if returnCode != 0:
    #     #     print("xrdfs ls returned nonzero with output/error above; try to guess direct URL")
    #     # else:
    #     #     lfnToUse = "root://" + host + "/" + prefixEnd
    #     #     print("OK, use {}".format(lfnToUse))
    #     #     return lfnToUse
    #     lfnToTest = prefix + lfn
    #     handle = ROOT.TFile.Open(lfnToTest)
    #     print("OK, use prefix  {}".format(prefix))
    #     return handle
    # except (FileNotFoundError, KeyError, RuntimeError, StopIteration, OSError) as ex:
    #     print("\tCaught exception trying to get and check prefix from storage.json on cvmfs:", ex, "; try to guess direct URL")
    if cmsSite is not None:
        lfnTest = "/store/test/xrootd/" + cmsSite
    else:
        print("cmsSite is None; use default XRootD redirector {}".format(lfnPrefixDefault))
        return ROOT.TFile.Open(lfnPrefixDefault + lfn)
    lfnToUse = lfnPrefixDefault + lfnTest
    oldLevel = ROOT.gErrorIgnoreLevel
    try:
        # we now suppress the error messages printed by TFile::Open for this test
        ROOT.gErrorIgnoreLevel = 6001
        handle = ROOT.TFile.Open(lfnToUse + lfn)
        print("OK, use {}".format(lfnToUse))
    except OSError as ex:
        print("\tCaught exception:", ex)
        print("ROOT.TFile.Open() could not open the file for some reason; use default XRootD redirector {}".format(lfnPrefixDefault))
        handle = ROOT.TFile.Open(lfnPrefixDefault + lfn)
    finally:
        ROOT.gErrorIgnoreLevel = oldLevel
    return handle


# if os.getenv("GLIDEIN_CMSSite").split("_")[1] == "US":
#     lfnPrefixDefault = "root://cmsxrootd.fnal.gov/"
#     hostDefault = "cmsxrootd.fnal.gov"
# else:
#     lfnPrefixDefault = "root://cms-xrd-global.cern.ch/"
#     hostDefault = "cms-xrd-global.cern.ch"
lfnPrefixDefault = "root://cms-xrd-global.cern.ch/"
hostDefault = "cms-xrd-global.cern.ch"
cmsSite = os.getenv("GLIDEIN_CMSSite")
host = hostDefault
print("found GLIDEIN_CMSSite={}".format(cmsSite))
fileHandles=[]
goFast=True
start = time.time()
host = ""
for idx, fn in enumerate(files):
    lfn = fn.replace("root://cms-xrd-global.cern.ch/", "")
    handle = GetTFile(lfn)
    fileHandles.append(handle)
    if fileHandles[-1].GetCompressionSettings() != fileHandles[0].GetCompressionSettings() :
        goFast=False
        print("Disabling fast merging as inputs have different compressions")
of=ROOT.TFile(ofname,"recreate")
if goFast :
        of.SetCompressionSettings(fileHandles[0].GetCompressionSettings())
of.cd()

for e in fileHandles[0].GetListOfKeys() :
        name=e.GetName()
        if name=='rootTupleTree':
            obj=fileHandles[0].Get("rootTupleTree/tree")
            obj.SetTitle("tree")
            obj.SetName("tree")
        else:
            obj=e.ReadObj()
        print("Merging" ,name)
        inputs=ROOT.TList()
        inputs.SetName("inputsList")
        dirInputs=dict()
        #isDir=ROOT.TClass.GetClass(e.GetClassName())=='TDirectoryFile'
        isDir=obj.IsA().InheritsFrom(ROOT.TDirectoryFile.Class())
        if isDir:
            for dirObj in obj.GetListOfKeys():
                thisObj = dirObj.ReadObj()
                dirInputs[thisObj.GetName()] = ROOT.TList()
                dirInputs[thisObj.GetName()].SetName(thisObj.GetName()+"list")
        isTree= obj.IsA().InheritsFrom(ROOT.TTree.Class())
        if isTree:
                branchNames=set([x.GetName() for x in obj.GetListOfBranches()])
                obj.SetBranchStatus("*", 1)  # need to explicitly activate branches which have no entries for CloneTree
                obj=obj.CloneTree(-1,"fast" if goFast else "")
        for fh in fileHandles[1:] :
                if name!='rootTupleTree':
                    otherObj=fh.GetListOfKeys().FindObject(name).ReadObj()
                else:
                    otherObj=fh.Get("rootTupleTree/tree")
                inputs.Add(otherObj)
                if isTree and (obj.GetName()=='Events' or obj.GetName()=='tree'): 
                        otherObj.SetAutoFlush(0)
                        otherBranches=set([ x.GetName() for x in otherObj.GetListOfBranches() ])
                        missingBranches=list(branchNames-otherBranches)
                        additionalBranches=list(otherBranches-branchNames)
                        print("missing:",missingBranches,"\n Additional:",additionalBranches)
                        for br in missingBranches :
                                #fill "Other"
                                zeroFill(otherObj,br,obj.GetListOfBranches().FindObject(br))
                        for br in additionalBranches :
                                #fill main
                                branchNames.add(br)
                                zeroFill(obj,br,otherObj.GetListOfBranches().FindObject(br))
                        #merge immediately for trees
                if isTree and obj.GetName()=='Runs':
                        otherObj.SetAutoFlush(0)
                        otherBranches=set([ x.GetName() for x in otherObj.GetListOfBranches() ])
                        missingBranches=list(branchNames-otherBranches)
                        additionalBranches=list(otherBranches-branchNames)
                        print("missing:",missingBranches,"\n Additional:",additionalBranches)
                        for br in missingBranches :
                                #fill "Other"
                                zeroFill(otherObj,br,obj.GetListOfBranches().FindObject(br),allowNonBool=True)
                        for br in additionalBranches :
                                #fill main
                                branchNames.add(br)
                                zeroFill(obj,br,otherObj.GetListOfBranches().FindObject(br),allowNonBool=True)
                        #merge immediately for trees
                if isTree:
                        obj.Merge(inputs,"fast" if goFast else "")
                        inputs.Clear()
                if isDir:
                    for objName in dirInputs.keys():
                        otherDirObj = otherObj.GetListOfKeys().FindObject(objName).ReadObj()
                        dirInputs[objName].Add(otherDirObj)
        
        if isTree and obj.GetTitle()=="tree":
            of.mkdir("rootTupleTree")
            of.cd("rootTupleTree")
            obj.Write()
            of.cd()
        elif isTree:
            obj.Write()
        elif isDir:
            of.mkdir(e.GetName())
            of.cd(e.GetName())
            for dirObj in obj.GetListOfKeys():
                thisObj = dirObj.ReadObj()
                thisObj.Merge(dirInputs[thisObj.GetName()])
                thisObj.Write()
            of.cd()
        elif obj.IsA().InheritsFrom(ROOT.TH1.Class()) :         
                obj.Merge(inputs)
                obj.Write()
                inputs.Clear()
        elif obj.IsA().InheritsFrom(ROOT.TObjString.Class()) :  
                for st in inputs:
                        if  st.GetString()!=obj.GetString():
                                print("Strings are not matching")
                obj.Write()
                inputs.Clear()
        else:
                print("Cannot handle ", obj.IsA().GetName())
        
print("Time used for hadding files: %.2f s" % (time.time() - start))
