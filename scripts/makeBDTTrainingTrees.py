#!/usr/bin/env python3
#
import os
import shutil
import sys
import math
import multiprocessing
import traceback
import tempfile
from rich.progress import Progress
from rich.console import Console
from time import sleep

from tmvaBDT import LoadChainFromTxtFile, GetTotalEventsHist

import ROOT
from ROOT import TMVA, TFile, TString, TCut, TChain, TFileCollection, gROOT, gDirectory, gInterpreter, TEntryList, TH1D, TProfile, RDataFrame, TLorentzVector

gROOT.SetBatch()

manager = multiprocessing.Manager()
result_list = manager.list()

logString = "INFO: running {} parallel jobs for {} separate LQ masses requested..."

@ROOT.Numba.Declare(["float","float","float","float"], "float")
def CalcMAsym(M_e1j1, M_e2j2, M_e1j2, M_e2j1):
    if math.fabs(M_e1j1-M_e2j2) < math.fabs(M_e1j2-M_e2j1):
        if M_e1j1 < M_e2j2:
           M_ej_max = M_e2j2
           M_ej_min = M_e1j1
        else:
           M_ej_max = M_e1j1
           M_ej_min = M_e2j2
    else:
        if M_e1j2 < M_e2j1:
           M_ej_max = M_e2j1
           M_ej_min = M_e1j2
        else:
           M_ej_max = M_e1j2
           M_ej_min = M_e2j1
    masym = (M_ej_max-M_ej_min)/(M_ej_max+M_ej_min)
    return masym


@ROOT.Numba.Declare(["float","float","float","float"], "float")
def CalcMejMax(M_e1j1, M_e2j2, M_e1j2, M_e2j1):
    if math.fabs(M_e1j1-M_e2j2) < math.fabs(M_e1j2-M_e2j1):
        if M_e1j1 < M_e2j2:
           M_ej_max = M_e2j2
        else:
           M_ej_max = M_e1j1
    else:
        if M_e1j2 < M_e2j1:
           M_ej_max = M_e2j1
        else:
           M_ej_max = M_e1j2
    return M_ej_max


@ROOT.Numba.Declare(["float","float","float","float"], "float")
def CalcMejMin(M_e1j1, M_e2j2, M_e1j2, M_e2j1):
    if math.fabs(M_e1j1-M_e2j2) < math.fabs(M_e1j2-M_e2j1):
        if M_e1j1 < M_e2j2:
           M_ej_min = M_e1j1
        else:
           M_ej_min = M_e2j2
    else:
        if M_e1j2 < M_e2j1:
           M_ej_min = M_e1j2
        else:
           M_ej_min = M_e2j1
    return M_ej_min


# imported from BDT script
#@ROOT.Numba.Declare(["int"], "float")
#def GetMassFloat(mass):
#    return float(mass)


ROOT.gInterpreter.Declare('''
       float CalcMeejj(float Ele1_Pt, float Ele1_Eta, float Ele1_Phi, float Ele2_Pt, float Ele2_Eta, float Ele2_Phi,
           float Jet1_Pt, float Jet1_Eta, float Jet1_Phi, float Jet2_Pt, float Jet2_Eta, float Jet2_Phi) {
         TLorentzVector e1, j1, e2, j2, eejj;
         e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
         e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
         j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
         j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
         eejj = e1 + e2 + j1 + j2; 
         return eejj.M();
       }
''')

ROOT.gInterpreter.Declare('''
        float CalcMeejjj(float Ele1_Pt, float Ele1_Eta, float Ele1_Phi, float Ele2_Pt, float Ele2_Eta, float Ele2_Phi,
        float Jet1_Pt, float Jet1_Eta, float Jet1_Phi, float Jet2_Pt, float Jet2_Eta, float Jet2_Phi, float Jet3_Pt, float Jet3_Eta, float Jet3_Phi) {
        TLorentzVector e1, j1, e2, j2, j3, eejjj;
        e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
        e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
        j1.SetPtEtaPhiM ( Jet1_Pt, Jet1_Eta, Jet1_Phi, 0.0 );
        j2.SetPtEtaPhiM ( Jet2_Pt, Jet2_Eta, Jet2_Phi, 0.0 );
        j3.SetPtEtaPhiM ( Jet3_Pt, Jet3_Eta, Jet3_Phi, 0.0 );
        eejjj = e1 + e2 + j1 + j2 + j3;
        return eejjj.M();
        }
''')


def log_result(result):
    # This is called whenever foo_pool(i) returns a result.
    # result_list is modified only by the main process, not the pool workers.
    result_list.append(result)
    #sys.stdout.write("\r"+logString.format(jobCount, len(massList)))
    #sys.stdout.write("\t"+str(len(result_list))+" jobs done")
    #sys.stdout.flush()


def handle_error(error):
    print("ERROR: caught exception in job: {}".format(error), flush=True)


def ProcessDataset(args):
    txtFile, mass, branchesToSave, year, qcdLabel, makeWeightsNegative = args
    #print("got args = {}".format(args))
    try:
        tchainBkg = LoadChainFromTxtFile(txtFile)
        if tchainBkg is None:
            return True
        nEntriesBeforeCut = tchainBkg.GetEntries()
        maxTrials = 5
        trials = 0
        redo = True
        df = RDataFrame(tchainBkg)
        df = df.Filter(mycut.GetTitle())  # will work for expressions valid in C++
        df = df.Define("Masym", "Numba::CalcMAsym(M_e1j1, M_e2j2, M_e1j2, M_e2j1)")
        df = df.Define("MejMax", "Numba::CalcMejMax(M_e1j1, M_e2j2, M_e1j2, M_e2j1)")
        df = df.Define("MejMin", "Numba::CalcMejMin(M_e1j1, M_e2j2, M_e1j2, M_e2j1)")
        df = df.Define("Meejj", "CalcMeejj(Ele1_Pt, Ele1_Eta, Ele1_Phi, Ele2_Pt, Ele2_Eta, Ele2_Phi,Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet2_Pt, Jet2_Eta, Jet2_Phi)")
        df = df.Define("Meejjj", "CalcMeejjj(Ele1_Pt, Ele1_Eta, Ele1_Phi, Ele2_Pt, Ele2_Eta, Ele2_Phi,Jet1_Pt, Jet1_Eta, Jet1_Phi, Jet2_Pt, Jet2_Eta, Jet2_Phi, Jet3_Pt, Jet3_Eta, Jet3_Phi)")
        df = df.Define("LQCandidateMassInt", str(mass))
        df = df.Define("LQCandidateMass", "Numba::GetMassFloat(LQCandidateMassInt)")
        if makeWeightsNegative:
            df = df.Redefine("EventWeight", "-1.0*EventWeight")
        # elif "powhegMiNNLO" in txtFile:
        #     print("INFO: Using sign of Weight for this powhegMiNNLO sample: {}".format(txtFile))
        #     df = df.Redefine("Weight", "float(TMath::Sign(1, Weight))")
        expectedEvents = df.Count().GetValue()
        datasetName = os.path.basename(txtFile).replace(".txt", "")
        tempDir = tempfile.mkdtemp()
        if len(qcdLabel):
            tfilepath = "{}/{}.root".format(tempDir,datasetName)
            outputOnEos = outputTFileDir+"/{}/{}".format(qcdLabel, mass)
            #if not os.path.isdir(outputTFileDir+"/{}/{}".format(qcdLabel, mass)):
             #   os.makedirs(outputTFileDir+"/{}/{}".format(qcdLabel, mass))
        else:
            tfilepath = "{}/{}.root".format(tempDir,datasetName)
            outputOnEos = outputTFileDir+"/{}".format(mass)
            #if not os.path.isdir(outputTFileDir+"/{}".format(mass)):
             #   os.mkdir(outputTFileDir+"/{}".format(mass))
        if "root://" not in outputOnEos:
            os.makedirs(outputOnEos, exist_ok=True)
            print("proccess dataset ",datasetName)
        if signalNameTemplate.format(mass) in tfilepath:
            eventCounterHist = GetTotalEventsHist(mass, year, allSignalDatasetsDict, signalNameTemplate)
        else:
            eventCounterHist = GetTotalEventsHist(mass, year, {datasetName: [txtFile]}, datasetName)
        eventCounterHistEntries = eventCounterHist.GetEntries()
        while redo and trials < maxTrials:
            print("INFO: Writing snapshot to:", tfilepath)
        #    if "TTToHadronic" in tfilepath:
        #        print("columns for TTToHadronic sample")
        #        df.Describe().Print()
            #if expectedEvents > 0:
            df.Snapshot("rootTupleTree/tree", tfilepath, branchesToSave)
            #sleep(3)    
            #if signalNameTemplate.format(mass) in tfilepath:
            #    eventCounterHist = GetTotalEventsHist(mass, allSignalDatasetsDict, signalNameTemplate)
            #    tfile = TFile(tfilepath, "update")
            #    tfile.cd()
            #    eventCounterHist.Write()
            #    tfile.Close()
            tfile = TFile.Open(tfilepath, "update")
            tfile.cd()
            eventCounterHist.Write()
            tfile.Close()
         #   else:
        #    print("INFO: expectedEvents={} for txtFile={}, mass={}. Events before cut={}".format(expectedEvents,txtFile,mass,nEntriesBeforeCut))
            # check -- sometimes, for some reason, the snapshot or hist writing fails, so this should catch that
            # expectedEvents = expectedEvents.GetValue()
            if expectedEvents > 0:
                tfile = TFile.Open(tfilepath)
                tree = tfile.Get("rootTupleTree/tree")
                if not tree or tree == None:
                    print("WARN: didn't get a tree from the file {}; redoing training tree for txtFile={}, mass={}".format(tfilepath, txtFile, mass), flush=True)
                elif tree.ClassName() != "TTree":
                    print("WARN: read tree of class {} instead of TTree; redoing training tree for txtFile={}, mass={}".format(tree.ClassName(), txtFile, mass), flush=True)
                else:
                    if tree.GetEntries() != expectedEvents:
                        print("WARN: didn't get expected number of events; redoing training tree for txtFile={}, mass={}".format(txtFile, mass), flush=True)
                    else:
                        redo = False
                        break
                # check for hist in written file
                hist = tfile.Get("EventCounter")
                if not hist or hist is None or hist.GetEntries() != eventCounterHistEntries:
                    print("WARN: didn't get proper EventCounter hist; redoing training tree/hist for txtFile={}, mass={}".format(txtFile, mass), flush=True)
                tfile.Close()
            else:
                redo = False
                break
            trials+=1
        if redo:
            #print("ERROR: After {} trials, didn't write the tree properly into the root file; txtFile={}, mass={}".format(maxTrials, txtFile, mass), flush=True)
            raise RuntimeError("After {} trials, didn't write the tree or hist properly into the root file; txtFile={}, mass={}".format(maxTrials, txtFile, mass))
        else:
            #if len(qcdLabel):
             #   shutil.copy(tfilepath, outputOnEos+"/{}.root".format(datasetName))
            #else:
            shutil.copy(tfilepath, outputOnEos)
            print("copied to ", outputOnEos+"/{}.root".format(datasetName))
            #clean out /tmp/
            for f in os.listdir(tempDir):
                os.remove(tempDir+"/"+f)
            if os.path.isdir(tempDir):
                os.rmdir(tempDir)
    except Exception as e:
        #print("ERROR: exception in ProcessDataset for txtFile={}, mass={}".format(txtFile, mass), flush=True)
        traceback.print_exc()
        # raise e
        raise RuntimeError("Caught exception in ProcessDataset for txtFile={}, mass={}".format(txtFile, mass))
    return True


def GetBackgroundDatasetsDict(inputListBkgBase):
    if not inputListBkgBase.endswith("/"):
        inputListBkgBase += "/"
    backgroundDatasetsDict = {
            "ZJet_amcatnlo_ptBinned" : [
                # inclusive stitched  empty
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                ],
 #           "ZJet_powhegminnlo" : [
 #               inputListBkgBase+"DYJetsToEE_M-1000to1500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-100to200_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-1500to2000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-2000toInf_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-200to400_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-400to500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-500to700_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-700to800_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               inputListBkgBase+"DYJetsToEE_M-800to1000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
 #               ],
            "ZJet_HT_LO"  :  [
                 inputListBkgBase+"DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_HT-70to100_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8.txt",
                 inputListBkgBase+"DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
                 ],
            "TTbar_powheg" : [
                inputListBkgBase+"TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.txt",
                #inputListBkgBase+"TTToHadronic_TuneCP5_13TeV-powheg-pythia8.txt",
                #inputListBkgBase+"TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8.txt",
                ],
            "DIBOSON_nlo" : [
                # inputListBkgBase+"WWTo4Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"ZZTo4L_TuneCP5_13TeV_powheg_pythia8.txt",
                # inputListBkgBase+"ZZTo2Nu2Q_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8.txt",
                inputListBkgBase+"WZTo1L3Nu_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                #inputListBkgBase+"WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WZTo3LNu_mllmin4p0_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                inputListBkgBase+"WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
                ],
            # "TRIBOSON" : [
            #     inputListBkgBase+"WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"WZZ_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     inputListBkgBase+"ZZZ_TuneCP5_13TeV-amcatnlo-pythia8.txt",
            #     ],
            # "TTW" : [
            #     inputListBkgBase+"TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.txt",
            #     inputListBkgBase+"TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8.txt",
            #     ],
            # "TTZ" : [
            #     inputListBkgBase+"ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8.txt"
            #     ],
            "SingleTop" : [
                inputListBkgBase+"ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ST_t-channel_top_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ST_t-channel_antitop_5f_InclusiveDecays_TuneCP5_13TeV-powheg-pythia8.txt",
                inputListBkgBase+"ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8.txt",
                ],
            # "WJet_amcatnlo_jetBinned" : [
            #     inputListBkgBase+"WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            #     ],
            # "WJet_amcatnlo_ptBinned" : [
            #     inputListBkgBase+"WJetsToLNu_Wpt-0To50_ext1_amcatnloFXFX_pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_Wpt-50To100_ext1_amcatnloFXFX_pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_Pt-100To250_amcatnloFXFX_pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_Pt-250To400_amcatnloFXFX_pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_Pt-400To600_amcatnloFXFX_pythia8.txt",
            #     inputListBkgBase+"WJetsToLNu_Pt-600ToInf_amcatnloFXFX_pythia8.txt",
            #     ],
            # "PhotonJets_Madgraph" : [
            #     #inputListBkgBase+"GJets_HT-40To100_madgraphMLM.txt",
            #     inputListBkgBase+"GJets_DR-0p4_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #     inputListBkgBase+"GJets_DR-0p4_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #     inputListBkgBase+"GJets_DR-0p4_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #     inputListBkgBase+"GJets_DR-0p4_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8.txt",
            #     ],
    }
    return backgroundDatasetsDict

def GetQCDDatasetList(year):
    year = str(year)
    qcdFakes = {
        "QCDFakes_DYJ" : [
            inputListQCD1FRBase+"DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-0To50_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            inputListQCD1FRBase+"DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8.txt",
            ],
       # "QCDFakes_DYJ" : [
       #     inputListQCD1FRBase+"DYJetsToEE_M-1000to1500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-100to200_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-10to50_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-1500to2000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-2000toInf_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-200to400_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-400to500_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-500to700_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-50_massWgtFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-700to800_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     inputListQCD1FRBase+"DYJetsToEE_M-800to1000_H2ErratumFix_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos.txt",
       #     ]
    }
    if year == "2016preVFP":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"SinglePhoton_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v4.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
            inputListQCD2FRBase+"SinglePhoton_Run2016B-ver2_HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016C-HIPM_UL2016_MiniAODv2_NanoAODv9-v4.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016D-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016E-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            inputListQCD2FRBase+"SinglePhoton_Run2016F-HIPM_UL2016_MiniAODv2_NanoAODv9-v2.txt",
            ]
        qcdFakes["QCDFakes_DYJ"] = [txtFile.replace(".txt", "_APV.txt") for txtFile in qcdFakes["QCDFakes_DYJ"]]
    elif year == "2016postVFP":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"SinglePhoton_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
                    inputListQCD2FRBase+"SinglePhoton_Run2016H-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2016G-UL2016_MiniAODv2_NanoAODv9-v2.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2016F-UL2016_MiniAODv2_NanoAODv9-v1.txt",
                    ]
    elif year == "2017":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"SinglePhoton_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"SinglePhoton_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
                    inputListQCD2FRBase+"SinglePhoton_Run2017B-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2017C-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2017D-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2017E-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"SinglePhoton_Run2017F-UL2017_MiniAODv2_NanoAODv9-v1.txt",
                    ]
    elif year == "2018":
        qcdFakes["QCDFakes_DATA"] = [
                    inputListQCD1FRBase+"EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD1FRBase+"EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3.txt",
                    ]
        qcdFakes["QCDFakes_DATA_2FR"] = [
                    inputListQCD2FRBase+"EGamma_Run2018A-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"EGamma_Run2018B-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"EGamma_Run2018C-UL2018_MiniAODv2_NanoAODv9-v1.txt",
                    inputListQCD2FRBase+"EGamma_Run2018D-UL2018_MiniAODv2_NanoAODv9-v3.txt",
                    ]
    else:
        print("ERROR: did not find QCD fakes datasets for year {}".format(year))

    return qcdFakes


mycut = TCut("M_e1e2 > 220 && sT_eejj > 400 && Pt_e1e2 > 70 && Ele1_Pt > 50 && Ele2_Pt > 50 && Jet1_Pt > 50 && Jet2_Pt > 50 && PassTrigger==1")


branchesToSave = [
        "run",
        "ls",
        "event",
        "Weight",
        "PrefireWeight",
        "puWeight",
        "EventTriggerScaleFactor",
        # "ZVtxSF",
        # "Ele1_RecoSF",
        # "Ele2_RecoSF",
        # "Ele1_TrigSF",
        # "Ele2_TrigSF",
        # "Ele1_EGMLooseIDSF",
        # "Ele2_EGMLooseIDSF",
        # "Ele1_HEEPSF",
        # "Ele2_HEEPSF",
        # "FakeRateEffective",
        # "MinPrescale",
        "PassTrigger",
        "EventWeight",
        "sT_eejj",
        "Pt_e1e2",
        "Ele1_Pt",
        "Ele2_Pt",
        "Ele1_Eta",
        "Ele2_Eta",
        "Ele1_Phi",
        "Ele2_Phi",
        "Jet1_Pt",
        "Jet2_Pt",
        "Jet3_Pt",
        "Jet1_Eta",
        "Jet2_Eta",
        "Jet3_Eta",
        "Jet1_Phi",
        "Jet2_Phi",
        "Jet3_Phi",
        "DR_Ele1Jet1",
        "DR_Ele1Jet2",
        "DR_Ele2Jet1",
        "DR_Ele2Jet2",
        "DR_Jet1Jet2",
        "PFMET_Type1_Pt",
        "PFMET_Type1_Phi",
        "M_e1e2",
        "M_e1j1",
        "M_e1j2",
        "M_e2j1",
        "M_e2j2",
        "Masym",
        "MejMin",
        "MejMax",
        "Meejj",
        "Meejjj",
        "LQCandidateMass"
]
qcdSamplesToSub = ["QCDFakes_DATA_2FR", "QCDFakes_DYJ"]
####################################################################################################
# Run
####################################################################################################
ROOT.DisableImplicitMT()
parallelize = True
date = "9oct2023"
includeQCD = True
year = sys.argv[1]
inputListBkgBase = "$LQANA/config/myDatasets/BDT/"+str(year)+"/16sep/trainingTreeInputs/preselOnly"
inputListQCD1FRBase = "$LQANA/config/myDatasets/BDT/"+str(year)+"/16sep/trainingTreeInputs/singleFR/"
inputListQCD2FRBase = "$LQANA/config/myDatasets/BDT/"+str(year)+"/16sep/trainingTreeInputs/doubleFR/"
inputListSignalBase = inputListBkgBase
outputTFileDir = os.getenv("LQDATAEOS")+"/BDTTrainingTrees/LQToBEle/"+str(year)+"/16SepSkims"
#signalNameTemplate = "LQToDEle_M-{}_pair_bMassZero_TuneCP2_13TeV-madgraph-pythia8"
signalNameTemplate = "LQToBEle_M-{}_pair_TuneCP2_13TeV-madgraph-pythia8"
if __name__ == "__main__":
    print("INFO: Using year = {}".format(year))
    print("INFO: Using signal name template: {}".format(signalNameTemplate))
    print("INFO: Using inputListBkgBase: {}".format(inputListBkgBase))
    if includeQCD:
        print("INFO: Using inputListQCD1FRBase: {}".format(inputListQCD1FRBase))
        print("INFO: Using inputListQCD2FRBase: {}".format(inputListQCD2FRBase))
    print("INFO: Using inputListSignalBase: {}".format(inputListSignalBase))
    print("INFO: Saving training trees to: {}".format(outputTFileDir))
    sys.stdout.flush()
    console = Console()
    allSignalDatasetsDict = {}
    backgroundDatasetsDict = GetBackgroundDatasetsDict(inputListBkgBase)
    if year == "2016preVFP":
        for dataset in backgroundDatasetsDict.keys():
            #if "QCDFakes_DATA" in dataset:
            #    continue
            backgroundDatasetsDict[dataset] = [txtFile.replace(".txt", "_APV.txt") for txtFile in backgroundDatasetsDict[dataset]]
        #backgroundDatasetsDict.update(qcdFakes2016pre)
        signalNameTemplate+="_APV"
    #elif year == "2016postVFP":
    #    backgroundDatasetsDict.update(qcdFakes2016post)
    if includeQCD:
        backgroundDatasetsDict.update(GetQCDDatasetList(year))
    massList = list(range(300, 3100, 100))
    massList.extend([3500, 4000])
    #massList = [1200]
    for mass in massList:
        signalName = signalNameTemplate.format(mass)
        if not inputListSignalBase.endswith("/"):
            inputListSignalBase += "/"
        allSignalDatasetsDict[signalName] = [inputListSignalBase+signalName+".txt"]

    if parallelize:
        jobCount = 0
        #for signalSample in allSignalDatasetsDict.keys():
        #    for idx, txtFile in enumerate(allSignalDatasetsDict[signalSample]):
        #        jobCount += 1
        #        for bkgSample in backgroundDatasetsDict.keys():
        #            for bkgTxtFile in backgroundDatasetsDict[bkgSample]:
        #                jobCount += 1
        # ncores = multiprocessing.cpu_count()
        ncores = 4  # only use 4 parallel jobs to be nice
        pool = multiprocessing.Pool(ncores)
        for signalSample in allSignalDatasetsDict.keys():
            for idx, txtFile in enumerate(allSignalDatasetsDict[signalSample]):
                mass = int(signalSample[signalSample.find("M")+2:signalSample.find("_", signalSample.find("M"))])
                #tfilepath = ProcessDataset(txtFile, mass, branchesToSave)
                try:
                    pool.apply_async(ProcessDataset, [[txtFile, mass, branchesToSave, year, "", False]], callback=log_result, error_callback=handle_error)
                    #pool.apply_async(ProcessDataset, [[txtFile, mass, branchesToSave, console]], callback=lambda x: progress.advance(task_id))
                    jobCount += 1
                except KeyboardInterrupt:
                    print("\n\nCtrl-C detected: Bailing.")
                    pool.terminate()
                    exit(-1)
                except Exception as e:
                    print("ERROR: caught exception in job for LQ mass: {}; exiting".format(mass))
                    traceback.print_exc()
                    exit(-2)
                for bkgSample in backgroundDatasetsDict.keys():
                    for bkgTxtFile in backgroundDatasetsDict[bkgSample]:
                        try:
                            makeWeightsNegative = True if bkgSample in qcdSamplesToSub else False
                            qcdLabel = bkgSample if "qcd" in bkgSample.lower() else ""
                            pool.apply_async(ProcessDataset, [[bkgTxtFile, mass, branchesToSave, year, qcdLabel, makeWeightsNegative]], callback=log_result, error_callback=handle_error)
                            #pool.apply_async(ProcessDataset, [[bkgTxtFile, mass, branchesToSave, console]], callback=lambda x: progress.advance(task_id))
                            jobCount += 1
                        except KeyboardInterrupt:
                            print("\n\nCtrl-C detected: Bailing.")
                            pool.terminate()
                            exit(-1)
                        except Exception as e:
                            print("ERROR: caught exception in job for LQ mass: {}; exiting".format(mass))
                            traceback.print_exc()
                            exit(-2)
                        # ProcessDataset(bkgTxtFile, mass, branchesToSave)
        # now close the pool and wait for jobs to finish
        pool.close()
        #sys.stdout.write(logString.format(jobCount, len(massList)))
        #sys.stdout.write("\t"+str(len(result_list))+" jobs done")
        #sys.stdout.flush()
        pool.join()
        # check results?
        if len(result_list) < jobCount:
            print("ERROR: {} jobs had errors. Exiting.".format(jobCount-len(result_list)))
            exit(-2)
