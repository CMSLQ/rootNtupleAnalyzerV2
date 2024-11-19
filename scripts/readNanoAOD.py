#!/usr/bin/env python3
import ROOT
import copy


#@ROOT.Numba.Declare(["unsigned int","unsigned int","unsigned long"], "bool")
def IsMatchingEvent(run, ls, event):
    for selEvt in selectedEvents:
        selEvtSplit = selEvt.split(":")
        selRun = int(selEvtSplit[0])
        selLs = int(selEvtSplit[1])
        selEvt = int(selEvtSplit[2])
        if run!=selRun:
            continue
        if ls!=selLs:
            continue
        if event!=selEvt:
            continue
        return True
    return False


ROOT.gInterpreter.Declare('''
       float CalcMee(float Ele1_Pt, float Ele1_Eta, float Ele1_Phi, float Ele2_Pt, float Ele2_Eta, float Ele2_Phi) {
         TLorentzVector e1, e2, ee;
         e1.SetPtEtaPhiM ( Ele1_Pt, Ele1_Eta, Ele1_Phi, 0.0 );
         e2.SetPtEtaPhiM ( Ele2_Pt, Ele2_Eta, Ele2_Phi, 0.0 );
         ee = e1 + e2; 
         return ee.M();
       }
''')

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

lorentz = ROOT.Math.LorentzVector("ROOT::Math::PtEtaPhiM4D<double>")
maxEvents = -1
#selectedEvents = ["275658:197:387100603", "275658:197:386939841", "275658:197:387721840", "275658:197:387386104", "275658:197:386009644", "275658:197:386328252", "275658:197:387394030", "275658:197:387609551", "275658:197:386583344", ]
selectedEvents = []

if maxEvents > 0:
    print("Run over", maxEvents, "events only.")

filesInclusive = [
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/555795A5-1818-AF43-8312-2FB955755C10.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/40000/78C9AA75-CE75-5345-9511-6B399F5B350A.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/47F26205-1E71-5148-A980-30676F8262D6.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/51C95A21-3DBA-2145-B06F-3EB33AF34CEF.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/3ACA7CA4-5327-154C-928D-F866F3C7F7CF.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/932B2BF1-08E2-3D46-A804-E988B75727A7.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/8B0AC3DF-B6A9-4147-B4EA-BE63E001F264.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/25C63D54-F9E0-124A-8B29-8D942975DD07.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/3824692C-DDD6-1946-8D61-20A2C1463819.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/3316AD85-024D-984C-92A0-B37E2F1B8199.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/BA416FBA-E485-334C-86C2-05ECC297B048.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/FEA553C1-A192-144B-A511-A26C47C6819F.root",  # contains run 275658, ls=197
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/E031FEFA-EAD3-BB40-960E-2172F21F1538.root",
    # "root://xrootd-cms.infn.it//store/data/Run2016C/SingleElectron/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v2/2520000/45BCD4F7-7E67-3F42-85EE-7FF4CFF4671D.root",
    "root://xrootd-cms.infn.it//store/data/Run2016C/SinglePhoton/NANOAOD/HIPM_UL2016_MiniAODv2_NanoAODv9-v4/2530000/574BF391-DB30-3744-AC33-EA56ECEB518B.root"
]
tree = ROOT.TChain("Events")
for tfile in filesInclusive:
    tree.Add(tfile)
hltBranchNames = [branch.GetName() for branch in tree.GetListOfBranches() if "HLT" in branch.GetName()]
nevents = tree.GetEntries()
for iev, event in enumerate(tree):
    if maxEvents > 0 and iev >= maxEvents:
        break
    if iev < 10 or iev % 10000 == 0:
        print("event {}/{}".format(iev, nevents))
    if len(selectedEvents) > 0:
        if not IsMatchingEvent(event.run, event.luminosityBlock, event.event):
            continue
    # selections
    if event.nElectron <= 0 or event.Electron_pt[0] < 35 or event.HLT_Photon175:
        continue
    if not event.HLT_Photon120:
        continue
    if event.Electron_pt[0] > 120:
        continue
    #
    print("Run = {}, ls = {}, event = {}".format(event.run, event.luminosityBlock, event.event))
    nIDEles = 0
    eleTotalPt = 0
    p4 = lorentz()
    for i in range(event.nElectron):
        #if event.Electron_cutBased[i] > 1 and nIDEles < 2:
        #    p4 += lorentz(event.Electron_pt[i], event.Electron_eta[i], event.Electron_phi[i], event.Electron_mass[i])
        #    print("\tEle{} Pt = {}, Eta = {}, Phi = {}".format(i, event.Electron_pt[i], event.Electron_eta[i], event.Electron_phi[i]))
        #    eleTotalPt+=event.Electron_pt[i]
        #    nIDEles+=1
        print("\tEle{} Pt = {}, Eta = {}, Phi = {}".format(i, event.Electron_pt[i], event.Electron_eta[i], event.Electron_phi[i]))
    for i in range(event.nPhoton):
        print("\tPhoton{} Pt = {}, Eta = {}, Phi = {}".format(i, event.Photon_pt[i], event.Photon_eta[i], event.Photon_phi[i]))
    jetTotalPt = 0
    nIDJets = 0
    #for i in range(event.nJet):
    #    if event.Jet_jetId > 0 and nIDJets < 2:
    #        print("\tJet{} Pt = {}, Eta = {}, Phi = {}".format(i, event.Jet_pt[i], event.Jet_eta[i], event.Jet_phi[i]))
    #        jetTotalPt+=event.Jet_pt[i]
    #        nIDJets+=1
    #print("\tMee = {}".format(p4.M()))
    #print("\tsT = {}".format(eleTotalPt+jetTotalPt))
    for i in range(event.nTrigObj):
        print("\tTrigObj{} Pt = {}, Eta = {}, Phi = {}, id = {}, filterBits = {}".format(i, event.TrigObj_pt[i], event.TrigObj_eta[i], event.TrigObj_phi[i], event.TrigObj_id[i], event.TrigObj_filterBits[i]))
    firedTriggers = []
    for i in hltBranchNames:
        if tree.GetBranch(i).GetLeaf(i).GetValue():
            firedTriggers.append(i)
    print("Fired triggers: {}".format(firedTriggers))

# df = RDataFrame(tree)
# df = df.Define("selectedEvent", "Numba::IsMatchingEvent(run, luminosityBlock, event)")
# df = df.Filter("selectedEvent==True")
# df = df.Define("ElectronPassLooseID", "Electron_cutBased > 1")
# df = df.Define("Mee", "CalcMee(ElectronPassLooseID_Pt, ElectronPassLooseID_Eta, ElectronPassLooseID_Phi, ElectronPassLooseID_Pt, ElectronPassLooseID_Eta, ElectronPassLooseID_Phi)")
# cols = ROOT.vector('string')()
# cols.push_back("run")
# cols.push_back("luminosityBlock")
# cols.push_back("event")
# cols.extend([""])
# d2 = myDf.Display(cols)
# d2.Print()
