{
  TChain ch("rootTupleTree/tree");
  TFileCollection fc("dum","","/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/nanoV7_2016_pskQCDEEJJ_egLoose_24mar2022_comb/SinglePhoton_Run2016B-02Apr2020_ver2-v1.txt");
  ch.AddFileInfoList(fc.GetList());

  TMVA::Experimental::RReader BDT("/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/versionsOfAnalysis/2016/nanoV7/eejj/mar24_2022_egLoose_BDT_qcd/train_14-17/dataset/weights/TMVAClassification_LQToDEle_M-1400_pair_BDTG.weights.xml");
  computeBDT = TMVA::Experimental::Compute<28, float>(BDT);

  float cutValForIntegral = 0.986;
  ROOT::RDataFrame df(ch);
  TCut mycutb = TCut("M_e1e2 > 200 && sT_eejj > 400");
  auto df1 = df.Filter(mycutb.GetTitle());
  auto df2 = df1.Define("BDTv", computeBDT, BDT.GetVariableNames());
  auto df3 = df2.Define("BDT", "BDTv[0]");
  auto df4 = df3.Define("eventWeight", "Weight");
  std::string histName = "BDTVal";
  TH1D hbkg (histName.c_str(), histName.c_str(), 10000, -1, 1);
  histBkg = df4.Histo1D(ROOT::RDF::TH1DModel(hbkg), "BDT", "eventWeight");
  auto bkgIntegralOverCut = df4.Filter("BDT > 0.986").Sum("eventWeight").GetValue();
  auto bkgEntriesOverCut = df4.Filter("BDT > 0.986").Count().GetValue();

  std::cout << "entries with BDT > 0.986 = " << bkgEntriesOverCut << ", integral unweighted = " << bkgIntegralOverCut << std::endl;
  std::vector<std::string> cols;
  cols.push_back("run");
  cols.push_back("ls");
  cols.push_back("event");
  cols.push_back("BDT");
  cols.push_back("eventWeight");
  cols.push_back("sT_eejj");
  cols.push_back("PFMET_Type1_Pt");
  cols.push_back("PFMET_Type1_Phi");
  cols.push_back("M_e1e2");
  cols.push_back("M_e1j1");
  cols.push_back("M_e1j2");
  cols.push_back("M_e2j1");
  cols.push_back("M_e2j2");
  cols.push_back("Ele1_Pt");
  cols.push_back("Ele2_Pt");
  cols.push_back("Ele1_Eta");
  cols.push_back("Ele2_Eta");
  cols.push_back("Ele1_Phi");
  cols.push_back("Ele2_Phi");
  cols.push_back("Jet1_Pt");
  cols.push_back("Jet2_Pt");
  cols.push_back("Jet3_Pt");
  cols.push_back("Jet1_Eta");
  cols.push_back("Jet2_Eta");
  cols.push_back("Jet3_Eta");
  cols.push_back("Jet1_Phi");
  cols.push_back("Jet2_Phi");
  cols.push_back("Jet3_Phi");
  cols.push_back("DR_Ele1Jet1");
  cols.push_back("DR_Ele1Jet2");
  cols.push_back("DR_Ele2Jet1");
  cols.push_back("DR_Ele2Jet2");
  cols.push_back("DR_Jet1Jet2");
  auto display = df4.Filter("BDT > 0.986").Display(cols, 20);
  //display->Print();
  std::cout << display->AsString() << std::endl;
}
