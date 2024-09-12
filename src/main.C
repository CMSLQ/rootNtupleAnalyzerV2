#include "analysisClass.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h> 
#include <iomanip>

#include <TTree.h>
#include <TChain.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TVector2.h>
#include <TVector3.h>

using namespace std;

int main(int argc, char* argv[])
{
  const int Nparam=4;   // NUMBER OF PARAMETERS

  if(argc<Nparam+1 || argc>Nparam+2)
    {
      cout << "main() : arcg = " << argc << ", while we expected " << Nparam+1 << " or " << Nparam+2 << " arguments. Exiting." <<endl;
      cout << "Usage  : ./main inputList cutFile treeName outputRootFileWithoutExtension [outputEfficiencyFileWithoutExtension]" << endl;
      cout << "Example: ./main config/inputListExample.txt config/cutFileExample.txt RootTupleMaker data/output/rootFile data/output/cutEfficiencyFile" << endl;
      cout << "Example: ./main config/inputListExample.txt config/cutFileExample.txt treeCreator/RootTupleMakerPAT data/output/rootFile data/output/cutEfficiencyFile" << endl;
      cout << "Example: ./main config/inputListExample.txt config/cutFileExample.txt rootTupleTree/tree  data/output/rootFile data/output/cutEfficiencyFile" << endl;
      exit (1);
    };

  std::unique_ptr<string> inputList      = std::make_unique<string>(argv[1]);
  std::unique_ptr<string> cutFile        = std::make_unique<string>(argv[2]);
  std::unique_ptr<string> treeName       = std::make_unique<string>(argv[3]);
  std::unique_ptr<string> outputFileName = std::make_unique<string>(argv[4]);
  std::unique_ptr<string> cutEfficFile = std::make_unique<string>(*outputFileName.get());
  if(argc==6)
    cutEfficFile   = std::make_unique<string>(argv[5]);

  analysisClass analysisClass_(inputList.get(), cutFile.get(), treeName.get(), outputFileName.get(), cutEfficFile.get());
  analysisClass_.Loop();

}

