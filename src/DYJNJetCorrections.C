#include "DYJNJetCorrections.h"
#include <vector>
#include <stdexcept>
#include <TH1D.h>

namespace DYJNJetCorrections {

     void initHist(const std::string& year) {
       dyjnjetHist.reset(new TH1D("dyjnjet", "dyjnjet", dyjnjetBinEdges.size()-1, dyjnjetBinEdges.data()));
      dyjnjetHist->SetDirectory(0);
      for(int xbin=1; xbin <= dyjnjetHist->GetNbinsX(); ++xbin) {
       if(year == "2016pre")
        dyjnjetHist->SetBinContent(xbin, dyjnjetBinContent2016pre[xbin-1]);
      }
    }

     double LookupNJetCorrection(const std::string& year, const int nJets) {
      if(dyjnjetHist==nullptr)
        initHist(year);
      int jetBin = dyjnjetHist->FindBin(nJets);
      if(jetBin > dyjnjetHist->GetNbinsX())
        jetBin = dyjnjetHist->GetNbinsX();
      return dyjnjetHist->GetBinContent(jetBin);
    }

}
