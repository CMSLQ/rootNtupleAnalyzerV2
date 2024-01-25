#include "EWKCorrections.h"
#include <vector>
#include <stdexcept>
#include <TH1D.h>

namespace EWKCorrections {

     void initHist() {
      ewkHist.reset(new TH1D("ewk", "ewk", ewkBinEdges.size()-1, ewkBinEdges.data()));
      ewkHist->SetDirectory(0);
      for(int xbin=1; xbin <= ewkHist->GetNbinsX(); ++xbin)
        ewkHist->SetBinContent(xbin, ewkBinContent[xbin-1]);
    }

     double LookupEWKCorrection(const double GenZPt) {
      if(ewkHist==nullptr)
        initHist();
      int ptBin = ewkHist->FindBin(GenZPt);
      if(ptBin > ewkHist->GetNbinsX())
        throw std::runtime_error("Found overflow bin while getting EWKCorrection for GenZPt=" + std::to_string(GenZPt));
      return ewkHist->GetBinContent(ptBin);
    }

}
