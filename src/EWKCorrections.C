#include "EWKCorrections.h"
#include <vector>
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
      return ewkHist->GetBinContent(ewkHist->FindBin(GenZPt));
    }

}
