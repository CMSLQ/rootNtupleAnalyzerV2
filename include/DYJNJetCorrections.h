#include <vector>
#include <TH1D.h>

namespace DYJNJetCorrections {

     const std::vector<float> dyjnjetBinEdges = {1.5, 2.5, 3.5, 4.5, 5.5, 100.0};
     const std::vector<float> dyjnjetBinContent2016pre = {0.878685, 1.135846, 1.454879, 2.026025, 3.597194};
     const std::vector<float> dyjnjetBinContent2016post = {0.937645, 1.150977, 1.626356, 2.105785, 2.748908};
     const std::vector<float> dyjnjetBinContent2017 = {0.971482, 1.203601, 1.586480, 2.066689, 2.751360};
     const std::vector<float> dyjnjetBinContent2018 = {0.934070, 1.147625, 1.560774, 2.124003, 2.780355};
     inline std::unique_ptr<TH1D> dyjnjetHist = nullptr;

     void initHist(const std::string& year);
     double LookupNJetCorrection(const std::string& year, const int nJets);
}
