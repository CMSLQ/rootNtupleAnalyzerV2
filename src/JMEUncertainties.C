#include "Collection.h"
#include "PFJet.h"
#include "JMESystematicsCalculators.h"
#include <variant>
#include <filesystem>

using namespace ROOT;
using METSystematicsResult = std::map<std::string, double>;

template<typename S>
class JMEUncertainties {

  public:
    using JetMETArgType = std::variant<RVecF, RVecI, float, int>;
    using JetVariationsResult = JetVariationsCalculator::result_t;
    using METVariationsResult = Type1METVariationsCalculator::result_t;

    JMEUncertainties() = default;

    template<typename T = JetVariationsCalculator>
    JMEUncertainties(analysisClass* d,
        const std::string& jecJSON, const std::string& jetAlgo, const std::string& jecTag,
        const std::string& jecLevel,
        const std::vector<std::string>& jesUncertainties, bool addHEM2018Issue,
        const std::string& jerTag, const std::string& jsonFileSmearingTool,
        const std::string& smearingToolName,
        bool splitJER=false, bool doGenMatch=true,
        float genMatchMaxDR=0.2, float genMatchMaxDPt = 3.0) :
      //requires(std::is_same_v<S, JetVariationsCalculator>)  : 
      m_data(d),
      m_splitJER (splitJER), // smearing
      m_doGenMatch (doGenMatch), // smearing
      m_genMatchMaxDR (genMatchMaxDR), // smearing
      m_genMatchMaxDPt (genMatchMaxDPt) // smearing
  {
    m_calc = JetVariationsCalculator::create(jecJSON, jetAlgo, jecTag, jecLevel, jesUncertainties, addHEM2018Issue, jerTag, jsonFileSmearingTool, smearingToolName, splitJER, doGenMatch, genMatchMaxDR, genMatchMaxDPt);
  }

    template<typename T = Type1METVariationsCalculator>
      JMEUncertainties(analysisClass* d,
          const std::string& jecJSON, const std::string& jetAlgo, const std::string& jecTag,
          const std::string& jecLevel,
          const std::string& l1Jec,
          float unclEnThreshold,
          float emEnFracThreshold,
          const std::vector<std::string>& jesUncertainties, bool addHEM2018Issue,
          bool isT1SmearedMET,
          const std::string& jerTag, const std::string& jsonFileSmearingTool,
          const std::string& smearingToolName,
          bool splitJER=false, bool doGenMatch=true,
          float genMatchMaxDR=0.2, float genMatchMaxDPt = 3.0) :
        m_data(d),
        m_splitJER (splitJER), // smearing
        m_doGenMatch (doGenMatch), // smearing
        m_genMatchMaxDR (genMatchMaxDR), // smearing
        m_genMatchMaxDPt (genMatchMaxDPt) // smearing
  {
    m_calc = Type1METVariationsCalculator::create(jecJSON, jetAlgo, jecTag, jecLevel, l1Jec, unclEnThreshold, emEnFracThreshold, jesUncertainties, addHEM2018Issue, isT1SmearedMET, jerTag, jsonFileSmearingTool, smearingToolName, splitJER, doGenMatch, genMatchMaxDR, genMatchMaxDPt);
  }

    void ComputeJetVariations(CollectionPtr jetCollection, CollectionPtr genJetCollection = nullptr,
        bool isMC=false, bool verbose=false) {
      JetVariationsResult result = produceJetVariations(jetCollection, genJetCollection, isMC);
      std::vector<std::string> variations = available();
      std::vector<RVecF> systematics;
      std::vector<RVecF> massSystematics;
      std::vector<std::string> systematicsNames;
      std::vector<unsigned short> rawIndices = jetCollection->GetRawIndices();
      sort(rawIndices.begin(), rawIndices.end());
      for(unsigned short varIdx = 0; varIdx < result.size(); ++varIdx) {
        RVecF calculatedPtVariations = result.pt(varIdx);
        RVecF calculatedMassVariations = result.mass(varIdx);
        RVecF ptVariationByJet;
        RVecF massVariationByJet;
        unsigned short vecIdx = 0;
        if(!rawIndices.size())
          break;
        unsigned short lastRawIdxIdx = 0;
        for(unsigned short idx = 0; idx <= rawIndices[rawIndices.size()-1]; ++idx) {
          bool foundElement = idx == rawIndices[lastRawIdxIdx];
          float pt = foundElement ? calculatedPtVariations[lastRawIdxIdx] : 0.0;
          ptVariationByJet.emplace_back(pt);
          float mass = foundElement ? calculatedMassVariations[lastRawIdxIdx] : 0.0;
          massVariationByJet.emplace_back(mass);
          if(foundElement)
            lastRawIdxIdx++;
        }
        systematics.push_back(ptVariationByJet);
        massSystematics.push_back(massVariationByJet);
      }
      systematics.insert(systematics.end(), massSystematics.begin(), massSystematics.end());
      for(unsigned short i = 0; i < result.size(); ++i) {
        systematicsNames.push_back("Pt_"+variations[i]);
      }
      for(unsigned short i = 0; i < result.size(); ++i) {
        systematicsNames.push_back("mass_"+variations[i]);
      }
      jetCollection->SetSystematics(std::move(systematicsNames), std::move(systematics));
    }
    
    METSystematicsResult ComputeType1METVariations(CollectionPtr jetCollection, CollectionPtr genJetCollection = nullptr,
        bool isMC=false, bool verbose=false) {
      METVariationsResult result = produceType1METVariations(jetCollection, genJetCollection, isMC);
      std::vector<std::string> variations = available();
      METSystematicsResult returnVal;
      for(unsigned short i = 0; i < result.size(); ++i)
        returnVal["Pt_"+variations[i]] = result.pt(i);
      for(unsigned short i = 0; i < result.size(); ++i) 
        returnVal["Phi_"+variations[i]] = result.phi(i);
      return returnVal;
    }

    std::vector<std::string> available() {
      return std::visit(
          [](auto& calc)->std::vector<std::string> {
          return calc.available();
          }, m_calc);
    }

    void setAddHEM2018Issue(bool addHEM2018) {
      std::visit(
          [&addHEM2018](auto& calc) {
          calc.setAddHEM2018Issue(addHEM2018);
          }, m_calc);
    }

  private:
    std::string getValueTypeName() const {
      return std::visit( [](auto&&x)->decltype(auto){ return typeid(x).name(); }, m_calc);
    }
    template <typename T = float> T GetValueFromBranchName(const std::string& branchName) {
      return m_data->readerTools_->ReadValueBranch<T>(branchName);
    }
    template <typename T = float> RVec<T> GetRVecFromBranchName(const std::string& branchName) {
      return RVec<T>(static_cast<T*>(m_data->readerTools_->ReadArrayBranch<T>(branchName).GetAddress()),
          m_data->readerTools_->ReadArrayBranch<T>(branchName).GetSize());
    }
    template <typename T = float> RVec<T> GetRVecFromJetCollection(CollectionPtr jetCollection, const std::string& varName) {
      RVec<T> result;
      std::vector<unsigned short> indices = jetCollection->GetRawIndices();
      for(const auto& idx : indices) {
        result.emplace_back(jetCollection->ReadArrayBranch<T>(varName, idx));
      }
      return result;
    }
    JetVariationsResult produceJetVariations(CollectionPtr jetCollection, CollectionPtr genJetCollection, bool isMC=false) {
      bool forMET=false;
      bool addHEM2018Issue=false;
      std::vector<JetMETArgType> jetArgs = GetJetMETArgs(jetCollection, isMC, forMET, addHEM2018Issue, genJetCollection);
      return std::get<JetVariationsCalculator>(m_calc).produce(std::get<RVecF >(jetArgs[0]), std::get<RVecF >(jetArgs[1]), std::get<RVecF >(jetArgs[2]),
          std::get<RVecF >(jetArgs[3]),std::get<RVecF >(jetArgs[4]),std::get<RVecF >(jetArgs[5]),
          std::get<RVecI >(jetArgs[6]),
          std::get<float>(jetArgs[7]),
          std::get<RVecI >(jetArgs[8]),std::get<RVecI >(jetArgs[9]),
          std::get<int>(jetArgs[10]),
          std::get<RVecF >(jetArgs[11]), std::get<RVecF >(jetArgs[12]), std::get<RVecF >(jetArgs[13]), std::get<RVecF >(jetArgs[14])
          );
    }

    METVariationsResult produceType1METVariations(CollectionPtr jetCollection, CollectionPtr genJetCollection, bool isMC=false) {
      //std::cout << "INFO: produceType1METVariations! " << jetCollection->GetSize() << " jets, " << genJetCollection->GetSize() << " genJets." << std::endl;
      bool forMET=true;
      bool addHEM2018Issue=false;
      std::vector<JetMETArgType> metArgs = GetJetMETArgs(jetCollection, isMC, forMET, addHEM2018Issue, genJetCollection);
      return std::get<Type1METVariationsCalculator>(m_calc).produce(std::get<RVecF >(metArgs[0]), std::get<RVecF >(metArgs[1]), std::get<RVecF >(metArgs[2]),
          std::get<RVecF >(metArgs[3]),std::get<RVecF >(metArgs[4]),std::get<RVecF >(metArgs[5]),
          std::get<RVecF >(metArgs[6]),
          std::get<RVecF >(metArgs[7]),std::get<RVecF >(metArgs[8]),
          std::get<RVecI >(metArgs[9]),
          std::get<float>(metArgs[10]),
          std::get<RVecI >(metArgs[11]),
          std::get<RVecI >(metArgs[12]),
          std::get<int>(metArgs[13]),
          std::get<RVecF >(metArgs[14]), std::get<RVecF >(metArgs[15]), std::get<RVecF >(metArgs[16]), std::get<RVecF >(metArgs[17]),
          std::get<float>(metArgs[18]),std::get<float>(metArgs[19]),
          std::get<RVecF >(metArgs[20]), std::get<RVecF >(metArgs[21]), std::get<RVecF >(metArgs[22]), std::get<RVecF >(metArgs[23]),
          std::get<RVecF >(metArgs[24]), std::get<RVecF >(metArgs[25]), std::get<RVecF >(metArgs[26]), std::get<float>(metArgs[27]),std::get<float>(metArgs[28])
          );
    }

  private:
    // reimplemented in C++ from https://gitlab.cern.ch/cp3-cms/CMSJMECalculators/-/blob/main/python/CMSJMECalculators/utils.py#L33
    std::vector<JetMETArgType> GetJetMETArgs(CollectionPtr jetCollection,
        bool isMC=true, bool forMET=false, bool addHEM2018Issue=false,
        CollectionPtr genJetCollection = nullptr) {
      std::vector<JetMETArgType> args;
      UInt_t nJet = jetCollection->GetSize();
      args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_pt"));
      RVecF jetEta = GetRVecFromJetCollection(jetCollection, "Jet_eta");
      args.push_back(jetEta);
      args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_phi"));
      args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_mass"));
      args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_rawFactor"));
      args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_area"));
      if(forMET) {
        args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_muonSubtrFactor"));
        args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_neEmEF"));
        args.push_back(GetRVecFromJetCollection(jetCollection, "Jet_chEmEF"));
      }
      if(addHEM2018Issue)
        args.push_back(GetRVecFromJetCollection<int>(jetCollection, "Jet_jetId"));
      else
        args.push_back(RVecI());
      args.push_back(GetValueFromBranchName("fixedGridRhoFastjetAll"));
      if(isMC) {
        args.push_back(GetRVecFromJetCollection<int>(jetCollection, "Jet_genJetIdx"));
        args.push_back(GetRVecFromJetCollection<int>(jetCollection, "Jet_partonFlavour"));
        //UInt_t run = GetValueFromBranchName<UInt_t>("run");
        //UInt_t luminosityBlock = GetValueFromBranchName<UInt_t>("luminosityBlock");
        ULong64_t event = GetValueFromBranchName<ULong64_t>("event");
        //args.push_back(static_cast<std::uint32_t>((run<<20) + (luminosityBlock<<10) + event + 1 + (nJet != 0 ? int(jetEta[0]/.01) : 0)));
        args.push_back(static_cast<int>(event));
        args.push_back(GetRVecFromJetCollection(genJetCollection, "GenJet_pt"));
        args.push_back(GetRVecFromJetCollection(genJetCollection, "GenJet_eta"));
        args.push_back(GetRVecFromJetCollection(genJetCollection, "GenJet_phi"));
        args.push_back(GetRVecFromJetCollection(genJetCollection, "GenJet_mass"));
      }
      else {
        args.push_back(RVecI());
        args.push_back(RVecI());
        args.push_back(0);
        args.push_back(RVecF());
        args.push_back(RVecF());
        args.push_back(RVecF());
        args.push_back(RVecF());
      }

      if(forMET) {
        args.push_back(GetValueFromBranchName("RawMET_phi"));
        args.push_back(GetValueFromBranchName("RawMET_pt"));
        args.push_back(GetRVecFromBranchName("CorrT1METJet_rawPt"));
        args.push_back(GetRVecFromBranchName("CorrT1METJet_eta"));
        args.push_back(GetRVecFromBranchName("CorrT1METJet_phi"));
        args.push_back(GetRVecFromBranchName("CorrT1METJet_area"));
        args.push_back(GetRVecFromBranchName("CorrT1METJet_muonSubtrFactor"));
        args.push_back(RVecF());
        args.push_back(RVecF());
        args.push_back(GetValueFromBranchName("MET_MetUnclustEnUpDeltaX"));
        args.push_back(GetValueFromBranchName("MET_MetUnclustEnUpDeltaY"));
      }
      return args;
    }

    void CheckFileExists(std::string& filePath) {
      if(!std::filesystem::exists(filePath)) {
        STDOUT("ERROR: JES/JER file " << filePath << " does not exist! Can't compute uncertainties. Exiting.");
        exit(-8);
      }
    }

    std::variant<JetVariationsCalculator, Type1METVariationsCalculator> m_calc;
    analysisClass* m_data;

    std::string m_jerTextFilePath, m_ptResTxtFile, m_sfTxtFile;
    bool m_splitJER;
    bool m_doGenMatch;
    float m_genMatchMaxDR;
    float m_genMatchMaxDPt;
};

namespace JMEUtils {


};
