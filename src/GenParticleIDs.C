#include <algorithm>
#include <cmath>
#include <iostream>

#include "GenParticle.h"
#include "Collection.h"

bool GenParticle::PassUserID (ID id, bool verbose){ 
  if      ( id == GEN_ELE_FROM_LQ         ) { return PassUserID_GenEleFromLQ     (verbose); }
  else if ( id == GEN_MUON_FROM_LQ        ) { return PassUserID_GenMuonFromLQ    (verbose); }
  else if ( id == GEN_TAU_FROM_LQ         ) { return PassUserID_GenTauFromLQ     (verbose); }
  else if ( id == GEN_ELE_HARD_SCATTER    ) { return PassUserID_GenEleHardScatter(verbose); }
  else if ( id == GEN_NU_HARD_SCATTER     ) { return PassUserID_GenNuHardScatter (verbose); }

  else if ( id == GEN_MU_HARD_SCATTER     ) { return PassUserID_GenMuHardScatter    (verbose); }
  else if ( id == GEN_QUARK_HARD_SCATTER  ) { return PassUserID_GenQuarkHardScatter (verbose); }
  else if ( id == GEN_QUARK_HARD_PROCESS  ) { return PassUserID_GenQuarkHardProcess (verbose); }
  else if ( id == GEN_Z_FROMHARDPROCESS_LASTCOPY ) { return PassUserID_GenZFromHardProcessLastCopy(verbose); }
  //else if ( id == GEN_W_HARD_SCATTER      ) { return PassUserID_GenWHardScatter     (verbose); }
  else if ( id == GEN_NU_FROM_W  	        ) { return PassUserID_GenNuFromW          (verbose); }
  //else if ( id == GEN_ELE_FROM_W  	      ) { return PassUserID_GenEleFromW         (verbose); }
  else if ( id == GEN_FROM_LQ  	          ) { return PassUserID_FromLQ          (verbose); }
  else if ( id == GEN_ELE_HARDPROCESS_FINALSTATE  	      ) { return PassUserID_GenEleHardProcessFinalState(verbose); }

  else if ( id == GEN_ELE_FIDUCIAL        ) { return PassUserID_ECALFiducial        (verbose); } 

  else if ( id == GEN_LQ                  ) { return PassUserID_GenLQ               (verbose); } 
  else if ( id == GEN_TOP                 ) { return PassUserID_GenTop              (verbose); }

  else if ( id == GEN_STATUS62            ) { return PassUserID_Status62            (verbose); }
  else if ( id == GEN_IS_LAST_COPY        ) { return PassUserID_GenIsLastCopy       (verbose); }

  else return false;
}

bool GenParticle::PassUserID_GenEleHardScatter(bool verbose){ 
  // pythia 8: outgoing hard electron status is 23 (still intermediate)
  // technically, this is not the really final status=1 particle
  //if ( !(Status() == 3 || Status() == 23) ) return false;  
  if ( IsHardProcess() && abs(PdgId()) == 11 ) return true;
  return false;
}

bool GenParticle::PassUserID_GenNuHardScatter(bool verbose){ 
  if ( IsHardProcess() && abs(PdgId()) == 12 ) return true;
  return false;
}

bool GenParticle::PassUserID_GenMuHardScatter(bool verbose){ 
  if ( IsHardProcess() && abs(PdgId()) == 13 ) return true;
  return false;
}

bool GenParticle::PassUserID_GenQuarkHardScatter(bool verbose) {
  if ( IsHardProcess() && abs(PdgId()) < 9 && abs(PdgId()) > 0) return true;
  return false;
}

bool GenParticle::PassUserID_GenQuarkHardProcess(bool verbose) {
  if ( IsFromHardProcess() && abs(PdgId()) < 9 && abs(PdgId()) > 0) return true;
  return false;
}

bool GenParticle::PassUserID_FromLQ(bool verbose){
  // pythia 8: outgoing hard electron status is 23 (still intermediate)
  // technically, this is not the really final status=1 particle
  //XXX replace with IsHardProcess() ? testing needed
  //if ( !(Status() == 3 || Status() == 23) ) return false;  
  //
  //std::cout << Name() << " " << ": "
	// << "PDG = "    << PdgId () << ", "
	// << "Status = " << Status () << ", "
	// << "Pt = "     << Pt ()    << ", "
	// << "Eta = "    << Eta()    << ", "
	// << "Phi = "    << Phi();// << std::endl;
  //std::cout << ", mother index = " << MotherIndex() << std::endl;
  //
  //int mpdgId = m_collection->GetData()->GenParticlePdgId->at(MotherIndex());
  //int mstatus = m_collection->GetData()->GenParticleStatus->at(MotherIndex());
  //float mpt = m_collection->GetData()->GenParticlePt->at(MotherIndex());
  //float meta = m_collection->GetData()->GenParticleEta->at(MotherIndex());
  //float mphi = m_collection->GetData()->GenParticlePhi->at(MotherIndex());
  //std::cout << "MOTHER: PDG = "    << mpdgId << ", "
	// << "Status = " << mstatus << ", "
	// << "Pt = "     << mpt    << ", "
	// << "Eta = "    << meta    << ", "
	// << "Phi = "    << mphi << std::endl;
  //
  int motherIndex = MotherIndex();
  while(motherIndex > 0) {
    int mother_pdg_id = m_collection->ReadArrayBranch<Int_t>("GenPart_pdgId", motherIndex);
    if ( abs(mother_pdg_id) == 42 || abs(mother_pdg_id) == 9000007) {
      m_MotherLQIndex = motherIndex;
      return true;
    }
    motherIndex = m_collection->ReadArrayBranch<Int_t>("GenPart_genPartIdxMother", motherIndex);
  }
  return false;
}

bool GenParticle::PassUserID_FromDY(bool verbose){
  if ( !IsHardProcess() && !(Status()==3) ) return false;
  if ( MotherIndex()<0) return false; // can't tell if it's from W if mother index not set
  int mother_pdg_id = MotherIndex() >= 0 ? m_collection->ReadArrayBranch<Int_t>("GenPart_pdgId",MotherIndex()) : -1;
  if ( abs(mother_pdg_id) != 22 && 
       abs(mother_pdg_id) != 23 ) return false;
  return true;
}

bool GenParticle::PassUserID_FromW(bool verbose){
  if ( !IsHardProcess() && !(Status()==3) ) return false;
  if ( MotherIndex()<0) return false; // can't tell if it's from W if mother index not set
  if(verbose) {
    std::cout << "GenParticle::PassUserID_FromW" << Name() << " " << ": "
      << "PDG = "    << PdgId () << ", "
      << "Status = " << Status () << ", "
      << "Pt = "     << Pt ()    << ", "
      << "Eta = "    << Eta()    << ", "
      << "Phi = "    << Phi();// << std::endl;
    std::cout << ", mother index = " << MotherIndex() << std::endl;
    std::cout << "GenParticle::PassUserID_FromW() MotherIndex=" << MotherIndex() << std::endl;
  }
  int mother_pdg_id = MotherIndex() >= 0 ? m_collection->ReadArrayBranch<Int_t>("GenPart_pdgId",MotherIndex()) : -1;
  if ( abs(mother_pdg_id) != 24 ) return false;
  return true;
}

bool GenParticle::PassUserID_GenEleFromLQ (bool verbose){
  if ( abs(PdgId())  != 11         ) return false;
  if ( !PassUserID_FromLQ(verbose) ) return false;
  return true;
}

bool GenParticle::PassUserID_GenMuonFromLQ(bool verbose){
  if ( abs(PdgId())  != 13         ) return false;
  if ( !PassUserID_FromLQ(verbose) ) return false;
  return true;
}

bool GenParticle::PassUserID_GenTauFromLQ (bool verbose){
  if ( abs(PdgId())  != 15         ) return false;
  if ( !PassUserID_FromLQ(verbose) ) return false;
  return true;
}

//bool GenParticle::PassUserID_GenZGammaHardScatter(bool verbose){
//  if ( abs(PdgId()) != 22 && 
//       abs(PdgId()) != 23 ) return false;
//  if(verbose) {
//    std::cout << "GenParticle::PassUserID_GenZGammaHardScatter" << Name() << " " << ": "
//      << "PDG = "    << PdgId () << ", "
//      << "Status = " << Status () << ", "
//      << "Pt = "     << Pt ()    << ", "
//      << "Eta = "    << Eta()    << ", "
//      << "Phi = "    << Phi() << std::endl;
//  }
//  if ( IsHardProcess() || Status() == 3) return true;
//  
//  //if ( Status()>20 && Status()<30) return true; //pythia8
//  //if (Status() == 62) return true;  
//  return false;
//}

 
//bool GenParticle::PassUserID_GenWHardScatter     (bool verbose){
//  if ( abs(PdgId()) != 24 ) return false;
//  if(verbose) {
//    std::cout << "GenParticle::PassUserID_GenWHardScatter" << Name() << " " << ": "
//      << "PDG = "    << PdgId () << ", "
//      << "Status = " << Status () << ", "
//      << "Pt = "     << Pt ()    << ", "
//      << "Eta = "    << Eta()    << ", "
//      << "Phi = "    << Phi() << std::endl;
//  }
//  // now that we have a W
//  if ( IsHardProcess() || Status() == 3) return true;
//  return false;
//} 

bool GenParticle::PassUserID_GenZFromHardProcessLastCopy(bool verbose){
  if ( abs(PdgId()) != 23 ) return false;
  if(verbose) {
    std::cout << "GenParticle::PassUserID_GenZFromHardProcessLastCopy" << Name() << " " << ": "
      << "PDG = "    << PdgId () << ", "
      << "Status = " << Status () << ", "
      << "Pt = "     << Pt ()    << ", "
      << "Eta = "    << Eta()    << ", "
      << "Phi = "    << Phi()    << ", "
      << "Mass = "   << Mass()    << ", "
      << "IsFromHardProcess = " << IsFromHardProcess() << ", "
      << "IsLastCopy = " << IsLastCopy() << std::endl;
  }
  if (IsFromHardProcess() && IsLastCopy()) return true;
  return false;
}


bool GenParticle::PassUserID_GenNuFromW          (bool verbose){
  if ( abs(PdgId()) != 12         ) return false;
  if ( !PassUserID_FromW(verbose) ) return false;
  return true;
} 

bool GenParticle::PassUserID_GenEleFromW         (bool verbose){
  if ( abs(PdgId()) != 11         ) return false;
  if ( !PassUserID_FromW(verbose) ) return false;
  return true;
} 

bool GenParticle::PassUserID_GenEleFromDY        (bool verbose){
  if ( abs(PdgId())  != 11         ) return false;
  if ( !PassUserID_FromDY(verbose) ) return false;
  return true;
} 

bool GenParticle::PassUserID_ECALFiducial(bool verbose){
  if ( IsGenElectronFiducial() ) return true;
  else return false;
}

bool GenParticle::PassUserID_MuonFiducial(bool verbose){
  if ( IsMuonFiducial() ) return true;
  else return false;
}

bool GenParticle::PassUserID_GenLQ(bool verbose)
{
  if ( abs(PdgId())  == 42         ) return true;
  if ( abs(PdgId())  == 9000007    ) return true; // for LQ toolbox [https://arxiv.org/abs/1801.07641]; ID depends on model. this is S3m43
  return false;
}

bool GenParticle::PassUserID_GenTop(bool verbose)
{
  if ( abs(PdgId())  != 6         ) return false;
  return true;
}

bool GenParticle::PassUserID_Status62(bool verbose)
{
  if ( Status()  != 62         ) return false;
  return true;
}

bool GenParticle::PassUserID_GenEleHardProcessFinalState(bool verbose)
{
  if ( abs(PdgId()) != 11         ) return false;
  if ( !IsFromHardProcessFinalState() ) return false;
  return true;
}

bool GenParticle::PassUserID_GenIsLastCopy(bool verbose)
{
  if ( IsLastCopy() ) return true;
  return false;
}
