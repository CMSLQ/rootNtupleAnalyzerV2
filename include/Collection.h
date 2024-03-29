#ifndef COLLECTION_H
#define COLLECTION_H

#include <memory>
#include <vector>
#include <algorithm>
#include <iostream>
#include <ostream> 
#include <typeinfo>
#include <TRandom3.h>
#include <TLorentzVector.h>
#include <ROOT/RVec.hxx>

#include "IDTypes.h"
#include "analysisClass.h"

#include "correction.h"

template<class Object1>
  void examineVec ( std::vector<Object1> objVec, const char * name, ID id = NULL_ID, bool verbose = false ){
  int n_constituents = objVec.size();
  std::cout << "N(" << name << ") = " << n_constituents << std::endl;
  for (int i = 0; i < n_constituents; ++i ){ 
    Object1 constituent = objVec.at(i);
    std::cout << "\t" << "Constituent" << "\t#" << i << ":" << "\t" << constituent << std::flush;
    if      ( id == NULL_ID  ) std::cout << std::endl;
    else {
      if ( verbose        ) { 
        constituent.PassUserID ( id, verbose );
      }
      else {
        std::cout << ", ID = " << constituent.PassUserID ( id ) << std::endl;
      }
    }
  }
}


class Collection;
using CollectionPtr = std::shared_ptr<Collection>;


class Collection {
 public:
  
  //-------------------------------------------------------------
  // Constructors and destructors
  //-------------------------------------------------------------
  
  Collection ( std::shared_ptr<TTreeReaderTools> tools);
  Collection ( std::shared_ptr<TTreeReaderTools> tools, unsigned short size );
  Collection ( Collection & c );
  
  //-------------------------------------------------------------
  // Set functions
  //-------------------------------------------------------------
  
  void SetTriggerObjectIndex ( short i ) { m_trigObj_index = i; } 
  void SetRawIndices( std::vector<unsigned short> i ) { m_raw_indices.swap(i); }
  void SetLeadNConstituents ( unsigned short n ) {
    m_raw_indices.clear();
    for (unsigned short i = 0; i < n; ++i) 
      m_raw_indices.push_back ( i );
  }
  
  //-------------------------------------------------------------
  // Get functions
  //-------------------------------------------------------------
  
  template<class Object1> Object1 GetConstituent(unsigned short i) { 
    if ( m_trigObj_index > 0 )    return Object1 (*this, m_raw_indices[i], m_trigObj_index); 
    else                          return Object1 (*this, m_raw_indices[i]);
  }

  void RemoveConstituent(unsigned short i) {
    m_raw_indices.erase(m_raw_indices.begin()+i);
  }
  template<class Object1> void RemoveConstituent(Object1& obj) { 
    auto it = std::find ( m_raw_indices.begin(), m_raw_indices.end(), obj.GetRawIndex() ) ;
    if(it != m_raw_indices.end())
      RemoveConstituent(std::distance(m_raw_indices.begin(), it));
  }

  
  std::vector<unsigned short>   GetRawIndices() { return  m_raw_indices; } 
  unsigned short                GetSize()       { return  m_raw_indices.size();  } 
  template <typename T> T ReadValueBranch(const std::string& branchName) {
    return m_readerTools->ReadValueBranch<T>(branchName);
  }
  template <typename T> T ReadArrayBranch(const std::string& branchName, unsigned int idx) {
    return m_readerTools->ReadArrayBranch<T>(branchName, idx);
  }

  bool HasBranch(const std::string& branchName) {
    return (m_readerTools->GetTree()->GetBranch(branchName.c_str())) != 0;
  }
  
  //-------------------------------------------------------------
  // Modify collection indices
  //-------------------------------------------------------------
  
  void Clear() { m_raw_indices.clear(); }
  void Append ( unsigned short i ) { m_raw_indices.push_back ( i ); } 

  template <class Object1> 
    bool Has ( const Object1 & o ) { 
    std::vector<unsigned short>::iterator it = std::find ( m_raw_indices.begin(), m_raw_indices.end(), o.GetRawIndex() ) ;
    bool not_found = ( it == m_raw_indices.end() );
    return (!not_found);
  }

  template <class Object1, class Object2> 
    int HasHowMany ( const CollectionPtr other_collection ){
    std::vector<unsigned short> other_collection_raw_indices = other_collection -> GetRawIndices();
    std::vector<unsigned short> common_raw_indices;
    // std::sort ( m_raw_indices.begin(), m_raw_indices.end() );
    // std::sort ( other_collection_raw_indices.begin(), other_collection_raw_indices.end() );
    std::set_intersection ( other_collection_raw_indices.begin(), other_collection_raw_indices.end(),
			    m_raw_indices.begin()               , m_raw_indices.end(),
			    std::back_inserter ( common_raw_indices ) );
    return common_raw_indices.size();
  }

  float GetSystematicValue(const unsigned short rawIdx, const std::string& systName) {
    unsigned short variationIdx = FindSystematicVariationIndex(systName);
    return m_systematicVariations[variationIdx][rawIdx];
  }

  unsigned short FindSystematicVariationIndex(const std::string& systName) {
    auto variationIdxItr = std::find(m_systematicVariationNames.begin(), m_systematicVariationNames.end(), systName);
    if(variationIdxItr == m_systematicVariationNames.end()) {
      STDOUT("ERROR: specified systName = " << systName << " does not exist in the stored systematic variations:");
      for(auto syst : m_systematicVariationNames)
        std::cout << syst << ", ";
      std::cout << endl;
      exit(-8);
    }
    return std::distance(m_systematicVariationNames.begin(), variationIdxItr);
  }

  void SetSystematics(std::vector<std::string>&& systNames, std::vector<ROOT::VecOps::RVec<float> >&& systematics) {
    m_systematicVariationNames = std::move(systNames);
    m_systematicVariations = std::move(systematics);
  }
  void SetSystematics(std::vector<std::string>& systNames, std::vector<ROOT::VecOps::RVec<float> >& systematics) {
    m_systematicVariationNames = systNames;
    m_systematicVariations = systematics;
  }

  std::vector<std::string> GetSystematicsNames() {
    return m_systematicVariationNames;
  }


  std::vector<ROOT::VecOps::RVec<float> > GetSystematics() {
    return m_systematicVariations;
  }
  //------------------------------------------------------------------------------------------
  // For skimming
  //------------------------------------------------------------------------------------------

  // For skimming by ID

  template<class Object1>
    CollectionPtr SkimByID( ID id, bool verbose = false ) { 
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    new_collection -> Clear();
    unsigned short size = GetSize();
    for (unsigned short i = 0; i < size ; ++i){
      Object1 constituent = GetConstituent<Object1> (i);
      if ( constituent.PassUserID (id, verbose) ) 
        new_collection -> Append ( constituent.GetRawIndex() );
    }
    return new_collection;
  }

  // For skimming by minimum pt 
  
  template<class Object1>
    CollectionPtr SkimByMinPt ( double min_pt ) { 
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    unsigned short size = GetSize();
    for (unsigned short i = 0; i < size ; ++i){
      Object1 constituent = GetConstituent<Object1> (i);
      if ( constituent.Pt() >= min_pt ) 
        new_collection -> Append ( constituent.GetRawIndex() );
    }
    return new_collection;
  }

  // For skimming by minimum ptHeep, for electrons
  
  template<class Object1>
    CollectionPtr SkimByMinPtHeep ( double min_pt ) { 
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    unsigned short size = GetSize();
    for (unsigned short i = 0; i < size ; ++i){
      Object1 constituent = GetConstituent<Object1> (i);
      if ( constituent.PtHeep() >= min_pt ) 
        new_collection -> Append ( constituent.GetRawIndex() );
    }
    return new_collection;
  }

  // For skimming by eta range
  
  template<class Object1>
    CollectionPtr SkimByEtaRange ( double min_eta, double max_eta ) { 
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    unsigned short size = GetSize();
    for (unsigned short i = 0; i < size ; ++i){
      Object1 constituent = GetConstituent<Object1> (i);
      if ( constituent.Eta() >= min_eta && constituent.Eta() <= max_eta ) 
        new_collection -> Append ( constituent.GetRawIndex() );
    }
    return new_collection;
  }
  
  //------------------------------------------------------------------------------------------
  // Skim by minimum delta R from objects in another collection
  // 
  // If all objects in another collection have deltaR from this 
  //    object greater than some minimum: keep the object
  // If any one object in another collection has deltaR from this 
  //    object less than some minimum: throw out the object
  // 
  // Use this to protect against overlaps
  //------------------------------------------------------------------------------------------
  
  template < class Object1, class Object2 > 
    CollectionPtr SkimByVetoDRMatch ( const CollectionPtr other_collection, double min_dr ){
    unsigned short this_collection_size = GetSize();
    unsigned short other_collection_size = other_collection -> GetSize();
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      double tmp_min_dr = 999.0;
      Object1 this_collection_constituent  = GetConstituent<Object1>(i);
      for ( unsigned short j = 0; j < other_collection_size; ++j ){
        Object2 other_collection_constituent = other_collection -> GetConstituent<Object2> (j);
        double dr = this_collection_constituent.DeltaR ( &other_collection_constituent );
        if ( dr < tmp_min_dr ) tmp_min_dr = dr;
      }
      if ( tmp_min_dr >= min_dr ) { 
        new_collection -> Append ( this_collection_constituent.GetRawIndex() );
      }
      //else
      //  std::cout << "SkimByVetoDRMatch(): too close to an object in other_collection; remove object1: " << this_collection_constituent << std::endl;
    }
    return new_collection;
  }
  
  template < class Object1, class Object2 > 
    CollectionPtr SkimByVetoDRMatch ( Object2 & other_object, double min_dr ){
    unsigned short this_collection_size = GetSize();
    CollectionPtr new_collection ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      Object1 this_collection_constituent  = GetConstituent<Object1>(i);
      double dr = this_collection_constituent.DeltaR ( & other_object );
      if ( dr >= min_dr ) new_collection -> Append ( this_collection_constituent.GetRawIndex() );
    }
    return new_collection;
  }
  
  //------------------------------------------------------------------------------------------
  // Skim by maximum delta R from objects in another collection
  // 
  // If all objects in another collection have deltaR from this object
  //    greater than some maximum: throw out the object
  // If any one object in another collection has deltaR from this object
  //    less than some maximum: keep the object
  //
  // Use this for deltaR matching one group of objects to another group
  //------------------------------------------------------------------------------------------
  
  template <class Object1, class Object2> 
    CollectionPtr SkimByRequireDRMatch ( const CollectionPtr other_collection, double max_dr ){
    unsigned short this_collection_size = GetSize();
    unsigned short other_collection_size = other_collection -> GetSize();
    CollectionPtr new_collection  ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      double tmp_min_dr = 999.0;
      Object1 this_collection_constituent  = GetConstituent<Object1>(i);
      for ( unsigned short j = 0; j < other_collection_size; ++j ){
        Object2 other_collection_constituent = other_collection -> GetConstituent<Object2>(j);
        double dr = this_collection_constituent.DeltaR ( & other_collection_constituent );
        if ( dr < tmp_min_dr ) tmp_min_dr = dr;
      }
      if ( tmp_min_dr <= max_dr ) new_collection -> Append ( this_collection_constituent.GetRawIndex() );
    }
    return new_collection;
  }

  //------------------------------------------------------------------------------------------
  // Get closest object in DR
  //------------------------------------------------------------------------------------------
  
  template <class Object1, class Object2 >
    Object1 GetClosestInDR ( Object2 & other_object ){ 
    Object1 retval;
    double tmp_min_dr = 999.0;
    unsigned short this_collection_size = GetSize();
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      Object1 this_collection_constituent  = GetConstituent<Object1>(i);
      double dr = this_collection_constituent.DeltaR ( & other_object );
      if ( dr < tmp_min_dr ) {
        tmp_min_dr = dr;
        retval = this_collection_constituent;
      }
    }
    return retval;
  }

  //------------------------------------------------------------------------------------------
  // Remove duplicates (DR = 0)
  //------------------------------------------------------------------------------------------

  template <class Object1> 
    CollectionPtr RemoveDuplicates () { 
    CollectionPtr new_collection  ( new Collection(m_readerTools, 0));
    new_collection -> SetSystematics(GetSystematicsNames(), GetSystematics());
    new_collection -> SetTriggerObjectIndex ( m_trigObj_index );
    unsigned short this_collection_size = GetSize();
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      Object1 this_collection_constituent  = GetConstituent<Object1>(i);
      bool is_duplicate = false;
      for (unsigned short j = 0; j < i; ++j){
	Object1 previous_collection_constituent  = GetConstituent<Object1>(j);
	double dr = this_collection_constituent.DeltaR ( & previous_collection_constituent );
	if ( dr < 1e-4 ) is_duplicate = true;
      }
      if ( !is_duplicate ) new_collection -> Append ( this_collection_constituent.GetRawIndex() );
    }
    return new_collection;
  }

  //------------------------------------------------------------------------------------------
  // Match and smear energy (e.g. for JER/EER studies)
  //------------------------------------------------------------------------------------------

  template <class Object1, class Object2>
    std::vector<Object1> MatchAndSmearEnergy ( const CollectionPtr matching_collection, double max_dr, TRandom3 * engine, TLorentzVector & v_delta_met,
        const correction::Correction* correctionPtRes = nullptr, const correction::Correction* correctionJERSF = nullptr, std::string variation = "nom",
        bool verbose=false) {
    unsigned short this_collection_size = GetSize();
    std::vector<Object1> smearedObjVec;
    smearedObjVec.reserve(this_collection_size);
    TLorentzVector v_old, v_new, v_delta;
    double smearFactor = 1;
    for (unsigned short i = 0; i < this_collection_size ; ++i) {
      Object1 this_collection_constituent = GetConstituent<Object1>(i);
      if(verbose)
        std::cout << "MatchAndSmearEnergy(): Obj " << this_collection_constituent.Name() << " constituent #" << i << "/" << this_collection_size << ":" << std::endl;
      Object2 matched_object;
      double resolution = correctionPtRes ? this_collection_constituent.EnergyResFromCorrection(correctionPtRes) : this_collection_constituent.EnergyRes();
      double res_scale_factor   = correctionJERSF ? this_collection_constituent.EnergyResScaleFactorFromCorrection(correctionJERSF, variation) : this_collection_constituent.EnergyResScaleFactor();
      double old_pt = this_collection_constituent.Pt();
      bool matched = this_collection_constituent.template MatchByDRAndDPt < Object2 > ( matching_collection, matched_object, max_dr, 3*resolution*old_pt );
      double new_pt = -1.0;
      if ( matched ) { 
        double matched_pt           = matched_object.Pt();

        //double scale_error          = this_collection_constituent.EnergyResScaleError ();
        //double smearing             = engine -> Gaus ( 1.0 , scale_error );
        //double smeared_scale_factor = scale_factor * smearing;
        //double delta_pt             = smeared_scale_factor * ( old_pt - matched_pt );
        //new_pt               = std::max ( double (0.0), matched_pt + delta_pt ) ;
        double delta_pt = old_pt - matched_pt;
        smearFactor = 1 + (res_scale_factor - 1.) * delta_pt / old_pt;
        
        //if(std::string(this_collection_constituent.Name())=="PFJet") {
        if(verbose) {
          std::cout << "\tMatched obj " << this_collection_constituent.Name() << " constituent #" << i << " with pt = " << old_pt << " GeV and Eta= " <<
            this_collection_constituent.Eta() << " to " << matched_object.Name() << " with pt = " << matched_pt << " GeV" << std::endl;
          std::cout << "\t" << "old RECO pt          = " << old_pt               << std::endl;
          std::cout << "\t" << "GEN pt               = " << matched_pt           << std::endl;
          std::cout << "\t" << "res. scale factor    = " << res_scale_factor     << std::endl;
          //std::cout << "\t" << "smearing             = " << smearing             << std::endl;
          std::cout << "\t" << "smeared scale factor = " << smearFactor          << std::endl;
          std::cout << "\t" << "delta pt             = " << delta_pt             << std::endl;
          std::cout << "\t" << "max DPt for matching = " << 3*resolution*old_pt  << std::endl;
          //std::cout << "\t" << "new RECO pt          = " << new_pt               << std::endl;
        }
        
      }
      else {
        // not well-matched to GenParticle
        double scale_factor   = correctionJERSF ? this_collection_constituent.EnergyResScaleFactorFromCorrection(correctionJERSF, variation) : this_collection_constituent.EnergyResScaleFactor();
        if(scale_factor < 0)
          scale_factor = 0;
        double sigma = resolution;
        smearFactor = 1. + engine->Gaus(0.0, sigma) * std::sqrt(std::max(scale_factor * scale_factor - 1, 0.0));
        //std::cout << "Not well-matched jet" << std::endl;
        //std::cout << "\t" << "old RECO pt          = " << old_pt               << std::endl;
        //std::cout << "\t" << "scale factor         = " << scale_factor         << std::endl;
        //std::cout << "\t" << "sigma                = " << sigma             << std::endl;
        //std::cout << "\t" << "smeared scale factor = " << smearFactor << std::endl;
        //v_new.SetPtEtaPhiM( old_pt, this_collection_constituent.Eta(), this_collection_constituent.Phi(), 0.0 );
        //v_new*=smearFactor;
        ////std::cout << "\t" << "delta pt             = " << delta_pt             << std::endl;
        //std::cout << "\t" << "new RECO pt          = " << v_new.Pt()               << std::endl;
      }

      if( smearFactor < 0) smearFactor = 0;
      v_old.SetPtEtaPhiM( old_pt, this_collection_constituent.Eta(), this_collection_constituent.Phi(), 0.0 );
      v_new.SetPtEtaPhiM( old_pt, this_collection_constituent.Eta(), this_collection_constituent.Phi(), 0.0 );
      v_new*=smearFactor;
      new_pt = v_new.Pt();
      v_delta = v_old - v_new;
      v_delta_met = v_delta_met + v_delta;

      this_collection_constituent.SetPt(v_new.Pt());
      this_collection_constituent.SetEta(v_new.Eta());
      this_collection_constituent.SetPhi(v_new.Phi());
      if(verbose) {
        std::cout << "\tsmeared object has Pt=" << this_collection_constituent.Pt() << ", Eta=" << this_collection_constituent.Eta() <<
          ", Phi=" << this_collection_constituent.Phi() << std::endl;
      }
      if ( new_pt >= 1e-6 )
        smearedObjVec.push_back(this_collection_constituent);

      if(verbose) {
        if(!matched) {
          std::cout << "\t" << "no match found! collection looks like:" << std::endl;
          matching_collection->examine<Object2>("matching collection");
        }
        std::cout << "\t" << "smeared obj vec size = " << smearedObjVec.size() << std::endl;
      }
    }

    return smearedObjVec;
  }

  //------------------------------------------------------------------------------------------
  // Scale energy 
  //------------------------------------------------------------------------------------------

  template <class Object1> 
    std::vector<Object1> ScaleEnergy ( int scale_sign, TLorentzVector & v_delta_met ){
    unsigned short this_collection_size = GetSize();
    std::vector<Object1> scaledObjVec;
    scaledObjVec.reserve(this_collection_size);
    TLorentzVector v_old, v_new, v_delta;
    for (unsigned short i = 0; i < this_collection_size ; ++i) {    
      Object1 this_collection_constituent = GetConstituent<Object1>(i);
      
      double old_pt       = this_collection_constituent.Pt();
      double scale        = this_collection_constituent.EnergyScaleFactor();
      double double_sign  = double (scale_sign);
      double scale_factor = 1.0 + ( double_sign * scale );
      double new_pt       = std::max ( double (0.0), old_pt * scale_factor);

      v_old.SetPtEtaPhiM( old_pt, this_collection_constituent.Eta(), this_collection_constituent.Phi(), 0.0 );
      v_new.SetPtEtaPhiM( new_pt, this_collection_constituent.Eta(), this_collection_constituent.Phi(), 0.0 );
      v_delta = v_old - v_new;
      v_delta_met = v_delta_met + v_delta;
      
      //std::cout << "Old pt = "       << old_pt       << ", "
      //		<< "scale = "        << scale        << ", "
      //		<< "sign = "         << scale_sign   << ", "
      //		<< "scale factor = " << scale_factor << ", "
      //		<< "new pt = "       << new_pt       << ", " << std::endl;
      //std::cout << "\tOld: "    << old_pt       << ", " << this_collection_constituent.Eta() << ", " << this_collection_constituent.Phi() << std::endl;
      //std::cout << "\tNew: "    << new_pt       << ", " << this_collection_constituent.Eta() << ", " << this_collection_constituent.Phi() << std::endl;
      //std::cout << "\tOld - New: " << v_delta.Pt() << ", " << v_delta.Eta()                     << ", " << v_delta.Phi() << std::endl;
      //std::cout << "\tDelta(MET) = " << v_delta_met.Pt() << ", " << v_delta_met.Eta()                     << ", " << v_delta_met.Phi() << std::endl;

      this_collection_constituent.SetPt(new_pt);
      if ( new_pt >= 1e-6 )
        scaledObjVec.push_back(this_collection_constituent);

    }
    
    return scaledObjVec;
  }

  //------------------------------------------------------------------------------------------
  // For debugging
  //------------------------------------------------------------------------------------------
  
  template<class Object1>
    void examine ( const char * name, ID id = NULL_ID, bool verbose = false ){
      unsigned short int n_constituents = GetSize();
      std::cout << "N(" << name << ") = " << n_constituents << std::endl;
      for (int i = 0; i < n_constituents; ++i ){ 
        Object1 constituent = GetConstituent<Object1>(i);
        std::cout << "\t" << "Constituent" << "\t#" << i << ":" << "\t" << constituent << std::flush;
        if ( id != NULL_ID  ) {
          if ( verbose        ) { 
            constituent.PassUserID ( id, verbose );
          }
          else {
            std::cout << ", ID = " << constituent.PassUserID ( id );
          }
        }
        std::cout << ", rawIndex=" << constituent.GetRawIndex() << std::endl;
      }
    }

  //------------------------------------------------------------------------------------------
  // Member variables
  //------------------------------------------------------------------------------------------
  
 protected:
  short m_trigObj_index;
  std::vector<unsigned short> m_raw_indices; 
  std::shared_ptr<TTreeReaderTools> m_readerTools;
  std::vector<ROOT::VecOps::RVec<float> > m_systematicVariations;
  std::vector<std::string> m_systematicVariationNames;
  
};

#endif 

