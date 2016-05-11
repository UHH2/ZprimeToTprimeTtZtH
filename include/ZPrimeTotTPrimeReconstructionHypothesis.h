#pragma once

#include "UHH2/core/include/Particle.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Jet.h"
#include <map>


class ZPrimeTotTPrimeReconstructionHypothesis {
 public:
  explicit ZPrimeTotTPrimeReconstructionHypothesis(){ m_tophad_topjet_ptr = 0; }
  LorentzVector toplep_v4() const{return m_toplep_v4;}
  LorentzVector tophad_v4() const{return m_tophad_v4;} 
  LorentzVector HZW_v4() const{return m_HZW_v4;} 
  LorentzVector neutrino_v4() const{return m_neutrino_v4;} 
  Particle lepton() const{return m_lepton;}
  // TopJet HZW() const{return m_HZW;}
  const std::vector<Jet>& HZW_subjets()const{return m_subjets;}
  const std::vector<Jet>& toplep_jets() const{return m_toplep_jets;}
  const std::vector<Jet>& tophad_jets() const{return m_tophad_jets;}
  //const std::vector<Jet>& HZW_jets() const{return m_HZW_jets;}
  // const std::vector<TopJet>& TOPtoplep_jets() const{return m_TOPtoplep_jets;}
  // const std::vector<TopJet>&  tophad_jets() const{return m_tophad_jets;}
  const std::vector<TopJet>&  HZW_jets() const{return m_HZW_jets;}

  const TopJet* tophad_topjet_ptr() const{return m_tophad_topjet_ptr;}
  LorentzVector top_v4() const{ return m_lepton.charge() > 0 ? m_toplep_v4 : m_tophad_v4;}
  LorentzVector antitop_v4() const{ return m_lepton.charge() < 0 ? m_toplep_v4 : m_tophad_v4;}
  LorentzVector wlep_v4() const{ return m_neutrino_v4+m_lepton.v4();}
  LorentzVector blep_v4() const{return m_blep_v4;}

  /// get the discriminator value for this hypothesis; thows a runtime_error if it does not exist.
  float discriminator(const std::string & l) const {
    auto it = m_discriminators.find(l);
    if(it == m_discriminators.end()){
      throw std::runtime_error("ReconstructionHypothesis::discriminator: discriminator with label '" + l + "' not set");
    }
    return it->second;
  }
  
  /// test if a discriminator value with a certian label has already been added
  bool has_discriminator(const std::string & label) const {
    return m_discriminators.find(label) != m_discriminators.end();
  }
  
  void set_blep_v4(LorentzVector v4){m_blep_v4=v4;}
  void set_toplep_v4(LorentzVector v4){m_toplep_v4=v4;}
  void set_tophad_v4(LorentzVector v4){m_tophad_v4=v4;} 
  void set_HZW_v4(LorentzVector v4){m_HZW_v4=v4;} 
  void set_neutrino_v4(LorentzVector v4){m_neutrino_v4=v4;}
  void add_toplep_jet(const Jet& j){m_toplep_jets.push_back(j);}
   void add_tophad_jet(const Jet& j){m_tophad_jets.push_back(j);}
  //void add_HZW_jet(const Jet& j){m_HZW_jets.push_back(j);}
  // void add_TOPtoplep_jet(const TopJet& j){m_TOPtoplep_jets.push_back(j);}
  // void add_tophad_jet(const TopJet& j){m_tophad_jets.push_back(j);}
  //void add_HZW_jet(const TopJet& j){m_HZW_jets.push_back(j);}
  void clear_toplep_jet(){m_toplep_jets.clear();}
  void clear_tophad_jet(){m_tophad_jets.clear();}
  void clear_subjets(){m_subjets.clear();}
  void set_tophad_topjet_ptr(const TopJet* const tjp){m_tophad_topjet_ptr = tjp;}
  void set_lepton(const Particle & l){m_lepton = l;}
  //  void set_HZW(const TopJet & l){m_HZW = l;}
  void set_subjets(const Jet& subjet){m_subjets.push_back(subjet);}
  void set_discriminator(const std::string & label, float discr){
    m_discriminators[label] = discr;
  }
  
 private:
  LorentzVector m_blep_v4;
  LorentzVector m_toplep_v4;
  LorentzVector m_tophad_v4;
  LorentzVector m_HZW_v4;
  LorentzVector m_neutrino_v4;
  std::vector<Jet>m_subjets;

  std::vector<Jet> m_toplep_jets;
  std::vector<Jet> m_tophad_jets;
  //  std::vector<Jet> m_HZW_jets;
  //  std::vector<TopJet> m_TOPtoplep_jets;
  //std::vector<TopJet> m_tophad_jets;
  std::vector<TopJet> m_HZW_jets;

  const TopJet* m_tophad_topjet_ptr;
  Particle m_lepton;
  // TopJet m_HZW;

  std::map<std::string, float> m_discriminators;
};
