#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TTbarGen.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"


  class EffHists: public uhh2::Hists {
  public:
   
    EffHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name);

    virtual void fill(const uhh2::Event & ev) override;

  private:
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    uhh2::Event::Handle<std::vector<Jet>> h_btag_medium; 
    uhh2::Event::Handle<std::vector<Jet>> h_btag_loose; 
    uhh2::Event::Handle<std::vector<Jet>> h_btag_tight;
    uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>> h_hyps;

    std::string m_discriminator_name;
  };


