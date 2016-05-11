#pragma once

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h>

const ZPrimeTotTPrimeReconstructionHypothesis * get_best_hypothesis(const std::vector<ZPrimeTotTPrimeReconstructionHypothesis> & hyps, const std::string & label);


class ZPrimeTotTPrimeChi2Discriminator: public uhh2::AnalysisModule {
public:
    struct cfg {
        std::string discriminator_label;
        cfg(): discriminator_label("Chi2"){}
    };
    
    ZPrimeTotTPrimeChi2Discriminator(uhh2::Context & ctx, const std::string & rechyps_name, const cfg & config = cfg());
    virtual bool process(uhh2::Event & event) override;
    
private:
    uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>> h_hyps;
    cfg config;
};
