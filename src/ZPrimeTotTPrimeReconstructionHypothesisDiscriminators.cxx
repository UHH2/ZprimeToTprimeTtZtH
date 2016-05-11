#include "UHH2/core/include/Utils.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"

#include <set>

using namespace uhh2;
using namespace std;

namespace {
    
// invariant mass of a lorentzVector, but save for timelike / spacelike vectors
float inv_mass(const LorentzVector & p4){
    if(p4.isTimelike()){
            return p4.mass();
    }
    else{
        return -sqrt(-p4.mass2());
    }
}

}

const ZPrimeTotTPrimeReconstructionHypothesis * get_best_hypothesis(const std::vector<ZPrimeTotTPrimeReconstructionHypothesis> & hyps, const std::string & label){
    const ZPrimeTotTPrimeReconstructionHypothesis * best = nullptr;
    float current_best_disc = numeric_limits<float>::infinity();
    for(const auto & hyp : hyps){
      //   std::cout <<"in for " << std::endl;
        if(!hyp.has_discriminator(label)) continue;
	//		std::cout << "hypothesis if" << std::endl;
        auto disc = hyp.discriminator(label);
        if(disc < current_best_disc){
            best = &hyp;
            current_best_disc = disc;
        }
    }
    //       std::cout << "current best" << current_best_disc << std::endl;
    if(std::isfinite(current_best_disc)){
      //          std::cout << "is finite"<< std::endl;
        return best;
    }
    else{
        return nullptr;
    }
}

////


ZPrimeTotTPrimeChi2Discriminator::ZPrimeTotTPrimeChi2Discriminator(Context & ctx, const std::string & rechyps_name, const cfg & config_): config(config_){
    h_hyps = ctx.get_handle<vector<ZPrimeTotTPrimeReconstructionHypothesis>>(rechyps_name);
}


bool ZPrimeTotTPrimeChi2Discriminator::process(uhh2::Event & event){
    auto & hyps = event.get(h_hyps);
    const double mass_thad = 181;
    const double mass_thad_sigma = 15;
    const double mass_tlep = 174;
    const double mass_tlep_sigma = 18;
    const double mass_Higgs=125;
    const double mass_Higgs_sigma = 18;
    const double mass_Z=91;
    const double mass_Z_sigma = 21;

    for(auto & hyp: hyps){
        double mass_thad_rec = inv_mass(hyp.tophad_v4());
        double mass_tlep_rec = inv_mass(hyp.toplep_v4());
	double mass_Higgs_rec = inv_mass(hyp.HZW_v4());
	double chi2_thad = pow((mass_thad_rec - mass_thad) / mass_thad_sigma, 2);
        double chi2_tlep = pow((mass_tlep_rec - mass_tlep) / mass_tlep_sigma, 2);
	double chi2_H = pow((mass_Higgs_rec - mass_Higgs)/mass_Higgs_sigma,2);
	double chi2_Z = pow((mass_Higgs_rec - mass_Z)/mass_Z_sigma,2);
	double bestchi2=0;

	// if((chi2_thad + chi2_tlep + chi2_H) > (chi2_thad + chi2_tlep + chi2_Z)){
	//   std::cout << "Falsch" <<endl;
	// }
	// if((chi2_thad + chi2_tlep + chi2_H) > (chi2_thad + chi2_tlep + chi2_Z) ) {bestchi2 = chi2_thad + chi2_tlep + chi2_H ;}else{bestchi2 = chi2_thad + chi2_tlep + chi2_Z;}
	bestchi2 = chi2_thad + chi2_tlep;
	hyp.set_discriminator(config.discriminator_label, bestchi2 );
        hyp.set_discriminator(config.discriminator_label + "H", chi2_thad + chi2_tlep + chi2_H );
	hyp.set_discriminator(config.discriminator_label + "Z", chi2_thad + chi2_tlep + chi2_Z );
        hyp.set_discriminator(config.discriminator_label + "_tlep", chi2_tlep);
        hyp.set_discriminator(config.discriminator_label + "_thad", chi2_thad);
	hyp.set_discriminator(config.discriminator_label + "_H", chi2_H);
	hyp.set_discriminator(config.discriminator_label + "_Z", chi2_Z);
  }
  return true;
}
////
