#include "UHH2/ZPrimeTotTPrime/include/BackgroundGenSelection.h"
#include "UHH2/core/include/LorentzVector.h"
#include <stdexcept>
#include <cmath>



using namespace uhh2;
using namespace std;

BackgroundGen::BackgroundGen(const vector<GenParticle> & genparticles){

  for(unsigned int i=0; i<genparticles.size(); i++){
    const GenParticle & genp = genparticles[i];
    if(abs(genp.pdgId()) == 24){

      if(abs(genp.daughter(&genparticles,1)->pdgId()) <=6 ){
	Quark1 = *genp.daughter(&genparticles,1);
	Quark2 = *genp.daughter(&genparticles,2);
	}
    }else if(abs(genp.pdgId()) == 23){

      if(abs(genp.daughter(&genparticles,1)->pdgId()) <=6 ){
	Quark1 = *genp.daughter(&genparticles,1);
	Quark2 = *genp.daughter(&genparticles,2);
	}   
    }
}
}


///////////////////////////////////////////   Variables //////////////////////////////////////////////////////


LorentzVector BackgroundGen::Quark1_v4()const{return Quark1.v4();};
LorentzVector BackgroundGen::Quark2_v4()const{return Quark2.v4();};




///////////////////////////////////////////////////////////////////////////////////////////////////////////////
BackgroundGenProducer::BackgroundGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
   h_background = ctx.get_handle<BackgroundGen>(name);
}


bool BackgroundGenProducer::process(Event & event){
    event.set(h_background, BackgroundGen(*event.genparticles));
    return true;
}
