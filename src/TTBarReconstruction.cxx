#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"
#include <cassert>
#include "UHH2/core/include/Event.h"
#include <vector>
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"

#include "UHH2/ZPrimeTotTPrime/include/TTBarReconstruction.h"

using namespace uhh2;
using namespace std;

LepTopReconstruction::LepTopReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label): m_neutrinofunction(neutrinofunction) {
    h_recohyps = ctx.declare_event_output<vector<ZPrimeTotTPrimeReconstructionHypothesis>>(label);
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
    h_btag = ctx.get_handle<std::vector<Jet> >("BTag_loose");
    berror =false;
}

LepTopReconstruction::~LepTopReconstruction() {}



bool LepTopReconstruction::process(uhh2::Event & event) {
  if(berror)  std::cout<<"TTBarReconstruction::LepTop Beginn"<<std::endl;
 
  assert(event.jets);
  assert(event.met);

 std::vector<TopJet>* topjets =event.topjets;
 if(berror) std::cout<<"Size of topjets am Anfang "<< topjets->size()<<std::endl;
  const Particle & lepton = event.get(h_primlep);
  std::vector<LorentzVector> neutrinos = m_neutrinofunction( lepton.v4(), event.met->v4());
  ZPrimeTotTPrimeReconstructionHypothesis hyp;

  if(!event.is_valid(h_btag))return false;
  auto & btag = event.get(h_btag);

  unsigned int n_jets = btag.size();
  if(berror)std::cout<<"Size der gebtaggedten Jets " << n_jets << std::endl;

  //LEP. TOP RECONSTRUCTION
  if(n_jets>8) n_jets=8;
    
  //nearest btagged jet to muon is for leptonic top 
  double dRmin = deltaR(lepton.v4(),btag.at(0));

  if(berror) std::cout<<"Kleinster Abstand erster btag zu muon "<<dRmin<<std::endl;

  auto blep =  btag.at(0);
  std::unique_ptr< std::vector<Jet> > other_btags(new std::vector<Jet> (*event.jets));
  other_btags->clear();
  other_btags->reserve(btag.size());
  std::unique_ptr< std::vector<Jet> > had_btags(new std::vector<Jet> (*event.jets));
  had_btags->clear();
  had_btags->reserve(btag.size());
  std::unique_ptr< std::vector<Jet> > bhad(new std::vector<Jet> (*event.jets));
  bhad->clear();
  bhad->reserve(btag.size());
  double n_blep=0; 

  for(unsigned int k=1; k<n_jets; k++){

    double delta_lep = deltaR(lepton.v4(),btag.at(k));
   if(berror) std::cout<<"DeltaR um naechsten btagged jet zu finden " << delta_lep << std::endl;
   if(delta_lep < dRmin){other_btags -> push_back(blep);blep =  btag.at(k);dRmin = delta_lep;n_blep++;}
    else other_btags -> push_back(btag.at(k));
  }//nearest btagged jet to muon is for leptonic top 

  if(berror) std::cout<<"Kleinster Abstand vom blep zu muon "<<dRmin<<std::endl;
  
  //if(!n_blep)return false;

  if(berror)std::cout <<"Size of other bttaged jets (not blep) "<<other_btags->size()<<std::endl;

  //loop over all neutrino solutions
  for(const auto & neutrino_p4 : neutrinos){
    
    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
    LorentzVector toplep_v4 = wlep_v4;
    int lepjets=0;
          
    hyp.set_lepton(lepton);
    hyp.set_neutrino_v4(neutrino_p4);
    hyp.clear_toplep_jet();
    
    toplep_v4 = toplep_v4 +blep.v4();
    hyp.add_toplep_jet(blep);
    lepjets++;
    //leptop reconstructed

    //reconstruct had top with condition: W has to be AK8, b jet AK4, W should not overlap with both b jets (lep. top and had. top), W & b has to be back to back to lep. top

    double n_bhad =0;
    //b jet is back to back
    for(unsigned int k=0; k<other_btags->size(); k++){
      double delta_had = deltaPhi(toplep_v4,other_btags->at(k));
      if(delta_had>(M_Pi/2)){ bhad ->push_back( other_btags->at(k));n_bhad++;}
   if(berror)   std::cout<<"deltaR(bhad,blep) "<<delta_had << std::endl;
    }

    if(!n_bhad)return false;
    if(berror)std::cout<<"Size of had b cand. "<< bhad -> size()<< "   Abstand zum mu "<<  deltaPhi(toplep_v4,bhad->at(0)) <<std::endl;

    // over all had b cand
    for(unsigned int l=0; l<bhad->size(); l++){
      //W cand should have no overlap with both b jets and has to be back to back to lep top
      std::unique_ptr< std::vector<TopJet> > w_cand(new std::vector<TopJet> (*event.topjets));
      w_cand->clear();
      w_cand->reserve(topjets->size());

      // if(berror)std::cout<<"Size W_cand nach reserve "<< w_cand->at(0).v4() <<std::endl;
      for(unsigned int a = 0; a < topjets->size(); a++){
	double deltaR_WB = deltaR(bhad->at(l),topjets->at(a));
	if(berror)	std::cout<<"deltaRWB ueber alle topjets "<< deltaR_WB << std::endl;
	double deltaR_WBB = deltaR(blep,topjets->at(a));
	double deltaPhi_WB = deltaPhi(toplep_v4,topjets->at(a));
	if(deltaR_WB > 0.8 && deltaR_WBB > 0.8 && deltaPhi_WB >(M_Pi/2)){ w_cand -> push_back(topjets->at(a));}
      }

   
      if(!w_cand->size())return false;
      if(berror)std::cout<<"Size W_cand "<< w_cand->size() <<std::endl;
      

      for(unsigned int s = 0; s < w_cand->size();s++){
	if(berror)	std::cout<<"Abstand W had zu blep "<< deltaR(blep,w_cand->at(s).v4())<<"    Abstand Whad zu bhad"<<deltaR(bhad->at(l),w_cand->at(s).v4())<<std::endl;

	hyp.set_W(w_cand->at(s));
	hyp.set_W_v4(w_cand->at(s).v4());
	
	//EXPERIMENT
	hyp.clear_Wsubjets();
	for(const Jet m:w_cand->at(s).subjets()) hyp.set_Wsubjets(m);
	//ENDE

	hyp.set_HZW_v4(w_cand->at(s).v4());
	hyp.add_tophad_jet(bhad->at(l));
	hyp.set_bhad_v4(bhad->at(l).v4());
	hyp.set_tophad_v4(w_cand->at(s).v4() + bhad->at(l).v4() );
	//search jet with highest pt assigned to leptonic top
	int blep_idx(-1);
	float maxpt(-1.);
	for(unsigned int i=0; i<hyp.toplep_jets().size(); ++i){
	  if(maxpt< hyp.toplep_jets().at(i).pt()){
	    maxpt = hyp.toplep_jets().at(i).pt();
	    blep_idx = i;
	  }
	}
	if(blep_idx != -1) hyp.set_blep_v4(hyp.toplep_jets().at(blep_idx).v4());
	if(berror) std::cout<<"for dem setzen der Hyp "<< lepjets <<std::endl;
	if( lepjets>0 ) {
	  hyp.set_toplep_v4(toplep_v4);
	  recoHyps.emplace_back(move(hyp));
	  if(berror)	  std::cout<<"in if"<<recoHyps.size() <<std::endl;
	}
      }//w_cand
    }// bhad 
  } //neutrinos
  

  event.set(h_recohyps, move(recoHyps));
 
  return true;
}
