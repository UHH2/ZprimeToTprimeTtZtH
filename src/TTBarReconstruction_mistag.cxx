#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"
#include <cassert>
#include "UHH2/core/include/Event.h"
#include <vector>
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"

#include "UHH2/ZprimeToTprimeTtZtH/include/TTBarReconstruction_mistag.h"

using namespace uhh2;
using namespace std;

LepTopReconstruction_mistag::LepTopReconstruction_mistag(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label): m_neutrinofunction(neutrinofunction) {
  h_recohyps = ctx.declare_event_output<vector<ZPrimeTotTPrimeReconstructionHypothesis>>(label);
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
  h_btag = ctx.get_handle<std::vector<Jet> >("BTag_medium");
  berror =false;
  distance = string2double( ctx.get("distance"));
}

LepTopReconstruction_mistag::~LepTopReconstruction_mistag() {}



bool LepTopReconstruction_mistag::process(uhh2::Event & event) {
  if(berror)  std::cout<<"///// TTBarReconstruction::LepTop Beginn //////"<<std::endl;
 
  assert(event.jets);
  assert(event.met);

  std::vector<TopJet>* topjets =event.topjets;
  if(berror) std::cout<<"TTBarReconstruction:: Size of topjets at the beginning "<< topjets->size()<<std::endl;
  const Particle & lepton = event.get(h_primlep);
  std::vector<LorentzVector> neutrinos = m_neutrinofunction( lepton.v4(), event.met->v4());
  ZPrimeTotTPrimeReconstructionHypothesis hyp;

  if(!event.is_valid(h_btag))return false;
  auto & btag = event.get(h_btag);

  unsigned int n_jets = btag.size();
  if(berror)std::cout<<"TTBarReconstruction:: Size of btagged jets " << n_jets << std::endl;

  //LEP. TOP RECONSTRUCTION
  if(n_jets>8) n_jets=8;
    
  //nearest btagged jet to muon is for leptonic top 
  double dRmin = deltaR(lepton.v4(),btag.at(0));

  if(berror) std::cout<<"TTBarReconstruction:: distance between first btag and muon "<<dRmin<<std::endl;

  auto blep =  btag.at(0);
  std::unique_ptr< std::vector<Jet> > other_btags(new std::vector<Jet> (*event.jets));
  other_btags->clear();
  other_btags->reserve(btag.size());
  std::unique_ptr< std::vector<Jet> > had_btags(new std::vector<Jet> (*event.jets));
  had_btags->clear();
  had_btags->reserve(btag.size());
  double n_blep=0; 
  for(unsigned int k=1; k<n_jets; k++){
    double delta_lep = deltaR(lepton.v4(),btag.at(k));
    if(delta_lep < dRmin){other_btags -> push_back(blep);blep =  btag.at(k);dRmin = delta_lep;n_blep++;}
    else other_btags -> push_back(btag.at(k));
  }//nearest btagged jet to muon is for leptonic top 

  if(berror) std::cout<<"TTBarReconstruction:: smallest distance of btag to muon "<<dRmin<<std::endl;
  if(berror)std::cout <<"Size of other btaged jets (not blep) "<<other_btags->size()<<std::endl;

  //loop over all neutrino solutions
  for(const auto & neutrino_p4 : neutrinos){
    if(berror)std::cout<<"///////////////////  In neutrion solution"<<std::endl;
    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
    LorentzVector toplep_v4 = wlep_v4;
    int lepjets=0;
    
    toplep_v4 = toplep_v4 +blep.v4();
    lepjets++;

    //reconstruct had top with condition: W has to be AK8, b jet AK4, W should not overlap with both b jets (lep. top and had. top), W & b has to be back to back to lep. top
    //find all AK4 jets that are not overlapping with blep
    std::unique_ptr< std::vector<Jet> > other_ak4jets(new std::vector<Jet> (*event.jets));
    other_ak4jets->clear();
    other_ak4jets->reserve(event.jets->size());
    for(const auto &ak4jet:*event.jets){
      double deltar = deltaR(ak4jet,blep.v4());
      if(deltar>0.8) other_ak4jets->push_back(ak4jet);
    }

    
   
   
    //all other AK8 jets on the other hemisphere
    double n_ak8had =0;
    std::unique_ptr< std::vector<TopJet> > ak8had(new std::vector<TopJet> (*event.topjets));
    ak8had->clear();
    ak8had->reserve(topjets->size());
   
    //ak8 jet is back to back
  
    for(unsigned int k=0; k<topjets->size(); k++){
      double delta_blep = deltaR(blep,topjets->at(k));
      double delta_had = deltaPhi(toplep_v4,topjets->at(k));
    //Hier
      //    if(delta_had>(M_Pi/2)&&delta_blep>1.2){ak8had->push_back( topjets->at(k));n_ak8had++;}
      if(delta_blep>1.2){ak8had->push_back( topjets->at(k));n_ak8had++;}
   
      if(berror)   std::cout<<"deltaR(bhad,blep) "<<delta_had << std::endl;
    }

    if(!n_ak8had)return false;
    if(berror)std::cout<<"Size of had b cand. "<< ak8had -> size()<< " deltaR to mu "<<  deltaPhi(toplep_v4,ak8had->at(0)) <<std::endl;
   
    //for all AK8 jets find not overlapping ak4 jets
    int help_i=0;
    for(const auto &ak8jet:*ak8had){
      help_i++;
      if(berror) std::cout<<"over all ak8 jets  "<< help_i <<" / "<< ak8had ->size()<<std::endl;
    
      double n_ak4had =0;
      std::vector<Jet>* ak4had(new std::vector<Jet> (*event.jets));
      ak4had->clear();
      ak4had->reserve(other_ak4jets->size());
      
      //ak4 jet is back to back and does not overlap with current AK8 jet
      for(unsigned int k=0; k<other_ak4jets->size(); k++){
	double delta_had = deltaPhi(toplep_v4,other_ak4jets->at(k));
	double delta_wcand = deltaR(other_ak4jets->at(k),ak8jet);
	if(delta_had>(M_Pi/2)&&delta_wcand>1.2){ak4had  ->push_back( other_ak4jets->at(k));n_ak4had++;}
	if(berror)   std::cout<<"deltaR(bhad,blep) "<<delta_had << std::endl;
      }
    
      if(!n_ak4had)return false;
      if(berror)std::cout<<"Size of had b cand. "<< ak4had -> size()<< " deltaR to mu "<<  deltaPhi(toplep_v4,ak4had->at(0)) <<std::endl;


      // // over all had b cand
      // for(unsigned int l=0; l<bhad->size(); l++){
      // if(berror)  std::cout<<"///////////////////////  In bhad cand  "<<l<<std::endl;
 
      // //W cand should have no overlap with both b jets and has to be back to back to lep top
      // std::unique_ptr< std::vector<TopJet> > w_cand(new std::vector<TopJet> (*event.topjets));
      // w_cand->clear();
      // w_cand->reserve(topjets->size());
     
      // for(unsigned int a = 0; a < topjets->size(); a++){
      // 	if(berror)std::cout<<"//////////////////////// In w cand"<<std::endl;
      // 	double deltaR_WB = deltaR(bhad->at(l),topjets->at(a));
      // 	double deltaR_WBB = deltaR(blep,topjets->at(a));
      // 	double deltaPhi_WB = deltaPhi(toplep_v4,topjets->at(a));
      // 	if(deltaR_WB > distance && deltaR_WBB > distance && deltaPhi_WB >(M_Pi/2)){ w_cand -> push_back(topjets->at(a));}
      // }
    
     
      // if(!w_cand->size())return false;
      // if(berror)std::cout<<"Size W_cand "<< w_cand->size() <<std::endl;
      

      for(unsigned int s = 0; s <ak4had ->size();s++){
	if(berror) std::cout<<"over all ak4 jets  "<< s+1 <<" / "<< ak4had ->size()<<std::endl;
	if(berror)	std::cout<<"deltaR W had zu blep "<< deltaR(blep,ak8jet.v4())<<"deltaR Whad zu bhad"<<deltaR(ak4had->at(s),ak8jet.v4())<<std::endl;

	if(berror) std::cout << "Mass of W cand .M() "<<ak8jet.v4().M()<<std::endl;

	//EXPERIMENT
	hyp.clear_Wsubjets();
	if(berror)	std::cout <<"number of subjets "<<hyp.W_subjets().size() <<"und viererimpuls "<<std::endl;
	for(const Jet m:ak8jet.subjets()) hyp.set_Wsubjets(m);
	//ENDE

	hyp.set_W(ak8jet);
	hyp.set_W_v4(ak8jet.v4());
	hyp.set_lepton(lepton);
	hyp.set_neutrino_v4(neutrino_p4);
	hyp.clear_toplep_jet();
	hyp.add_toplep_jet(blep);
       	hyp.clear_tophad_jet();
	hyp.add_tophad_jet(ak4had->at(s));
	 hyp.set_bhad_v4(ak4had->at(0).v4());
	hyp.set_tophad_v4(ak8jet.v4() + ak4had->at(s).v4() );

	if(berror)	std::cout<<"with _v4() v4 blep " <<hyp.blep_v4() << " v4 bhad "<<hyp.bhad_v4() << "v4 W "<< hyp.W_v4()<<std::endl;


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
	if(berror) std::cout<<"before set Hyp number of lepjets >0  "<< lepjets <<std::endl;
	if( lepjets>0   ) {
	  hyp.set_toplep_v4(toplep_v4);
	  recoHyps.emplace_back(move(hyp));
	  if(berror)   std::cout<<"in if  "<< recoHyps.size() <<std::endl;
	}
	}//AK4 jets
	//  if(berror) std::cout<<"Before LepHadTop_mistag"<<std::endl;
      // LepHadTop_mistag(ak4had,neutrinos,lepton,hyp, ak8jet);

    }// AK8 jets 
  } //neutrinos
  
  event.set(h_recohyps, move(recoHyps));
  return true;
}




void LepTopReconstruction_mistag::LepHadTop_mistag(std::vector<Jet>* ak4,  std::vector<LorentzVector> neutrinos, const Particle & lepton , ZPrimeTotTPrimeReconstructionHypothesis hyp, TopJet Tag){

  unsigned int n_jets = ak4->size();
  //LEP. TOP RECONSTRUCTION
  if(n_jets>6) n_jets=6;
  const unsigned int max_j = pow(3, n_jets);// (0)leptonisches Top, (1)had Top, (2) nichts
  
  for(const auto & neutrino_p4 : neutrinos){
     if(berror) std::cout<<"neutrino solution  "<<std::endl;
    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
    for (unsigned int j=0; j < max_j; j++) {
      if(berror)std::cout<<"different hyps "<<std::endl;
      //EXPERIMENT
      hyp.clear_subjets(); 
      for(const Jet s:Tag.subjets()) hyp.set_subjets(s);
      //ENDE

      LorentzVector toplep_v4 = wlep_v4;
      LorentzVector tophad_v4;
      int lepjets=0;
      int hadjets=0;
      int num = j;
      
      hyp.set_lepton(lepton);
      hyp.set_neutrino_v4(neutrino_p4);
      // clear hyp.add
      // hyp.clear_toplep_jet();
      hyp.clear_tophad_jet();
      for (unsigned int k=0; k<n_jets; k++) {


	if(num%4==0) {
	  tophad_v4 = tophad_v4 +ak4->at(k).v4();
	  hyp.add_tophad_jet(ak4->at(k));
	  hadjets++;
	}
	num /= 2;
      }//assign lep 

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
      
      if(  hadjets>0) {
       	hyp.set_toplep_v4(toplep_v4);
	hyp.set_tophad_v4(tophad_v4);
	recoHyps.emplace_back(move(hyp));
	if(berror) std::cout<<"number of hyps  "<< recoHyps.size() <<std::endl;
      }
    } //combinations
  }
}
