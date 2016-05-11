#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"
#include <cassert>
#include "UHH2/core/include/Event.h"
#include <vector>
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"


#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeSidebandReconstruction2.h"

using namespace uhh2;
using namespace std;

ZPrimeTotTPrimeSidebandReconstruction2::ZPrimeTotTPrimeSidebandReconstruction2(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label): m_neutrinofunction(neutrinofunction) {
    h_recohyps = ctx.declare_event_output<vector<ZPrimeTotTPrimeReconstructionHypothesis>>(label);
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
    h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
    //  h_toptag = ctx.get_handle<std::vector<TopJet> >("TopTag");
    //   h_higgstag = ctx.get_handle<std::vector<TopJet> >("HiggsTag");
    // h_ZWtag =ctx.get_handle<std::vector<TopJet> >("ZWTag");
    berror=false;

}

ZPrimeTotTPrimeSidebandReconstruction2::~ZPrimeTotTPrimeSidebandReconstruction2() {}



bool ZPrimeTotTPrimeSidebandReconstruction2::process(uhh2::Event & event) {

  if(berror)  std::cout<<"############################################### SidebandReconstruction Beginn ###########################################################"<<std::endl;
 
  assert(event.jets);
  assert(event.met);

  const Particle & lepton = event.get(h_primlep);
  std::vector<LorentzVector> neutrinos = m_neutrinofunction( lepton.v4(), event.met->v4());

  ZPrimeTotTPrimeReconstructionHypothesis hyp;

  //TTbar Reconstruction with ONLY AK4 jets
  if(!LepHadTop(neutrinos,lepton,hyp, event))return false;
       
  event.set(h_recohyps, move(recoHyps));
  return true;
}




// std::vector<TopJet>* ZPrimeTotTPrimeSidebandReconstruction2::AK8cleaning(std::vector<Jet>* usedAK4,uhh2::Event & event ){
//   if(berror)   std::cout << "SidebandReconstruction: L 53 Begin AK8cleaning "<< std::endl;


//   std::vector<TopJet>* AK8Jets(new std::vector<TopJet> (*event.topjets));
//   AK8Jets->clear();
//   AK8Jets->reserve(event.topjets->size());
//   if(berror) std::cout<< "SidebandReconstruction L 59 Size aller AK8 jets vor cleaning "<< event.topjets->size()<<std::endl;
//   for(const TopJet ak8:*event.topjets){
// if(berror)   std::cout << "SidebandReconstruction: L 61 all AK8 Jets "<< std::endl;
//   bool bdeltaR=true;
//     for(const Jet ak4:*usedAK4){
//     double deltar = deltaR(ak4,ak8);
//     if(berror)   std::cout << "SidebandReconstruction: L 64 all AK4 Jets (deltaR(ak4,ak8)<1.2)"<<deltar<< std::endl;
//     if(deltar < 1.2)bdeltaR=false;
//     }
//     if(bdeltaR)AK8Jets ->push_back(ak8);
//   }
//   sort_by_pt<TopJet>(*AK8Jets);
//   if(berror)   std::cout << "SidebandReconstruction: L 71 cleaned ak8 jets (size)"<<AK8Jets->size()<< std::endl;
//   return AK8Jets;
// }



bool ZPrimeTotTPrimeSidebandReconstruction2::LepHadTop(std::vector<LorentzVector> neutrinos, const Particle & lepton , ZPrimeTotTPrimeReconstructionHypothesis hyp,uhh2::Event & event ){
  if(berror)   std::cout << "SidebandReconstruction: L 74 Begin LepHadTop"<< std::endl;

  std::vector<Jet>* jets =event.jets;
  unsigned int n_jets = jets->size();


  //LEP. TOP RECONSTRUCTION
  if(n_jets>8) n_jets=8;
  const unsigned int max_j = pow(3, n_jets);// (0)leptonisches Top, (1)had Top, (2) nichts
     
  //bool bhyp=false;   
  for(const auto & neutrino_p4 : neutrinos){
    if(berror)   std::cout << "SidebandReconstruction: L 86 in NeutrinoSolutions"<< std::endl;

    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
   
    for (unsigned int j=0; j < max_j; j++) {
      if(berror)   std::cout << "SidebandReconstruction: L 91 in all combinations"<< std::endl;
      LorentzVector toplep_v4 = wlep_v4;
      LorentzVector tophad_v4;
      int lepjets=0;
      int hadjets=0;
      int num = j;
  
      hyp.set_lepton(lepton);
      hyp.set_neutrino_v4(neutrino_p4);
      hyp.clear_toplep_jet();
      hyp.clear_tophad_jet();

      if(berror)   std::cout << "SidebandReconstruction: L 105 in asign to lep or had top (njets)"<<n_jets<< std::endl;
      for (unsigned int k=0; k<n_jets; k++) {
	if(num%4==0) {
	  toplep_v4 = toplep_v4 +jets->at(k).v4();
	  hyp.add_toplep_jet(jets->at(k));
	  lepjets++;
	}
	if(num%4==1) {
	  tophad_v4 = tophad_v4 +jets->at(k).v4();
	  hyp.add_tophad_jet(jets->at(k));
	  hadjets++;
	}
	num /= 3;
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
      if(berror)   std::cout << "SidebandReconstruction: L 131 if(lepjets>0 && hadjets>0)"<<lepjets << "had"<<hadjets<< std::endl;
      if( lepjets>0 && hadjets>0) {
       	hyp.set_toplep_v4(toplep_v4);
	hyp.set_tophad_v4(tophad_v4);

	// //cleaning all AK8 overlapping with the used AK4
	// std::vector<Jet>* usedak4(new std::vector<Jet> (*event.jets));
	// usedak4->clear();
	// usedak4->reserve(event.jets->size());
	// for(const Jet jet: hyp.toplep_jets()) usedak4->push_back(jet);
	// for(const Jet jet: hyp.tophad_jets()) usedak4->push_back(jet);

	// if(berror)   std::cout << "SidebandReconstruction: L 142 defined usedak4 (size should be greater 0)"<<usedak4->size()<< std::endl;

	// std::vector<TopJet>* topjets;

	// topjets = AK8cleaning(usedak4, event);
	 if(berror) std::cout << "SidebandREconstruction> L 147  topjets size "<< event.topjets->size()<<std::endl;
	// //	if(!(topjets->size()))break;
	//	if(event.topjets->size()){
	  hyp.set_HZW_v4(event.topjets->at(0).v4()); 
	  recoHyps.emplace_back(move(hyp));
	  //   bhyp=true;
	  //  	} 

      }
    } //combinations
  }//neutrino

   //if(bhyp) return true;
  //else return false;
   return true;
}

