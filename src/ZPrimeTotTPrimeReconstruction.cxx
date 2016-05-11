#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/core/include/Utils.h"
#include <cassert>
#include "UHH2/core/include/Event.h"
#include <vector>
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"


#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstruction.h"

using namespace uhh2;
using namespace std;

ZPrimeTotTPrimeReconstruction::ZPrimeTotTPrimeReconstruction(Context & ctx, const NeutrinoReconstructionMethod & neutrinofunction, const string & label): m_neutrinofunction(neutrinofunction) {
    h_recohyps = ctx.declare_event_output<vector<ZPrimeTotTPrimeReconstructionHypothesis>>(label);
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
    h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
    h_toptag = ctx.get_handle<std::vector<TopJet> >("TopTag");
    h_higgstag = ctx.get_handle<std::vector<TopJet> >("HiggsTag");
    h_ZWtag =ctx.get_handle<std::vector<TopJet> >("ZWTag");
    berror =false;
}

ZPrimeTotTPrimeReconstruction::~ZPrimeTotTPrimeReconstruction() {}



bool ZPrimeTotTPrimeReconstruction::process(uhh2::Event & event) {
  if(berror)  std::cout<<"Reconstruction Beginn"<<std::endl;
 
  assert(event.jets);
  assert(event.met);

  const Particle & lepton = event.get(h_primlep);
 
 
  std::vector<LorentzVector> neutrinos = m_neutrinofunction( lepton.v4(), event.met->v4());


  std::vector<Jet>* jets =event.jets;
  ZPrimeTotTPrimeReconstructionHypothesis hyp;

 

  if((event.is_valid(h_higgstag) && event.is_valid(h_toptag)) ||(event.is_valid(h_ZWtag) && event.is_valid(h_toptag)) ){
    const auto & toptag = event.get(h_toptag);
    if(toptag.size()>0){
 

      if(event.is_valid(h_higgstag)){
	std::vector<TopJet> higgstag = event.get(h_higgstag);
 if(berror) 	std::cout<<"in higgstag mit toptag"<<higgstag.size()<<std::endl;

	if(higgstag.size()>0){
	  hyp.set_HZW_v4(higgstag.at(0).v4());
	  //  hyp.set_HZW(higgstag.at(0));
	  // //	  hyp.clear_subjets();
	  // for(const Jet s:higgstag.at(0).subjets()) hyp.set_subjets(s);
	  // std::cout<<"higgscleaning"<<std::endl;
	  //cleaning AK8 higgstag from AK4
	  std::vector<Jet>* higgsak4jets;
	  higgsak4jets = AK4cleaning(jets,higgstag.at(0) );

	 
	  hyp.set_tophad_v4(toptag.at(0).v4());
	  std::cout<<"topcleaning"<<std::endl;
	  //cleaning AK8 Toptag from AK4
	  std::vector<Jet>* topak4jets;
	  topak4jets = AK4cleaning(higgsak4jets,toptag.at(0) );

	  if(topak4jets->size()== 0) return false;
	  //all neutrino solutions and top lep
	  LepTop(topak4jets,neutrinos,lepton,hyp, higgstag);
	}else if(event.is_valid(h_ZWtag)){
	std::vector<TopJet> ZWtag = event.get(h_ZWtag);
 if(berror) 	std::cout<<"in ZWtag mit toptag  "<< ZWtag.size()<<std::endl;
	if(ZWtag.size()>0){
	  hyp.set_HZW_v4(ZWtag.at(0).v4());
	  //  hyp.set_HZW(ZWtag.at(0));
	  //	  hyp.clear_subjets();
 // std::cout <<"in REconstrucktion Size der subjets zwtag "<< ZWtag.at(0).subjets().size()<<" v4 des ersten " <<ZWtag.at(0).subjets().at(0).v4() << std::endl;
 // 	  for(const Jet s:ZWtag.at(0).subjets()) hyp.set_subjets(s);
 // 	  std::cout << "Test in REconstruction " << hyp.HZW_subjets().size()<< std::endl;
	  //cleaning AK8 ZWtag from AK4
	  std::vector<Jet>* ZWak4jets;
 if(berror) 	   std::cout<<"zwtagcleaning"<<std::endl;
	  ZWak4jets = AK4cleaning(jets,ZWtag.at(0) );
	
	 
	    hyp.set_tophad_v4(toptag.at(0).v4());
 if(berror) 	    std::cout<<"topcleaning"<<std::endl;
	    //cleaning AK8 Toptag from AK4
	    std::vector<Jet>* topak4jets;
	    topak4jets = AK4cleaning(ZWak4jets,toptag.at(0) );
	  
	    if(topak4jets->size()== 0) return false; 
	    //all neutrino solutions and top lep
	    LepTop(topak4jets,neutrinos,lepton,hyp, ZWtag);

	 
	}//ZWtag >0
      }//ZW valid
      }
      event.set(h_recohyps, move(recoHyps));
   
 
    }// at least one toptag
    else if(event.is_valid(h_higgstag) ||event.is_valid(h_ZWtag)){
      if(event.is_valid(h_higgstag)){
	std::vector<TopJet> higgstag = event.get(h_higgstag);
 if(berror) 	std::cout<<"Reconstruction L:110 in Higgstag und kein Toptag"<<higgstag.size()<<std::endl;
      
	if(higgstag.size()>0){
	  hyp.set_HZW_v4(higgstag.at(0).v4());
	  // hyp.set_HZW(higgstag.at(0));
	  //	  hyp.clear_subjets();
	  // for(const Jet s:higgstag.at(0).subjets()) hyp.set_subjets(s);
	  //cleaning AK8 higgstag from AK4
	  std::vector<Jet>*  higgsak4jets;
	  higgsak4jets = AK4cleaning(jets,higgstag.at(0) );
	  if(higgsak4jets->size() <=1) return false;
	//all neutrino solutions and top lep
	  LepHadTop(higgsak4jets,neutrinos,lepton,hyp, higgstag);

    }else if(event.is_valid(h_ZWtag)){
      std::vector<TopJet> ZWtag = event.get(h_ZWtag);
  if(berror)      std::cout<<"Reconstruction L:127 in zwstag und kein Toptag"<<ZWtag.size() << std::endl;
     
      if(ZWtag.size()>0){
	hyp.set_HZW_v4(ZWtag.at(0).v4());
	//	hyp.set_HZW(ZWtag.at(0));
 // 	//	hyp.clear_subjets();
 // std::cout <<"in REconstrucktion Size der subjets zwtag "<< ZWtag.at(0).subjets().size()<<" v4 des ersten " <<ZWtag.at(0).subjets().at(0).v4() << std::endl;
 // 	for(const Jet s:ZWtag.at(0).subjets()) hyp.set_subjets(s);
 //  std::cout << "Test in REconstruction " << hyp.HZW_subjets().size()<< std::endl;
	//cleaning AK8 ZWtag from AK4
	std::vector<Jet>* ZWak4jets;
	ZWak4jets = AK4cleaning(jets,ZWtag.at(0) );
	//all neutrino solutions and top lep
	if(ZWak4jets->size() <= 1)return false;
	LepHadTop(ZWak4jets,neutrinos,lepton,hyp, ZWtag);

    }//all ZWtag
  }//ZW valid

      }
  event.set(h_recohyps, move(recoHyps));
  }//higgs||ZW
 }// higss&&Top || ZW&&Top
 
  return true;
}


std::vector<LorentzVector> NeutrinoReconstruction(const LorentzVector & lepton, const LorentzVector & met) {
    TVector3 lepton_pT = toVector(lepton);
    lepton_pT.SetZ(0);
    TVector3 neutrino_pT = toVector(met);
    neutrino_pT.SetZ(0);
    constexpr float mass_w = 80.399f;
    float mu = mass_w * mass_w / 2 + lepton_pT * neutrino_pT;
    float A = - (lepton_pT * lepton_pT);
    float B = mu * lepton.pz();
    float C = mu * mu - lepton.e() * lepton.e() * (neutrino_pT * neutrino_pT);
    float discriminant = B * B - A * C;
    std::vector<LorentzVector> solutions;
    if (0 >= discriminant) {
        // Take only real part of the solution for pz:
        LorentzVectorXYZE solution (met.Px(),met.Py(),-B / A,0);
        solution.SetE(solution.P());
        solutions.emplace_back(toPtEtaPhi(solution));
    }
    else {
        discriminant = sqrt(discriminant);
        LorentzVectorXYZE solution (met.Px(),met.Py(),(-B - discriminant) / A,0);
        solution.SetE(solution.P());
        solutions.emplace_back(toPtEtaPhi(solution));

        LorentzVectorXYZE solution2 (met.Px(),met.Py(),(-B + discriminant) / A,0);
        solution2.SetE(solution2.P());
        solutions.emplace_back(toPtEtaPhi(solution2));
    }
    return solutions;
}


std::vector<Jet>* ZPrimeTotTPrimeReconstruction::AK4cleaning(std::vector<Jet>* ak4, TopJet ak8){
  //  std::cout << "size ak4 jets "<<ak4->size()<<std::endl; 
  std::vector<Jet>* ak4jets(new std::vector<Jet> (*ak4)); 
  ak4jets->clear();
  ak4jets->reserve(ak4->size());
  
  for(unsigned int j=0; j<ak4->size(); j++){
    double deltar = deltaR(ak4->at(j),ak8);
    //  std::cout<<"deltaR " <<deltar<<std::endl;
    if(deltar > 1.2){
      //  std::cout<< "setze ak4"<< std::endl;
      ak4jets->push_back(ak4->at(j));
    }
  }
 
  return ak4jets;
}



void ZPrimeTotTPrimeReconstruction::LepTop(std::vector<Jet>* ak4,  std::vector<LorentzVector> neutrinos, const Particle & lepton , ZPrimeTotTPrimeReconstructionHypothesis hyp,std::vector<TopJet> Tag ){

  unsigned int n_jets = ak4->size();
  //  std::cout << "n_jets "<<n_jets <<std::endl;

  //LEP. TOP RECONSTRUCTION
  if(n_jets>8) n_jets=8;
  const unsigned int max_j = pow(2, n_jets);// (0)leptonisches Top, (1)nichts
    
  for(const auto & neutrino_p4 : neutrinos){

    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
    for (unsigned int j=0; j < max_j; j++) {
      
      //EXPERIMENT
    hyp.clear_subjets();
    for(const Jet s:Tag.at(0).subjets()) hyp.set_subjets(s);
    //ENDE

      LorentzVector toplep_v4 = wlep_v4;
      int lepjets=0;
      int num = j;

      
      hyp.set_lepton(lepton);
      hyp.set_neutrino_v4(neutrino_p4);
      // clear hyp.add
      hyp.clear_toplep_jet();
      for (unsigned int k=0; k<n_jets; k++) {

	if(num%4==0) {
	  toplep_v4 = toplep_v4 +ak4->at(k).v4();
	  hyp.add_toplep_jet(ak4->at(k));
	  lepjets++;
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
      //   std::cout<<"lepjets "<<lepjets<<std::endl;
      if( lepjets>0 ) {
	hyp.set_toplep_v4(toplep_v4);
	//	std::cout<<"Setze Hyp LepTop"<<std::endl;
	recoHyps.emplace_back(move(hyp));
      }

    } //combinations
  }
}

void ZPrimeTotTPrimeReconstruction::LepHadTop(std::vector<Jet>* ak4,  std::vector<LorentzVector> neutrinos, const Particle & lepton , ZPrimeTotTPrimeReconstructionHypothesis hyp, std::vector<TopJet> Tag){

  unsigned int n_jets = ak4->size();


  //LEP. TOP RECONSTRUCTION
  if(n_jets>8) n_jets=8;
  const unsigned int max_j = pow(3, n_jets);// (0)leptonisches Top, (1)had Top, (2) nichts
    
  for(const auto & neutrino_p4 : neutrinos){
    
    const LorentzVector wlep_v4 = lepton.v4() + neutrino_p4;
    for (unsigned int j=0; j < max_j; j++) {
      
      //EXPERIMENT
      hyp.clear_subjets(); 
      for(const Jet s:Tag.at(0).subjets()) hyp.set_subjets(s);
      //ENDE


      LorentzVector toplep_v4 = wlep_v4;
      LorentzVector tophad_v4;
      int lepjets=0;
      int hadjets=0;
      int num = j;

      
      hyp.set_lepton(lepton);
      hyp.set_neutrino_v4(neutrino_p4);
      // clear hyp.add
      hyp.clear_toplep_jet();
      hyp.clear_tophad_jet();
      for (unsigned int k=0; k<n_jets; k++) {

	if(num%4==0) {
	  toplep_v4 = toplep_v4 +ak4->at(k).v4();
	  hyp.add_toplep_jet(ak4->at(k));
	  lepjets++;
	}
	if(num%4==1) {
	  tophad_v4 = tophad_v4 +ak4->at(k).v4();
	  hyp.add_tophad_jet(ak4->at(k));
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
     
      if( lepjets>0 && hadjets>0) {
	//	std::cout<<"Setze Hyp LepHadTop"<<std::endl;
       	hyp.set_toplep_v4(toplep_v4);
	hyp.set_tophad_v4(tophad_v4);
	//	std::cout << "in REconstruction hyp setzen in LEp Had" << hyp.HZW_subjets().size()<< std::endl;

	recoHyps.emplace_back(move(hyp));
 
      }

    } //combinations
  }
 
}
