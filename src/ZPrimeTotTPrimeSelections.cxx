#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h"
#include <limits>
#include <cmath>
#include "TH2.h"
#include "TFile.h"
#include "TH1F.h"
#include <stdexcept>

using namespace uhh2examples;
using namespace uhh2;


ZPrimeTotTPrimeDijetSelection::ZPrimeTotTPrimeDijetSelection(float dphi_min_, float third_frac_max_): dphi_min(dphi_min_), third_frac_max(third_frac_max_){}
    
bool ZPrimeTotTPrimeDijetSelection::passes(const Event & event){
    assert(event.jets); // if this fails, it probably means jets are not read in
    if(event.jets->size() < 2) return false;
    const auto & jet0 = event.jets->at(0);
    const auto & jet1 = event.jets->at(1);
    auto dphi = deltaPhi(jet0, jet1);
    if(dphi < dphi_min) return false;
    if(event.jets->size() == 2) return true;
    const auto & jet2 = event.jets->at(2);
    auto third_jet_frac = jet2.pt() / (0.5 * (jet0.pt() + jet1.pt()));
    return third_jet_frac < third_frac_max;
}


// MuonptSlection:: MuonptSlection(float ptmin){}

// bool  MuonptSlection::passes(const Event & event){
//   assert(event.muons);
//   if(event.muons)
// }

////////////////////////////////////////////////////////

ZPrimeTotTPrimeNJetCut::ZPrimeTotTPrimeNJetCut(int nmin_, int nmax_, float ptmin1_,float ptmin2_, float etamax_):
  nmin(nmin_), nmax(nmax_), ptmin1(ptmin1_),ptmin2(ptmin2_), etamax(etamax_) {}

bool ZPrimeTotTPrimeNJetCut::passes(const uhh2::Event& event){
  assert(event.jets);
  std::vector<Jet>* jets = event.jets;

  int njet(0);
  if(jets->at(0).pt()>= ptmin1){
    for(auto & jet : *jets){
      if(fabs(jet.eta()) < etamax) ++njet;
    }
  }
  if(njet >= 2){
    if(jets->at(1).pt()>= ptmin2){return false;}
  }

  return (njet >= nmin) && (njet <= nmax || nmax < 0);
}


////////////////////////////////////////////////////////


////////////////////////////////////////////////////////

ZPrimeTotTPrimeNTopJetCut::ZPrimeTotTPrimeNTopJetCut(int nmin_, int nmax_, float ptmin1_,float ptmin2_, float etamax_):
  nmin(nmin_), nmax(nmax_), ptmin1(ptmin1_),ptmin2(ptmin2_), etamax(etamax_) {}

bool ZPrimeTotTPrimeNTopJetCut::passes(const uhh2::Event& event){

  std::vector<TopJet>* jets = event.topjets;

  int njet(0);
  if(jets->at(0).pt()<= ptmin1){
    for(auto & jet : *jets){
      if(fabs(jet.eta()) < etamax) ++njet;
    }
    if(njet>=2){
      if(!(jets->at(1).pt()>= ptmin2)){return false;}
    }
   
  }
  return (njet >= nmin) && (njet <= nmax || nmax < 0);
}


////////////////////////////////////////////////////////

ZPrimeTotTPrimeMassCut::ZPrimeTotTPrimeMassCut(float mmin_, float mmax_): mmin(mmin_), mmax(mmax_) {}

bool ZPrimeTotTPrimeMassCut::passes(const uhh2::Event& event){
  std::vector<Jet>* jets = event.jets;
  if(jets->at(0).v4().M() <= mmax && jets->at(0).v4().M() >= mmin ){
    return true;  
  }else{
    return false;
  }
}
////////////////////////////////////////////////////////

ZPrimeTotTPrimeTopMassCut::ZPrimeTotTPrimeTopMassCut(float mmin_, float mmax_): mmin(mmin_), mmax(mmax_) {}

bool ZPrimeTotTPrimeTopMassCut::passes(const uhh2::Event& event){
  std::vector<TopJet>* jets = event.topjets;
  if(jets->at(0).v4().M() <= mmax && jets->at(0).v4().M() >= mmin ){
    return true;  
  }else{
    return false;
  }
}
////////////////////////////////////////////////////////




ZPrimeTotTPrimeMETCut::ZPrimeTotTPrimeMETCut(float min_MET_, float max_MET_):min_MET(min_MET_),max_MET(max_MET_){}

bool ZPrimeTotTPrimeMETCut::passes(const Event & event){
  assert(event.met);
  float MET = event.met->pt();
  //  return(MET > min_MET && (MET < max_MET || max_MET <0));
return(MET < min_MET );
}

////////////////////////////////////////////////////////
ZPrimeTotTPrimePtJetCut::ZPrimeTotTPrimePtJetCut(float min_pt_, float max_pt_, unsigned int jetnumber_):min_pt(min_pt_),max_pt(max_pt_), jetnumber(jetnumber_){}

bool ZPrimeTotTPrimePtJetCut::passes(const Event & event){
 std::vector<Jet>* jets = event.jets;
  assert(jets);
  if (jets->size()>jetnumber){
    const auto & jet0 = jets->at(jetnumber);
    float PT=jet0.pt() ;
    return(PT > min_pt && PT < max_pt);
  }
  else return false;
}

////////////////////////////////////////////////////////

ZPrimeTotTPrimePtTopJetCut::ZPrimeTotTPrimePtTopJetCut(float min_pt_, float max_pt_, unsigned int jetnumber_):min_pt(min_pt_),max_pt(max_pt_), jetnumber(jetnumber_){}

bool ZPrimeTotTPrimePtTopJetCut::passes(const Event & event){
 std::vector<TopJet>* jets = event.topjets;
  assert(jets);
  if (jets->size()>jetnumber){
    const auto & jet0 = jets->at(jetnumber);
    float PT=jet0.pt() ;
    return(PT > min_pt && PT < max_pt);
  }
  else return false;
}

////////////////////////////////////////////////////////



// see https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation53XReReco
float btagging::csv_threshold(const csv_wp & wp){
    using namespace btagging;
    switch(wp){
        case csv_wp::loose: return 0.244f;
        case csv_wp::medium: return 0.679f;
        case csv_wp::tight: return 0.898f;
        default : return 0.679f;
    }
    // This should never happen; even if, the coompiler should warn in the switch.
    // But to avoid a compiler warning that no value is returned, include this line:
    // throw invalid_argument("unknown working point given to btagging::csv_threshold");
}

NBTagSelection::NBTagSelection(int nmin_, int nmax_, btagging::csv_wp wp): nmin(nmin_), nmax(nmax_), min_csv(btagging::csv_threshold(wp)){}

bool NBTagSelection::passes(const Event & event){
 std::vector<Jet>* jets = event.jets;
    int nbtag = 0;
    for(const Jet & j : *jets){
        if(j.btag_combinedSecondaryVertex() >= min_csv) ++nbtag;
    }
    return nbtag >= nmin && (nmax < 0 || nbtag <= nmax);
}



HtSelection::HtSelection(double ht_min_, double ht_max_):ht_min(ht_min_), ht_max(ht_max_){}
bool HtSelection::passes(const Event & event){
  auto met = event.met->pt();
 std::vector<Jet>* jets = event.jets;

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;

  for(const auto & jet : *jets){
    ht_jets += jet.pt();
    }
  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
    }*/

  ht = ht_lep + ht_jets + met;

  bool pass = false;
  // pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  pass = ht > ht_min;
  return pass;
}


////////////////////////////////////////////////////////

uhh2::GenMttbarCut::GenMttbarCut(uhh2::Context& ctx, const float mtt_min, const float mtt_max, const std::string& ttgen_name):
  mttbar_min_(mtt_min), mttbar_max_(mtt_max), h_ttbargen_(ctx.get_handle<TTbarGen>(ttgen_name)) {}

bool uhh2::GenMttbarCut::passes(const uhh2::Event& event){

  const TTbarGen& ttbargen = event.get(h_ttbargen_);

  if(ttbargen.DecayChannel() == TTbarGen::e_notfound)
    throw std::runtime_error("GenMttbarCut::passes -- undefined decay-channel for TTbarGen object");

  const float mttbargen = (ttbargen.Top().v4() + ttbargen.Antitop().v4()).M();

  return (mttbar_min_ < mttbargen) && (mttbargen < mttbar_max_);
}
////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////
ZPrimeTotTPrimeDeltaRCut::ZPrimeTotTPrimeDeltaRCut(uhh2::Context& ctx ,float deltar_,const std::string & discriminator_name ):deltar(deltar_),name(discriminator_name),h_hyps(ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>("ZPrimeTotTPrimeReconstruction")),h_zprimegen( ctx.get_handle<ZPrimeGen>("zprimegen")){}

bool ZPrimeTotTPrimeDeltaRCut::passes(const Event & event){
  std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = event.get(h_hyps);
  const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, name );
  const auto & zprimegen = event.get(h_zprimegen);
  if(deltaR(zprimegen.Higgs(), hyp->HZW_v4() )<=0.6){return true;}
  else return false;


	 }

// // //////////////////////////////////////////////////////// 
// ZPrimeTotTPrimeAK4cleaner::ZPrimeTotTPrimeAK4cleaner(float deltar_  ):deltar(deltar_){}

// bool ZPrimeTotTPrimeAK4cleaner::operator()(const Particle & p,const Event & event)const{
//  std::vector<TopJet>* Topjets = event.topjets;
   
//  float delta = deltaR(Topjets ->at(0), p);
//  return delta > deltar;
 
// }
// // //////////////////////////////////////////////////////// 

// ////////////////////////////////////////////////////////
ZPrimeTotTPrimenumbersub::ZPrimeTotTPrimenumbersub(unsigned int num_):num(num_){}

bool ZPrimeTotTPrimenumbersub::passes(const Event & event){

  if(event.topjets->at(0).subjets().size()< num)return false;

  return true;
}

// // //////////////////////////////////////////////////////// 
// ////////////////////////////////////////////////////////
ZPrimeTotTPrimedrmin::ZPrimeTotTPrimedrmin(float dr_):dr(dr_){}

bool ZPrimeTotTPrimedrmin::passes(const Event & event){
  const auto & jets = *event.topjets;
  for(unsigned int i = 0; i <jets.size(); i++){
    const auto & jet = jets[i];
 
    auto next_jet = closestParticle(jet, jets);
    drmin = next_jet ? deltaR(jet, *next_jet) : std::numeric_limits<float>::infinity();

  }
  if(drmin < dr)return true;
  else return false;     
 
}

// // //////////////////////////////////////////////////////// 
bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert(event.muons && event.jets);
  if(event.muons->size() < 1){
    std::cout << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  float drmin, ptrel;  
  if(event.muons->size()) std::tie(drmin, ptrel) = drmin_pTrel(event.muons->at(0), *event.jets);
 

  return (drmin > min_deltaR_) || (ptrel > min_pTrel_);
}
 // //////////////////////////////////////////////////////// 
ZPrimeTotTPrimeChiCut::ZPrimeTotTPrimeChiCut(uhh2::Context& ctx,float chi_,const std::string & hyps_name,const std::string & discriminator_name):chi(chi_){
  m_discriminator_name  = discriminator_name;
  h_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(hyps_name);
}

bool ZPrimeTotTPrimeChiCut::passes(const Event & event){
 std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = event.get(h_hyps);
  const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
  double Chi_hyp = hyp->discriminator(m_discriminator_name);
  return (Chi_hyp < chi);
 


}
// //////////////////////////////////////////////////////// 
ZPrimeTotTPrimeRelIsoCut::ZPrimeTotTPrimeRelIsoCut(double relmin_):relisomin(relmin_){}

bool ZPrimeTotTPrimeRelIsoCut::passes(const Event & event){
  assert(event.muons);
  double reliso =event.muons->at(0).relIso();
 
  return (reliso < relisomin);

}
///////////////////////////////////////////////////////// 
bool TopJetLeptonDeltaRCleaner::process(uhh2::Event& event){

  assert(event.topjets);
  std::vector<TopJet> cleaned_topjets;

  for(const auto & tjet : *event.topjets){
    bool skip_tjet(false);

    if(event.muons){
      for(const auto & muo : *event.muons)
        if(uhh2::deltaR(tjet, muo) < minDR_) skip_tjet = true;
    }

    if(skip_tjet) continue;

    if(!skip_tjet) cleaned_topjets.push_back(tjet);
  }

  event.topjets->clear();
  event.topjets->reserve(cleaned_topjets.size());
  for(const auto& j : cleaned_topjets) event.topjets->push_back(j);

  return true;
}
// //////////////////////////////////////////////////////// 
ZPrimeTotTPrimeMuonPT::ZPrimeTotTPrimeMuonPT(double ptmax_):ptmax(ptmax_){}

bool ZPrimeTotTPrimeMuonPT::passes(const Event & event){
  assert(event.muons);
  if(event.muons->size()){
  double pt =event.muons->at(0).v4().Pt();
 
  return (pt < ptmax);
  }
  else return false;

}
///////////////////////////////////////////////////////// 
bool Tau32_inverted::operator()(const TopJet & topjet, const uhh2::Event &) const {
    auto tau2 = topjet.tau2();
    if(!std::isfinite(tau2) || tau2 == 0.0) return false;
    auto tau3 = topjet.tau3();
    if(!std::isfinite(tau3)) return false;
    return tau3 / tau2 > threshold;
}

///////////////////////////////////////////////////////////////////
bool UnHiggsTag::operator()(TopJet const & topjet, uhh2::Event const & event) const {
    auto subjets = topjet.subjets();
    if(subjets.size() < 2) return false;
    clean_collection(subjets, event, btagid_);
    if (subjets.size() < 2) return false;
    sort_by_pt(subjets);

    LorentzVector firsttwosubjets = subjets[0].v4() + subjets[1].v4();
    if(!firsttwosubjets.isTimelike()) {
        return false;
    }
    auto mjet = firsttwosubjets.M();
    if(mjet > minmass_) return false;
    if(mjet < maxmass_) return false;
    return true;
}
///////////////////////////////////////////////////////////////////
bool UnType2TopTag::operator()(const TopJet & topjet, const uhh2::Event & event) const {
   float mjet = 0.0;
   switch(m_typeOfMass)
    {
    case MassType::ungroomed : mjet=topjet.v4().M(); break;
    case MassType::groomed : mjet=m_groomed(topjet); break;
    default: std::cout<<"Type2TopTag Mass type not valid"<<std::endl;
    }
   if(m_SubjetId){
     bool pass_SubjetId = false;
     for(auto subjet : topjet.subjets()){
       if((*m_SubjetId)(subjet,event)){
	 pass_SubjetId = true;
	 break;
       }
     }
     if(!pass_SubjetId) return false;
   }
   if(mjet > m_mjetLower) return false;
   if(mjet < m_mjetUpper) return false;
   return true;
}

UnType2TopTag::UnType2TopTag(double mjetLower, double mjetUpper, MassType typeOfMass, boost::optional<JetId> SubjetId): m_mjetLower(mjetLower), m_mjetUpper(mjetUpper),  m_typeOfMass(typeOfMass), m_SubjetId(SubjetId){}

UnType2TopTag::UnType2TopTag(MassType typeOfMass): m_mjetLower(60.), m_mjetUpper(100.), m_typeOfMass(typeOfMass), m_SubjetId(boost::none){}

//////////////////////////////////////////////////////////// Scale factors ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ZPrimeTotTPrimeEff::ZPrimeTotTPrimeEff(uhh2::Context & ctx){
//   TString unc_name = ctx.get("unc_name");
//   TString unc_folder = "/nfs/dust/cms/user/abenecke/ZPrimeTotTPrime/25ns/rootfile/eff/hists/";
//   TFile* file = new TFile(unc_folder +"Eff_"+unc_name+".root", "READ");
//   TH1F* sf_hist = (TH1F*) file->Get("tot_eff_h");
//   SF = sf_hist->GetBinContent(1);

     
// }
    
// bool ZPrimeTotTPrimeEff::process(const Event & event){
//   cout << "zwtag "<<pass_zwtag << " higgs tag "<< pass_higgstag << endl;
//   cout<<"SF" << SF << " weight " << event.weight<<endl;
//   event.weight = event.weight * SF;
//   cout <<"weight "<<event.weight<<endl;
// }


bool ZPrimeTotTPrimeHiggsTag::operator()(TopJet const & topjet, uhh2::Event const & event) const {
    auto subjets = topjet.subjets();
    if(subjets.size() < 2) return false;
    clean_collection(subjets, event, btagid_);
    if (!(subjets.size() ==1)) return false;
    sort_by_pt(subjets);

    LorentzVector firsttwosubjets = subjets[0].v4() + subjets[1].v4();
    if(!firsttwosubjets.isTimelike()) {
        return false;
    }
    auto mjet = firsttwosubjets.M();
    if(mjet < minmass_) return false;
    if(mjet > maxmass_) return false;
    return true;
}
