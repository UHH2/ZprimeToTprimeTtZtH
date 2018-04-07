#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h"
#include <limits>
#include <cmath>
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include <stdexcept>
#include <iostream>
#include <memory>
#include <fstream>
#include <sstream>

#include "UHH2/core/include/Hists.h"
using namespace uhh2examples;
using namespace uhh2;
using namespace std;

////////////////////////////////////////////////////////                                                                                                                             

ElectronTriggerSF::ElectronTriggerSF(Context & ctx){
  sysdirection=ctx.get("Systematic_EleTrigger");
}

bool   ElectronTriggerSF::process(uhh2::Event & event){
  Particle elec = event.electrons->at(0);
  double elec_pt=elec.pt();
  double elec_eta =abs(elec.eta());
  double add_unc=0.01;

  if(elec_pt<=140){
    sf=0.9802;
    if(sysdirection == "up")sf=sqrt(pow(0.9834,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.977,2)-pow(add_unc,2));

  }else if(elec_pt<=160){
    sf=0.9873;
    if(sysdirection == "up")sf=sqrt(pow(0.9909,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9835,2)-pow(add_unc,2));

  }else if(elec_pt<=180){
    sf=0.9867;
    if(sysdirection == "up")sf=sqrt(pow(0.9911,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9818,2)-pow(add_unc,2));

  }else if(elec_pt<=200){
    sf =0.9905;
    if(sysdirection == "up")sf=sqrt(pow(0.9959,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9845,2)-pow(add_unc,2));
  }else if(elec_pt<=250){
    sf =0.979421;
    if(sysdirection == "up")sf=sqrt(pow(0.9852,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9731,2)-pow(add_unc,2));
  }else if(elec_pt<=400){
    sf=0.996993;
    if(sysdirection == "up")sf=sqrt(pow(0.9989,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9939,2)-pow(add_unc,2));

  }else{
    sf=1;
    if(sysdirection == "up")sf=sqrt(pow(1.001,2) + pow(add_unc,2));
    if(sysdirection == "down")sf=sqrt(pow(0.9837,2)-pow(add_unc,2));

  }
  double weight= event.weight *sf;
  event.weight = weight;


  return true;
}

////////////////////////////////////////////////////////                                                                                                                            
                                                                                                                                                                                     

MuonTriggerSF::MuonTriggerSF(Context & ctx){
  sysdirection=ctx.get("Systematic_MuonTrigger");
}

bool   MuonTriggerSF::process(uhh2::Event & event){
  Particle muon = event.muons->at(0);
  double muon_pt=muon.pt();
  double muon_eta =abs(muon.eta());

  //pT eta dependency
  if(muon_pt<=100){
    if(muon_eta<=0.9){
      sf = 0.972;
      if(sysdirection == "up")sf+= 0.002;
      if(sysdirection == "down")sf-= 0.002;
    }else if(muon_eta<=1.2){
      sf = 0.947;
      if(sysdirection == "up")sf+= 0.003;
      if(sysdirection == "down")sf-= 0.003;
    }else if(muon_eta<=1.7){
      sf = 1.002;
      if(sysdirection == "up")sf+= 0.004;
      if(sysdirection == "down")sf-= 0.004;
    }else if(muon_eta<=2.4){
      sf = 0.957;
      if(sysdirection == "up")sf+= 0.005;
      if(sysdirection == "down")sf-= 0.005;
    }
  }else if(muon_pt<=180){
    if(muon_eta<=0.9){
      sf = 0.975;
      if(sysdirection == "up")sf+= 0.003;
      if(sysdirection == "down")sf-= 0.003;
    }else if(muon_eta<=1.2){
      sf = 0.935;
      if(sysdirection == "up")sf+= 0.004;
      if(sysdirection == "down")sf-= 0.004;
    }else if(muon_eta<=1.7){
      sf = 0.982;
      if(sysdirection == "up")sf+= 0.006;
      if(sysdirection == "down")sf-= 0.006;
    }else if(muon_eta<=2.4){
      sf = 0.945;
      if(sysdirection == "up")sf+= 0.007;
      if(sysdirection == "down")sf-= 0.007;
    }
  }else if(muon_pt<=210){
    if(muon_eta<=0.9){
      sf = 0.937;
      if(sysdirection == "up")sf+= 0.011;
      if(sysdirection == "down")sf-= 0.011;
    }else if(muon_eta<=1.2){
      sf = 0.947;
      if(sysdirection == "up")sf+= 0.024;
      if(sysdirection == "down")sf-= 0.024;
    }else if(muon_eta<=1.7){
      sf = 0.954;
      if(sysdirection == "up")sf+= 0.022;
      if(sysdirection == "down")sf-= 0.022;
    }else if(muon_eta<=2.4){
      sf = 0.949;
      if(sysdirection == "up")sf+= 0.033;
      if(sysdirection == "down")sf-= 0.033;
    }
  }else if(muon_pt<=300){
    if(muon_eta<=0.9){
      sf = 0.951;
      if(sysdirection == "up")sf+= 0.017;
      if(sysdirection == "down")sf-= 0.017;
    }else if(muon_eta<=1.2){
      sf = 0.898;
      if(sysdirection == "up")sf+= 0.032;
      if(sysdirection == "down")sf-= 0.032;
    }else if(muon_eta<=1.7){
      sf = 0.906;
      if(sysdirection == "up")sf+= 0.031;
      if(sysdirection == "down")sf-= 0.031;
    }else if(muon_eta<=2.4){
      sf = 0.833;
      if(sysdirection == "up")sf+= 0.048;
      if(sysdirection == "down")sf-= 0.048;
    }
  }else if(muon_pt>300){
    if(muon_eta<=0.9){
      sf = 0.981;
      if(sysdirection == "up")sf+= 0.045;
      if(sysdirection == "down")sf-= 0.045;
    }else if(muon_eta<=1.2){
      sf = 0.884;
      if(sysdirection == "up")sf+= 0.048;
      if(sysdirection == "down")sf-= 0.048;
    }else if(muon_eta<=1.7){
      sf = 0.984;
      if(sysdirection == "up")sf+= 0.082;
      if(sysdirection == "down")sf-= 0.082;
    }else if(muon_eta<=2.4){
      sf = 0.775;
      if(sysdirection == "up")sf+= 0.086;
      if(sysdirection == "down")sf-= 0.086;
    }
  }


  double weight= event.weight *sf;
  event.weight = weight;


  return true;
}




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
ZPrimeTotTPrimePartonW::ZPrimeTotTPrimePartonW(){}

bool ZPrimeTotTPrimePartonW::passes(const uhh2::Event& event){
  std::vector<GenParticle> & genparticles = *event.genparticles;
  double Wpt=10000;

 
 for(const auto & gp : genparticles){
   if(abs(gp.pdgId())==24) Wpt = gp.pt();
 }
 // cout<<Wpt<<endl;
 return Wpt<100;
}
////////////////////////////////////////////////////////
ZPrimeTotTPrimePartonWCut::ZPrimeTotTPrimePartonWCut(TString filename_){
  filename = filename_;
}

bool ZPrimeTotTPrimePartonWCut::passes(const uhh2::Event& event){
  std::vector<GenParticle> & genparticles = *event.genparticles;
  double Wpt=10000;

 
 for(const auto & gp : genparticles){
   if(abs(gp.pdgId())==24) Wpt = gp.pt();
 }


 if(filename.Contains("100") && filename.Contains("250")) return 100<Wpt && Wpt<250;
 else if(filename.Contains("250") && filename.Contains("400")) return 250<Wpt && Wpt<400;
 else if(filename.Contains("400") && filename.Contains("600")) return 400<Wpt && Wpt<600;
 else if(filename.Contains("600")) return 600<Wpt;
 else if(filename.Contains("50") && filename.Contains("100")) return 50<Wpt && Wpt<100;
 else return true;

 return true;
}

////////////////////////////////////////////////////////

ElectronTriggerWeights::ElectronTriggerWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: ElectronTriggerWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  
  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  Eff_lowpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_TTbar_eff"));					 
  Eff_highpt_MC.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_TTbar_eff"));					 
  Eff_lowpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_lowpt_eta_DATA_eff"));					 
  Eff_highpt_DATA.reset((TGraphAsymmErrors*)file->Get("gr_highpt_eta_DATA_eff"));					 
  
  if(SysDirection != "nominal" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
  
}

bool ElectronTriggerWeights::process(Event & event){

  if(event.isRealData) return true;

  //cout << "Weight before SF: " << event.weight << endl;
  double prob_notrig_mc = 1, prob_notrig_data = 1;
  for(const auto & ele : *event.electrons){
    double eta = ele.eta();
    if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Ele-|eta| > 2.4 is not supported at the moment.");


    //find right bin in eta
    int idx = 0;
    bool lowpt = false;
    if(30 <= ele.pt() && ele.pt() < 120){
      lowpt = true;
      //lowpt trigger
      bool keep_going = true;
      while(keep_going){
	double x,y;
	Eff_lowpt_MC->GetPoint(idx,x,y);
	keep_going = eta > x + Eff_lowpt_MC->GetErrorXhigh(idx);
	if(keep_going) idx++;
      }
    }
    else if(ele.pt() >= 120){
     //highpt trigger
      bool keep_going = true;
      while(keep_going){
	double x,y;
	Eff_highpt_MC->GetPoint(idx,x,y);
	keep_going = eta > x + Eff_highpt_MC->GetErrorXhigh(idx);
	if(keep_going) idx++;
      }
    }
    else throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Electron has pt<30. Clean electron collection before applying weights.");

    //access efficiencies for MC and DATA, possibly accout for systematics = statistical + add. 2% up/down
    double eff_data = -1, eff_mc = -1, dummy_x;
    double stat_data = -1, stat_mc = -1, tp = 0.02, total_syst_data = -1, total_syst_mc = -1;
    if(lowpt){
      Eff_lowpt_MC->GetPoint(idx,dummy_x,eff_mc);
      Eff_lowpt_DATA->GetPoint(idx,dummy_x,eff_data);

      if(SysDirection == "up"){		
	stat_mc = Eff_lowpt_MC->GetErrorYlow(idx);	
	stat_data = Eff_lowpt_DATA->GetErrorYhigh(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));

      	eff_mc -= total_syst_mc;    	
      	eff_data += total_syst_data;	
      }							
      else if(SysDirection == "down"){
	stat_mc = Eff_lowpt_MC->GetErrorYhigh(idx);	
	stat_data = Eff_lowpt_DATA->GetErrorYlow(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
      	eff_mc += Eff_lowpt_MC->GetErrorYhigh(idx);    	
      	eff_data -= Eff_lowpt_DATA->GetErrorYlow(idx);	
      }                                                 
    }
    else{
      Eff_highpt_MC->GetPoint(idx,dummy_x,eff_mc);
      Eff_highpt_DATA->GetPoint(idx,dummy_x,eff_data);

      if(SysDirection == "up"){	
	stat_mc = Eff_highpt_MC->GetErrorYlow(idx);	
	stat_data = Eff_highpt_DATA->GetErrorYhigh(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
			
	eff_mc -= Eff_highpt_MC->GetErrorYlow(idx);    	
	eff_data += Eff_highpt_DATA->GetErrorYhigh(idx);	
      }							
      else if(SysDirection == "down"){	
	stat_mc = Eff_highpt_MC->GetErrorYhigh(idx);	
	stat_data = Eff_highpt_DATA->GetErrorYlow(idx);
	total_syst_mc = sqrt(pow(stat_mc,2) + pow(tp,2));
	total_syst_data = sqrt(pow(stat_data,2) + pow(tp,2));
					
	eff_mc += Eff_highpt_MC->GetErrorYhigh(idx);    	
	eff_data -= Eff_highpt_DATA->GetErrorYlow(idx);	
      }                                                 
    }

    //multiply to the efficiency for not triggering
    prob_notrig_mc *= 1-eff_mc;
    prob_notrig_data *= 1-eff_data;

    //cout << "Efficiency for this ele -- MC: " << eff_mc << ", DATA" << eff_data << endl;
    //cout << "prob for not triggering -- MC: " << prob_notrig_mc << ", DATA: " << prob_notrig_data << endl;

  }

  //Scale weight by (1-prob_notrig_data) / (1-prob_notrig_mc)
  double SF = (1-prob_notrig_data)/(1-prob_notrig_mc);
  event.weight *= SF;

  //cout << "Weight with SF: " << event.weight << endl; 
  return true;
}


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

ZPrimeTotTPrimeDRele::ZPrimeTotTPrimeDRele(float mmin_, float mmax_): mmin(mmin_), mmax(mmax_) {}

bool ZPrimeTotTPrimeDRele::passes(const uhh2::Event& event){
  std::vector<Electron>* electrons = event.electrons;
  std::vector<float> drmin_buf;
  for(const auto & ele : *event.electrons){
        
    if(event.jets){
      auto nj = nextJet(ele, *event.jets);
      auto drmin_val = nj ? deltaR(ele, *nj) : numeric_limits<float>::infinity();
      drmin_buf.push_back(drmin_val);
        
    }
  }



  return ( drmin_buf.back() >= mmin && (mmax<0||drmin_buf.back() <= mmax));
   
}

////////////////////////////////////////////////////////




ZPrimeTotTPrimeMETCut::ZPrimeTotTPrimeMETCut(float min_MET_, float max_MET_):min_MET(min_MET_),max_MET(max_MET_){}

bool ZPrimeTotTPrimeMETCut::passes(const Event & event){
  assert(event.met);
  float MET = event.met->pt();
  //  return(MET > min_MET && (MET < max_MET || max_MET <0));
return(MET > min_MET );
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

  // for(const auto & jet : *jets){
  //   ht_jets += jet.pt();
  //   }
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
  pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  // pass = ht > ht_min;
  return pass;
}

HtJetsSelection::HtJetsSelection(double ht_min_, double ht_max_):ht_min(ht_min_), ht_max(ht_max_){}
bool HtJetsSelection::passes(const Event & event){
  auto met = event.met->pt();
 std::vector<Jet>* jets = event.jets;

  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;

  for(const auto & jet : *jets){
    ht_jets += jet.pt();
    }
  // for(const auto & electron : *event.electrons){
  //   ht_lep += electron.pt();
  // }
  // for(const auto & muon : *event.muons){
  //   ht_lep += muon.pt();
  // }
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
    }*/

  ht = ht_lep + ht_jets;

  bool pass = false;
  pass = ht > ht_min && (ht_max < 0 || ht < ht_max);
  // pass = ht > ht_min;
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

const Particle* leading_lepton(const uhh2::Event& event){

  const Particle* lep(0);

  float ptL_max(0.);
  if(event.muons)    { for(const auto& mu : *event.muons)    { if(mu.pt() > ptL_max){ ptL_max = mu.pt(); lep = &mu; } } }
  if(event.electrons){ for(const auto& el : *event.electrons){ if(el.pt() > ptL_max){ ptL_max = el.pt(); lep = &el; } } }

  if(!lep) throw std::runtime_error("leading_lepton -- pt-leading lepton not found");

  return lep;
}

// // //////////////////////////////////////////////////////// 
bool uhh2::TwoDCut::passes(const uhh2::Event& event){

  assert((event.muons||event.electrons) && event.jets);
  // std::cout<<"number of electrons "<<event.electrons->size() << "  number of muons "<< event.muons->size()<<"  number of AK4 jets  "<< event.jets->size()<<std::endl;
  if((event.muons->size() < 1) && (event.electrons->size() < 1) ){
    std::cout << "\n @@@ WARNING -- TwoDCut::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  const Particle* lepton = leading_lepton(event);

  float drmin, ptrel;
  std::tie(drmin, ptrel) = drmin_pTrel(*lepton, *event.jets);

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
    if(event.electrons){
      for(const auto & muo : *event.electrons)
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
    subjets = topjet.subjets();
    // std::cout<<"subjets size == 2 "<< subjets.size() <<std::endl;
    sort_by_pt(subjets);

    LorentzVector firsttwosubjets = subjets[0].v4() + subjets[1].v4();
    if(!firsttwosubjets.isTimelike()) {
        return false;
    }
    auto mjet = firsttwosubjets.M();
    if(mjet < minmass_) return false;
    if(mjet > maxmass_) return false;
    // std::cout<<"ZPrimeHiggstag returnes true, mass jet "<< mjet<<std::endl;
    return true;
}


uhh2::TriangularCuts::TriangularCuts(float a, float b): a_(a), b_(b) {

  if(!b_) std::runtime_error("TriangularCuts -- incorrect initialization (parameter 'b' is null)");
}

bool uhh2::TriangularCuts::passes(const uhh2::Event& event){

  assert(event.muons || event.electrons);
  assert(event.jets && event.met);

  if((event.muons->size()+event.electrons->size()) != 1){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of muons+electrons in the event (!=1). returning 'false'\n";
    return false;
  }

  if(!event.jets->size()){
    std::cout << "\n @@@ WARNING -- TriangularCuts::passes -- unexpected number of jets in the event (==0). returning 'false'\n";
    return false;
  }

  // pt-leading charged lepton
  const Particle* lep1 = leading_lepton(event);

  // 1st entry in jet collection (should be the pt-leading jet)
  const Particle* jet1 = &event.jets->at(0);

  // MET-lepton triangular cut
  bool pass_tc_lep = fabs(fabs(deltaPhi(*event.met, *lep1)) - a_) < a_/b_ * event.met->pt();

  // MET-jet triangular cut
  bool pass_tc_jet = fabs(fabs(deltaPhi(*event.met, *jet1)) - a_) < a_/b_ * event.met->pt();

  return pass_tc_lep && pass_tc_jet;
}



MuonTrkWeights::MuonTrkWeights(Context & ctx, TString path_, TString SysDirection_): path(path_), SysDirection(SysDirection_){

  auto dataset_type = ctx.get("dataset_type");
  bool is_mc = dataset_type == "MC";
  if(!is_mc){
    cout << "Warning: MuonTrkWeights will not have an effect on this non-MC sample (dataset_type = '" + dataset_type + "')" << endl;
    return;
  }
  
  unique_ptr<TFile> file;												 
  file.reset(new TFile(path,"READ"));											 
  
  Trk_SF.reset((TGraphAsymmErrors*)file->Get("ratio_eff_eta3_dr030e030_corr"));					 
  
  if(SysDirection != "none" && SysDirection != "up" && SysDirection != "down") throw runtime_error("In LQToTopMuModules.cxx, ElectronTriggerWeights.process(): Invalid SysDirection specified.");
  
}

bool MuonTrkWeights::process(Event & event){

  if(event.isRealData) return true;

  double SF = 1.0;
  for(const auto & mu : *event.muons){
    double eta = mu.eta();
    if(fabs(eta) > 2.4) throw runtime_error("In LQToTopMuModules.cxx, MuonTrkWeights::process(): Mu-|eta| > 2.4 is not supported at the moment.");
    
    //find right bin in eta
    int idx = 0;
    bool keep_going = true;
    while(keep_going){
      double x,y;
      Trk_SF->GetPoint(idx,x,y);
      keep_going = eta > x + Trk_SF->GetErrorXhigh(idx);
      if(keep_going) idx++;
    }
    
    double eff_fnl = 1., dummy_x;
    Trk_SF->GetPoint(idx,dummy_x,eff_fnl);
    
    SF *= eff_fnl;
    if(SysDirection == "up"){
      SF *= 1.005;
    }
    if(SysDirection == "down"){
      SF *= 0.995;
    }
  
  }

  event.weight *= SF;

  return true;

}

MistagRateSF::MistagRateSF(Context & ctx,TString file_,TString SysDirection_ ):SysDirection(SysDirection_){

  h_toptag = ctx.get_handle<std::vector<TopJet>>("TopTag");
  h_higgstag = ctx.get_handle<std::vector<TopJet>>("HiggsTag");
  h_higgstag_1b = ctx.get_handle<std::vector<TopJet>>("HiggsTag_one_btag");
  h_ZWtag = ctx.get_handle<std::vector<TopJet>>("ZWTag");
  //define vetor with the different functions that are used to calculate the error
  //[[func1,params],[func2,params],...]
  ifstream test(file_);
  if (test.is_open()) {
    while (!test.eof() ) {
      std::string line;
      getline (test,line);
      std::istringstream iss(line);
      std::vector<std::string> params;
      do {
	std::string sub;
  	iss >> sub; 
	params.push_back(sub);
      } 
      while (iss);
      input.push_back(params);
    }
  }

  test.close();
  berror =false;

}

bool MistagRateSF::process(Event & event, TString tagger){

  TopJet topjet;
  double error;
  if(tagger == "higgs"){
    const auto & higgstagevt = event.get(h_higgstag);
    topjet = higgstagevt.at(0);
  }else if(tagger=="higgs_one"){
    const auto & higgstagevt = event.get(h_higgstag_1b);
    topjet = higgstagevt.at(0);
  }else if(tagger=="zw"){
    const auto & higgstagevt = event.get(h_ZWtag);
    topjet = higgstagevt.at(0);
  }else if(tagger=="top"){
    const auto & higgstagevt = event.get(h_toptag);
   topjet = higgstagevt.at(0);
 }else throw runtime_error("in ZPrimeTotTPrimeSelections tagger not spezified. Only 'higgs', 'higgs_one' , 'zw', 'top' are spezified.");
  
  if(berror)std::cout<<"============================================= Ganz am Anfang"<<std::endl;
  double pt = topjet.pt();
  if(berror)std::cout<<"pt  "<<pt<<std::endl;
  double resulting_error=0;
  Double_t rangemin = stod(input[0][1]);
  Double_t rangemax = stod(input[0][2]);
  TString function = input[0][0];

  //define first default function
  TF1* func_default = new TF1("Ratio",function ,rangemin,  rangemax);
  //Hier
  func_default->FixParameter(0,stod(input[0][3]));

  resulting_error = pow(stod(input[0][4]),2);
  if(berror)std::cout<< "error_default  " << stod(input[0][4])<<std::endl;

  //define one function for each other function than the default in vector
  for(unsigned int i=1;i<input.size()-1;i++){
    Double_t rangemin = stod(input[i][1]);
    Double_t rangemax = stod(input[i][2]);
    TString function =input[i][0] ;
    TF1* func = new TF1("Ratio",function ,rangemin,  rangemax);
    for(unsigned int k=3;k<input[i].size()-1;k++)  func->FixParameter(k-3,stod(input[i][k]));
    
    if(berror)std::cout<<"Eval func_Default  "<<func_default->Eval(pt)<<"  Eval func  "<<func->Eval(pt)<<std::endl;
    if(berror)std::cout<<"substract values such that |func1.eval-func2.eval|  "<<func_default->Eval(pt)-func->Eval(pt)<<std::endl;

    //Hier
    //substract values such that |func1.eval-func2.eval|^2 
    //When pt of topjet out of range use func->Eval(rangemax)
    if(pt>rangemax) error = pow(func_default->Eval(rangemax)-func->Eval(rangemax),2);
    else error = pow(func_default->Eval(pt)-func->Eval(pt),2);

    // add difference to error
    resulting_error += error;
  }
  if(berror)std::cout<<"resulting error with all errors  "<<resulting_error <<std::endl;
  //take sqrt for getting error
  resulting_error = sqrt(resulting_error);
  if(berror)std::cout<<"Final resulting error  "<<resulting_error <<std::endl;

  //apply one SF or SF plus error to event.weight
  double sf = stod(input[0][3]);
  if(berror)std::cout<<"SF  "<<sf<<std::endl;
  if(SysDirection == "up") sf+=resulting_error;
  if(berror)std::cout<<"SF up  "<<sf<<std::endl;
  if(SysDirection == "down") sf-=resulting_error;
  if(berror)std::cout<<"SF down "<<sf<<std::endl;


  event.weight *= sf;
  return true;
}


IsMistag::IsMistag(Context & ctx){
  h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
  h_background = ctx.get_handle<BackgroundGen>("background");
  h_ttbar = ctx.get_handle<TTbarGen>("ttbar");

  h_toptag = ctx.get_handle<std::vector<TopJet>>("TopTag");
  h_higgstag = ctx.get_handle<std::vector<TopJet>>("HiggsTag");
  h_higgstag_1b = ctx.get_handle<std::vector<TopJet>>("HiggsTag_one_btag");
  h_ZWtag = ctx.get_handle<std::vector<TopJet>>("ZWTag");
}

//signal samples
bool IsMistag::process(uhh2::Event & event, TString tagger){

  //if(!event.is_valid(h_zprimegen) && !event.is_valid(h_background)&& !event.is_valid(h_ttbar))throw runtime_error("in ZPrimeTotTPrimeSelections GenSelection not spezified.");
  
 if(!(tagger=="higgs"||tagger=="higgs_one" || tagger== "zw" ||  tagger== "top"))throw runtime_error("in ZPrimeTotTPrimeSelections tagger not spezified. Only 'higgs', 'higgs_one' , 'zw', 'top' are spezified.");

 TopJet topjet;
 if(tagger == "higgs"){
   const auto & higgstagevt = event.get(h_higgstag);
   topjet = higgstagevt.at(0);
 }else if(tagger=="higgs_one"){
   const auto & higgstagevt = event.get(h_higgstag_1b);
   topjet = higgstagevt.at(0);
 }else if(tagger=="zw"){
   const auto & higgstagevt = event.get(h_ZWtag);
   topjet = higgstagevt.at(0);
 }else if(tagger=="top"){
   const auto & higgstagevt = event.get(h_toptag);
   topjet = higgstagevt.at(0);
 }else throw runtime_error("in ZPrimeTotTPrimeSelections tagger not spezified. Only 'higgs', 'higgs_one' , 'zw', 'top' are spezified.");

  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);
    double deltar_higgs_Q1 = deltaR(topjet,zprimegen.HQ1());
    double deltar_higgs_Q2 = deltaR(topjet,zprimegen.HQ2());
    double deltar_higgs = deltaR(topjet,zprimegen.Higgs());
    double deltar_z_Q1 = deltaR(topjet,zprimegen.ZQ1());
    double deltar_z_Q2 = deltaR(topjet,zprimegen.ZQ2());
    double deltar_w_Q1 = deltaR(topjet,zprimegen.WQ1());
    double deltar_w_Q2 = deltaR(topjet,zprimegen.WQ2());
    double deltar_top_bhad = deltaR(topjet,zprimegen.BHad());
    double deltar_top_WQ1 = deltaR(topjet,zprimegen.WHadQ1());
    double deltar_top_WQ2 = deltaR(topjet,zprimegen.WHadQ2());

   


    if((tagger == "higgs") && deltar_higgs_Q1<=0.8 && deltar_higgs_Q2<=0.8) return false;
    else if((tagger== "zw") &&( (deltar_z_Q1<=0.8 &&deltar_z_Q2<=0.8) || (deltar_w_Q1<=0.8 &&deltar_w_Q2<=0.8) ) ) return false;
    else if((tagger== "top") && deltar_top_bhad <=0.8 && deltar_top_WQ1<=0.8 && deltar_top_WQ2<=0.8) return false;
    else return true;
    
  }// valid h_zprimegen
  else if (event.is_valid(h_background)){
    const auto & background = event.get(h_background);
    double deltar_z_Q1 = deltaR(topjet,background.Quark1_v4());
    double deltar_z_Q2 = deltaR(topjet,background.Quark2_v4());
 
    if((tagger== "zw") && (deltar_z_Q1<=0.8 &&deltar_z_Q2<=0.8)  ) return false;
    else return true;
  }//valid h_backgrounds
  else if(event.is_valid(h_ttbar)){
    const auto & ttbar = event.get(h_ttbar);
    double deltar_w_Q1 = deltaR(topjet,ttbar.Q1());
    double deltar_w_Q2 = deltaR(topjet,ttbar.Q2());
    double deltar_top_bhad = deltaR(topjet,ttbar.BHad());
    if((tagger== "zw") && deltar_w_Q1<=0.8 && deltar_w_Q2 <=0.8)return false;
    else if(tagger == "top" && deltar_w_Q1<=0.8 && deltar_w_Q2 <=0.8 && deltar_top_bhad <=0.8) return false;
    else return true;

  }//valid h_ttbar
  
   return true;
}
