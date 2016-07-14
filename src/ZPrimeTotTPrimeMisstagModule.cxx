#include <iostream>
#include <memory>

#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>
#include "UHH2/common/include/PrintingModules.h"
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/TriggerSelection.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/TopJetIds.h>
#include <UHH2/common/include/Utils.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/EventHists.h>
#include "UHH2/common/include/CollectionProducer.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include <UHH2/common/include/TriggerSelection.h>
#include "UHH2/common/include/LuminosityHists.h"
#include <fstream>

#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHists.h>
//#include <UHH2/ZPrimeTotTPrime/include/MistagHists.h>
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h>


using namespace uhh2examples;
using namespace uhh2;
using namespace std;

class ZPrimeTotTPrimeMisstagModule: public uhh2::AnalysisModule {
  
public:
  explicit ZPrimeTotTPrimeMisstagModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  
private:
  enum lepton { muon, elec };
  lepton channel_;
  std::unique_ptr<AnalysisModule> printer;
  
  //cleaner
  std::unique_ptr<TopJetCleaner>      topjet_cleaner;
  std::unique_ptr<TopJetCleaner>      mass_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner; 
  
  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;
  
  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::unique_ptr<uhh2::Selection> lumi_sel;

  //Selections
  std::unique_ptr<AndSelection>  metfilters_selection;
  std::unique_ptr<uhh2::Selection> topjet2_sel; //two AK8 pT>400GeV eta<2.4 deltaPhi>2.1
  std::unique_ptr<uhh2::Selection> tagger_sel;


  //Handle
  //before tagger set a handle on the different topjets: pt sorted first topjet,seconde topjet...
  uhh2::Event::Handle< std::vector<TopJet> > h_first_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_seconde_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_third_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_fourth_topjets;
  //after tagger set a handle on the different topjets on pT and pass tagger or not
  uhh2::Event::Handle< std::vector<TopJet> > h_first_pass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_seconde_pass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_third_pass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_fourth_pass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_first_nopass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_seconde_nopass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_third_nopass_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_fourth_nopass_topjets;
  // all topjets passed the tagger
  uhh2::Event::Handle< std::vector<TopJet> > h_zwtag_topjets;
  uhh2::Event::Handle< std::vector<TopJet> > h_nozwtag_topjets;


  //Hists
  //input
  std::unique_ptr<uhh2::Hists> input_h_event;
  std::unique_ptr<uhh2::Hists> input_h_muo;
  std::unique_ptr<uhh2::Hists> input_h_ele;
  std::unique_ptr<uhh2::Hists> input_h_jet;
  std::unique_ptr<Hists> input_h_eff;
  std::unique_ptr<Hists>input_h_topjet;
  std::unique_ptr<Hists> input_h_mistag;
 
  //Hist Topjet2
  std::unique_ptr<Hists> eff_topjet2_h;
  std::unique_ptr<Hists> jet_topjet2_h;
  std::unique_ptr<Hists> muon_topjet2_h;
  std::unique_ptr<Hists> event_topjet2_h;
  std::unique_ptr<Hists> topjet_topjet2_h;
  std::unique_ptr<Hists> tagger_topjet2_h;

  //before tagger fill hist on the different topjets: pt sorted first topjet,seconde topjet...
  std::unique_ptr<Hists> topjet_first_topjet2_h;
  std::unique_ptr<Hists> topjet_seconde_topjet2_h;
  std::unique_ptr<Hists> topjet_third_topjet2_h;
  std::unique_ptr<Hists> topjet_fourth_topjet2_h;


  //output
  std::unique_ptr<uhh2::Hists> output_h_event;
  std::unique_ptr<uhh2::Hists> output_h_muo;
  std::unique_ptr<uhh2::Hists> output_h_ele;
  std::unique_ptr<uhh2::Hists> output_h_jet;
  std::unique_ptr<Hists> output_h_eff;
  std::unique_ptr<Hists>output_h_topjet;
  std::unique_ptr<Hists> output_h_mistag;
  //after tagger fill hist on the different topjets on pT and zwtag
  std::unique_ptr<Hists> topjet_first_pass_zwtag_h;
  std::unique_ptr<Hists> topjet_seconde_pass_zwtag_h;
  std::unique_ptr<Hists> topjet_third_pass_zwtag_h;
  std::unique_ptr<Hists> topjet_fourth_pass_zwtag_h;
  std::unique_ptr<Hists> topjet_first_nopass_zwtag_h;
  std::unique_ptr<Hists> topjet_seconde_nopass_zwtag_h;
  std::unique_ptr<Hists> topjet_third_nopass_zwtag_h;
  std::unique_ptr<Hists> topjet_fourth_nopass_zwtag_h;
  // all topjets passed the tagger
  std::unique_ptr<Hists> topjet_pass_zwtag_h;
  std::unique_ptr<Hists> topjet_nopass_zwtag_h;

  //General
  std::string filename;
  // Reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  bool b_error;
};

ZPrimeTotTPrimeMisstagModule::ZPrimeTotTPrimeMisstagModule(uhh2::Context& ctx){ 

  //GenParticleprinter
  printer.reset(new GenParticlesPrinter(ctx));

  //corrector
  std::vector<std::string> JEC_AK4, JEC_AK8;
  const bool isMC = (ctx.get("dataset_type") == "MC");
   
  //// Data/MC scale
  if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx));
  } else     lumi_sel.reset(new LumiSelection(ctx));


  PrimaryVertexId pvid=StandardPrimaryVertexId();
  metfilters.emplace_back(new PrimaryVertexCleaner(pvid));
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("CSCTightHalo2015Filter", "Flag_CSCTightHalo2015Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  metfilters_selection->add<TriggerSelection>("chargedHadronTrackResolutionFilter", "Flag_chargedHadronTrackResolutionFilter"); 
  metfilters_selection->add<TriggerSelection>("muonBadTrackFilter", "Flag_muonBadTrackFilter");
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid);


  filename =  ctx.get("dataset_version");

  //calculator
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));

  //cleaner
  topjet_cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(200., 2.5))));
  const TopJetId massID = Type2TopTag(24,182,Type2TopTag::MassType::groomed);
  mass_cleaner.reset(new TopJetCleaner(ctx,massID));

  //Selections
  topjet2_sel.reset(new NTopJetSelection(2,-1,TopJetId(PtEtaCut( 250., 2.4)))); 
  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
  tagger_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
  
  //Genral
  const std::string ttbar_gen_label ("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else   genmttbar_sel.reset(new uhh2::AndSelection(ctx));

  //Handle
  //before tagger set a handle on the different topjets: pt sorted first topjet,seconde topjet...
  h_first_topjets= ctx.declare_event_output< std::vector<TopJet> > ("first_topjets");
  h_seconde_topjets= ctx.declare_event_output< std::vector<TopJet> > ("seconde_topjets");
  h_third_topjets= ctx.declare_event_output< std::vector<TopJet> > ("third_topjets");
  h_fourth_topjets= ctx.declare_event_output< std::vector<TopJet> > ("fourth_topjets");
  //after tagger set a handle on the different topjets on pT and pass tagger or not
  h_first_pass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("first_pass_topjets");
  h_seconde_pass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("seconde_pass_topjets");
  h_third_pass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("third_pass_topjets");
  h_fourth_pass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("fourth_pass_topjets");
  h_first_nopass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("first_nopass_topjets");
  h_seconde_nopass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("seconde_nopass_topjets");
  h_third_nopass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("third_nopass_topjets");
  h_fourth_nopass_topjets= ctx.declare_event_output< std::vector<TopJet> > ("fourth_nopass_topjets");
  // all topjets passed the tagger
  h_zwtag_topjets= ctx.declare_event_output< std::vector<TopJet> > ("zwtag_topjets");
  h_nozwtag_topjets= ctx.declare_event_output< std::vector<TopJet> > ("nozwtag_topjets");


  //Hist
  input_h_event.reset(new EventHists   (ctx, "input_Event"));
  input_h_muo.reset(new MuonHists (ctx, "input_Muons"));
  input_h_ele.reset(new ElectronHists(ctx, "input_Electrons"));
  input_h_jet.reset(new JetHists     (ctx, "input_Jets"));
  input_h_eff.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
  input_h_topjet.reset(new TopJetHists     (ctx, "input_TopJets")); 

  //Hists topjet2
  topjet_topjet2_h.reset(new TopJetHists(ctx, "topjet_topjet2"));
  eff_topjet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_topjet2"));
  jet_topjet2_h.reset(new JetHists(ctx, "jet_topjet2"));
  muon_topjet2_h.reset(new MuonHists(ctx, "muon_topjet2"));
  event_topjet2_h.reset(new EventHists(ctx, "event_topjet2"));
 

  //before tagger fill hist on the different topjets: pt sorted first topjet,seconde topjet...
  topjet_first_topjet2_h.reset(new TopJetHists(ctx, "topjet_first_topjet2",4,"first_topjets"));
  topjet_seconde_topjet2_h.reset(new TopJetHists(ctx, "topjet_seconde_topjet2",4,"seconde_topjets"));
  topjet_third_topjet2_h.reset(new TopJetHists(ctx, "topjet_third_topjet2",4,"third_topjets"));
  topjet_fourth_topjet2_h.reset(new TopJetHists(ctx, "topjet_fourth_topjet2",4,"fourth_topjets"));

  //Output hists
  output_h_event.reset(new EventHists   (ctx, "output_Event"));
  output_h_muo.reset(new MuonHists (ctx, "output_Muons"));
  output_h_ele.reset(new ElectronHists(ctx, "output_Electrons"));
  output_h_jet.reset(new JetHists     (ctx, "output_Jets"));
  output_h_eff.reset(new ZPrimeTotTPrimeHists(ctx, "output_eff"));
  output_h_topjet.reset(new TopJetHists     (ctx, "output_TopJets"));
  
  //after tagger fill hist on the different topjets on pT and zwtag
   topjet_first_pass_zwtag_h.reset(new TopJetHists(ctx, "topjet_first_pass_zwtag",4,"first_pass_topjets"));
   topjet_seconde_pass_zwtag_h.reset(new TopJetHists(ctx, "topjet_seconde_pass_zwtag",4,"seconde_pass_topjets"));
   topjet_third_pass_zwtag_h.reset(new TopJetHists(ctx, "topjet_third_pass_zwtag",4,"third_pass_topjets"));
   topjet_fourth_pass_zwtag_h.reset(new TopJetHists(ctx, "topjet_fourth_pass_zwtag",4,"fourth_pass_topjets"));
   topjet_first_nopass_zwtag_h.reset(new TopJetHists(ctx, "topjet_first_nopass_zwtag",4,"first_nopass_topjets"));
   topjet_seconde_nopass_zwtag_h.reset(new TopJetHists(ctx, "topjet_seconde_nopass_zwtag",4,"seconde_nopass_topjets"));
   topjet_third_nopass_zwtag_h.reset(new TopJetHists(ctx, "topjet_third_nopass_zwtag",4,"third_nopass_topjets"));
   topjet_fourth_nopass_zwtag_h.reset(new TopJetHists(ctx, "topjet_fourth_nopass_zwtag",4,"fourth_nopass_topjets"));
  // all topjets passed the tagger
   topjet_pass_zwtag_h.reset(new TopJetHists(ctx, "topjet_pass_zwtag",4,"zwtag_topjets"));
   topjet_nopass_zwtag_h.reset(new TopJetHists(ctx, "topjet_nopass_zwtag",4,"nozwtag_topjets"));

   b_error=false;
}

bool ZPrimeTotTPrimeMisstagModule::process(Event & event) {
  if(b_error) std::cout<<"In MisstagModule: Am Anfang"<<std::endl;
  for (auto & mod : htcalc) {
    mod->process(event);
  }

  if(!metfilters_selection->passes(event)) return false;
  for(auto & m : metfilters){
    m->process(event);
  }


  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* pileup SF */
  if(!event.isRealData){ pileup_SF->process(event);lumiweight->process(event);}

  //Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(filename  == "TTbarAll"){
    ttgenprod->process(event);
    if(!genmttbar_sel->passes(event)) return false;
  }

  uhh2::Event::TriggerIndex ti_HT;
  ti_HT=event.get_trigger_index("HLT_PFHT800_v*");
  bool HT_trigger = event.passes_trigger(ti_HT);
  if(!HT_trigger) return false;

  topjet_cleaner->process(event);
  mass_cleaner->process(event);

  //  std::cout<<"In MisstagModule: vor Input"<<std::endl;
  input_h_event ->fill(event);
  input_h_muo ->fill(event);
  input_h_ele   ->fill(event);
  input_h_jet   ->fill(event);
  input_h_eff ->fill(event);
  input_h_topjet   ->fill(event);
  
  /////////////////////////////////////////////////////////////////// Topjet2 Selection ///////////////////////////////////////////////////////////////////////////////////////////////////////
  if(b_error) std::cout<<"In MisstagModule: vor topjet2"<<std::endl;
  bool pass_topjet = topjet2_sel->passes(event);
  if(!pass_topjet)return false;
  muon_topjet2_h->fill(event);
  eff_topjet2_h ->fill(event);
  event_topjet2_h->fill(event);
  topjet_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);

  //Handle
  //before tagger set a handle on the different topjets: pt sorted first topjet,seconde topjet...
  //first_topjets
  //sort the topjets to be sure they are arrange with increasing PT
  sort_by_pt<TopJet>(*event.topjets);
  std::vector<TopJet>* topjet_first(new std::vector<TopJet> (*event.topjets));
  topjet_first->clear();
  topjet_first->reserve(event.topjets->size());
  topjet_first->push_back(event.topjets->at(0));
  event.set(h_first_topjets,*topjet_first);
  //  seconde_topjets
  std::vector<TopJet>* topjet_seconde(new std::vector<TopJet> (*event.topjets));
  topjet_seconde->clear();
  topjet_seconde->reserve(event.topjets->size());
  topjet_seconde->push_back(event.topjets->at(1));
  event.set(h_seconde_topjets,*topjet_seconde);
 
  // third_topjets
  std::vector<TopJet>* topjet_third(new std::vector<TopJet> (*event.topjets));
  topjet_third->clear();
  topjet_third->reserve(event.topjets->size());
  if(event.topjets->size()>=3){
    topjet_third->push_back(event.topjets->at(2));
  }  
  event.set(h_third_topjets,*topjet_third);
  topjet_third_topjet2_h->fill(event);
  
  
  //  fourth_topjets
  std::vector<TopJet>* topjet_fourth(new std::vector<TopJet> (*event.topjets));
  topjet_fourth->clear();
  topjet_fourth->reserve(event.topjets->size());
  if(event.topjets->size()>=4){
    topjet_fourth->push_back(event.topjets->at(3));
  }  
  event.set(h_fourth_topjets,*topjet_fourth);
  topjet_fourth_topjet2_h->fill(event);
  
  //before tagger fill hist on the different topjets: pt sorted first topjet,seconde topjet...
  topjet_first_topjet2_h->fill(event);
  topjet_seconde_topjet2_h->fill(event);

  /////////////////////////////////////////////////////////////////// tagger Selection ///////////////////////////////////////////////////////////////////////////////////////////////////////
  if(b_error) std::cout<<"In MisstagModule: vor tagger"<<std::endl;
  bool pass_zwtag = tagger_sel->passes(event);
  if(!pass_zwtag)return false;

  output_h_event ->fill(event);
  output_h_muo ->fill(event);
  output_h_ele   ->fill(event);
  output_h_jet   ->fill(event);
  output_h_eff ->fill(event);
  output_h_topjet   ->fill(event);

  if(b_error) std::cout <<"----------------------------------------  Anfang ---------------------------------------------------------------------------------------------------------------------"<<std::endl;
  // change topjet collection and apply tagger again to see which topjet (pt sorted is tagged)
  //first
  std::vector<TopJet>* AK8Jets(new std::vector<TopJet> (*event.topjets));
  std::vector<TopJet>* tagged(new std::vector<TopJet> (*event.topjets));
  tagged-> clear();
  tagged->reserve(event.topjets->size());
  std::vector<TopJet>* untagged(new std::vector<TopJet> (*event.topjets));
  untagged-> clear();
  untagged->reserve(event.topjets->size());
  event.topjets->clear();
  event.topjets->reserve(topjet_first->size());
  for(const auto & j : *topjet_first) event.topjets->push_back(j); 
  if(b_error) std::cout <<"Size of topjets should be one  " <<event.topjets->size()<<std::endl;
  pass_zwtag = tagger_sel->passes(event);
  if(b_error) std::cout<<"pass_zwtagt  "<<pass_zwtag<<std::endl;
  if(pass_zwtag){
    if(b_error) std::cout<<"in pass zwtag first topjet"<<std::endl;
    event.set(h_first_pass_topjets, *event.topjets);
    tagged->push_back(event.topjets->at(0));
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    event.set(h_first_nopass_topjets,*event.topjets);
  }else {
    if(b_error) std::cout<<"in else zwtag first topjet"<<std::endl;
    event.set(h_first_nopass_topjets, *event.topjets);
    untagged->push_back(event.topjets->at(0));
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    event.set(h_first_pass_topjets,*event.topjets);
  }
  event.topjets->clear();
  event.topjets->reserve(AK8Jets->size());
  for(const auto & j : *AK8Jets) event.topjets->push_back(j);

  //seconde
  event.topjets->clear();
  event.topjets->reserve(topjet_seconde->size());
  for(const auto & j : *topjet_seconde) event.topjets->push_back(j); 
  if(b_error) std::cout <<"Size of topjets should be one  " <<event.topjets->size()<<std::endl;
  pass_zwtag = tagger_sel->passes(event);
  if(b_error) std::cout<<"pass_zwtagt  "<<pass_zwtag<<std::endl;
  if(pass_zwtag){
    if(b_error) std::cout<<"in pass zwtag seconde topjet"<<std::endl;
    event.set(h_seconde_pass_topjets, *event.topjets);
    tagged->push_back(event.topjets->at(0));
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    event.set(h_seconde_nopass_topjets,*event.topjets);
  }else{
    if(b_error) std::cout<<"in else zwtag seconde topjet"<<std::endl;
    event.set(h_seconde_nopass_topjets, *event.topjets);
    untagged->push_back(event.topjets->at(0));
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    event.set(h_seconde_pass_topjets,*event.topjets);
  }
  event.topjets->clear();
  event.topjets->reserve(AK8Jets->size());
  for(const auto & j : *AK8Jets) event.topjets->push_back(j);

 
    //third
    event.topjets->clear();
    event.topjets->reserve(topjet_third->size());
    for(const auto & j : *topjet_third) event.topjets->push_back(j); 
    if(b_error) std::cout <<"Size of topjets should be one  " <<event.topjets->size()<<std::endl;
    pass_zwtag = tagger_sel->passes(event);
    if(b_error) std::cout<<"pass_zwtagt  "<<pass_zwtag<<std::endl;
    if(pass_zwtag){
      if(b_error) std::cout<<"in pass zwtag third topjet"<<std::endl;
      event.set(h_third_pass_topjets, *event.topjets);
      if(event.topjets->size())tagged->push_back(event.topjets->at(0));
      event.topjets->clear();
      event.topjets->reserve(AK8Jets->size());
      event.set(h_third_nopass_topjets,*event.topjets);
    } else{
      if(b_error) std::cout<<"in else zwtag third topjet"<<std::endl;
      event.set(h_third_nopass_topjets, *event.topjets);
      if(event.topjets->size()) untagged->push_back(event.topjets->at(0));
      event.topjets->clear();
      event.topjets->reserve(AK8Jets->size());
      event.set(h_third_pass_topjets,*event.topjets);
    }
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    for(const auto & j : *AK8Jets) event.topjets->push_back(j);
    topjet_third_pass_zwtag_h->fill(event);
    topjet_third_nopass_zwtag_h->fill(event);
  
    if(b_error)cout<<"Size of tagged"<<tagged->size() <<endl;
    if(b_error)cout<<"Size of untagged"<<untagged->size() <<endl;
    if(b_error)cout<<"Size of topjets"<<event.topjets->size() <<endl;

    //foruth
    event.topjets->clear();
    event.topjets->reserve(topjet_fourth->size());
    for(const auto & j : *topjet_fourth) event.topjets->push_back(j); 
    if(b_error) std::cout <<"Size of topjets should be one  " <<event.topjets->size()<<std::endl;
    pass_zwtag = tagger_sel->passes(event);
    if(b_error) std::cout<<"pass_zwtagt  "<<pass_zwtag<<std::endl;
    if(pass_zwtag){
      if(b_error) std::cout<<"in pass zwtag fourth topjet"<<std::endl;
      event.set(h_fourth_pass_topjets, *event.topjets);
      if(event.topjets->size()>=4) tagged->push_back(event.topjets->at(0));
      event.topjets->clear();
      event.topjets->reserve(AK8Jets->size());
      event.set(h_fourth_nopass_topjets,*event.topjets);
    }else{
      if(b_error) std::cout<<"in else zwtag fourth topjet"<<std::endl;
      event.set(h_fourth_nopass_topjets, *event.topjets);
      if(event.topjets->size()>=4) untagged->push_back(event.topjets->at(0));
      event.topjets->clear();
      event.topjets->reserve(AK8Jets->size());
      event.set(h_fourth_pass_topjets,*event.topjets);
    }
    event.topjets->clear();
    event.topjets->reserve(AK8Jets->size());
    for(const auto & j : *AK8Jets) event.topjets->push_back(j);
    topjet_fourth_pass_zwtag_h->fill(event);
    topjet_fourth_nopass_zwtag_h->fill(event);
  
if(b_error) std::cout <<"----------------------------------------  Ende ---------------------------------------------------------------------------------------------------------------------"<<std::endl;

   // all topjets passed the tagger
  event.set(h_zwtag_topjets, *tagged);
  event.set(h_nozwtag_topjets, *untagged);


  //after tagger fill hist on the different topjets on pT and zwtag
  topjet_first_pass_zwtag_h->fill(event);
  topjet_seconde_pass_zwtag_h->fill(event);
  topjet_first_nopass_zwtag_h->fill(event);
  topjet_seconde_nopass_zwtag_h->fill(event);
  // all topjets passed the tagger
  topjet_pass_zwtag_h->fill(event);
  topjet_nopass_zwtag_h->fill(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeMisstagModule)


