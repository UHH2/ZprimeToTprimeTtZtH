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

#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHists.h>
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h>


#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeSidebandReconstruction.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstruction.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHypothesisHists.h" 

// #include <UHH2/common/include/HypothesisHists.h>
// #include <UHH2/common/include/TTbarReconstruction.h>
// #include <UHH2/common/include/ReconstructionHypothesis.h>
// #include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

using namespace uhh2examples;
using namespace uhh2;

class ZPrimeTotTPrimeSidebandModule : public uhh2::AnalysisModule {

 public:
  explicit ZPrimeTotTPrimeSidebandModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners
  std::unique_ptr<MuonCleaner>     muo_cleaner;
  std::unique_ptr<ElectronCleaner> ele_cleaner;
  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCleaner>      jet_cleaner2;
  std::unique_ptr<JetCleaner>      jet_cleaner1;
  std::unique_ptr<JetCleaner>      topjet_IDcleaner;
  std::unique_ptr<TopJetCleaner>   topjet_cleaner;
  std::unique_ptr<TopJetCleaner>   toptag_cleaner;
  std::unique_ptr<TopJetCleaner>   higgstag_cleaner;
  std::unique_ptr<TopJetCleaner>   ZWtag_cleaner;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
  std::unique_ptr<JetCleaner>  ak4_cleaner;
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;
  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<AndSelection>  metfilters_selection;

  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;

  //Selections
  std::unique_ptr<uhh2::Selection> lumi_sel;
  std::unique_ptr<uhh2::AndSelection> lep1_sel; //  one lepton (comment out )
  std::unique_ptr<uhh2::Selection> TOPjet2_sel; // pT,1 > 250 GeV, pT,2 > 100 GeV
  std::unique_ptr<uhh2::Selection> TOPjet1_sel; // at least 2 jets
  std::unique_ptr<uhh2::Selection> chi2cut_sel; // chi2min <50
  std::unique_ptr<uhh2::Selection> jet2_sel;
  std::unique_ptr<uhh2::Selection> btag1_sel;
 std::unique_ptr<uhh2::Selection> btag0_sel;

  std::unique_ptr<uhh2::Selection> twodcut_sel;// pt 20 rel 0.4
  std::unique_ptr<uhh2::Selection> reliso_sel; //relIso < 0.35
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> ht_sel;
  std::unique_ptr<uhh2::Selection> muonpt_sel;

  //TOP TAGGER
  std::unique_ptr<Selection> toptag_sel; // at least one toptag
   uhh2::Event::Handle< std::vector<TopJet> > h_toptag;

  //HIGGS TAGGER
  std::unique_ptr<Selection> higgstag_sel;
   uhh2::Event::Handle< std::vector<TopJet> > h_higgstag;

  //Z/W TAGGER
  std::unique_ptr<Selection> ZWtag_sel;
   uhh2::Event::Handle< std::vector<TopJet> > h_ZWtag;

  // reconstruction ZPrime for Signal
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrimeprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_reco;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_selection_reco;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_chi;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_chi;

  // reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

  
  ////////////////////HISTS///////////////////////////////////////
  //Input Hists
  std::unique_ptr<Hists> input_event_h;
  std::unique_ptr<Hists> input_lep_h;
  std::unique_ptr<Hists> input_eff_h;
  std::unique_ptr<Hists> input_jet_h; 
  std::unique_ptr<Hists> input_topjet_h;
  
  //Hists lep1
  std::unique_ptr<Hists> event_lep1_h;
  std::unique_ptr<Hists> topjet_lep1_h;
  std::unique_ptr<Hists> jet_lep1_h;
  std::unique_ptr<Hists> muon_lep1_h;
  std::unique_ptr<Hists> eff_lep1_h;
 
  //Hist Topjet2
  std::unique_ptr<Hists> eff_topjet2_h;
  std::unique_ptr<Hists> jet_topjet2_h;
  std::unique_ptr<Hists> muon_topjet2_h;
  std::unique_ptr<Hists> event_topjet2_h;
  std::unique_ptr<Hists> topjet_topjet2_h;
 //Hist Btag1
  std::unique_ptr<Hists> eff_btag1_h;
  std::unique_ptr<Hists> jet_btag1_h;
  std::unique_ptr<Hists> muon_btag1_h;
  std::unique_ptr<Hists> event_btag1_h;
  std::unique_ptr<Hists> topjet_btag1_h;
 std::unique_ptr<Hists> chi2min_btag1_h;


 //Hist Btag0
  std::unique_ptr<Hists> eff_btag0_h;
  std::unique_ptr<Hists> jet_btag0_h;
  std::unique_ptr<Hists> muon_btag0_h;
  std::unique_ptr<Hists> event_btag0_h;
  std::unique_ptr<Hists> topjet_btag0_h;
 std::unique_ptr<Hists> chi2min_btag0_h;


 //  //reliso Hists
 //  std::unique_ptr<Hists> eff_reliso_h;
 //  std::unique_ptr<Hists> jet_reliso_h;
 //  std::unique_ptr<Hists> topjet_reliso_h;
 //  std::unique_ptr<Hists> muon_reliso_h;
 // std::unique_ptr<Hists> event_reliso_h;

  //2DCut Hists
  std::unique_ptr<Hists> eff_twodcut_h;
  std::unique_ptr<Hists> jet_twodcut_h;
  std::unique_ptr<Hists> topjet_twodcut_h;
  std::unique_ptr<Hists> muon_twodcut_h;
 std::unique_ptr<Hists> event_twodcut_h;
   
//Higgstag Hists
  std::unique_ptr<TopJetHists> input_higgstag_h;
  std::unique_ptr<TopJetHists> output_higgstag_h;

  std::unique_ptr<Hists> topjet_higgstag_h;
  std::unique_ptr<Hists> eff_higgstag_h;
  std::unique_ptr<Hists> jet_higgstag_h;
  std::unique_ptr<Hists> muon_higgstag_h;
  std::unique_ptr<Hists> event_higgstag_h;

  //ZWTagg Hists
  std::unique_ptr<TopJetHists> input_zwtag_h;
  std::unique_ptr<TopJetHists> output_zwtag_h;

  std::unique_ptr<Hists> topjet_zwtag_h;
  std::unique_ptr<Hists> eff_zwtag_h;
  std::unique_ptr<Hists> jet_zwtag_h;
  std::unique_ptr<Hists> muon_zwtag_h;
  std::unique_ptr<Hists> event_zwtag_h;

  //Higgstag || ZWtag Hists
  std::unique_ptr<Hists> output_combined_h;

  //Toptag Hists
  std::unique_ptr<TopJetHists> input_toptag_h;
  std::unique_ptr<TopJetHists> output_toptag_h;

  std::unique_ptr<Hists> topjet_toptag_h;
  std::unique_ptr<Hists> eff_toptag_h;
  std::unique_ptr<Hists> muon_toptag_h;
  std::unique_ptr<Hists> jet_toptag_h;
  std::unique_ptr<Hists> event_toptag_h;

  std::unique_ptr<TopJetHists> higgs_top_topjet_toptag_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_toptag_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_toptag_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_toptag_h;

  //Hist Tagger
  std::unique_ptr<Hists> topjet_tagger_h;
  std::unique_ptr<Hists> eff_tagger_h;
  std::unique_ptr<Hists> jet_tagger_h;
  std::unique_ptr<Hists> muon_tagger_h;
  std::unique_ptr<Hists> event_tagger_h;

  std::unique_ptr<TopJetHists> higgs_top_topjet_tagger_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_tagger_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_tagger_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_tagger_h;

  std::unique_ptr<JetHists> higgs_top_jet_tagger_h;
  std::unique_ptr<JetHists> zw_top_jet_tagger_h;
  std::unique_ptr<JetHists> higgs_notop_jet_tagger_h;
  std::unique_ptr<JetHists> zw_notop_jet_tagger_h;

  std::unique_ptr<MuonHists> higgs_top_muon_tagger_h;
  std::unique_ptr<MuonHists> zw_top_muon_tagger_h;
  std::unique_ptr<MuonHists> higgs_notop_muon_tagger_h;
  std::unique_ptr<MuonHists> zw_notop_muon_tagger_h;

  std::unique_ptr<EventHists> higgs_top_event_tagger_h;
  std::unique_ptr<EventHists> zw_top_event_tagger_h;
  std::unique_ptr<EventHists> higgs_notop_event_tagger_h;
  std::unique_ptr<EventHists> zw_notop_event_tagger_h;

// //Hist Clean
//   std::unique_ptr<Hists> chi2min_clean_h;
//   std::unique_ptr<Hists> eff_clean_h;
//   std::unique_ptr<Hists> jet_clean_h;
//   std::unique_ptr<Hists> muon_clean_h;
//   std::unique_ptr<Hists> topjet_clean_h;
//   std::unique_ptr<Hists> event_clean_h;

  // //Hist Reco
  // std::unique_ptr<Hists> chi2min_reco_h;
  // std::unique_ptr<Hists> eff_reco_h;
  // std::unique_ptr<Hists> jet_reco_h;
  // std::unique_ptr<Hists> muon_reco_h;
  // std::unique_ptr<Hists> topjet_reco_h;
  // std::unique_ptr<Hists> event_reco_h;

  // std::unique_ptr<TopJetHists> higgs_top_topjet_reco_h;
  // std::unique_ptr<TopJetHists> zw_top_topjet_reco_h;
  // std::unique_ptr<TopJetHists> higgs_notop_topjet_reco_h;
  // std::unique_ptr<TopJetHists> zw_notop_topjet_reco_h;
  // std::unique_ptr<Hists> higgs_top_chi2min_reco_h;
  // std::unique_ptr<Hists> zw_top_chi2min_reco_h;
  // std::unique_ptr<Hists> higgs_notop_chi2min_reco_h;
  // std::unique_ptr<Hists> zw_notop_chi2min_reco_h;
 std::unique_ptr<Hists> chi2min_reco_h;
std::unique_ptr<Hists> top_chi2min_reco_h;
std::unique_ptr<Hists> notop_chi2min_reco_h;

  //Hist Chi2cut
  std::unique_ptr<Hists> chi2min_chi2cut_h;
  std::unique_ptr<Hists> eff_chi2cut_h;
  std::unique_ptr<Hists> jet_chi2cut_h;
  std::unique_ptr<Hists> muon_chi2cut_h;
  std::unique_ptr<Hists> topjet_chi2cut_h;
  std::unique_ptr<Hists> event_chi2cut_h;

 std::unique_ptr<TopJetHists> higgs_top_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_chi2cut_h;

  std::unique_ptr<Hists> higgs_top_chi2min_chi2cut_h;
  std::unique_ptr<Hists> zw_top_chi2min_chi2cut_h;
  std::unique_ptr<Hists> higgs_notop_chi2min_chi2cut_h;
  std::unique_ptr<Hists> zw_notop_chi2min_chi2cut_h;

 std::unique_ptr<JetHists> higgs_top_jet_chi2cut_h;
  std::unique_ptr<JetHists> zw_top_jet_chi2cut_h;
  std::unique_ptr<JetHists> higgs_notop_jet_chi2cut_h;
  std::unique_ptr<JetHists> zw_notop_jet_chi2cut_h;
  
  std::unique_ptr<MuonHists> higgs_top_muon_chi2cut_h;
  std::unique_ptr<MuonHists> zw_top_muon_chi2cut_h;
  std::unique_ptr<MuonHists> higgs_notop_muon_chi2cut_h;
  std::unique_ptr<MuonHists> zw_notop_muon_chi2cut_h;

  std::unique_ptr<EventHists> higgs_top_event_chi2cut_h;
  std::unique_ptr<EventHists> zw_top_event_chi2cut_h;
  std::unique_ptr<EventHists> higgs_notop_event_chi2cut_h;
  std::unique_ptr<EventHists> zw_notop_event_chi2cut_h;

  //counting signal region events
  std::unique_ptr<Hists> eff_selection_h;
  std::unique_ptr<Hists> chi2min_selection_h;
  std::unique_ptr<Hists> jet_selection_h;
  std::unique_ptr<Hists> muon_selection_h;
  std::unique_ptr<Hists> topjet_selection_h;
  std::unique_ptr<Hists> event_selection_h;

  // std::unique_ptr<TopJetHists> higgs_top_topjet_chi2cut_h;
 //  std::unique_ptr<TopJetHists> zw_top_topjet_chi2cut_h;
 //  std::unique_ptr<TopJetHists> higgs_notop_topjet_chi2cut_h;
 //  std::unique_ptr<TopJetHists> zw_notop_topjet_chi2cut_h;

 //  std::unique_ptr<Hists> higgs_top_chi2min_chi2cut_h;
 //  std::unique_ptr<Hists> zw_top_chi2min_chi2cut_h;
 //  std::unique_ptr<Hists> higgs_notop_chi2min_chi2cut_h;
 //  std::unique_ptr<Hists> zw_notop_chi2min_chi2cut_h;

 // std::unique_ptr<JetHists> higgs_top_jet_chi2cut_h;
 //  std::unique_ptr<JetHists> zw_top_jet_chi2cut_h;
 //  std::unique_ptr<JetHists> higgs_notop_jet_chi2cut_h;
 //  std::unique_ptr<JetHists> zw_notop_jet_chi2cut_h;
  
 //  std::unique_ptr<MuonHists> higgs_top_muon_chi2cut_h;
 //  std::unique_ptr<MuonHists> zw_top_muon_chi2cut_h;
 //  std::unique_ptr<MuonHists> higgs_notop_muon_chi2cut_h;
 //  std::unique_ptr<MuonHists> zw_notop_muon_chi2cut_h;

 //  std::unique_ptr<EventHists> higgs_top_event_chi2cut_h;
 //  std::unique_ptr<EventHists> zw_top_event_chi2cut_h;
 //  std::unique_ptr<EventHists> higgs_notop_event_chi2cut_h;
 //  std::unique_ptr<EventHists> zw_notop_event_chi2cut_h;

  // //met
  // std::unique_ptr<Hists> topjet_met_h;
  // std::unique_ptr<Hists>   eff_met_h;
  // std::unique_ptr<Hists>   jet_met_h;
  // std::unique_ptr<Hists>   muon_met_h;
  // std::unique_ptr<Hists>   chi2min_met_h;
  // std::unique_ptr<Hists> event_met_h;

  // std::unique_ptr<TopJetHists> higgs_top_topjet_met_h;
  // std::unique_ptr<TopJetHists> zw_top_topjet_met_h;
  // std::unique_ptr<TopJetHists> higgs_notop_topjet_met_h;
  // std::unique_ptr<TopJetHists> zw_notop_topjet_met_h;

  // std::unique_ptr<Hists> higgs_top_chi2min_met_h;
  // std::unique_ptr<Hists> zw_top_chi2min_met_h;
  // std::unique_ptr<Hists> higgs_notop_chi2min_met_h;
  // std::unique_ptr<Hists> zw_notop_chi2min_met_h;

  // std::unique_ptr<JetHists> higgs_top_jet_met_h;
  // std::unique_ptr<JetHists> zw_top_jet_met_h;
  // std::unique_ptr<JetHists> higgs_notop_jet_met_h;
  // std::unique_ptr<JetHists> zw_notop_jet_met_h;
  
  // std::unique_ptr<MuonHists> higgs_top_muon_met_h;
  // std::unique_ptr<MuonHists> zw_top_muon_met_h;
  // std::unique_ptr<MuonHists> higgs_notop_muon_met_h;
  // std::unique_ptr<MuonHists> zw_notop_muon_met_h;

  // std::unique_ptr<EventHists> higgs_top_event_met_h;
  // std::unique_ptr<EventHists> zw_top_event_met_h;
  // std::unique_ptr<EventHists> higgs_notop_event_met_h;
  // std::unique_ptr<EventHists> zw_notop_event_met_h;


 // //ht
 //  std::unique_ptr<Hists> topjet_ht_h;
 //  std::unique_ptr<Hists>   eff_ht_h;
 //  std::unique_ptr<Hists>   jet_ht_h;
 //  std::unique_ptr<Hists>   muon_ht_h;
 //  std::unique_ptr<Hists>   chi2min_ht_h;
 //  std::unique_ptr<Hists> event_ht_h;

  // std::unique_ptr<TopJetHists> higgs_top_topjet_ht_h;
  // std::unique_ptr<TopJetHists> zw_top_topjet_ht_h;
  // std::unique_ptr<TopJetHists> higgs_notop_topjet_ht_h;
  // std::unique_ptr<TopJetHists> zw_notop_topjet_ht_h;

  // std::unique_ptr<Hists> higgs_top_chi2min_ht_h;
  // std::unique_ptr<Hists> zw_top_chi2min_ht_h;
  // std::unique_ptr<Hists> higgs_notop_chi2min_ht_h;
  // std::unique_ptr<Hists> zw_notop_chi2min_ht_h;

  // std::unique_ptr<JetHists> higgs_top_jet_ht_h;
  // std::unique_ptr<JetHists> zw_top_jet_ht_h;
  // std::unique_ptr<JetHists> higgs_notop_jet_ht_h;
  // std::unique_ptr<JetHists> zw_notop_jet_ht_h;
  
  // std::unique_ptr<MuonHists> higgs_top_muon_ht_h;
  // std::unique_ptr<MuonHists> zw_top_muon_ht_h;
  // std::unique_ptr<MuonHists> higgs_notop_muon_ht_h;
  // std::unique_ptr<MuonHists> zw_notop_muon_ht_h;

  // std::unique_ptr<EventHists> higgs_top_event_ht_h;
  // std::unique_ptr<EventHists> zw_top_event_ht_h;
  // std::unique_ptr<EventHists> higgs_notop_event_ht_h;
  // std::unique_ptr<EventHists> zw_notop_event_ht_h;

  // std::unique_ptr<Hists> topjet_top_masswindow_h;
  //std::unique_ptr<Hists> topjet_higgs_masswindow_h;
  //std::unique_ptr<Hists> topjet_zw_masswindow_h;
 std::unique_ptr<Hists> lumi_h;
  //general
  std::string filename;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ZprimeTotTPrime_hyps;
  uhh2::Event::Handle<double> h_ht;
  uhh2::Event::Handle< std::vector<TopJet> > h_AK8;
  uhh2::Event::Handle< std::vector<Jet> > h_AK4;

  bool berror=false;
};


ZPrimeTotTPrimeSidebandModule::ZPrimeTotTPrimeSidebandModule(uhh2::Context& ctx){


  //choose channel
  const std::string& channel = ctx.get("channel", "");
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else throw std::runtime_error("ZprimeSelectionModule -- undefined argument for 'channel' key in xml file (must be 'muon' or 'elec'): "+channel);

  const bool isMC = (ctx.get("dataset_type") == "MC");
  //// COMMON MODULES
  if(isMC){ pileup_SF.reset(new MCPileupReweight(ctx)); lumiweight.reset(new MCLumiWeight(ctx));}
  else     lumi_sel.reset(new LumiSelection(ctx));

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

  std::vector<std::string> JEC_AK4, JEC_AK8;
  if(isMC){

    JEC_AK4 = JERFiles::Summer15_25ns_L123_AK4PFchs_MC;
    JEC_AK8 = JERFiles::Summer15_25ns_L123_AK8PFchs_MC;
  }
  else {

    JEC_AK4 = JERFiles::Summer15_25ns_L123_AK4PFchs_DATA;
    JEC_AK8 = JERFiles::Summer15_25ns_L123_AK8PFchs_DATA;
  }

  //// OBJ CLEANING
  muo_cleaner.reset(new MuonCleaner    (AndId<Muon>    (PtEtaCut  (50., 2.1), MuonIDMedium())));
  ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaSCCut(50., 2.4), ElectronID_MVAnotrig_Spring15_25ns_loose)));

  const JetId jetID(JetPFID(JetPFID::WP_LOOSE));
  jet_IDcleaner.reset(new JetCleaner(ctx,jetID));
  jet_cleaner2.reset(new JetCleaner(ctx,15., 2.4));
  jet_cleaner1.reset(new JetCleaner(ctx,30., 2.4));
  jetlepton_cleaner.reset(new JetLeptonCleaner(ctx, JEC_AK4));
  jetlepton_cleaner->set_drmax(.4);
  ak4_cleaner.reset(new JetCleaner(ctx,JetId(ZPrimeTotTPrimeAK4cleaner(1.2))));
  // ak4_cleaner.reset(new TopJetCleaner(ctx,TopJetId(ZPrimeTotTPrimeAK4cleaner(1.2))));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));
  topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));

  // SELECTIONS
  //jet1_sel.reset(new NJetSelection(1,-1, JetId(PtEtaCut(0.,2.5))));
  
  // TOPJET SELECTIONS
  TOPjet1_sel.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut( 250., 2.4))));
  //  TOPjet2_sel.reset(new ZPrimeTotTPrimeNTopJetCut(1,-1,400,100, 2.4));
  jet2_sel.reset(new NJetSelection(2, -1));

  //BTAG
  btag1_sel.reset(new NBTagSelection(1, -1));
  btag0_sel.reset(new NBTagSelection(0, 0));

  //2D Cut
   twodcut_sel.reset(new TwoDCut(.4, 25.));
  // reliso_sel.reset(new ZPrimeTotTPrimeRelIsoCut(0.15));
   //  muonpt_sel.reset(new ZPrimeTotTPrimeMuonPT(170));
  //TOP TAGGER
  const TopJetId topjetID = AndId<TopJet>(Type2TopTag(150,240,Type2TopTag::MassType::groomed), Tau32(0.6));
  toptag_sel.reset(new NTopJetSelection(1, -1, topjetID));
  toptag_cleaner.reset(new TopJetCleaner(ctx,topjetID));


  //Higgs TAGGER
   const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150), Tau21(1));
  higgstag_sel.reset(new NTopJetSelection(1, -1, higgsjetID));
  higgstag_cleaner.reset(new TopJetCleaner(ctx,higgsjetID));

  //W/Z TAGGER
  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
  ZWtag_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
  ZWtag_cleaner.reset(new TopJetCleaner(ctx,ZWjetID));


  //LEPTONSELECTION
  lep1_sel.reset(new uhh2::AndSelection(ctx));

  if(channel_ == muon){

    lep1_sel->add<NMuonSelection>    ("muoN == 1", 1, 1);
    lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);
  }
  else if(channel_ == elec){

    lep1_sel->add<NMuonSelection>    ("muoN == 0", 0, 0);
    lep1_sel->add<NElectronSelection>("eleN == 1", 1, 1);
  }


  /* KINEMATICAL RECO */
  const std::string ZprimeTotTPrime_gen_label ("zprimegen");
  const std::string ZprimeTotTPrime_hyps_label("ZPrimeTotTPrimeSidebandReconstruction");
  const std::string ZprimeTotTPrime_selection_hyps_label("ZPrimeTotTPrimeReconstruction");
  const std::string ttbar_hyps_label("TTbarReconstruction");
  const std::string ZprimeTotTPrime_chi2_label("Chi2");
  const std::string ttbar_chi2_label("Chi2");
  const std::string ttbar_gen_label ("ttbargen");

  reco_primlep.reset(new PrimaryLepton(ctx));

  ZprimeTotTPrime_reco.reset(new ZPrimeTotTPrimeSidebandReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrime_selection_reco.reset(new ZPrimeTotTPrimeReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_selection_hyps_label));
  ZprimeTotTPrime_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrimeprod.reset(new ZPrimeGenProducer(ctx, ZprimeTotTPrime_gen_label, false));
  h_ZprimeTotTPrime_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ZprimeTotTPrime_hyps_label);
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));


  chi2cut_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));


  // //MetCut
  // met_sel.reset(new ZPrimeTotTPrimeMETCut(250,-1));
  // ht_sel.reset(new HtSelection(400,-1));

  // Hists Input
  input_event_h.reset(new EventHists(ctx, "input"));
  input_lep_h.reset(new MuonHists(ctx, "input_Lep"));
  input_jet_h.reset(new JetHists (ctx, "input_Jet"));
  input_topjet_h.reset(new TopJetHists (ctx, "input_TopJet"));
  input_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
lumi_h.reset(new LuminosityHists(ctx,"lumi"));
  // Hists lep1
  topjet_lep1_h.reset(new TopJetHists(ctx, "topjet_lep1"));
  event_lep1_h.reset(new EventHists(ctx, "event_lep1"));
  muon_lep1_h.reset(new MuonHists(ctx, "muon_lep1"));
  jet_lep1_h.reset(new JetHists(ctx, "jet_lep1"));
  eff_lep1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_lep1"));

  //Hists topjet2
  topjet_topjet2_h.reset(new TopJetHists(ctx, "topjet_topjet2"));
  eff_topjet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_topjet2"));
  jet_topjet2_h.reset(new JetHists(ctx, "jet_topjet2"));
  muon_topjet2_h.reset(new MuonHists(ctx, "muon_topjet2"));
  event_topjet2_h.reset(new EventHists(ctx, "event_topjet2"));

 //Hists btag0
  topjet_btag0_h.reset(new TopJetHists(ctx, "topjet_btag0"));
  eff_btag0_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag0"));
  jet_btag0_h.reset(new JetHists(ctx, "jet_btag0"));
  muon_btag0_h.reset(new MuonHists(ctx, "muon_btag0"));
  event_btag0_h.reset(new EventHists(ctx, "event_btag0"));
chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

//Hists btag1
  topjet_btag1_h.reset(new TopJetHists(ctx, "topjet_btag1"));
  eff_btag1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag1"));
  jet_btag1_h.reset(new JetHists(ctx, "jet_btag1"));
  muon_btag1_h.reset(new MuonHists(ctx, "muon_btag1"));
  event_btag1_h.reset(new EventHists(ctx, "event_btag1"));
chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));
  // //Hist reliso
  // topjet_reliso_h.reset(new TopJetHists(ctx, "topjet_reliso"));
  // eff_reliso_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reliso"));
  // jet_reliso_h.reset(new JetHists(ctx, "jet_reliso"));
  // muon_reliso_h.reset(new MuonHists(ctx, "muon_reliso"));
  // event_reliso_h.reset(new EventHists(ctx, "event_reliso"));


  //Hist twodcut
  topjet_twodcut_h.reset(new TopJetHists(ctx, "topjet_twodcut"));
  eff_twodcut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_twodcut"));
  jet_twodcut_h.reset(new JetHists(ctx, "jet_twodcut"));
  muon_twodcut_h.reset(new MuonHists(ctx, "muon_twodcut"));
  event_twodcut_h.reset(new EventHists(ctx, "event_twodcut"));

  //HIGGSTAG HISTS
  input_higgstag_h.reset(new TopJetHists(ctx, "input_higgstag"));
  output_higgstag_h.reset(new TopJetHists(ctx, "output_higgstag"));
  output_higgstag_h->set_TopJetId(higgsjetID);

  eff_higgstag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_higgstag"));
  jet_higgstag_h.reset(new JetHists(ctx, "jet_higgstag"));
  muon_higgstag_h.reset(new MuonHists(ctx, "muon_higgstag"));
  event_higgstag_h.reset(new EventHists(ctx, "event_higgstag"));
  topjet_higgstag_h.reset(new TopJetHists(ctx, "topjet_higgstag"));
 
  //Z/W TAG HISTS
  input_zwtag_h.reset(new TopJetHists(ctx, "input_zwtag"));
  output_zwtag_h.reset(new TopJetHists(ctx, "output_zwtag"));
  output_zwtag_h->set_TopJetId(ZWjetID);

  eff_zwtag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_zwtag"));
  jet_zwtag_h.reset(new JetHists(ctx, "jet_zwtag"));
  muon_zwtag_h.reset(new MuonHists(ctx, "muon_zwtag"));
  event_zwtag_h.reset(new EventHists(ctx, "event_zwtag"));
  topjet_zwtag_h.reset(new TopJetHists(ctx, "topjet_zwtag"));

  // Hists Toptag
  input_toptag_h.reset(new TopJetHists(ctx, "input_toptag"));
  output_toptag_h.reset(new TopJetHists(ctx, "output_toptag"));
  output_toptag_h->set_TopJetId(topjetID);

  eff_toptag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_toptag"));
  jet_toptag_h.reset(new JetHists(ctx, "jet_toptag"));
  muon_toptag_h.reset(new MuonHists(ctx, "muon_toptag"));
  event_toptag_h.reset(new EventHists(ctx, "event_toptag"));
  topjet_toptag_h.reset(new TopJetHists(ctx, "topjet_toptag"));

  higgs_top_topjet_toptag_h.reset(new TopJetHists(ctx, "higgs_top_topjet_toptag"));
  higgs_notop_topjet_toptag_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_toptag"));
  zw_top_topjet_toptag_h.reset(new TopJetHists(ctx, "zw_top_topjet_toptag"));
  zw_notop_topjet_toptag_h.reset(new TopJetHists(ctx, "zw_notop_topjet_toptag"));

  //Hists tagger
  eff_tagger_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_tagger"));
  jet_tagger_h.reset(new JetHists(ctx, "jet_tagger"));
  muon_tagger_h.reset(new MuonHists(ctx, "muon_tagger"));
  event_tagger_h.reset(new EventHists(ctx, "event_tagger"));
  topjet_tagger_h.reset(new TopJetHists(ctx, "topjet_tagger"));

  higgs_top_topjet_tagger_h.reset(new TopJetHists(ctx, "higgs_top_topjet_tagger"));
  higgs_notop_topjet_tagger_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_tagger"));
  zw_top_topjet_tagger_h.reset(new TopJetHists(ctx, "zw_top_topjet_tagger"));
  zw_notop_topjet_tagger_h.reset(new TopJetHists(ctx, "zw_notop_topjet_tagger"));

  higgs_top_jet_tagger_h.reset(new JetHists(ctx, "higgs_top_jet_tagger"));
  higgs_notop_jet_tagger_h.reset(new JetHists(ctx, "higgs_notop_jet_tagger"));
  zw_top_jet_tagger_h.reset(new JetHists(ctx, "zw_top_jet_tagger"));
  zw_notop_jet_tagger_h.reset(new JetHists(ctx, "zw_notop_jet_tagger"));

  higgs_top_muon_tagger_h.reset(new MuonHists(ctx, "higgs_top_muon_tagger"));
  higgs_notop_muon_tagger_h.reset(new MuonHists(ctx, "higgs_notop_muon_tagger"));
  zw_top_muon_tagger_h.reset(new MuonHists(ctx, "zw_top_muon_tagger"));
  zw_notop_muon_tagger_h.reset(new MuonHists(ctx, "zw_notop_muon_tagger"));

  higgs_top_event_tagger_h.reset(new EventHists(ctx, "higgs_top_event_tagger"));
  higgs_notop_event_tagger_h.reset(new EventHists(ctx, "higgs_notop_event_tagger"));
  zw_top_event_tagger_h.reset(new EventHists(ctx, "zw_top_event_tagger"));
  zw_notop_event_tagger_h.reset(new EventHists(ctx, "zw_notop_event_tagger"));
 
 // //Hists Clean
 //  topjet_clean_h.reset(new TopJetHists(ctx, "topjet_clean"));
 //  eff_clean_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_clean"));
 //  jet_clean_h.reset(new JetHists(ctx, "jet_clean"));
 //  muon_clean_h.reset(new MuonHists(ctx, "muon_clean"));
 //  event_clean_h.reset(new EventHists(ctx, "event_clean"));
 //  chi2min_clean_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_clean",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

 //  //Hists Reco
 //  topjet_reco_h.reset(new TopJetHists(ctx, "topjet_reco"));
 //  eff_reco_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reco"));
 //  jet_reco_h.reset(new JetHists(ctx, "jet_reco"));
 //  muon_reco_h.reset(new MuonHists(ctx, "muon_reco"));
 //  event_reco_h.reset(new EventHists(ctx, "event_reco"));
 //  chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

  // higgs_top_topjet_reco_h.reset(new TopJetHists(ctx, "higgs_top_topjet_reco"));
  // higgs_notop_topjet_reco_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_reco"));
  // zw_top_topjet_reco_h.reset(new TopJetHists(ctx, "zw_top_topjet_reco"));
  // zw_notop_topjet_reco_h.reset(new TopJetHists(ctx, "zw_notop_topjet_reco"));

  // higgs_top_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // higgs_notop_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_top_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_notop_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));
 top_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "top_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));
 notop_chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "notop_chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

  //Hists Chi2cut
  topjet_chi2cut_h.reset(new TopJetHists(ctx, "topjet_chi2cut"));
  eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut"));
  jet_chi2cut_h.reset(new JetHists(ctx, "jet_chi2cut"));
  muon_chi2cut_h.reset(new MuonHists(ctx, "muon_chi2cut"));
  event_chi2cut_h.reset(new EventHists(ctx, "event_chi2cut"));
  chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));


  higgs_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_top_topjet_chi2cut"));
  higgs_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_chi2cut"));
  zw_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_top_topjet_chi2cut"));
  zw_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_notop_topjet_chi2cut"));

  higgs_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  higgs_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  higgs_top_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_top_jet_chi2cut"));
  higgs_notop_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_notop_jet_chi2cut"));
  zw_top_jet_chi2cut_h.reset(new JetHists(ctx, "zw_top_jet_chi2cut"));
  zw_notop_jet_chi2cut_h.reset(new JetHists(ctx, "zw_notop_jet_chi2cut"));

  higgs_top_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_top_muon_chi2cut"));
  higgs_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_notop_muon_chi2cut"));
  zw_top_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_top_muon_chi2cut"));
  zw_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_notop_muon_chi2cut"));

  higgs_top_event_chi2cut_h.reset(new EventHists(ctx, "higgs_top_event_chi2cut"));
  higgs_notop_event_chi2cut_h.reset(new EventHists(ctx, "higgs_notop_event_chi2cut"));
  zw_top_event_chi2cut_h.reset(new EventHists(ctx, "zw_top_event_chi2cut"));
  zw_notop_event_chi2cut_h.reset(new EventHists(ctx, "zw_notop_event_chi2cut"));


  eff_selection_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_selection"));
  topjet_selection_h.reset(new TopJetHists(ctx, "topjet_selection"));
  jet_selection_h.reset(new JetHists(ctx, "jet_selection"));
  muon_selection_h.reset(new MuonHists(ctx, "muon_selection"));
  event_selection_h.reset(new EventHists(ctx, "event_selection"));
  chi2min_selection_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_selection",ZprimeTotTPrime_selection_hyps_label,ZprimeTotTPrime_chi2_label ));

  // higgs_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_top_topjet_chi2cut"));
  // higgs_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_chi2cut"));
  // zw_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_top_topjet_chi2cut"));
  // zw_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_notop_topjet_chi2cut"));

  // higgs_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // higgs_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  // higgs_top_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_top_jet_chi2cut"));
  // higgs_notop_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_notop_jet_chi2cut"));
  // zw_top_jet_chi2cut_h.reset(new JetHists(ctx, "zw_top_jet_chi2cut"));
  // zw_notop_jet_chi2cut_h.reset(new JetHists(ctx, "zw_notop_jet_chi2cut"));

  // higgs_top_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_top_muon_chi2cut"));
  // higgs_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_notop_muon_chi2cut"));
  // zw_top_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_top_muon_chi2cut"));
  // zw_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_notop_muon_chi2cut"));

  // higgs_top_event_chi2cut_h.reset(new EventHists(ctx, "higgs_top_event_chi2cut"));
  // higgs_notop_event_chi2cut_h.reset(new EventHists(ctx, "higgs_notop_event_chi2cut"));
  // zw_top_event_chi2cut_h.reset(new EventHists(ctx, "zw_top_event_chi2cut"));
  // zw_notop_event_chi2cut_h.reset(new EventHists(ctx, "zw_notop_event_chi2cut"));

  // //met
  // topjet_met_h.reset(new TopJetHists(ctx, "topjet_met"));
  // eff_met_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_met"));
  // jet_met_h.reset(new JetHists(ctx, "jet_met"));
  // muon_met_h.reset(new MuonHists(ctx, "muon_met"));
  // event_met_h.reset(new EventHists(ctx, "event_met"));
  // chi2min_met_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_met",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

  // higgs_top_topjet_met_h.reset(new TopJetHists(ctx, "higgs_top_topjet_met"));
  // higgs_notop_topjet_met_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_met"));
  // zw_top_topjet_met_h.reset(new TopJetHists(ctx, "zw_top_topjet_met"));
  // zw_notop_topjet_met_h.reset(new TopJetHists(ctx, "zw_notop_topjet_met"));

  // higgs_top_chi2min_met_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_met",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // higgs_notop_chi2min_met_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_met",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_top_chi2min_met_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_met",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  // zw_notop_chi2min_met_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_met",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  // higgs_top_jet_met_h.reset(new JetHists(ctx, "higgs_top_jet_met"));
  // higgs_notop_jet_met_h.reset(new JetHists(ctx, "higgs_notop_jet_met"));
  // zw_top_jet_met_h.reset(new JetHists(ctx, "zw_top_jet_met"));
  // zw_notop_jet_met_h.reset(new JetHists(ctx, "zw_notop_jet_met"));

  // higgs_top_muon_met_h.reset(new MuonHists(ctx, "higgs_top_muon_met"));
  // higgs_notop_muon_met_h.reset(new MuonHists(ctx, "higgs_notop_muon_met"));
  // zw_top_muon_met_h.reset(new MuonHists(ctx, "zw_top_muon_met"));
  // zw_notop_muon_met_h.reset(new MuonHists(ctx, "zw_notop_muon_met"));

  // higgs_top_event_met_h.reset(new EventHists(ctx, "higgs_top_event_met"));
  // higgs_notop_event_met_h.reset(new EventHists(ctx, "higgs_notop_event_met"));
  // zw_top_event_met_h.reset(new EventHists(ctx, "zw_top_event_met"));
  // zw_notop_event_met_h.reset(new EventHists(ctx, "zw_notop_event_met"));

 // //ht
 //  topjet_ht_h.reset(new TopJetHists(ctx, "topjet_ht"));
 //  eff_ht_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_ht"));
 //  jet_ht_h.reset(new JetHists(ctx, "jet_ht"));
 //  muon_ht_h.reset(new MuonHists(ctx, "muon_ht"));
 //  event_ht_h.reset(new EventHists(ctx, "event_ht"));
 //  chi2min_ht_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_ht",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

 //  higgs_top_topjet_ht_h.reset(new TopJetHists(ctx, "higgs_top_topjet_ht"));
 //  higgs_notop_topjet_ht_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_ht"));
 //  zw_top_topjet_ht_h.reset(new TopJetHists(ctx, "zw_top_topjet_ht"));
 //  zw_notop_topjet_ht_h.reset(new TopJetHists(ctx, "zw_notop_topjet_ht"));

 //  higgs_top_chi2min_ht_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_ht",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
 //  higgs_notop_chi2min_ht_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_ht",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
 //  zw_top_chi2min_ht_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_ht",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
 //  zw_notop_chi2min_ht_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_ht",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

 //  higgs_top_jet_ht_h.reset(new JetHists(ctx, "higgs_top_jet_ht"));
 //  higgs_notop_jet_ht_h.reset(new JetHists(ctx, "higgs_notop_jet_ht"));
 //  zw_top_jet_ht_h.reset(new JetHists(ctx, "zw_top_jet_ht"));
 //  zw_notop_jet_ht_h.reset(new JetHists(ctx, "zw_notop_jet_ht"));

 //  higgs_top_muon_ht_h.reset(new MuonHists(ctx, "higgs_top_muon_ht"));
 //  higgs_notop_muon_ht_h.reset(new MuonHists(ctx, "higgs_notop_muon_ht"));
 //  zw_top_muon_ht_h.reset(new MuonHists(ctx, "zw_top_muon_ht"));
 //  zw_notop_muon_ht_h.reset(new MuonHists(ctx, "zw_notop_muon_ht"));

 //  higgs_top_event_ht_h.reset(new EventHists(ctx, "higgs_top_event_ht"));
 //  higgs_notop_event_ht_h.reset(new EventHists(ctx, "higgs_notop_event_ht"));
 //  zw_top_event_ht_h.reset(new EventHists(ctx, "zw_top_event_ht"));
 //  zw_notop_event_ht_h.reset(new EventHists(ctx, "zw_notop_event_ht"));


 //  //Higgstag || ZWTag
 //  output_combined_h.reset(new ZPrimeTotTPrimeHists(ctx, "output_combined"));
 // // topjet_higgs_masswindow_h.reset(new TopJetHists(ctx, "masswindow_higgs"));
 // // topjet_zw_masswindow_h.reset(new TopJetHists(ctx, "masswindow_zw"));
 // // topjet_top_masswindow_h.reset(new TopJetHists(ctx, "masswindow_top"));


  //general
  filename =  ctx.get("dataset_version");
  h_toptag = ctx.declare_event_output< std::vector<TopJet> > ("TopTag");
  h_higgstag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag");
  h_ZWtag = ctx.declare_event_output< std::vector<TopJet> > ("ZWTag");
  h_ht = ctx.get_handle<double>("HT");
  h_AK8 = ctx.get_handle<std::vector<TopJet>>("AK8");
  h_AK4 = ctx.get_handle<std::vector<Jet>>("AK4");

  //Trigger
 const std::string& trigger = ctx.get("trigger", "NULL");
if(trigger != "NULL") trigger_sel = make_unique<TriggerSelection>(trigger);
    else                  trigger_sel = make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*");
  
}

bool ZPrimeTotTPrimeSidebandModule::process(uhh2::Event& event){
  
  if(berror)  std::cout<<"/////////////////////////////////////SelectionModule L:218 Am Anfang///////////////////////////////////////////////////////////////////////////////"<<std::endl;

  for (auto & mod : htcalc) {
      mod->process(event);
    }

  if(!metfilters_selection->passes(event)) return false;
  for(auto & m : metfilters){
        m->process(event);
  }
  //Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(filename  == "TTbarAll"){
    ttgenprod->process(event);
    if(!genmttbar_sel->passes(event)) return false;
  }

  //ZPrime Genrator Level
  if(filename.find("MC_ZPrime")!=std::string::npos){
    ZprimeTotTPrimeprod->process(event); 
  }

 if(berror)  std::cout<<"SelectionModule L:232 vor Input Histogrammen"<<std::endl;
  
  /////////////////////////////////////////////////////////// Common Modules   ///////////////////////////////////////////////////////////////////////////////


  //common Modules
  /* luminosity sections from CMS golden-JSON file */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* pileup SF */
  if(!event.isRealData){ pileup_SF->process(event);lumiweight->process(event);}
  ////

  // OBJ CLEANING
  muo_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);

  ele_cleaner->process(event);
  sort_by_pt<Electron>(*event.electrons);

  jet_IDcleaner->process(event);
  jetlepton_cleaner->process(event);
  jet_cleaner2->process(event); 
  
  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

///////Trigger///////
  const bool pass_trigger = trigger_sel->passes(event);
  if(!pass_trigger) return false;

/////////////////////////////////////////////////////////// Input Histogramme ///////////////////////////////////////////////////////////////////////////////
 
  input_eff_h ->fill(event);
  input_event_h->fill(event);
  input_topjet_h->fill(event);
  input_jet_h->fill(event);
  input_lep_h->fill(event);



 if(berror) std::cout<<"SelectionModule L:268 vor LeptonSelection"<<std::endl;
  /////////////////////////////////////////////////////////// LEPTON selection ///////////////////////////////////////////////////////////////////////////////
  const bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  muon_lep1_h->fill(event);
  eff_lep1_h ->fill(event);
  event_lep1_h->fill(event);
  topjet_lep1_h->fill(event);
  jet_lep1_h->fill(event);

  


   topjetlepton_cleaner->process(event);

 if(berror)  std::cout<<"SelectionModule L:268 vor TopJetSelection"<<std::endl;
  ////////////////////////////////////////////////////////// TopJET selection//////////////////////////////////////////////////////////////////////////////////

  const bool pass_Topjet1 = TOPjet1_sel->passes(event);
  if(!pass_Topjet1) return false;

 
  // PT Cuts 
  // const bool pass_Topjet2 = TOPjet2_sel->passes(event);
  // if(!pass_Topjet2) return false;
  topjet_topjet2_h->fill(event);
  eff_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);
  muon_topjet2_h->fill(event);
  event_topjet2_h->fill(event);

 if(berror)  std::cout<<"SelectionModule L:268 vor BTag"<<std::endl;
  ////////////////////////////////////////////////////////// Btag//////////////////////////////////////////////////////////////////////////////////

  // const bool pass_btag = btag_sel->passes(event);
  // if(!pass_btag) return false;

  // topjet_btag_h->fill(event);
  // eff_btag_h->fill(event);
  // jet_btag_h->fill(event);
  // muon_btag_h->fill(event);
  // event_btag_h->fill(event);




  // std::cout<<"SelectionModule L:294 vor RelIso"<<std::endl;  
  // ////////////////////////////////////////////////////////// Reliso //////////////////////////////////////////////////////////////////////////////////
  // bool pass_reliso= reliso_sel->passes(event);
  // if(pass_reliso){
  //   topjet_reliso_h->fill(event);
  //   eff_reliso_h->fill(event);
  //   jet_reliso_h->fill(event);
  //   muon_reliso_h->fill(event);
  //   event_reliso_h->fill(event);
  // }

 if(berror)  std::cout<<"SelectionModule L:294 vor TwoDCut"<<std::endl;  
  ////////////////////////////////////////////////////////// TwoDCut //////////////////////////////////////////////////////////////////////////////////
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  topjet_twodcut_h->fill(event);
  jet_twodcut_h->fill(event);
 eff_twodcut_h->fill(event);
 muon_twodcut_h->fill(event);

  jet_cleaner1->process(event);
  event_twodcut_h->fill(event);

  //save uncleaned AK8 and AK4 jets
  std::unique_ptr< std::vector<TopJet> > uncleaned_AK8 (new std::vector<TopJet> (*event.topjets));
  event.set(h_AK8,*uncleaned_AK8 );
  std::unique_ptr< std::vector<Jet> > uncleaned_AK4 (new std::vector<Jet> (*event.jets));
  event.set(h_AK4,* uncleaned_AK4);


const bool pass_jet2 = jet2_sel->passes(event);
  if(!pass_jet2) return false;

 if(berror) std::cout<<"SelectionModule L:338 vor SidebandReco"<<std::endl;
  /////////////////////////////////////////////////////////// Sideband Reco //////////////////////////////////////////////////////////////////////////////////
 reco_primlep->process(event);

 if(! ZprimeTotTPrime_reco->process(event))return false;
 ZprimeTotTPrime_chi->process(event); 
 // save only the chi2-best ttbar hypothesis in output sub-ntuple
 std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps = event.get(h_ZprimeTotTPrime_hyps);
 const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis(hyps, "Chi2");
 if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
 const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj(*hyp);

// chi2min_clean_h->fill(event);
//  topjet_clean_h->fill(event);
//  eff_clean_h->fill(event);
//  jet_clean_h->fill(event);
//  muon_clean_h->fill(event);
//  event_clean_h->fill(event);
  ////////////////////////////////////////////////////////////////
 std::vector<Jet>* usedak4(new std::vector<Jet> (*event.jets));
 usedak4->clear();
 usedak4->reserve(event.jets->size());
 for(const Jet jet: hyp->toplep_jets()) usedak4->push_back(jet);
 for(const Jet jet: hyp->tophad_jets()) usedak4->push_back(jet);

 if(berror) std::cout<<"SelectionModule L:858 Size used ak4 "<<usedak4->size() <<std::endl;

 //Cleaning AK8 by overlap of AK4
 std::vector<TopJet>* AK8Jets(new std::vector<TopJet> (*event.topjets));
 AK8Jets->clear();
 AK8Jets->reserve(event.topjets->size());

 if(berror) std::cout<<"SelectionModule L:858 Size AK8 before cleaning "<<event.topjets->size() <<std::endl;
 for(const TopJet ak8:*event.topjets){
 bool bdeltaR=true;
   for(const Jet ak4:*usedak4){
     double deltar = deltaR(ak4,ak8);
if(berror) std::cout<<"SelectionModule L:858 DeltaR(ak4, ak8)<1.2 "<<deltar <<std::endl;
     if(deltar < 1.2) bdeltaR=false;
if(berror) std::cout<<"SelectionModule L:858 bdeltaR  "<<bdeltaR <<std::endl;
    }
    if(bdeltaR)AK8Jets ->push_back(ak8);
  }
if(berror) std::cout<<"SelectionModule L:858 Size ak8 "<<AK8Jets->size() <<std::endl;
  sort_by_pt<TopJet>(*AK8Jets);
  if(!(AK8Jets->size()) )return false;
  //  if((60 < AK8Jets->at(0).v4().M())&&( AK8Jets->at(0).v4().M()<115))return false;
  event.topjets->clear();
  event.topjets->reserve(AK8Jets->size());
 for(const auto & j : *AK8Jets) event.topjets->push_back(j); 
 sort_by_pt<TopJet>(*event.topjets);
if(berror) std::cout<<"SelectionModule L:858 Size topjets Collection "<<event.topjets->size() <<std::endl;
 ak4_cleaner->process(event);
  ////////////////////////////////////////////////////////////////



  chi2min_reco_h->fill(event);
 // topjet_reco_h->fill(event);
 // eff_reco_h->fill(event);
 // jet_reco_h->fill(event);
 // muon_reco_h->fill(event);
 // event_reco_h->fill(event);

 hyps.clear();
 hyps.push_back(hyp_obj);

 chi2min_reco_h->fill(event);




 //////////////////////////////////////////////////////////  CHi2 Cut  ////////////////////////////////////////////////////////

 bool pass_chi2cut = chi2cut_sel->passes(event);
 if(pass_chi2cut) return false;
 chi2min_chi2cut_h->fill(event);
 topjet_chi2cut_h->fill(event);
 eff_chi2cut_h->fill(event);
 jet_chi2cut_h->fill(event);
 muon_chi2cut_h->fill(event);
 event_chi2cut_h->fill(event);

/////////////////////////////////////////////////////////  btag  ////////////////////////////////////////////////////////
bool pass_btag1 = btag1_sel->passes(event);
 bool pass_btag0 = btag0_sel->passes(event);

 if(pass_btag1){
   topjet_btag1_h->fill(event);
   eff_btag1_h->fill(event);
   jet_btag1_h->fill(event);
   muon_btag1_h->fill(event);
   event_btag1_h->fill(event); 
   chi2min_btag1_h->fill(event);
 }

 if(pass_btag0){
   topjet_btag0_h->fill(event);
   eff_btag0_h->fill(event);
   jet_btag0_h->fill(event);
   muon_btag0_h->fill(event);
   event_btag0_h->fill(event);
   chi2min_btag0_h->fill(event);
 }


 ///////////      counting signal events      /////////////////////////
 /////AK8 und AK4 auf urspruenglich zurueck setzten

 event.topjets->clear();
 event.topjets->reserve(uncleaned_AK8->size());
 for(const auto & j : *uncleaned_AK8) event.topjets->push_back(j); 
 sort_by_pt<TopJet>(*event.topjets);


 const bool pass_higgstag = higgstag_sel->passes(event);
 bool pass_zwtag=false;
 if(!pass_higgstag){
   pass_zwtag = ZWtag_sel->passes(event);
   if(!pass_zwtag) return false;
 }
 const bool pass_toptag = toptag_sel->passes(event);

std::unique_ptr< std::vector<TopJet> > higgsjets_all (new std::vector<TopJet> (*event.topjets));
  //cleanen der topjets um sie in einen neuen Vektor zu speichern
  higgstag_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);
  ak4_cleaner->process(event);
  std::unique_ptr< std::vector<TopJet>>  topjets_higgstag(new std::vector<TopJet> (*event.topjets));
  if(topjets_higgstag->size()) { std::cout <<"Size der subjets Higgstag"<< topjets_higgstag->at(0).subjets().size()<< std::endl;}
  // zurueck speichern der all topjets
  event.topjets->clear();
  event.topjets->reserve(higgsjets_all->size());
  for(const auto & j : *higgsjets_all) event.topjets->push_back(j); 
  sort_by_pt<TopJet>(*event.topjets);
  //handle auf toptag Topjets
  event.set(h_higgstag,*topjets_higgstag );

   std::unique_ptr< std::vector<TopJet> > ZWjets_all(new std::vector<TopJet> (*event.topjets));
    //cleanen der topjets um sie in einen neuen Vektor zu speichern
    ZWtag_cleaner->process(event);
    sort_by_pt<TopJet>(*event.topjets);
  if(pass_zwtag){
    ak4_cleaner->process(event);
  }
    std::unique_ptr< std::vector<TopJet>>  topjets_ZWtag(new std::vector<TopJet> (*event.topjets));
    // zurueck speichern der all topjets
    event.topjets->clear();
    event.topjets->reserve(ZWjets_all->size());
    for(const auto & j : *ZWjets_all) event.topjets->push_back(j); 
    sort_by_pt<TopJet>(*event.topjets);
    //handle auf toptag Topjets
    event.set(h_ZWtag,*topjets_ZWtag );

 //alle topjets
  std::unique_ptr< std::vector<TopJet> > topjets_all (new std::vector<TopJet> (*event.topjets));
  //cleanen der topjets um sie in einen neuen Vektor zu speichern
  toptag_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);
  if(pass_toptag){
    ak4_cleaner->process(event);
  }
  std::unique_ptr< std::vector<TopJet>>  topjets_toptag(new std::vector<TopJet> (*event.topjets));
  // zurueck speichern der all topjets
  event.topjets->clear();
  event.topjets->reserve(topjets_all->size());
  for(const auto & j : *topjets_all) event.topjets->push_back(j); 
  sort_by_pt<TopJet>(*event.topjets);
  //handle auf toptag Topjets
  event.set(h_toptag,*topjets_toptag );


 reco_primlep->process(event);
 if(!(ZprimeTotTPrime_selection_reco->process(event)))return false;
 ZprimeTotTPrime_chi->process(event); 
 // save only the chi2-best ttbar hypothesis in output sub-ntuple
 std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps2 = event.get(h_ZprimeTotTPrime_hyps);
 const ZPrimeTotTPrimeReconstructionHypothesis* hyp2 = get_best_hypothesis(hyps2, "Chi2");
 if(!hyp2) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
 const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj2(*hyp2);
 hyps2.clear();
 hyps2.push_back(hyp_obj2);
  pass_chi2cut = chi2cut_sel->passes(event);
 if(!pass_chi2cut) return false;

 eff_selection_h->fill(event);
 chi2min_selection_h->fill(event);
 topjet_selection_h->fill(event);
 jet_selection_h->fill(event);
 muon_selection_h->fill(event);
 event_selection_h->fill(event);

 // if(berror)  std::cout<<"SelectionModule L:338 vor HiggsTAGGER"<<std::endl;
 // /////////////////////////////////////////////////////////// Higgs TAGGER //////////////////////////////////////////////////////////////////////////////////
 
  // std::unique_ptr< std::vector<TopJet> > higgsjets_all (new std::vector<TopJet> (*event.topjets));
  // //cleanen der topjets um sie in einen neuen Vektor zu speichern
  // higgstag_cleaner->process(event);
  // sort_by_pt<TopJet>(*event.topjets);
  // ak4_cleaner->process(event);
  // std::unique_ptr< std::vector<TopJet>>  topjets_higgstag(new std::vector<TopJet> (*event.topjets));
  // if(topjets_higgstag->size()) { std::cout <<"Size der subjets Higgstag"<< topjets_higgstag->at(0).subjets().size()<< std::endl;}
  // // zurueck speichern der all topjets
  // event.topjets->clear();
  // event.topjets->reserve(higgsjets_all->size());
  // for(const auto & j : *higgsjets_all) event.topjets->push_back(j); 
  // sort_by_pt<TopJet>(*event.topjets);
  // //handle auf toptag Topjets
  // event.set(h_higgstag,*topjets_higgstag );

  // // bool b_higgs =false;
  // // bool b_zw =false;
  // // bool b_top =false;

  // // for(const auto j: *event.topjets){
  // //   LorentzVector subjet_sum;
  // //   for (const auto s : j.subjets()) {
  // //     subjet_sum += s.v4();
  // //   }
  // //   double  mHiggs_rec=subjet_sum.M();
  
  // //   if(mHiggs_rec < 150 && mHiggs_rec > 100 ) b_higgs = true;
  // //   if(mHiggs_rec < 115 && mHiggs_rec > 60 ) b_zw = true;
  // //   if(mHiggs_rec < 240 && mHiggs_rec > 150 )  b_top = true;
  // // }
  // // if(b_higgs)topjet_higgs_masswindow_h->fill(event);
  // // if(b_zw)topjet_zw_masswindow_h->fill(event);
  // // if(b_top)topjet_top_masswindow_h->fill(event);


  // //Higgstag
  //  const bool pass_higgstag = higgstag_sel->passes(event);
  //  // const bool pass_higgstag =false;
  // input_zwtag_h->fill(event);
  // input_higgstag_h->fill(event);
  // if(pass_higgstag){
  //   output_higgstag_h->fill(event);
  //   topjet_higgstag_h->fill(event);
  //   eff_higgstag_h->fill(event);
  //   jet_higgstag_h->fill(event);
  //   muon_higgstag_h->fill(event);
  //   event_higgstag_h->fill(event);
  // }


 // ///////////////////////////////////////////////////////// HIGGSTAGGER || ZW TAGGER //////////////////////////////////////////////////////////////////////////////////
 //  bool pass_zwtag=false;
 //  if(!pass_higgstag){
 //    pass_zwtag = ZWtag_sel->passes(event);
 //    if(!pass_zwtag) return false;
    
 //    topjet_zwtag_h->fill(event);
 //    eff_zwtag_h->fill(event);
 //    jet_zwtag_h->fill(event);
 //    muon_zwtag_h->fill(event);
 //    event_zwtag_h->fill(event);
 //  }

 //  std::cout<<"SelectionModule L:338 vor ZWTAGGER"<<std::endl;  
 //    ///////////////////////////////////////////////////////// ZW TAGGER //////////////////////////////////////////////////////////////////////////////////

 //    std::unique_ptr< std::vector<TopJet> > ZWjets_all(new std::vector<TopJet> (*event.topjets));
 //    //cleanen der topjets um sie in einen neuen Vektor zu speichern
 //    ZWtag_cleaner->process(event);
 //    sort_by_pt<TopJet>(*event.topjets);
 //  if(pass_zwtag){
 //    ak4_cleaner->process(event);
 //  }
 //    std::unique_ptr< std::vector<TopJet>>  topjets_ZWtag(new std::vector<TopJet> (*event.topjets));
 //    // zurueck speichern der all topjets
 //    event.topjets->clear();
 //    event.topjets->reserve(ZWjets_all->size());
 //    for(const auto & j : *ZWjets_all) event.topjets->push_back(j); 
 //    sort_by_pt<TopJet>(*event.topjets);
 //    //handle auf toptag Topjets
 //    event.set(h_ZWtag,*topjets_ZWtag );
 //    ///////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////////
 //    if(pass_zwtag){  
 //      output_zwtag_h->fill(event);
 //    }

 //  output_combined_h->fill(event); 


 //  std::cout<<"SelectionModule L:338 vor TopTAGGER"<<std::endl;  
  ///////////////////////////////////////////////////////// TOP TAGGER //////////////////////////////////////////////////////////////////////////////////
 // const bool pass_toptag = toptag_sel->passes(event);

 //  //alle topjets
 //  std::unique_ptr< std::vector<TopJet> > topjets_all (new std::vector<TopJet> (*event.topjets));
 //  //cleanen der topjets um sie in einen neuen Vektor zu speichern
 //  toptag_cleaner->process(event);
 //  sort_by_pt<TopJet>(*event.topjets);
 //  if(pass_toptag){
 //    ak4_cleaner->process(event);
 //  }
 //  std::unique_ptr< std::vector<TopJet>>  topjets_toptag(new std::vector<TopJet> (*event.topjets));
 //  // zurueck speichern der all topjets
 //  event.topjets->clear();
 //  event.topjets->reserve(topjets_all->size());
 //  for(const auto & j : *topjets_all) event.topjets->push_back(j); 
 //  sort_by_pt<TopJet>(*event.topjets);
 //  //handle auf toptag Topjets
 //  event.set(h_toptag,*topjets_toptag );

 //  //Toptag
 
 //  output_toptag_h->fill(event);
 //  topjet_toptag_h->fill(event);
 //  eff_toptag_h->fill(event);
 //  jet_toptag_h->fill(event);
 //  muon_toptag_h->fill(event);
 //  event_toptag_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_topjet_toptag_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_toptag_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_toptag_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_toptag_h->fill(event);
 
 ///////////////////////////////////////////////////////////////// Tagger Hist ///////////////////////////////////////////////////////
  // topjet_tagger_h->fill(event);
  // eff_tagger_h->fill(event);
  // jet_tagger_h->fill(event);
  // muon_tagger_h->fill(event);
  // event_tagger_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_topjet_tagger_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_tagger_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_tagger_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_tagger_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_jet_tagger_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_jet_tagger_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_jet_tagger_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_jet_tagger_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_muon_tagger_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_muon_tagger_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_muon_tagger_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_muon_tagger_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_event_tagger_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_event_tagger_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_event_tagger_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_event_tagger_h->fill(event);

 // if(berror) std::cout<<"SelectionModule L:303 vor KinReco"<<std::endl;
  // //////////////////////////////////////////////////////////// KIN RECO///////////////////////////////////////////////////////////////////////////////////
  
  // reco_primlep->process(event);
  // if(!(ZprimeTotTPrime_reco->process(event)))return false;
  // ZprimeTotTPrime_chi->process(event); 
  // // save only the chi2-best ttbar hypothesis in output sub-ntuple
  // std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps = event.get(h_ZprimeTotTPrime_hyps);
  // const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis(hyps, "Chi2");
  // if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
  // const ZPrimeTotTPrimeReconstructionnHypothesis hyp_obj(*hyp);
  // // std::cout << "Test in Selection Module " << hyp->HZW_subjets().size() << " v4 des ersten " <<hyp->HZW_subjets().at(0).v4()<<std::endl;
  // // chi2min_reco_h->fill(event);
  // // topjet_reco_h->fill(event);
  // // eff_reco_h->fill(event);
  // // jet_reco_h->fill(event);
  // // muon_reco_h->fill(event);
  // // event_reco_h->fill(event);

  // // if(pass_higgstag && pass_toptag) higgs_top_topjet_reco_h->fill(event);
  // // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_reco_h->fill(event);
  // // if(pass_zwtag && pass_toptag)zw_top_topjet_reco_h->fill(event);
  // // if(pass_zwtag && !pass_toptag)zw_notop_topjet_reco_h->fill(event);

  // // if(pass_higgstag && pass_toptag) higgs_top_chi2min_reco_h->fill(event);
  // // if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_reco_h->fill(event);
  // // if(pass_zwtag && pass_toptag)zw_top_chi2min_reco_h->fill(event);
  // // if(pass_zwtag && !pass_toptag)zw_notop_chi2min_reco_h->fill(event);

  // hyps.clear();
  // hyps.push_back(hyp_obj);
  // std::cout << "SelectionModule L:715 vor Chi2Cut"<<std::endl;
  // //////////////////////////////////////////////////////////  CHi2 Cut  ////////////////////////////////////////////////////////
  // bool pass_chi2cut = chi2cut_sel->passes(event);
  // if(!pass_chi2cut) return false;

 // //////////////////////////////////////////////////////////// Sideband KIN RECO///////////////////////////////////////////////////////////////////////////////////
// reco_primlep->process(event);


//  //////////////////////////////////////////////////////////  CHi2 Cut  ////////////////////////////////////////////////////////
//  bool pass_chi2cut = chi2cut_sel->passes(event);
//    if(!pass_chi2cut) return false;
//   chi2min_chi2cut_h->fill(event);
//   topjet_chi2cut_h->fill(event);
//   eff_chi2cut_h->fill(event);
//   jet_chi2cut_h->fill(event);
//   muon_chi2cut_h->fill(event);
//   event_chi2cut_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_topjet_chi2cut_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_chi2cut_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_chi2cut_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_chi2cut_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_chi2min_chi2cut_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_chi2cut_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_chi2min_chi2cut_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_chi2min_chi2cut_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_jet_chi2cut_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_jet_chi2cut_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_jet_chi2cut_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_jet_chi2cut_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_muon_chi2cut_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_muon_chi2cut_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_muon_chi2cut_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_muon_chi2cut_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_event_chi2cut_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_event_chi2cut_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_event_chi2cut_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_event_chi2cut_h->fill(event);

 


 // if(berror)  std::cout << "SelectionModule L:877 vor TestSection"<<std::endl;
 //  ////////////////////////////////////////////////////////Met Selection////////////////////////////////////////////////////////////////////////////////////////
  
 //  bool pass_met = met_sel->passes(event);
 //  if(!pass_met) return false;
 //    chi2min_met_h->fill(event);
 //    topjet_met_h->fill(event);
 //    eff_met_h->fill(event);
 //    jet_met_h->fill(event);
 //    muon_met_h->fill(event);
 //    event_met_h->fill(event);
 
   //  if(pass_higgstag && pass_toptag) higgs_top_topjet_met_h->fill(event);
//     if(pass_higgstag && !pass_toptag) higgs_notop_topjet_met_h->fill(event);
//     if(pass_zwtag && pass_toptag)zw_top_topjet_met_h->fill(event);
//     if(pass_zwtag && !pass_toptag)zw_notop_topjet_met_h->fill(event);

//     if(pass_higgstag && pass_toptag) higgs_top_chi2min_met_h->fill(event);
//     if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_met_h->fill(event);
//     if(pass_zwtag && pass_toptag)zw_top_chi2min_met_h->fill(event);
//     if(pass_zwtag && !pass_toptag)zw_notop_chi2min_met_h->fill(event);

// if(pass_higgstag && pass_toptag) higgs_top_jet_met_h->fill(event);
//   if(pass_higgstag && !pass_toptag) higgs_notop_jet_met_h->fill(event);
//   if(pass_zwtag && pass_toptag)zw_top_jet_met_h->fill(event);
//   if(pass_zwtag && !pass_toptag)zw_notop_jet_met_h->fill(event);

//   if(pass_higgstag && pass_toptag) higgs_top_muon_met_h->fill(event);
//   if(pass_higgstag && !pass_toptag) higgs_notop_muon_met_h->fill(event);
//   if(pass_zwtag && pass_toptag)zw_top_muon_met_h->fill(event);
//   if(pass_zwtag && !pass_toptag)zw_notop_muon_met_h->fill(event);

//   if(pass_higgstag && pass_toptag) higgs_top_event_met_h->fill(event);
//   if(pass_higgstag && !pass_toptag) higgs_notop_event_met_h->fill(event);
//   if(pass_zwtag && pass_toptag)zw_top_event_met_h->fill(event);
//   if(pass_zwtag && !pass_toptag)zw_notop_event_met_h->fill(event);    
  


 // ////////////////////////////////////////////////////////HT Selection////////////////////////////////////////////////////////////////////////////////////////
 //  bool pass_ht = ht_sel->passes(event); 

 //  if(!pass_ht)return false;
 //  chi2min_ht_h->fill(event);
 //  topjet_ht_h->fill(event);
 //  eff_ht_h->fill(event);
 //  jet_ht_h->fill(event);
 //  muon_ht_h->fill(event);
 //  event_ht_h->fill(event);
 
  // if(pass_higgstag && pass_toptag) higgs_top_topjet_ht_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_ht_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_ht_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_ht_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_chi2min_ht_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_ht_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_chi2min_ht_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_chi2min_ht_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_jet_ht_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_jet_ht_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_jet_ht_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_jet_ht_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_muon_ht_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_muon_ht_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_muon_ht_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_muon_ht_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_event_ht_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_event_ht_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_event_ht_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_event_ht_h->fill(event);

   
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeSidebandModule)
