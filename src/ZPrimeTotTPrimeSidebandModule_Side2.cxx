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
#include "UHH2/common/include/TopPtReweight.h"

#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h>


#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstruction.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHypothesisHists.h" 
#include "UHH2/ZprimeToTprimeTtZtH/include/SR_PDFHists.h"

// #include <UHH2/common/include/HypothesisHists.h>
// #include <UHH2/common/include/TTbarReconstruction.h>
// #include <UHH2/common/include/ReconstructionHypothesis.h>
// #include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>

using namespace uhh2examples;
using namespace uhh2;

class ZPrimeTotTPrimeSidebandModule_Side2 : public uhh2::AnalysisModule {

 public:
  explicit ZPrimeTotTPrimeSidebandModule_Side2(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  enum lepton { muon, elec };
  lepton channel_;

  // cleaners
  std::unique_ptr<MuonCleaner>     muo_cleaner;
  std::unique_ptr<ElectronCleaner> ele_cleaner1; 
  std::unique_ptr<ElectronCleaner> ele_cleaner2; 
  std::unique_ptr<ElectronCleaner> photon_cleaner;

  std::unique_ptr<JetCorrector> jet_corrector;
  std::unique_ptr<JetCorrector> jet_corrector_BCD;
  std::unique_ptr<JetCorrector> jet_corrector_EF;
  std::unique_ptr<JetCorrector> jet_corrector_G;
  std::unique_ptr<JetCorrector> jet_corrector_H;

  std::unique_ptr<TopJetCorrector> topjet_corrector;
  std::unique_ptr<TopJetCorrector> topjet_corrector_BCD;
  std::unique_ptr<TopJetCorrector> topjet_corrector_EF;
  std::unique_ptr<TopJetCorrector> topjet_corrector_G;
  std::unique_ptr<TopJetCorrector> topjet_corrector_H;

  std::unique_ptr<SubJetCorrector> subjet_corrector;
  std::unique_ptr<SubJetCorrector> subjet_corrector_BCD;
  std::unique_ptr<SubJetCorrector> subjet_corrector_EF;
  std::unique_ptr<SubJetCorrector> subjet_corrector_G;
  std::unique_ptr<SubJetCorrector> subjet_corrector_H;

  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner_BCD;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner_EF;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner_G;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner_H;


  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCleaner>      jet_cleaner2;
  std::unique_ptr<JetCleaner>      jet_cleaner1;
  std::unique_ptr<TopJetCleaner>   topjet_cleaner;
  std::unique_ptr<TopJetCleaner>   topjet_mass_cleaner;
  std::unique_ptr<TopJetCleaner>   toptag_cleaner;
  std::unique_ptr<TopJetCleaner>   higgstag_cleaner;
  std::unique_ptr<TopJetCleaner>   ZWtag_cleaner;
  std::unique_ptr<TopJetCleaner>   toptag_uncleaner;
  std::unique_ptr<TopJetCleaner>   higgstag_uncleaner;
  std::unique_ptr<TopJetCleaner>   ZWtag_uncleaner;
  std::unique_ptr<JetCleaner>  ak4_cleaner;
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;
  std::unique_ptr<uhh2::Selection> triggerMu50_sel;
  std::unique_ptr<uhh2::Selection> triggerTrkMu50_sel;
  std::unique_ptr<uhh2::Selection> triggerElec_sel1;
  std::unique_ptr<uhh2::Selection> triggerElec_sel2;
  std::unique_ptr<uhh2::Selection> triggerPhoton_sel;
  std::unique_ptr<AndSelection>  metfilters_selection;

  std::unique_ptr<uhh2::Selection> genwpt_cut;  
 

  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::unique_ptr<uhh2::AnalysisModule> muonscale;
  std::unique_ptr<MuonTrkWeights> muontracker_sf;
  std::unique_ptr<uhh2::AnalysisModule> elecid2;
  std::unique_ptr<uhh2::AnalysisModule> elecid1;
  std::unique_ptr<uhh2::AnalysisModule> elecreco;
  std::unique_ptr<uhh2::AnalysisModule> electrigger;
  std::unique_ptr<uhh2::AnalysisModule> triggerscale;
  std::unique_ptr<uhh2::AnalysisModule> btagwAK4;
  std::unique_ptr<uhh2::AnalysisModule> topptreweighting_all;

  //Selections
  std::unique_ptr<uhh2::Selection> lumi_sel;
  std::unique_ptr<uhh2::AndSelection> lep1_sel; //  one lepton (comment out )
  std::unique_ptr<uhh2::Selection> TOPjet1_sel; // at least 2 jets
  std::unique_ptr<uhh2::Selection> chi2cut_sel; // chi2min <50
  std::unique_ptr<uhh2::Selection> btag1_sel;
  std::unique_ptr<uhh2::Selection> btag0_sel;
  std::unique_ptr<uhh2::Selection> met_sel7;

  std::unique_ptr<uhh2::Selection> twodcut_sel;// pt 20 rel 0.4
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> met_sel;
  std::unique_ptr<uhh2::Selection> ht_sel;
 

  //TOP TAGGER
  std::unique_ptr<Selection> toptag_sel; // at least one toptag
  uhh2::Event::Handle< std::vector<TopJet> > h_toptag;

  //HIGGS TAGGER
  std::unique_ptr<Selection> higgstag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_higgstag;

  // //Z/W TAGGER
  // std::unique_ptr<Selection> ZWtag_sel;
  // uhh2::Event::Handle< std::vector<TopJet> > h_ZWtag;

  // reconstruction ZPrime for Signal
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrimeprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_reco;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_chi;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_chi;

  // reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
  std::unique_ptr<uhh2::Selection> genwpt_sel;

  
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
 std::unique_ptr<Hists> pdf_btag1_h;

std::unique_ptr<Hists> higgs_top_chi2min_btag1_h;
  std::unique_ptr<Hists> zw_top_chi2min_btag1_h;
  std::unique_ptr<Hists> higgs_notop_chi2min_btag1_h;
  std::unique_ptr<Hists> zw_notop_chi2min_btag1_h;


//Hist Btag0
  std::unique_ptr<Hists> eff_btag0_h;
  std::unique_ptr<Hists> jet_btag0_h;
  std::unique_ptr<Hists> muon_btag0_h;
  std::unique_ptr<Hists> event_btag0_h;
  std::unique_ptr<Hists> topjet_btag0_h;
 std::unique_ptr<Hists> chi2min_btag0_h;
 std::unique_ptr<Hists> pdf_btag0_h;

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

  //Triangularcut Hists
  std::unique_ptr<Hists> eff_triangcut_h;
  std::unique_ptr<Hists> jet_triangcut_h;
  std::unique_ptr<Hists> topjet_triangcut_h;
  std::unique_ptr<Hists> muon_triangcut_h;
  std::unique_ptr<Hists> electron_triangcut_h;
  std::unique_ptr<Hists> event_triangcut_h;
  std::unique_ptr<Hists> lumi_triangcut_h;
   
//Higgstag Hists
  std::unique_ptr<TopJetHists> input_higgstag_h;
  std::unique_ptr<TopJetHists> output_higgstag_h;

  std::unique_ptr<Hists> topjet_higgstag_h;
  std::unique_ptr<Hists> eff_higgstag_h;
  std::unique_ptr<Hists> jet_higgstag_h;
  std::unique_ptr<Hists> muon_higgstag_h;
  std::unique_ptr<Hists> event_higgstag_h;

  // //ZWTagg Hists
  // std::unique_ptr<TopJetHists> input_zwtag_h;
  // std::unique_ptr<TopJetHists> output_zwtag_h;

  // std::unique_ptr<Hists> topjet_zwtag_h;
  // std::unique_ptr<Hists> eff_zwtag_h;
  // std::unique_ptr<Hists> jet_zwtag_h;
  // std::unique_ptr<Hists> muon_zwtag_h;
  // std::unique_ptr<Hists> event_zwtag_h;

  //Higgstag || ZWtag Hists
  // std::unique_ptr<Hists> output_combined_h;

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
  std::unique_ptr<TopJetHists> topjet_selection_h;
  std::unique_ptr<Hists> event_selection_h;


  //   std::unique_ptr<Hists> topjet_top_masswindow_h;
  //  std::unique_ptr<Hists> topjet_higgs_masswindow_h;
  // std::unique_ptr<Hists> topjet_zw_masswindow_h;

  std::unique_ptr<Hists> lumi_h;


  //general
  TString filename;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ZprimeTotTPrime_hyps;
  uhh2::Event::Handle<double> h_ht;
  bool berror;
  bool isMC;

  std::unique_ptr<AnalysisModule> syst_module;
  bool do_scale_variation;

  const int runnr_BCD = 276811;
  const int runnr_EF = 278802;
  const int runnr_G = 280385;

  TString elec_up_down;
  TString  toppt_sys;
};


ZPrimeTotTPrimeSidebandModule_Side2::ZPrimeTotTPrimeSidebandModule_Side2(uhh2::Context& ctx){


  //choose channel
  const std::string& channel = ctx.get("channel", "");
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else throw std::runtime_error("ZprimeSelectionModule -- undefined argument for 'channel' key in xml file (must be 'muon' or 'elec'): "+channel);

  isMC = (ctx.get("dataset_type") == "MC");
  //// Data/MC scale
  auto data_dir_path = ctx.get("data_dir_path");
  string sysAK4=ctx.get("btagging_sys");

  if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx,ctx.get("puVariation"))); 
    lumiweight.reset(new MCLumiWeight(ctx)); 
    btagwAK4.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM,"jets",sysAK4,"mujets","incl","MCBtagEfficienciesAK4","_AK4","BTagCalibration"));  
    if(channel_==muon){
      muonscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta", 1.,"MediumID2016", true, ctx.get("Systematic_MuonID")));
      muontracker_sf.reset(new MuonTrkWeights(ctx, data_dir_path +"Tracking_EfficienciesAndSF_BCDEFGH.root",  ctx.get("Systematic_MuonTrk")));
      // triggerscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins", 1., "trigger", true, ctx.get("Systematic_MuonTrigger") ));
      triggerscale.reset(new MuonTriggerSF(ctx));
    }else if(channel_==elec){
      // elecid2.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "", ctx.get("Systematic_EleID")));
      // electrigger.reset(new ElectronTriggerWeights(ctx,"/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/ZprimeToTprimeTtZtH/hist/ElectronEfficiencies.root",ctx.get("Systematic_EleTrigger")));
      electrigger.reset(new ElectronTriggerSF(ctx));
      elecreco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", ctx.get("Systematic_EleReco")));
      elecid1.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_MVA90_ID.root", 1, "", ctx.get("Systematic_EleID")));
    } 
    topptreweighting_all.reset(new TopPtReweight(ctx, 0.0615 , -0.0005 ,"ttbargen","weight_ttbar", true,1));
  }
  else     lumi_sel.reset(new LumiSelection(ctx));

  PrimaryVertexId pvid=StandardPrimaryVertexId();
  metfilters.emplace_back(new PrimaryVertexCleaner(pvid));
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");
  //  metfilters_selection->add<TriggerSelection>("chargedHadronTrackResolutionFilter", "Flag_chargedHadronTrackResolutionFilter"); 
  //  metfilters_selection->add<TriggerSelection>("muonBadTrackFilter", "Flag_muonBadTrackFilter");
  metfilters_selection->add<NPVSelection>("1 good PV",1,-1,pvid);

  
 //JEC
  std::vector<std::string> JEC_AK4, JEC_AK8,JEC_AK4_BCD,JEC_AK4_EF,JEC_AK4_G,JEC_AK4_H,JEC_AK8_BCD,JEC_AK8_EF,JEC_AK8_G,JEC_AK8_H;
  if(isMC){

    JEC_AK4 = JERFiles::Summer16_23Sep2016_V4_L123_AK4PFchs_MC;
    JEC_AK8 = JERFiles::Summer16_23Sep2016_V4_L123_AK8PFchs_MC;
  }
  else {

    JEC_AK4_BCD =  JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK4PFchs_DATA;
    JEC_AK4_EF = JERFiles::Summer16_23Sep2016_V4_EF_L123_AK4PFchs_DATA;
    JEC_AK4_G =  JERFiles::Summer16_23Sep2016_V4_G_L123_AK4PFchs_DATA;
    JEC_AK4_H =  JERFiles::Summer16_23Sep2016_V4_H_L123_AK4PFchs_DATA;
    
    JEC_AK8_BCD =  JERFiles::Summer16_23Sep2016_V4_BCD_L123_AK8PFchs_DATA;
    JEC_AK8_EF =  JERFiles::Summer16_23Sep2016_V4_EF_L123_AK8PFchs_DATA;
    JEC_AK8_G =  JERFiles::Summer16_23Sep2016_V4_G_L123_AK8PFchs_DATA;
    JEC_AK8_H =  JERFiles::Summer16_23Sep2016_V4_H_L123_AK8PFchs_DATA;
   
  }

 if(isMC){ 
    jet_corrector.reset(new JetCorrector(ctx, JEC_AK4));
    topjet_corrector.reset(new TopJetCorrector(ctx, JEC_AK8));
    subjet_corrector.reset(new SubJetCorrector(ctx,JEC_AK4));
    jetlepton_cleaner.reset(new JetLeptonCleaner(ctx,JEC_AK4));
    jetlepton_cleaner->set_drmax(.4);
  }
  else {
   
    jet_corrector_BCD.reset(new JetCorrector(ctx, JEC_AK4_BCD));
    jet_corrector_EF.reset(new JetCorrector(ctx, JEC_AK4_EF));
    jet_corrector_G.reset(new JetCorrector(ctx,JEC_AK4_G ));
    jet_corrector_H.reset(new JetCorrector(ctx,JEC_AK4_H ));

    topjet_corrector_BCD.reset(new TopJetCorrector(ctx, JEC_AK8_BCD));
    topjet_corrector_EF.reset(new TopJetCorrector(ctx, JEC_AK8_EF));
    topjet_corrector_G.reset(new TopJetCorrector(ctx,JEC_AK8_G ));
    topjet_corrector_H.reset(new TopJetCorrector(ctx,JEC_AK8_H ));

    subjet_corrector_BCD.reset(new SubJetCorrector(ctx, JEC_AK4_BCD));
    subjet_corrector_EF.reset(new SubJetCorrector(ctx, JEC_AK4_EF));
    subjet_corrector_G.reset(new SubJetCorrector(ctx,JEC_AK4_G ));
    subjet_corrector_H.reset(new SubJetCorrector(ctx,JEC_AK4_H ));

    jetlepton_cleaner_BCD.reset(new JetLeptonCleaner(ctx, JEC_AK4_BCD));
    jetlepton_cleaner_EF.reset(new JetLeptonCleaner(ctx, JEC_AK4_EF));
    jetlepton_cleaner_G.reset(new JetLeptonCleaner(ctx,JEC_AK4_G ));
    jetlepton_cleaner_H.reset(new JetLeptonCleaner(ctx,JEC_AK4_H ));

    jetlepton_cleaner_BCD->set_drmax(.4);
    jetlepton_cleaner_EF->set_drmax(.4);
    jetlepton_cleaner_G->set_drmax(.4);
    jetlepton_cleaner_H->set_drmax(.4);

  }



 //// OBJ CLEANING
 muo_cleaner.reset(new MuonCleaner    (AndId<Muon>    (PtEtaCut  (53., 2.4), MuonIDMedium())));
 // ele_cleaner2.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Spring16_tight)));
 ele_cleaner1.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(125, 2.4), ElectronID_MVAGeneralPurpose_Spring16_loose)));//Arne
 photon_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(250, 2.4),ElectronID_MVAGeneralPurpose_Spring16_loose )));
 topjet_mass_cleaner.reset(new TopJetCleaner(ctx,TopJetId(TopjetMassCleaner(30))));

 const JetId jetID(JetPFID(JetPFID::WP_LOOSE));
 jet_IDcleaner.reset(new JetCleaner(ctx,jetID));
 jet_cleaner2.reset(new JetCleaner(ctx,15., 2.4));
 jet_cleaner1.reset(new JetCleaner(ctx,30., 2.4));
 ak4_cleaner.reset(new JetCleaner(ctx,JetId(ZPrimeTotTPrimeAK4cleaner(1.2))));
 htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
 htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
 htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));
 topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));

  topjet_cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(200., 2.4))));


 // TOPJET SELECTIONS
 TOPjet1_sel.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut( 250., 2.4))));
 
 
 //BTAG
 btag1_sel.reset(new NBTagSelection(1, -1));
 btag0_sel.reset(new NBTagSelection(0, 0));


 //2D Cut
 twodcut_sel.reset(new TwoDCut(.4, 40.));

 if     (channel_ == elec) triangc_sel.reset(new uhh2::AndSelection(ctx));
 else if(channel_ == muon) triangc_sel.reset(new uhh2::AndSelection(ctx)); // always true (no triangular cuts for muon channel)


 //TOP TAGGER
 const TopJetId topjetID = AndId<TopJet>(Type2TopTag(150,240,Type2TopTag::MassType::groomed), Tau32(0.6));
 //  const TopJetId untopjetID =VetoId<TopJet>(topjetID);
 toptag_sel.reset(new NTopJetSelection(1, -1, topjetID));
 toptag_cleaner.reset(new TopJetCleaner(ctx,topjetID));
 //toptag_uncleaner.reset(new TopJetCleaner(ctx,untopjetID));


 //Higgs TAGGER
 // const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150,VetoId<Jet>(CSVBTag(CSVBTag::WP_LOOSE))), Tau21(1) );
 // const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150), Tau21(1) );
  const TopJetId higgsjetID = AndId<TopJet>(ZPrimeTotTPrimeUnHiggstag(60,150), Tau21(1)  );
 higgstag_sel.reset(new NTopJetSelection(1, -1, higgsjetID));
 higgstag_cleaner.reset(new TopJetCleaner(ctx,higgsjetID));
 // higgstag_uncleaner.reset(new TopJetCleaner(ctx,unhiggsjetID));

 //W/Z TAGGER
 // const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed),VetoId<TopJet>( Tau21(0.4)));
 //  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
 // const TopJetId ZWjetID = AndId<TopJet>(ZPrimeTotTPrimeUnHiggstag(60,150) , Tau21(1) );
 // const TopJetId unZWjetID =VetoId<TopJet>(ZWjetID);
 // ZWtag_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
 //ZWtag_cleaner.reset(new TopJetCleaner(ctx,ZWjetID));
 // ZWtag_uncleaner.reset(new TopJetCleaner(ctx,unZWjetID));

 // const TopJetId untagID = AndId<TopJet>(VetoId<TopJet>(higgsjetID),VetoId<TopJet>(ZWjetID));

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
  const std::string ZprimeTotTPrime_hyps_label("ZPrimeTotTPrimeReconstruction");
 const std::string ZprimeTotTPrime_selection_hyps_label("ZPrimeTotTPrimeReconstruction"); 
 const std::string ttbar_hyps_label("TTbarReconstruction");
  const std::string ZprimeTotTPrime_chi2_label("Chi2");
  const std::string ttbar_chi2_label("Chi2");
  const std::string ttbar_gen_label ("ttbargen");

  reco_primlep.reset(new PrimaryLepton(ctx));

  ZprimeTotTPrime_reco.reset(new ZPrimeTotTPrimeReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrime_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrimeprod.reset(new ZPrimeGenProducer(ctx, ZprimeTotTPrime_gen_label, false));
  h_ZprimeTotTPrime_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ZprimeTotTPrime_hyps_label);
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));


  chi2cut_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));
// Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "MC_WJetsAll"){ genwpt_sel.reset(new ZPrimeTotTPrimePartonW());}
  else                                          genwpt_sel.reset(new uhh2::AndSelection(ctx));


  // //MetCut
  met_sel.reset(new ZPrimeTotTPrimeMETCut(90,-1));
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

  // //Hist reliso
  // topjet_reliso_h.reset(new TopJetHists(ctx, "topjet_reliso"));
  // eff_reliso_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reliso"));
  // jet_reliso_h.reset(new JetHists(ctx, "jet_reliso"));
  // muon_reliso_h.reset(new MuonHists(ctx, "muon_reliso"));
  // event_reliso_h.reset(new EventHists(ctx, "event_reliso"));

  //Hists btag0
  topjet_btag0_h.reset(new TopJetHists(ctx, "topjet_btag0"));
  eff_btag0_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag0"));
  jet_btag0_h.reset(new JetHists(ctx, "jet_btag0"));
  muon_btag0_h.reset(new MuonHists(ctx, "muon_btag0"));
  event_btag0_h.reset(new EventHists(ctx, "event_btag0"));
  chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));
  if(isMC){
    pdf_btag0_h.reset(new SR_PDFHists(ctx,"pdf_btag0" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    pdf_btag1_h.reset(new SR_PDFHists(ctx,"pdf_btag1" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
  }

  //Hists btag1
  topjet_btag1_h.reset(new TopJetHists(ctx, "topjet_btag1"));
  eff_btag1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag1"));
  jet_btag1_h.reset(new JetHists(ctx, "jet_btag1"));
  muon_btag1_h.reset(new MuonHists(ctx, "muon_btag1"));
  event_btag1_h.reset(new EventHists(ctx, "event_btag1"));
  chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

  higgs_top_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  higgs_notop_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_top_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_notop_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));


  //Hist twodcut
  topjet_twodcut_h.reset(new TopJetHists(ctx, "topjet_twodcut"));
  eff_twodcut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_twodcut"));
  jet_twodcut_h.reset(new JetHists(ctx, "jet_twodcut"));
  muon_twodcut_h.reset(new MuonHists(ctx, "muon_twodcut"));
  event_twodcut_h.reset(new EventHists(ctx, "event_twodcut"));

  //Hist triangcut
  topjet_triangcut_h.reset(new TopJetHists(ctx, "topjet_triangcut"));
  eff_triangcut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_triangcut"));
  jet_triangcut_h.reset(new JetHists(ctx, "jet_triangcut"));
  muon_triangcut_h.reset(new MuonHists(ctx, "muon_triangcut"));
  electron_triangcut_h.reset(new ElectronHists(ctx, "electron_triangcut"));
  event_triangcut_h.reset(new EventHists(ctx, "event_triangcut"));
  lumi_triangcut_h.reset(new LuminosityHists(ctx,"lumi_triangcut"));

  //HIGGSTAG HISTS
  input_higgstag_h.reset(new TopJetHists(ctx, "input_higgstag"));
  output_higgstag_h.reset(new TopJetHists(ctx, "output_higgstag"));
  output_higgstag_h->set_TopJetId(higgsjetID);

  eff_higgstag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_higgstag"));
  jet_higgstag_h.reset(new JetHists(ctx, "jet_higgstag"));
  muon_higgstag_h.reset(new MuonHists(ctx, "muon_higgstag"));
  event_higgstag_h.reset(new EventHists(ctx, "event_higgstag"));
  topjet_higgstag_h.reset(new TopJetHists(ctx, "topjet_higgstag"));
 
  // //Z/W TAG HISTS
  // input_zwtag_h.reset(new TopJetHists(ctx, "input_zwtag"));
  // output_zwtag_h.reset(new TopJetHists(ctx, "output_zwtag"));
  // output_zwtag_h->set_TopJetId(ZWjetID);

  // eff_zwtag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_zwtag"));
  // jet_zwtag_h.reset(new JetHists(ctx, "jet_zwtag"));
  // muon_zwtag_h.reset(new MuonHists(ctx, "muon_zwtag"));
  // event_zwtag_h.reset(new EventHists(ctx, "event_zwtag"));
  // topjet_zwtag_h.reset(new TopJetHists(ctx, "topjet_zwtag"));

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
 
  // //Hists Reco
  // topjet_reco_h.reset(new TopJetHists(ctx, "topjet_reco"));
  // eff_reco_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reco"));
  // jet_reco_h.reset(new JetHists(ctx, "jet_reco"));
  // muon_reco_h.reset(new MuonHists(ctx, "muon_reco"));
  // event_reco_h.reset(new EventHists(ctx, "event_reco"));
  // chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_reco",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

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
  // topjet_selection_h->set_TopJetId(untagID);
  jet_selection_h.reset(new JetHists(ctx, "jet_selection"));
  muon_selection_h.reset(new MuonHists(ctx, "muon_selection"));
  event_selection_h.reset(new EventHists(ctx, "event_selection"));
  chi2min_selection_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_selection",ZprimeTotTPrime_selection_hyps_label,ZprimeTotTPrime_chi2_label ));


  //Higgstag || ZWTag
  //  output_combined_h.reset(new ZPrimeTotTPrimeHists(ctx, "output_combined"));
 // topjet_higgs_masswindow_h.reset(new TopJetHists(ctx, "masswindow_higgs"));
 // topjet_zw_masswindow_h.reset(new TopJetHists(ctx, "masswindow_zw"));
 // topjet_top_masswindow_h.reset(new TopJetHists(ctx, "masswindow_top"));


  //general
  filename =  ctx.get("dataset_version");

  // if(filename.Contains("WJets")) genwpt_cut.reset(new ZPrimeTotTPrimePartonWCut(filename));
  // else genwpt_cut.reset(new uhh2::AndSelection(ctx));


  h_toptag = ctx.declare_event_output< std::vector<TopJet> > ("TopTag");
  h_higgstag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag");
  //  h_ZWtag = ctx.declare_event_output< std::vector<TopJet> > ("ZWTag");
  h_ht = ctx.get_handle<double>("HT");

  //systematicen
  syst_module.reset(new MCScaleVariation(ctx));
  
  do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");


  //Trigger
  triggerMu50_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
  triggerTrkMu50_sel.reset(new TriggerSelection("HLT_TkMu50_v*"));
  triggerElec_sel1.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));//Arne
  //  triggerElec_sel2.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));//Iso
  triggerPhoton_sel.reset(new TriggerSelection("HLT_Photon175_v*"));

  elec_up_down = ctx.get("Systematic_EleTrigger");
  toppt_sys = ctx.get("Systematic_toppt");
  berror=false;
  
}

bool ZPrimeTotTPrimeSidebandModule_Side2::process(uhh2::Event& event){
  
  if(berror)  std::cout<<"/////////////////////////////////////SelectionModule L:218 Am Anfang///////////////////////////////////////////////////////////////////////////////"<<std::endl;

  for (auto & mod : htcalc) {
    mod->process(event);
  }

  if(!metfilters_selection->passes(event)) return false;
  for(auto & m : metfilters){
    m->process(event);
  }

  //Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(filename.Contains("TTbar")){
    ttgenprod->process(event);
      if(!genmttbar_sel->passes(event)) return false;
      if (toppt_sys=="none") topptreweighting_all->process(event);
  }
  
  if(filename.Contains("WJetsAll")){
     genwpt_sel->passes(event);
     if(!genwpt_sel->passes(event)) return false;
  }
 

  //ZPrime Genrator Level
  if(filename.Contains("MC_ZPrime")){
    ZprimeTotTPrimeprod->process(event); 
  }

  // if(filename.Contains("WJets")){
  //   bool pass_genwptcut=genwpt_cut->passes(event);
  //   if(!pass_genwptcut) return false;
  // }


  if(berror)   std::cout<<"SelectionModule L:232 vor Input Histogrammen"<<std::endl;
  
  /////////////////////////////////////////////////////////// Common Modules   ///////////////////////////////////////////////////////////////////////////////


  //common Modules
  /* luminosity sections from CMS golden-JSON file */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* pileup SF */
  if(!event.isRealData){ 
    pileup_SF->process(event);
    lumiweight->process(event);  
    if(channel_ == muon)    muonscale->process(event);
    if(channel_ == muon) triggerscale->process(event);
    if(channel_ == elec) elecreco->process(event);
  }
  ////

//systematicen
  if(do_scale_variation) syst_module->process(event); 

 ///correctors
  if(isMC){
    jetlepton_cleaner->process(event);
    jet_corrector->process(event);
    topjet_corrector->process(event);
    subjet_corrector->process(event);
    jet_corrector->correct_met(event);
  }else{
    if(event.run <= runnr_BCD)  {   
      jetlepton_cleaner_BCD->process(event);    
      jet_corrector_BCD->process(event);
      topjet_corrector_BCD->process(event);
      subjet_corrector_BCD->process(event);
      jet_corrector_BCD->correct_met(event);
    }
    else if(event.run < runnr_EF){  
      jetlepton_cleaner_EF->process(event);     
      jet_corrector_EF->process(event);
      topjet_corrector_EF->process(event);
      subjet_corrector_EF->process(event);
      jet_corrector_EF->correct_met(event);
    } 
    else if(event.run <= runnr_G) {  
      jetlepton_cleaner_G->process(event);     
      jet_corrector_G->process(event);
      topjet_corrector_G->process(event);
      subjet_corrector_G->process(event);
      jet_corrector_G->correct_met(event);
    } 
    else if(event.run > runnr_G) {   
      jetlepton_cleaner_H->process(event);    
      jet_corrector_H->process(event);
      topjet_corrector_H->process(event);
      subjet_corrector_H->process(event);
      jet_corrector_H->correct_met(event);
    } 
  }


  // OBJ CLEANING
  muo_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);


  jet_IDcleaner->process(event);
  jet_cleaner2->process(event); 
  
  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);


  // ///////Trigger///////
  if(channel_==muon){
    if (!isMC){
      bool pass_Mu50 = triggerMu50_sel->passes(event);
      if (event.run<274954) {
	if (!pass_Mu50) return false;
      }
      else {
	bool pass_TrkMu50 = triggerTrkMu50_sel->passes(event);
	if (!(pass_Mu50||pass_TrkMu50)) return false;
      }
    }else{
      bool pass_Mu50 = triggerMu50_sel->passes(event);
      bool pass_TrkMu50 = triggerTrkMu50_sel->passes(event);
      if (!(pass_Mu50||pass_TrkMu50)) return false;
    }
  }

  bool pass_elecTrigger = triggerElec_sel1 ->passes(event);
  bool pass_photonTrigger = triggerPhoton_sel->passes(event);
  // bool pass_elecTrigger2 = triggerElec_sel2 ->passes(event);
  if(channel_==elec){

    if(!(pass_elecTrigger||pass_photonTrigger)) return false;
    if(filename.Contains("photon") && event.electrons->at(0).pt()<250) return false;
    if(filename.Contains("rereco") && event.electrons->at(0).pt()>250) return false;

    if(pass_elecTrigger) ele_cleaner1->process(event);
    else if(pass_photonTrigger) photon_cleaner->process(event);
    sort_by_pt<Electron>(*event.electrons);
    if(isMC)  elecid1->process(event);
   

  }


  /////////////////////////////////////////////////////////// Input Histogramme ///////////////////////////////////////////////////////////////////////////////
 
  input_eff_h ->fill(event);
  input_event_h->fill(event);
  input_topjet_h->fill(event);
  input_jet_h->fill(event);
  input_lep_h->fill(event);



 if(berror)  std::cout<<"SelectionModule L:268 vor LeptonSelection"<<std::endl;
  /////////////////////////////////////////////////////////// LEPTON selection ///////////////////////////////////////////////////////////////////////////////
  const bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  if(channel_==elec )if(event.electrons->at(0).eta()<1.566 && event.electrons->at(0).eta()>1.4442) return false;

  // muon_lep1_h->fill(event);
  // eff_lep1_h ->fill(event);
  // event_lep1_h->fill(event);
  // topjet_lep1_h->fill(event);
  // jet_lep1_h->fill(event);

  if(!event.isRealData){
    if(channel_ == muon) muontracker_sf->process(event);
    if(channel_ == elec) electrigger->process(event);
  }

  


  topjetlepton_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  if(berror)  std::cout<<"SelectionModule L:268 vor TopJetSelection"<<std::endl;
  ////////////////////////////////////////////////////////// TopJET selection//////////////////////////////////////////////////////////////////////////////////

  topjet_mass_cleaner->process(event);
  topjet_cleaner->process(event);
  const bool pass_Topjet1 = TOPjet1_sel->passes(event);
  if(!pass_Topjet1) return false;

 


  topjet_topjet2_h->fill(event);
  eff_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);
  muon_topjet2_h->fill(event);
  event_topjet2_h->fill(event);

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

  if(berror)   std::cout<<"SelectionModule L:294 vor TwoDCut"<<std::endl;  
  ////////////////////////////////////////////////////////// TwoDCut //////////////////////////////////////////////////////////////////////////////////
  sort_by_pt<Jet>(*event.jets);
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  topjet_twodcut_h->fill(event);
  jet_twodcut_h->fill(event);
  eff_twodcut_h->fill(event);
  muon_twodcut_h->fill(event);

  jet_cleaner1->process(event);
  event_twodcut_h->fill(event);


 if(berror)   std::cout<<"SelectionModule L:294 vor Triangcut"<<std::endl;  
  ////////////////////////////////////////////////////////// Triangcut //////////////////////////////////////////////////////////////////////////////////
 // const bool pass_triangc = triangc_sel->passes(event);
  //  if(!pass_triangc) return false;
  topjet_triangcut_h->fill(event);
  jet_triangcut_h->fill(event);
  eff_triangcut_h->fill(event);
  muon_triangcut_h->fill(event);
  electron_triangcut_h->fill(event);
  lumi_triangcut_h->fill(event);
  event_triangcut_h->fill(event);

 
  if(berror)  std::cout<<"SelectionModule L:338 vor HiggsTAGGER"<<std::endl;
  /////////////////////////////////////////////////////////// Higgs TAGGER //////////////////////////////////////////////////////////////////////////////////
 
  std::unique_ptr< std::vector<TopJet> > higgsjets_all (new std::vector<TopJet> (*event.topjets));

  //cleanen der topjets um sie in einen neuen Vektor zu speichern
  higgstag_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);
  ak4_cleaner->process(event);
  std::unique_ptr< std::vector<TopJet>>  topjets_higgstag(new std::vector<TopJet> (*event.topjets));

  // zurueck speichern der all topjets
  event.topjets->clear();
  event.topjets->reserve(higgsjets_all->size());
  for(const auto & j : *higgsjets_all) event.topjets->push_back(j); 
  sort_by_pt<TopJet>(*event.topjets);
  //handle auf toptag Topjets
  event.set(h_higgstag,*topjets_higgstag );
 //  //cleanen der AK8 cleanen von tags
 //  higgstag_uncleaner->process(event);
 // std::unique_ptr< std::vector<TopJet> >unhiggstag (new std::vector<TopJet> (*event.topjets));
 // event.topjets->clear();
 //  event.topjets->reserve(higgsjets_all->size());
 //  for(const auto & j : *unhiggstag) event.topjets->push_back(j);
 //  for(unsigned int i=1; i<topjets_higgstag->size(); i++) event.topjets->push_back(topjets_higgstag->at(i));
 //  sort_by_pt<TopJet>(*event.topjets);

  // bool b_higgs =false;
  // bool b_zw =false;
  // bool b_top =false;

  // for(const auto j: *event.topjets){
  //   LorentzVector subjet_sum;
  //   for (const auto s : j.subjets()) {
  //     subjet_sum += s.v4();
  //   }
  //   double  mHiggs_rec=subjet_sum.M();
  
  //   if(mHiggs_rec < 150 && mHiggs_rec > 100 ) b_higgs = true;
  //   if(mHiggs_rec < 115 && mHiggs_rec > 60 ) b_zw = true;
  //   if(mHiggs_rec < 240 && mHiggs_rec > 150 )  b_top = true;
  // }
  // if(b_higgs)topjet_higgs_masswindow_h->fill(event);
  // if(b_zw)topjet_zw_masswindow_h->fill(event);
  // if(b_top)topjet_top_masswindow_h->fill(event);


  //Higgstag
   const bool pass_higgstag = higgstag_sel->passes(event);
   // const bool pass_higgstag =false;
   //  input_zwtag_h->fill(event);
  input_higgstag_h->fill(event);
  if(!pass_higgstag) return false;
  if(pass_higgstag){
    output_higgstag_h->fill(event);
    topjet_higgstag_h->fill(event);
    eff_higgstag_h->fill(event);
    jet_higgstag_h->fill(event);
    muon_higgstag_h->fill(event);
    event_higgstag_h->fill(event);
  }


 ///////////////////////////////////////////////////////// HIGGSTAGGER || ZW TAGGER //////////////////////////////////////////////////////////////////////////////////
  // bool pass_zwtag=false;
  // if(!pass_higgstag){
  //   pass_zwtag = ZWtag_sel->passes(event);
  //   if(!pass_zwtag) return false;
    
  //   topjet_zwtag_h->fill(event);
  //   eff_zwtag_h->fill(event);
  //   jet_zwtag_h->fill(event);
  //   muon_zwtag_h->fill(event);
  //   event_zwtag_h->fill(event);
  // }

  if(berror)  std::cout<<"SelectionModule L:338 vor ZWTAGGER"<<std::endl;  
    ///////////////////////////////////////////////////////// ZW TAGGER //////////////////////////////////////////////////////////////////////////////////

    // std::unique_ptr< std::vector<TopJet> > ZWjets_all(new std::vector<TopJet> (*event.topjets));
  //   //cleanen der topjets um sie in einen neuen Vektor zu speichern
  //   ZWtag_cleaner->process(event);
  //   sort_by_pt<TopJet>(*event.topjets);
  // if(pass_zwtag){
  //   ak4_cleaner->process(event);

  // }
  //   std::unique_ptr< std::vector<TopJet>>  topjets_ZWtag(new std::vector<TopJet> (*event.topjets));
  //   // zurueck speichern der all topjets
  //   event.topjets->clear();
  //   event.topjets->reserve(ZWjets_all->size());
  //   for(const auto & j : *ZWjets_all) event.topjets->push_back(j); 
  //   sort_by_pt<TopJet>(*event.topjets);
  //   //handle auf toptag Topjets
  //   event.set(h_ZWtag,*topjets_ZWtag );

  //   if(pass_zwtag){
//  //cleanen der AK8 cleanen von tags
//   ZWtag_uncleaner->process(event);
//  std::unique_ptr< std::vector<TopJet> >unZWtag (new std::vector<TopJet> (*event.topjets));
//  event.topjets->clear();
//   event.topjets->reserve(ZWjets_all->size());
//   for(const auto & j : *unZWtag) event.topjets->push_back(j);
//   for(unsigned int i=1; i<topjets_ZWtag->size(); i++) event.topjets->push_back(topjets_ZWtag->at(i));
//   sort_by_pt<TopJet>(*event.topjets);
// }
    ///////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////////
    // if(pass_zwtag){  
    //   output_zwtag_h->fill(event);
    // }

    // //output_combined_h->fill(event); 


  if(berror)  std::cout<<"SelectionModule L:338 vor TopTAGGER"<<std::endl;  
  ///////////////////////////////////////////////////////// TOP TAGGER //////////////////////////////////////////////////////////////////////////////////
 const bool pass_toptag = toptag_sel->passes(event);

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

 // //cleanen der AK8 cleanen von tags
 //  toptag_uncleaner->process(event);
 // std::unique_ptr< std::vector<TopJet> >untoptag (new std::vector<TopJet> (*event.topjets));
 // event.topjets->clear();
 //  event.topjets->reserve(topjets_all->size());
 //  for(const auto & j : *untoptag) event.topjets->push_back(j);
 //  for(unsigned int i=1; i<topjets_toptag->size(); i++) event.topjets->push_back(topjets_toptag->at(i));
 //  sort_by_pt<TopJet>(*event.topjets);


  //Toptag
 
   output_toptag_h->fill(event);
  topjet_toptag_h->fill(event);
  eff_toptag_h->fill(event);
  jet_toptag_h->fill(event);
  muon_toptag_h->fill(event);
  event_toptag_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_topjet_toptag_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_toptag_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_toptag_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_toptag_h->fill(event);
 
 ///////////////////////////////////////////////////////////////// Tagger Hist ///////////////////////////////////////////////////////
  topjet_tagger_h->fill(event);
  eff_tagger_h->fill(event);
  jet_tagger_h->fill(event);
  muon_tagger_h->fill(event);
  event_tagger_h->fill(event);

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

 if(berror)   std::cout<<"SelectionModule L:303 vor KinReco"<<std::endl;
  //////////////////////////////////////////////////////////// KIN RECO///////////////////////////////////////////////////////////////////////////////////
  
  reco_primlep->process(event);
  if(!(ZprimeTotTPrime_reco->process(event)))return false;
  ZprimeTotTPrime_chi->process(event); 
  // save only the chi2-best ttbar hypothesis in output sub-ntuple
  std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps = event.get(h_ZprimeTotTPrime_hyps);
  const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis(hyps, "Chi2");
  if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
  const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj(*hyp);
  // std::cout << "Test in Selection Module " << hyp->HZW_subjets().size() << " v4 des ersten " <<hyp->HZW_subjets().at(0).v4()<<std::endl;
  // chi2min_reco_h->fill(event);
  // topjet_reco_h->fill(event);
  // eff_reco_h->fill(event);
  // jet_reco_h->fill(event);
  // muon_reco_h->fill(event);
  // event_reco_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_topjet_reco_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_topjet_reco_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_topjet_reco_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_topjet_reco_h->fill(event);

  // if(pass_higgstag && pass_toptag) higgs_top_chi2min_reco_h->fill(event);
  // if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_reco_h->fill(event);
  // if(pass_zwtag && pass_toptag)zw_top_chi2min_reco_h->fill(event);
  // if(pass_zwtag && !pass_toptag)zw_notop_chi2min_reco_h->fill(event);

  ///////////////////////////////////// AK8 cleaning from used AK4 //////////////////////////////
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
 sort_by_pt<TopJet>(*AK8Jets);
 ////put cleaned AK8 jets in event.topjet
 event.topjets->clear();
  event.topjets->reserve(AK8Jets->size());
 for(const auto & j : *AK8Jets) event.topjets->push_back(j); 
 sort_by_pt<TopJet>(*event.topjets);
if(berror) std::cout<<"SelectionModule L:858 Size topjets Collection "<<event.topjets->size() <<std::endl;

 /////////////////AK8 cleaning end ////////////////


  hyps.clear();
  hyps.push_back(hyp_obj);
 if(berror)   std::cout << "SelectionModule L:715 vor Chi2Cut"<<std::endl;
  //////////////////////////////////////////////////////////  CHi2 Cut  ////////////////////////////////////////////////////////

  chi2min_reco_h->fill(event);
  if(pass_toptag)top_chi2min_reco_h->fill(event);
  if(!pass_toptag)notop_chi2min_reco_h->fill(event);

  bool pass_chi2cut = chi2cut_sel->passes(event);
  if(!pass_chi2cut) return false;

  bool pass_met;
  if(channel_==elec)pass_met = met_sel->passes(event);
  if(channel_==elec) if(!pass_met) return false;

  chi2min_chi2cut_h->fill(event);
  topjet_chi2cut_h->fill(event);
  eff_chi2cut_h->fill(event);
  jet_chi2cut_h->fill(event);
  muon_chi2cut_h->fill(event);
  event_chi2cut_h->fill(event);

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

// bool pass_btag1 = btag1_sel->passes(event);
// if(pass_btag1){
//    topjet_btag50_h->fill(event);
//    eff_btag50_h->fill(event);
//    jet_btag50_h->fill(event);
//    muon_btag50_h->fill(event);
//    event_btag50_h->fill(event); 
//    chi2min_btag50_h->fill(event);

//  // if(pass_higgstag && pass_toptag) higgs_top_chi2min_btag50_h->fill(event);
//  //  if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_btag50_h->fill(event);
//  //  if(pass_zwtag && pass_toptag)zw_top_chi2min_btag50_h->fill(event);
//  //  if(pass_zwtag && !pass_toptag)zw_notop_chi2min_btag50_h->fill(event);

//  }
 

/////////////////////////////////////////////////////////  btag  ////////////////////////////////////////////////////////
bool pass_btag1 = btag1_sel->passes(event);
bool pass_btag0 = btag0_sel->passes(event);

 if(pass_btag1){
   if(isMC)  btagwAK4->process(event);
   topjet_btag1_h->fill(event);
   eff_btag1_h->fill(event);
   jet_btag1_h->fill(event);
   muon_btag1_h->fill(event);
   event_btag1_h->fill(event); 
   chi2min_btag1_h->fill(event);
   if(isMC) pdf_btag1_h->fill(event);
 }

 if(pass_btag0){
 if(isMC)   btagwAK4->process(event);
   topjet_btag0_h->fill(event);
   eff_btag0_h->fill(event);
   jet_btag0_h->fill(event);
   muon_btag0_h->fill(event);
   event_btag0_h->fill(event);
   chi2min_btag0_h->fill(event);
   if(isMC) pdf_btag0_h->fill(event);
 }

 // //////////////////////////////////////////////////////////  Eventnumber  ////////////////////////////////////////////////////////
 //  // Outputfile of Eventnumber, Runnumber, Lumiblock to compare to Z'>tT'(Wb)
 //  if(filename.find("Data_C")!=std::string::npos){
 //  std::fstream g;
 //  g.open("eventnumber_side2_C.txt", ios::out);
 //  g << event.run<<" "<<event.luminosityBlock<<" "<<event.event  << std::endl;
 //  g.close();
 // }

 // if(filename.find("Data_D")!=std::string::npos){
 //  std::fstream f;
 //  f.open("eventnumber_side2_D.txt", ios::out);
 //  f << event.run<<" "<<event.luminosityBlock<<" "<<event.event  << std::endl;
 //  f.close();
 // }
 // ////////
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
 lumi_h->fill(event);

   
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeSidebandModule_Side2)
