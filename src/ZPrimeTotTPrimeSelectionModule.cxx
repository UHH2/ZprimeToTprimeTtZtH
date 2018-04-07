#include <iostream>
#include <memory>
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TROOT.h"
#include "TGraphAsymmErrors.h"

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
#include "UHH2/common/include/TopPtReweight.h"

#include <fstream>

#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/BackgroundGenSelection.h>

#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstruction.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHypothesisHists.h" 
#include "UHH2/ZprimeToTprimeTtZtH/include/SR_PDFHists.h"

using namespace uhh2examples;
using namespace uhh2;

class ZPrimeTotTPrimeSelectionModule : public uhh2::AnalysisModule {

public:
  explicit ZPrimeTotTPrimeSelectionModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  enum lepton { muon, elec };
  lepton channel_;
  std::unique_ptr<AnalysisModule> printer;
  // cleaners
  std::unique_ptr<MuonCleaner>     muo_cleaner;
  std::unique_ptr<ElectronCleaner> ele_cleaner1;
  std::unique_ptr<ElectronCleaner> photon_cleaner;
  std::unique_ptr<ElectronCleaner> ele_cleaner2;
  std::unique_ptr<TopJetCleaner> topjet_PTcleaner;

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

  std::unique_ptr<JetResolutionSmearer> jet_resolution_smearer;

  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCleaner>      jet_cleaner2;
  std::unique_ptr<JetCleaner>      jet_cleaner1;
  std::unique_ptr<TopJetCleaner>   topjet_cleaner;
  std::unique_ptr<JetCleaner>      topjet_IDcleaner;
  std::unique_ptr<TopJetCleaner>   topjet_mass_cleaner;
  std::unique_ptr<JetCleaner>  ak4_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;

  std::unique_ptr<TopJetCleaner>   higgstag_one_btag_cleaner;
  std::unique_ptr<TopJetCleaner>   toptag_cleaner;
  std::unique_ptr<TopJetCleaner>   higgstag_cleaner;
  std::unique_ptr<TopJetCleaner>   ZWtag_cleaner;
  std::unique_ptr<JetCleaner>   btag_cleaner_medium;

  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;

  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::unique_ptr<uhh2::AnalysisModule> muonscale;
  std::unique_ptr<MuonTrkWeights> muontracker_sf;
  std::unique_ptr<MistagRateSF> higgstagger_mistag_sf;
  std::unique_ptr<IsMistag> ismistag_higgs;
  std::unique_ptr<MistagRateSF> higgs1btagger_mistag_sf;
  std::unique_ptr<IsMistag> ismistag_higgs1b;
  std::unique_ptr<MistagRateSF> zwtagger_mistag_sf;
  std::unique_ptr<IsMistag> ismistag_zw;
  std::unique_ptr<uhh2::AnalysisModule> elecid2;
  std::unique_ptr<uhh2::AnalysisModule> elecid1;
  std::unique_ptr<uhh2::AnalysisModule> elecreco;
  std::unique_ptr<uhh2::AnalysisModule> electrigger;
  std::unique_ptr<uhh2::AnalysisModule> btagwAK4;
  std::unique_ptr<uhh2::AnalysisModule> btagwAK8;
  std::unique_ptr<uhh2::AnalysisModule> triggerscale;
  std::unique_ptr<uhh2::AnalysisModule> triggerscale_data;
  std::unique_ptr<uhh2::AnalysisModule> muonIDscale_data;
  std::unique_ptr<uhh2::AnalysisModule> topptreweighting_all;
  // std::unique_ptr<uhh2::AnalysisModule> eff_zw:

  //Selections
  std::unique_ptr<AndSelection>  metfilters_selection;
  std::unique_ptr<uhh2::Selection> triggerMu50_sel;
  std::unique_ptr<uhh2::Selection> triggerTrkMu50_sel;
  std::unique_ptr<uhh2::Selection> triggerElec_sel1;
  std::unique_ptr<uhh2::Selection> triggerElec_sel;
  std::unique_ptr<uhh2::Selection> triggerPhoton_sel;
  std::unique_ptr<uhh2::Selection> lumi_sel;
  std::unique_ptr<uhh2::AndSelection> lep1_sel;     //  exactly one lepton(muon) 
  std::unique_ptr<uhh2::Selection> twodcut_sel;     // pt 40 rel 0.4
  std::unique_ptr<uhh2::Selection> triangc_sel;
  std::unique_ptr<uhh2::Selection> drele_sel;
  std::unique_ptr<uhh2::Selection> topjet_sel;      // one AK8 pT>250
  std::unique_ptr<uhh2::Selection> chi2cut_sel;     // chi2min <50
  std::unique_ptr<uhh2::Selection> met_sel;  
  // std::unique_ptr<uhh2::Selection> ht_sel;
  std::unique_ptr<uhh2::Selection> btag1_sel;       //for limits
  // std::unique_ptr<uhh2::Selection> btag0_sel;

  std::unique_ptr<uhh2::Selection> jet2_sel; //for elec channel
  std::unique_ptr<uhh2::Selection> jet1_sel; //for elec channel

  //TOP TAGGER
  std::unique_ptr<Selection> toptag_sel; 
  uhh2::Event::Handle< std::vector<TopJet> > h_toptag;
  std::unique_ptr<Hists>topjet_top_masswindow_h;

  //HIGGS TAGGER
  std::unique_ptr<Selection> higgstag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_higgstag;
  std::unique_ptr<Hists>topjet_higgs_masswindow_h;

  //HIGGS TAGGER one btag
  std::unique_ptr<Selection> higgstag_one_btag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_higgstag_one_btag;

  //Z/W TAGGER
  std::unique_ptr<Selection> ZWtag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_ZWtag;
  std::unique_ptr<Hists>topjet_zw_masswindow_h;

  //Handles
  uhh2::Event::Handle< std::vector<TopJet> > h_tagged;
  uhh2::Event::Handle< std::vector<TopJet> > h_other;
  uhh2::Event::Handle< std::vector<Jet> > h_used;
  uhh2::Event::Handle< std::vector<Jet> > h_btag_medium;
  uhh2::Event::Handle< std::vector<Jet> > h_btag_used;

  // Reconstruction ZPrime for Signal
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrimeprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_reco;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_chi;
 

  // Reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;
std::unique_ptr<uhh2::AnalysisModule> backgroundprod;
  std::unique_ptr<uhh2::Selection> genwpt_sel;
  std::unique_ptr<uhh2::Selection> genwpt_cut;  

  ////////////////////HISTS///////////////////////////////////////
  std::unique_ptr<Hists> output_preselection_eff_h;

  //PU reweighting
  std::unique_ptr<Hists> pu_event_h;
  std::unique_ptr<Hists> pu_after_event_h;

  //Trigger
  std::unique_ptr<Hists> trigger_electron_h;
  std::unique_ptr<Hists> trigger_electron_after_h;
  std::unique_ptr<Hists> trigger_muon_h;
  std::unique_ptr<Hists> trigger_muon_after_h;
  std::unique_ptr<Hists> before_trigger_eff_h;
  std::unique_ptr<Hists> after_trigger_eff_h;

  //Input Hists
  std::unique_ptr<Hists> input_event_h;
  std::unique_ptr<Hists> input_elep_h;
  std::unique_ptr<Hists> input_mulep_h;
  std::unique_ptr<Hists> input_eff_h;
  std::unique_ptr<Hists> input_jet_h; 
  std::unique_ptr<Hists> input_topjet_h;
  std::unique_ptr<Hists> input_lumi_h;
  std::unique_ptr<Hists> input_gen_h;

  //Hists lep1
  std::unique_ptr<Hists> event_lep1_h;
  std::unique_ptr<Hists> topjet_lep1_h;
  std::unique_ptr<Hists> jet_lep1_h;
  std::unique_ptr<Hists> muon_lep1_h;
  std::unique_ptr<Hists> electron_lep1_h;
  std::unique_ptr<Hists> eff_lep1_h;
  std::unique_ptr<Hists> lumi_lep1_h;
 
  //Hist Topjet2
  std::unique_ptr<Hists> eff_topjet2_h;
  std::unique_ptr<Hists> jet_topjet2_h;
  std::unique_ptr<Hists> muon_topjet2_h;
  std::unique_ptr<Hists> electron_topjet2_h;
  std::unique_ptr<Hists> event_topjet2_h;
  std::unique_ptr<Hists> topjet_topjet2_h;
  std::unique_ptr<Hists> lumi_topjet2_h;

  //2DCut Hists
  std::unique_ptr<Hists> eff_twodcut_h;
  std::unique_ptr<Hists> jet_twodcut_h;
  std::unique_ptr<Hists> topjet_twodcut_h;
  std::unique_ptr<Hists> muon_twodcut_h;
  std::unique_ptr<Hists> electron_twodcut_h;
  std::unique_ptr<Hists> event_twodcut_h;
  std::unique_ptr<Hists> lumi_twodcut_h;
  std::unique_ptr<Hists> muon_pdf_twodcut_h;
  std::unique_ptr<Hists> electron_pdf_twodcut_h;
  std::unique_ptr<Hists> jet_pdf_twodcut_h;
  std::unique_ptr<Hists> topjet_pdf_twodcut_h;
  std::unique_ptr<Hists> event_pdf_twodcut_h;

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
  std::unique_ptr<Hists> electron_higgstag_h;
  std::unique_ptr<Hists> event_higgstag_h;

  //ZWTagg Hists
  std::unique_ptr<TopJetHists> input_zwtag_h;
  std::unique_ptr<TopJetHists> output_zwtag_h;

  std::unique_ptr<Hists> topjet_zwtag_h;
  std::unique_ptr<Hists> eff_zwtag_h;
  std::unique_ptr<Hists> jet_zwtag_h;
  std::unique_ptr<Hists> muon_zwtag_h;
  std::unique_ptr<Hists> electron_zwtag_h;
  std::unique_ptr<Hists> event_zwtag_h;


  //Toptag Hists
  std::unique_ptr<TopJetHists> input_toptag_h;
  std::unique_ptr<TopJetHists> output_toptag_h;

  std::unique_ptr<Hists> topjet_toptag_h;
  std::unique_ptr<Hists> eff_toptag_h;
  std::unique_ptr<Hists> muon_toptag_h;
  std::unique_ptr<Hists> electron_toptag_h;
  std::unique_ptr<Hists> jet_toptag_h;
  std::unique_ptr<Hists> event_toptag_h;

  std::unique_ptr<TopJetHists> higgs_top_topjet_toptag_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_toptag_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_toptag_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_toptag_h;

  //Hist Tagger
  std::unique_ptr<Hists> topjet_tagger_comb_h;
  std::unique_ptr<Hists> topjet_tagger_tagged_h;
  std::unique_ptr<Hists> topjet_tagger_other_h;
  std::unique_ptr<Hists> topjet_mistagger_h;
  std::unique_ptr<Hists> eff_tagger_h;
  std::unique_ptr<Hists> jet_tagger_h;
  std::unique_ptr<Hists> muon_tagger_h;
  std::unique_ptr<Hists> electron_tagger_h;
  std::unique_ptr<Hists> event_tagger_h;
  std::unique_ptr<Hists> lumi_tagger_h;

  std::unique_ptr<TopJetHists> higgs_top_topjet_tagger_h;
  std::unique_ptr<TopJetHists> one_higgs_top_topjet_tagger_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_tagger_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_tagger_h;
  std::unique_ptr<TopJetHists> one_higgs_notop_topjet_tagger_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_tagger_h;

  std::unique_ptr<JetHists> higgs_top_jet_tagger_h;
  std::unique_ptr<JetHists> zw_top_jet_tagger_h;
  std::unique_ptr<JetHists> higgs_notop_jet_tagger_h;
  std::unique_ptr<JetHists> zw_notop_jet_tagger_h;

  std::unique_ptr<MuonHists> higgs_top_muon_tagger_h;
  std::unique_ptr<MuonHists> zw_top_muon_tagger_h;
  std::unique_ptr<MuonHists> higgs_notop_muon_tagger_h;
  std::unique_ptr<MuonHists> zw_notop_muon_tagger_h;

  std::unique_ptr<ElectronHists> higgs_top_electron_tagger_h;
  std::unique_ptr<ElectronHists> zw_top_electron_tagger_h;
  std::unique_ptr<ElectronHists> higgs_notop_electron_tagger_h;
  std::unique_ptr<ElectronHists> zw_notop_electron_tagger_h;

  std::unique_ptr<EventHists> higgs_top_event_tagger_h;
  std::unique_ptr<EventHists> zw_top_event_tagger_h;
  std::unique_ptr<EventHists> higgs_notop_event_tagger_h;
  std::unique_ptr<EventHists> zw_notop_event_tagger_h;


  // //Hist Reco
  // std::unique_ptr<Hists> chi2min_reco_h;
  std::unique_ptr<Hists> eff_reco_h;
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
  std::unique_ptr<Hists> jet_chi2cut_comb_h;
  std::unique_ptr<Hists> jet_chi2cut_used_h;
  std::unique_ptr<Hists> muon_chi2cut_h;
  std::unique_ptr<Hists> electron_chi2cut_h;
  std::unique_ptr<Hists> topjet_chi2cut_comb_h;
  std::unique_ptr<Hists> topjet_chi2cut_tagged_h;
  std::unique_ptr<Hists> higgs_topjet_chi2cut_tagged_h;
  std::unique_ptr<Hists> one_higgs_topjet_chi2cut_tagged_h;
  std::unique_ptr<Hists> zw_topjet_chi2cut_tagged_h;
  std::unique_ptr<Hists> topjet_chi2cut_other_h;
  std::unique_ptr<Hists> event_chi2cut_h;
  std::unique_ptr<Hists> lumi_chi2cut_h;

  std::unique_ptr<TopJetHists> higgs_top_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> zw_top_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> higgs_notop_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> zw_notop_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> zw_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> higgs_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> one_higgs_top_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> one_higgs_notop_topjet_chi2cut_h;
  std::unique_ptr<TopJetHists> one_higgs_topjet_chi2cut_h;

  std::unique_ptr<Hists> higgs_top_chi2min_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_top_chi2min_chi2cut_h;
  std::unique_ptr<Hists> zw_top_chi2min_chi2cut_h;
  std::unique_ptr<Hists> higgs_notop_chi2min_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_notop_chi2min_chi2cut_h;
  std::unique_ptr<Hists> zw_notop_chi2min_chi2cut_h;

  std::unique_ptr<Hists> higgs_top_pdf_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_top_pdf_chi2cut_h;
  std::unique_ptr<Hists> zw_top_pdf_chi2cut_h;
  std::unique_ptr<Hists> higgs_notop_pdf_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_notop_pdf_chi2cut_h;
  std::unique_ptr<Hists> zw_notop_pdf_chi2cut_h;

  std::unique_ptr<Hists> zw_chi2min_chi2cut_h;
  std::unique_ptr<Hists> higgs_chi2min_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_chi2min_chi2cut_h;

  std::unique_ptr<JetHists> higgs_top_jet_chi2cut_h;
  std::unique_ptr<JetHists> zw_top_jet_chi2cut_h;
  std::unique_ptr<JetHists> higgs_notop_jet_chi2cut_h;
  std::unique_ptr<JetHists> zw_notop_jet_chi2cut_h;
  std::unique_ptr<JetHists> higgs_jet_chi2cut_h;
  std::unique_ptr<JetHists> zw_jet_chi2cut_h;
  std::unique_ptr<JetHists> one_higgs_top_jet_chi2cut_h;
  std::unique_ptr<JetHists> one_higgs_notop_jet_chi2cut_h;
  std::unique_ptr<JetHists> one_higgs_jet_chi2cut_h;
  
  std::unique_ptr<MuonHists> higgs_top_muon_chi2cut_h;
  std::unique_ptr<MuonHists> zw_top_muon_chi2cut_h;
  std::unique_ptr<MuonHists> higgs_notop_muon_chi2cut_h;
  std::unique_ptr<MuonHists> zw_notop_muon_chi2cut_h;
  std::unique_ptr<MuonHists> higgs_muon_chi2cut_h;
  std::unique_ptr<MuonHists> zw_muon_chi2cut_h;
  std::unique_ptr<MuonHists> one_higgs_top_muon_chi2cut_h;
  std::unique_ptr<MuonHists> one_higgs_notop_muon_chi2cut_h;
  std::unique_ptr<MuonHists> one_higgs_muon_chi2cut_h;

  std::unique_ptr<ElectronHists> higgs_top_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> zw_top_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> higgs_notop_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> zw_notop_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> higgs_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> zw_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> one_higgs_top_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> one_higgs_notop_electron_chi2cut_h;
  std::unique_ptr<ElectronHists> one_higgs_electron_chi2cut_h;

  std::unique_ptr<EventHists> higgs_top_event_chi2cut_h;
  std::unique_ptr<EventHists> zw_top_event_chi2cut_h;
  std::unique_ptr<EventHists> higgs_notop_event_chi2cut_h;
  std::unique_ptr<EventHists> zw_notop_event_chi2cut_h;
  std::unique_ptr<EventHists> higgs_event_chi2cut_h;
  std::unique_ptr<EventHists> zw_event_chi2cut_h;
  std::unique_ptr<EventHists> one_higgs_top_event_chi2cut_h;
  std::unique_ptr<EventHists> one_higgs_notop_event_chi2cut_h;
  std::unique_ptr<EventHists> one_higgs_event_chi2cut_h;


  std::unique_ptr<Hists> higgs_top_eff_chi2cut_h;
  std::unique_ptr<Hists> zw_top_eff_chi2cut_h;
  std::unique_ptr<Hists> higgs_notop_eff_chi2cut_h;
  std::unique_ptr<Hists> zw_notop_eff_chi2cut_h;
  std::unique_ptr<Hists> higgs_eff_chi2cut_h;
  std::unique_ptr<Hists> zw_eff_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_top_eff_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_notop_eff_chi2cut_h;
  std::unique_ptr<Hists> one_higgs_eff_chi2cut_h;

  //Hist MET
  std::unique_ptr<Hists> chi2min_met1_h;
  std::unique_ptr<Hists> eff_met1_h;
  std::unique_ptr<Hists> event_met1_h;
  std::unique_ptr<Hists> higgs_top_chi2min_met1_h;
  std::unique_ptr<Hists> one_higgs_top_chi2min_met1_h;
  std::unique_ptr<Hists> zw_top_chi2min_met1_h;
  std::unique_ptr<Hists> higgs_notop_chi2min_met1_h;
  std::unique_ptr<Hists> one_higgs_notop_chi2min_met1_h;
  std::unique_ptr<Hists> zw_notop_chi2min_met1_h;
  std::unique_ptr<Hists> zw_chi2min_met1_h;
  std::unique_ptr<Hists> higgs_chi2min_met1_h;
  std::unique_ptr<Hists> one_higgs_chi2min_met1_h;


  //counting signal region events
  std::unique_ptr<Hists> eff_selection_h;
  std::unique_ptr<Hists> chi2min_selection_h;
  std::unique_ptr<Hists> jet_selection_h;
  std::unique_ptr<Hists> muon_selection_h;
  std::unique_ptr<TopJetHists> topjet_selection_h;
  std::unique_ptr<Hists> event_selection_h;

  //general hists
  std::unique_ptr<Hists> btag_jet_h;
  std::unique_ptr<Hists> btag_topjet_h;
 //systematicen
  std::unique_ptr<AnalysisModule> syst_module;
  bool do_scale_variation;

  //general
  TString filename;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ZprimeTotTPrime_hyps;
  uhh2::Event::Handle<double> h_ht;
  uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
  uhh2::Event::Handle<BackgroundGen> h_background;
  uhh2::Event::Handle<TTbarGen> h_ttbargen;
  bool berror;
  TString unc_name;
  bool isMC;
  bool b_topjet_cleaner;
  TString elec_up_down;

  const int runnr_BCD = 276811;
  const int runnr_EF = 278802;
  const int runnr_G = 280385;

  TString  top_up_down;
  TString  top_up_down_miss;
  TString  higgs_up_down;
  TString  one_higgs_up_down;
  TString  zw_up_down;
  TString  zw_up_down_miss;
  TString  toppt_sys;

};


ZPrimeTotTPrimeSelectionModule::ZPrimeTotTPrimeSelectionModule(uhh2::Context& ctx){
 
  //GenParticleprinter
  printer.reset(new GenParticlesPrinter(ctx));

  //choose channel
  const std::string& channel = ctx.get("channel", "");
  if     (channel == "muon") channel_ = muon;
  else if(channel == "electron") channel_ = elec;
  else throw std::runtime_error("ZprimeSelectionModule -- undefined argument for 'channel' key in xml file (must be 'muon' or 'elec'): "+channel);
  isMC = (ctx.get("dataset_type") == "MC");
  //// Data/MC scale
  auto data_dir_path = ctx.get("data_dir_path");
  string sysAK4=ctx.get("btagging_sys", "central");
  if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx,ctx.get("puVariation"))); 
    lumiweight.reset(new MCLumiWeight(ctx));
    //   btagwAK4.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM,"jets",sysAK4,"mujets","incl","MCBtagEfficienciesAK4","_AK4","BTagCalibration")); 
    //Hier
    btagwAK8.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "topjets",ctx.get("subjetbtag_sys"),"lt","incl","MCBtagEfficienciesAK8","_AK8","BTagCalibrationAK8"));
    if(channel_==muon){
      muonscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta", 1.,"MediumID2016", true, ctx.get("Systematic_MuonID")));
      muontracker_sf.reset(new MuonTrkWeights(ctx, data_dir_path +"Tracking_EfficienciesAndSF_BCDEFGH.root",  ctx.get("Systematic_MuonTrk")));
      // triggerscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","IsoMu50_OR_IsoTkMu50_PtEtaBins", 1., "trigger", true, ctx.get("Systematic_MuonTrigger") ));
      triggerscale.reset(new MuonTriggerSF(ctx));
    }else if(channel_==elec){
      //  elecid2.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "", ctx.get("Systematic_EleID")));
      electrigger.reset(new ElectronTriggerSF(ctx));
      elecreco.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "", ctx.get("Systematic_EleReco")));
      elecid1.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_MVA90_ID.root", 1, "", ctx.get("Systematic_EleID")));

    }
    higgstagger_mistag_sf.reset(new MistagRateSF(ctx,"/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/ZprimeToTprimeTtZtH/hist/corrections/CorrectionFactors_higgs.txt",ctx.get("Systematic_higgs")));
    ismistag_higgs.reset(new IsMistag(ctx));
    higgs1btagger_mistag_sf.reset(new MistagRateSF(ctx,"/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/ZprimeToTprimeTtZtH/hist/corrections/CorrectionFactors_higgs_one.txt",ctx.get("Systematic_one_higgs")));
    ismistag_higgs1b.reset(new IsMistag(ctx));
    zwtagger_mistag_sf.reset(new MistagRateSF(ctx,"/nfs/dust/cms/user/abenecke/CMSSW_8_0_24_patch1/src/UHH2/ZprimeToTprimeTtZtH/hist/corrections/CorrectionFactors_zw.txt",ctx.get("Systematic_zw_miss")));
    ismistag_zw.reset(new IsMistag(ctx));
    //  eff_zw.reset(new ZPrimeTotTPrimeEff(ctx));
  }
  else  {   
    lumi_sel.reset(new LumiSelection(ctx)); 
    // triggerscale_data.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/abenecke/scripts/hists_weighted/trigger/Diff_MuonTrigger_EfficienciesAndSF_average_RunBtoH.root","Mu50_OR_TkMu50_PtEtaBins", 1.));
    //  muonIDscale_data.reset(new MCMuonScaleFactor(ctx,"/nfs/dust/cms/user/abenecke/scripts/hists_weighted/muonID/Diff_MuonID_EfficienciesAndSF_average_RunBtoH.root","MC_NUM_MediumID2016_DEN_genTracks_PAR_pt_eta", 1.));
  }
  
  PrimaryVertexId pvid=StandardPrimaryVertexId();
  metfilters.emplace_back(new PrimaryVertexCleaner(pvid));
  metfilters_selection.reset(new AndSelection(ctx, "metfilters"));
  metfilters_selection->add<TriggerSelection>("HBHENoiseFilter", "Flag_HBHENoiseFilter");
  metfilters_selection->add<TriggerSelection>("HBHENoiseIsoFilter", "Flag_HBHENoiseIsoFilter");
  metfilters_selection->add<TriggerSelection>("globalTightHalo2016Filter", "Flag_globalTightHalo2016Filter");
  metfilters_selection->add<TriggerSelection>("EcalDeadCellTriggerPrimitiveFilter", "Flag_EcalDeadCellTriggerPrimitiveFilter");
  if(!isMC) metfilters_selection->add<TriggerSelection>("eeBadScFilter", "Flag_eeBadScFilter");

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
    //    jetlepton_cleaner->set_drmax(.4);
    //  jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
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
 muo_cleaner.reset(new MuonCleaner    (AndId<Muon>    (PtEtaCut  (53., 2.4), MuonIDMedium_ICHEP())));
 //  ele_cleaner2.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(30., 2.4), ElectronID_Spring16_tight)));
 //  ele_cleaner1.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(125, 2.4), ElectronID_Spring16_tight_noIso)));//Arne
 ele_cleaner1.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(125, 2.4),ElectronID_MVAGeneralPurpose_Spring16_loose )));
 photon_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(250, 2.4),ElectronID_MVAGeneralPurpose_Spring16_loose )));
 // ele_cleaner1.reset(new ElectronCleaner(AndId<Electron>(PtEtaCut(125, 2.4),ElectronID_MVAHZZ_Spring16_loose )));
 topjet_mass_cleaner.reset(new TopJetCleaner(ctx,TopJetId(TopjetMassCleaner(30))));
  topjet_PTcleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(300,2.4))));

  const JetId jetID(JetPFID(JetPFID::WP_LOOSE));
  jet_IDcleaner.reset(new JetCleaner(ctx,jetID));

  jet_cleaner2.reset(new JetCleaner(ctx,15., 2.4));
  jet_cleaner1.reset(new JetCleaner(ctx,30., 2.4));
  
  topjet_cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(200., 2.4))));
 
  ak4_cleaner.reset(new JetCleaner(ctx,JetId(ZPrimeTotTPrimeAK4cleaner(1.2))));
  topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));
  const JetId bjetID_medium(CSVBTag(CSVBTag::WP_MEDIUM));
  btag_cleaner_medium.reset(new JetCleaner(ctx,bjetID_medium));

  //calculator
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));

 
 
  ////////////////////////////////////// Selections /////////////////////////////////////////////////

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

  // TOPJET 
  topjet_sel.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut( 250., 2.4))));
 
  // 2D Cut
  twodcut_sel.reset(new TwoDCut(.4, 40.));

  triangc_sel.reset(new uhh2::AndSelection(ctx));

  // if     (channel_ == elec) triangc_sel.reset(new TriangularCuts(1.5, 75.));
  // else if(channel_ == muon) triangc_sel.reset(new uhh2::AndSelection(ctx)); // always true (no triangular cuts for muon channel)

  drele_sel.reset(new ZPrimeTotTPrimeDRele(.4,-1));

  // // jet sel for elec channel
  // if(channel_==elec){
  //   jet2_sel.reset(new NJetSelection(2, -1, JetId(PtEtaCut( 50., 2.4))));
  //   jet1_sel.reset(new NJetSelection(1, -1, JetId(PtEtaCut( 200., 2.4))));
  // }
 
  // TOP TAGGER
  const TopJetId topjetID = AndId<TopJet>(Type2TopTag(150,220,Type2TopTag::MassType::groomed), Tau32(0.57));
  //  const TopJetId topjetID = AndId<TopJet>(Type2TopTag(140,220,Type2TopTag::MassType::groomed), Tau32(0.57));
  toptag_sel.reset(new NTopJetSelection(1, -1, topjetID));
  toptag_cleaner.reset(new TopJetCleaner(ctx,topjetID));

  // Higgs TAGGER
  //  const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150), Tau21(1) );
  //Hier
  //  const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150,DoubleBTag(DoubleBTag::WP_MEDIUM)), Tau21(1) );
  const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150, CSVBTag(CSVBTag::WP_MEDIUM)), Tau21(1) );
  higgstag_sel.reset(new NTopJetSelection(1, -1, higgsjetID));
  higgstag_cleaner.reset(new TopJetCleaner(ctx,higgsjetID));

  // Higgs TAGGER
  const TopJetId higgs_one_btag_jetID = AndId<TopJet>(ZPrimeTotTPrimeHiggsTag(100,150), Tau21(1) );
  // const TopJetId higgs_one_btag_jetID = AndId<TopJet>(ZPrimeTotTPrimeHiggsTag(100,150,CSVBTag(CSVBTag::WP_LOOSE)), Tau21(1) );
  higgstag_one_btag_sel.reset(new NTopJetSelection(1, -1, higgs_one_btag_jetID));
  higgstag_one_btag_cleaner.reset(new TopJetCleaner(ctx,higgs_one_btag_jetID));


  // W/Z TAGGER
  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
  ZWtag_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
  ZWtag_cleaner.reset(new TopJetCleaner(ctx,ZWjetID));


  /* KINEMATICAL RECO */
  const std::string ZprimeTotTPrime_gen_label ("zprimegen");
  const std::string Background_gen_label ("backgroundgen");
  const std::string ttbar_gen_label ("ttbargen");
  const std::string ZprimeTotTPrime_hyps_label("ZPrimeTotTPrimeReconstruction");
  const std::string ZprimeTotTPrime_selection_hyps_label("ZPrimeTotTPrimeReconstruction"); 
  const std::string ZprimeTotTPrime_chi2_label("Chi2");

  //GenParticles
  ZprimeTotTPrimeprod.reset(new ZPrimeGenProducer(ctx, ZprimeTotTPrime_gen_label, false));
  backgroundprod.reset(new BackgroundGenProducer(ctx, Background_gen_label, false));
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  reco_primlep.reset(new PrimaryLepton(ctx));
  ZprimeTotTPrime_reco.reset(new ZPrimeTotTPrimeReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrime_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ZprimeTotTPrime_hyps_label));
  h_ZprimeTotTPrime_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ZprimeTotTPrime_hyps_label);
 
  topptreweighting_all.reset(new TopPtReweight(ctx, 0.0615 , -0.0005, ttbar_gen_label,"weight_ttbar", true,1));
  
  chi2cut_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  met_sel.reset(new ZPrimeTotTPrimeMETCut(90,0));
  

  // //BTAG
  btag1_sel.reset(new NBTagSelection(1, -1));
  // btag0_sel.reset(new NBTagSelection(0, 0));

  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));

 // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "MC_WJetsAll"){ genwpt_sel.reset(new ZPrimeTotTPrimePartonW());}
  else                                          genwpt_sel.reset(new uhh2::AndSelection(ctx));
  

  //handels
  h_tagged = ctx.declare_event_output< std::vector<TopJet> > ("tagged");
  h_other = ctx.declare_event_output< std::vector<TopJet> > ("other");
  h_used = ctx.declare_event_output< std::vector<Jet> > ("used");
 
 
  ////////////////////////////////////// Hists  /////////////////////////////////////////////////
  output_preselection_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "output_preselection_eff"));
  //PU reweighting
  pu_event_h.reset(new EventHists(ctx, "pu_event"));
  pu_after_event_h.reset(new EventHists(ctx, "pu_after_event"));

  //Trigger
  trigger_electron_h.reset(new ElectronHists(ctx, "trigger_electron"));
  trigger_muon_h.reset(new MuonHists(ctx, "trigger_muon"));

  trigger_electron_after_h.reset(new ElectronHists(ctx, "trigger_electron_after"));
  trigger_muon_after_h.reset(new MuonHists(ctx, "trigger_muon_after"));

  before_trigger_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "before_trigger_eff"));
  after_trigger_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "after_trigger_eff"));

  // Hists Input
  input_event_h.reset(new EventHists(ctx, "input"));
  input_mulep_h.reset(new MuonHists(ctx, "input_muLep"));
  input_elep_h.reset(new ElectronHists(ctx, "input_eLep"));
  input_jet_h.reset(new JetHists (ctx, "input_Jet"));
  input_topjet_h.reset(new TopJetHists (ctx, "input_TopJet"));
  input_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
  input_lumi_h.reset(new LuminosityHists(ctx,"input_lumi"));
  input_gen_h.reset(new ZPrimeGenHists(ctx,"input_gen"));

  btag_jet_h.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiencyHistsAK4",CSVBTag::WP_MEDIUM,"jets" ));
  btag_topjet_h.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiencyHistsAK8",CSVBTag::WP_MEDIUM,"topjets" ));


  // Hists lep1
  topjet_lep1_h.reset(new TopJetHists(ctx, "topjet_lep1"));
  event_lep1_h.reset(new EventHists(ctx, "event_lep1"));
  muon_lep1_h.reset(new MuonHists(ctx, "muon_lep1"));
  electron_lep1_h.reset(new ElectronHists(ctx, "electron_lep1"));
  jet_lep1_h.reset(new JetHists(ctx, "jet_lep1"));
  eff_lep1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_lep1"));
  lumi_lep1_h.reset(new LuminosityHists(ctx,"lumi_lep1"));

  //Hists topjet2
  topjet_topjet2_h.reset(new TopJetHists(ctx, "topjet_topjet2"));
  eff_topjet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_topjet2"));
  jet_topjet2_h.reset(new JetHists(ctx, "jet_topjet2"));
  muon_topjet2_h.reset(new MuonHists(ctx, "muon_topjet2"));
  electron_topjet2_h.reset(new ElectronHists(ctx, "electron_topjet2"));
  event_topjet2_h.reset(new EventHists(ctx, "event_topjet2"));
  lumi_topjet2_h.reset(new LuminosityHists(ctx,"lumi_topjet2"));

  //Hist twodcut
  topjet_twodcut_h.reset(new TopJetHists(ctx, "topjet_twodcut"));
  eff_twodcut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_twodcut"));
  jet_twodcut_h.reset(new JetHists(ctx, "jet_twodcut"));
  muon_twodcut_h.reset(new MuonHists(ctx, "muon_twodcut"));
  electron_twodcut_h.reset(new ElectronHists(ctx, "electron_twodcut"));
  event_twodcut_h.reset(new EventHists(ctx, "event_twodcut"));
  lumi_twodcut_h.reset(new LuminosityHists(ctx,"lumi_twodcut"));
  

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
  electron_higgstag_h.reset(new ElectronHists(ctx, "electron_higgstag"));
  event_higgstag_h.reset(new EventHists(ctx, "event_higgstag"));
  topjet_higgstag_h.reset(new TopJetHists(ctx, "topjet_higgstag"));
 
  //Z/W TAG HISTS
  input_zwtag_h.reset(new TopJetHists(ctx, "input_zwtag"));
  output_zwtag_h.reset(new TopJetHists(ctx, "output_zwtag"));
  output_zwtag_h->set_TopJetId(ZWjetID);

  eff_zwtag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_zwtag"));
  jet_zwtag_h.reset(new JetHists(ctx, "jet_zwtag"));
  muon_zwtag_h.reset(new MuonHists(ctx, "muon_zwtag"));
  electron_zwtag_h.reset(new ElectronHists(ctx, "electron_zwtag"));
  event_zwtag_h.reset(new EventHists(ctx, "event_zwtag"));
  topjet_zwtag_h.reset(new TopJetHists(ctx, "topjet_zwtag"));

  // Hists Toptag
  input_toptag_h.reset(new TopJetHists(ctx, "input_toptag"));
  output_toptag_h.reset(new TopJetHists(ctx, "output_toptag"));
  output_toptag_h->set_TopJetId(topjetID);

  eff_toptag_h.reset(new ZPrimeTotTPrimeHists(ctx,"eff_toptag"));
  jet_toptag_h.reset(new JetHists(ctx, "jet_toptag"));
  muon_toptag_h.reset(new MuonHists(ctx, "muon_toptag"));
  electron_toptag_h.reset(new ElectronHists(ctx, "electron_toptag"));
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
  electron_tagger_h.reset(new ElectronHists(ctx, "electron_tagger"));
  event_tagger_h.reset(new EventHists(ctx, "event_tagger"));
  topjet_tagger_comb_h.reset(new TopJetHists(ctx, "topjet_tagger_comb"));
  topjet_mistagger_h.reset(new TopJetHists(ctx, "topjet_mistagger"));
  topjet_tagger_tagged_h.reset(new TopJetHists(ctx, "topjet_tagger_tagged",4,"tagged"));
  topjet_tagger_other_h.reset(new TopJetHists(ctx, "topjet_tagger_other",4,"other"));
  higgs_top_topjet_tagger_h.reset(new TopJetHists(ctx, "higgs_top_topjet_tagger"));
  one_higgs_top_topjet_tagger_h.reset(new TopJetHists(ctx, "one_higgs_top_topjet_tagger"));
  higgs_notop_topjet_tagger_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_tagger"));
  one_higgs_notop_topjet_tagger_h.reset(new TopJetHists(ctx, "one_higgs_notop_topjet_tagger"));
  zw_top_topjet_tagger_h.reset(new TopJetHists(ctx, "zw_top_topjet_tagger"));
  zw_notop_topjet_tagger_h.reset(new TopJetHists(ctx, "zw_notop_topjet_tagger"));
  lumi_tagger_h.reset(new LuminosityHists(ctx,"lumi_tagger"));

  higgs_top_jet_tagger_h.reset(new JetHists(ctx, "higgs_top_jet_tagger"));
  higgs_notop_jet_tagger_h.reset(new JetHists(ctx, "higgs_notop_jet_tagger"));
  zw_top_jet_tagger_h.reset(new JetHists(ctx, "zw_top_jet_tagger"));
  zw_notop_jet_tagger_h.reset(new JetHists(ctx, "zw_notop_jet_tagger"));

  higgs_top_muon_tagger_h.reset(new MuonHists(ctx, "higgs_top_muon_tagger"));
  higgs_notop_muon_tagger_h.reset(new MuonHists(ctx, "higgs_notop_muon_tagger"));
  zw_top_muon_tagger_h.reset(new MuonHists(ctx, "zw_top_muon_tagger"));
  zw_notop_muon_tagger_h.reset(new MuonHists(ctx, "zw_notop_muon_tagger"));

  higgs_top_electron_tagger_h.reset(new ElectronHists(ctx, "higgs_top_electron_tagger"));
  higgs_notop_electron_tagger_h.reset(new ElectronHists(ctx, "higgs_notop_electron_tagger"));
  zw_top_electron_tagger_h.reset(new ElectronHists(ctx, "zw_top_electron_tagger"));
  zw_notop_electron_tagger_h.reset(new ElectronHists(ctx, "zw_notop_electron_tagger"));

  higgs_top_event_tagger_h.reset(new EventHists(ctx, "higgs_top_event_tagger"));
  higgs_notop_event_tagger_h.reset(new EventHists(ctx, "higgs_notop_event_tagger"));
  zw_top_event_tagger_h.reset(new EventHists(ctx, "zw_top_event_tagger"));
  zw_notop_event_tagger_h.reset(new EventHists(ctx, "zw_notop_event_tagger"));
 
  // //Hists Reco
  // topjet_reco_h.reset(new TopJetHists(ctx, "topjet_reco"));
  eff_reco_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reco"));
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
  topjet_chi2cut_comb_h.reset(new TopJetHists(ctx, "topjet_chi2cut_comb"));
  topjet_chi2cut_tagged_h.reset(new TopJetHists(ctx, "topjet_chi2cut_tagged",4,"tagged"));
  higgs_topjet_chi2cut_tagged_h.reset(new TopJetHists(ctx, "higgs_topjet_chi2cut_tagged",4,"HiggsTag"));
  one_higgs_topjet_chi2cut_tagged_h.reset(new TopJetHists(ctx, "one_higgs_topjet_chi2cut_tagged",4,"HiggsTag_one_btag"));
  zw_topjet_chi2cut_tagged_h.reset(new TopJetHists(ctx, "zw_topjet_chi2cut_tagged",4,"ZWTag"));
  topjet_chi2cut_other_h.reset(new TopJetHists(ctx, "topjet_chi2cut_other",4,"other"));
  eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut"));
  jet_chi2cut_comb_h.reset(new JetHists(ctx, "jet_chi2cut_comb"));
  jet_chi2cut_used_h.reset(new JetHists(ctx, "jet_chi2cut_used",4,"used"));
  muon_chi2cut_h.reset(new MuonHists(ctx, "muon_chi2cut"));
  electron_chi2cut_h.reset(new ElectronHists(ctx, "electron_chi2cut"));
  event_chi2cut_h.reset(new EventHists(ctx, "event_chi2cut"));
  chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));
  lumi_chi2cut_h.reset(new LuminosityHists(ctx,"lumi_chi2cut"));

  higgs_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_top_topjet_chi2cut"));
  higgs_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_notop_topjet_chi2cut"));
  zw_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_top_topjet_chi2cut"));
  zw_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_notop_topjet_chi2cut"));
  higgs_topjet_chi2cut_h.reset(new TopJetHists(ctx, "higgs_topjet_chi2cut"));
  zw_topjet_chi2cut_h.reset(new TopJetHists(ctx, "zw_topjet_chi2cut"));
  one_higgs_top_topjet_chi2cut_h.reset(new TopJetHists(ctx, "one_higgs_top_topjet_chi2cut"));
  one_higgs_notop_topjet_chi2cut_h.reset(new TopJetHists(ctx, "one_higgs_notop_topjet_chi2cut"));
  one_higgs_topjet_chi2cut_h.reset(new TopJetHists(ctx, "one_higgs_topjet_chi2cut"));

  higgs_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  one_higgs_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "one_higgs_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  higgs_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  one_higgs_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "one_higgs_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_top_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_notop_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  if(isMC){
    muon_pdf_twodcut_h.reset(new SR_PDFHists(ctx,"muon_pdf_twodcut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    electron_pdf_twodcut_h.reset(new SR_PDFHists(ctx,"electron_pdf_twodcut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    jet_pdf_twodcut_h.reset(new SR_PDFHists(ctx,"jet_pdf_twodcut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    topjet_pdf_twodcut_h.reset(new SR_PDFHists(ctx,"topjet_pdf_twodcut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    event_pdf_twodcut_h.reset(new SR_PDFHists(ctx,"event_pdf" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));

    higgs_top_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"higgs_top_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    one_higgs_top_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"one_higgs_top_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    higgs_notop_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"higgs_notop_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    one_higgs_notop_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"one_higgs_notop_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    zw_top_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"zw_top_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
    zw_notop_pdf_chi2cut_h.reset(new SR_PDFHists(ctx,"zw_notop_pdf_chi2cut" , ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label,true));
  }
 
  one_higgs_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "one_higgs_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  higgs_chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_chi2min_chi2cut",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));


  higgs_top_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_top_jet_chi2cut"));
  higgs_notop_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_notop_jet_chi2cut"));
  zw_top_jet_chi2cut_h.reset(new JetHists(ctx, "zw_top_jet_chi2cut"));
  zw_notop_jet_chi2cut_h.reset(new JetHists(ctx, "zw_notop_jet_chi2cut"));
  higgs_jet_chi2cut_h.reset(new JetHists(ctx, "higgs_jet_chi2cut"));
  zw_jet_chi2cut_h.reset(new JetHists(ctx, "zw_jet_chi2cut"));
  one_higgs_top_jet_chi2cut_h.reset(new JetHists(ctx, "one_higgs_top_jet_chi2cut"));
  one_higgs_notop_jet_chi2cut_h.reset(new JetHists(ctx, "one_higgs_notop_jet_chi2cut"));
  one_higgs_jet_chi2cut_h.reset(new JetHists(ctx, "one_higgs_jet_chi2cut"));

  higgs_top_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_top_muon_chi2cut"));
  higgs_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_notop_muon_chi2cut"));
  zw_top_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_top_muon_chi2cut"));
  zw_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_notop_muon_chi2cut"));
  higgs_muon_chi2cut_h.reset(new MuonHists(ctx, "higgs_muon_chi2cut"));
  zw_muon_chi2cut_h.reset(new MuonHists(ctx, "zw_muon_chi2cut"));
  one_higgs_top_muon_chi2cut_h.reset(new MuonHists(ctx, "one_higgs_top_muon_chi2cut"));
  one_higgs_notop_muon_chi2cut_h.reset(new MuonHists(ctx, "one_higgs_notop_muon_chi2cut"));
  one_higgs_muon_chi2cut_h.reset(new MuonHists(ctx, "one_higgs_muon_chi2cut"));

  higgs_top_electron_chi2cut_h.reset(new ElectronHists(ctx, "higgs_top_electron_chi2cut"));
  higgs_notop_electron_chi2cut_h.reset(new ElectronHists(ctx, "higgs_notop_electron_chi2cut"));
  zw_top_electron_chi2cut_h.reset(new ElectronHists(ctx, "zw_top_electron_chi2cut"));
  zw_notop_electron_chi2cut_h.reset(new ElectronHists(ctx, "zw_notop_electron_chi2cut"));
  higgs_electron_chi2cut_h.reset(new ElectronHists(ctx, "higgs_electron_chi2cut"));
  zw_electron_chi2cut_h.reset(new ElectronHists(ctx, "zw_electron_chi2cut"));
  one_higgs_top_electron_chi2cut_h.reset(new ElectronHists(ctx, "one_higgs_top_electron_chi2cut"));
  one_higgs_notop_electron_chi2cut_h.reset(new ElectronHists(ctx, "one_higgs_notop_electron_chi2cut"));
  one_higgs_electron_chi2cut_h.reset(new ElectronHists(ctx, "one_higgs_electron_chi2cut"));

  higgs_top_event_chi2cut_h.reset(new EventHists(ctx, "higgs_top_event_chi2cut"));
  higgs_notop_event_chi2cut_h.reset(new EventHists(ctx, "higgs_notop_event_chi2cut"));
  zw_top_event_chi2cut_h.reset(new EventHists(ctx, "zw_top_event_chi2cut"));
  zw_notop_event_chi2cut_h.reset(new EventHists(ctx, "zw_notop_event_chi2cut"));
  higgs_event_chi2cut_h.reset(new EventHists(ctx, "higgs_event_chi2cut"));
  zw_event_chi2cut_h.reset(new EventHists(ctx, "zw_event_chi2cut"));
  one_higgs_top_event_chi2cut_h.reset(new EventHists(ctx, "one_higgs_top_event_chi2cut"));
  one_higgs_notop_event_chi2cut_h.reset(new EventHists(ctx, "one_higgs_notop_event_chi2cut"));
  one_higgs_event_chi2cut_h.reset(new EventHists(ctx, "one_higgs_event_chi2cut"));


  higgs_top_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "higgs_top_eff_chi2cut"));
  higgs_notop_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "higgs_notop_eff_chi2cut"));
  zw_top_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "zw_top_eff_chi2cut"));
  zw_notop_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "zw_notop_eff_chi2cut"));
  higgs_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "higgs_eff_chi2cut"));
  zw_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "zw_eff_chi2cut"));
  one_higgs_top_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "one_higgs_top_eff_chi2cut"));
  one_higgs_notop_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "one_higgs_notop_eff_chi2cut"));
  one_higgs_eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "one_higgs_eff_chi2cut"));

  eff_selection_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_selection"));
  topjet_selection_h.reset(new TopJetHists(ctx, "topjet_selection"));
  jet_selection_h.reset(new JetHists(ctx, "jet_selection"));
  muon_selection_h.reset(new MuonHists(ctx, "muon_selection"));
  event_selection_h.reset(new EventHists(ctx, "event_selection"));
  chi2min_selection_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_selection",ZprimeTotTPrime_selection_hyps_label,ZprimeTotTPrime_chi2_label ));

  //Higgstag || ZWTag
  topjet_higgs_masswindow_h.reset(new TopJetHists(ctx, "masswindow_higgs"));
  topjet_zw_masswindow_h.reset(new TopJetHists(ctx, "masswindow_zw"));
  topjet_top_masswindow_h.reset(new TopJetHists(ctx, "masswindow_top"));


 
  //general
  filename =  ctx.get("dataset_version");

  // if(filename.Contains("WJets")) genwpt_cut.reset(new ZPrimeTotTPrimePartonWCut(filename));
  // else genwpt_cut.reset(new uhh2::AndSelection(ctx));

  h_toptag = ctx.declare_event_output< std::vector<TopJet> > ("TopTag");
  h_higgstag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag");
  h_higgstag_one_btag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag_one_btag");
  h_ZWtag = ctx.declare_event_output< std::vector<TopJet> > ("ZWTag");
  h_ht = ctx.get_handle<double>("HT");
  h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
  h_background = ctx.get_handle<BackgroundGen>("backgroundgen");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
 
  //Trigger
  if(channel_==muon){
    triggerMu50_sel.reset(new TriggerSelection("HLT_Mu50_v*"));
    triggerTrkMu50_sel.reset(new TriggerSelection("HLT_TkMu50_v*"));
  }
  //triggerElec_sel1.reset(new TriggerSelection("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_v*"));//ZPrime > ttbar
  triggerElec_sel1.reset(new TriggerSelection("HLT_Ele115_CaloIdVT_GsfTrkIdT_v*"));//Arne
  triggerPhoton_sel.reset(new TriggerSelection("HLT_Photon175_v*"));
  //  triggerElec_sel2.reset(new TriggerSelection("HLT_Ele27_WPTight_Gsf_v*"));//Iso
  

  
  h_btag_medium = ctx.declare_event_output< std::vector<Jet> > ("BTag_medium");
  h_btag_used = ctx.declare_event_output< std::vector<Jet> > ("BTag_used");
  
  //systematicen
  syst_module.reset(new MCScaleVariation(ctx));

  do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");

  berror=false;
  b_topjet_cleaner = (ctx.get("b_topjet_cleaner")=="true");
  unc_name = ctx.get("unc_name");
  top_up_down = ctx.get("Systematic_toptagger");
  top_up_down_miss = ctx.get("Systematic_toptagger_miss");
  higgs_up_down =  ctx.get("Systematic_higgs");
  one_higgs_up_down =  ctx.get("Systematic_one_higgs");
  zw_up_down =  ctx.get("Systematic_zw");
  zw_up_down_miss =  ctx.get("Systematic_zw_miss");
  elec_up_down = ctx.get("Systematic_EleTrigger");
  toppt_sys = ctx.get("Systematic_toppt");
}

bool ZPrimeTotTPrimeSelectionModule::process(uhh2::Event& event){
  
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
    if (toppt_sys=="none")topptreweighting_all->process(event);
  }else if (!filename.Contains("MC_ZPrime") && !filename.Contains("Data")){
    backgroundprod->process(event);
  }
  
   if(filename.Contains("WJetsAll")){
     genwpt_sel->passes(event);
     if(!genwpt_sel->passes(event)) return false;
  }
 
   // printer->process(event);
  //ZPrime Genrator Level
  if(filename.Contains("MC_ZPrime")){
    ZprimeTotTPrimeprod->process(event); 
  }

  
  // //genptWparton selection to remove steps in final histograms
  // if(filename.Contains("WJets")){
  //   bool pass_genwptcut=genwpt_cut->passes(event);
  //   if(!pass_genwptcut) return false;
  // }

  if(berror)   std::cout<<"SelectionModule L:232 vor Input Histogrammen"<<std::endl;
  
  /////////////////////////////////////////////////////////// Common Modules   ///////////////////////////////////////////////////////////////////////////////


  //common Modules
  /* luminosity sections from CMS silver-JSON file */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  output_preselection_eff_h->fill(event);

  pu_event_h->fill(event);
  /* pileup SF */
  if(!event.isRealData){ 
    pileup_SF->process(event);
    lumiweight->process(event);
  }
  pu_after_event_h->fill(event);
  
  // if(event.isRealData&& event.run>=runnr_EF){ 
  //   triggerscale_data->process(event);
  //   muonIDscale_data->process(event);
  // }



  //
  ///correctors
  if(isMC){
    jetlepton_cleaner->process(event);
    jet_corrector->process(event);
    topjet_corrector->process(event);
    subjet_corrector->process(event);
    jet_corrector->correct_met(event);
    // jet_resolution_smearer->process(event);
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


  // systematicen
  if(do_scale_variation) syst_module->process(event); 


  //////////////////////////////////////////////////////////////////Trigger//////////////////////////////////////////////////////////////////////////////////////
  before_trigger_eff_h->fill(event);
  trigger_muon_h->fill(event);
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
  trigger_muon_after_h->fill(event);
  
  trigger_electron_h->fill(event);
  
  bool pass_elecTrigger = triggerElec_sel1 ->passes(event);
  bool pass_photonTrigger = triggerPhoton_sel->passes(event);
  if(channel_==elec){
    
    if(!(pass_elecTrigger||pass_photonTrigger)) return false;
    //    if(filename.Contains("photon")) if (pass_elecTrigger && pass_photonTrigger) return false;
    if(filename.Contains("photon") && event.electrons->at(0).pt()<250) return false;
    if(filename.Contains("rereco") && event.electrons->at(0).pt()>250) return false;
    trigger_electron_after_h->fill(event);

    if(pass_elecTrigger) ele_cleaner1->process(event);
    else if(pass_photonTrigger) photon_cleaner->process(event);
    sort_by_pt<Electron>(*event.electrons);
    if(isMC)  elecid1->process(event);
   
  }
 
  after_trigger_eff_h->fill(event);

 // OBJ CLEANING
  if(channel_ == muon) muo_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);

  

  jet_IDcleaner->process(event);
  jet_cleaner2->process(event); 
  
  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

  /////////////////////////////////////////////////////////// Input Histogramme ///////////////////////////////////////////////////////////////////////////////
 
  input_eff_h ->fill(event);
  input_event_h->fill(event);
  input_topjet_h->fill(event);
  input_jet_h->fill(event);
  input_mulep_h->fill(event);
  input_elep_h->fill(event);
  input_lumi_h->fill(event);
  input_gen_h->fill(event);


  if(berror)  std::cout<<"SelectionModule L:268 vor LeptonSelection"<<std::endl;
  /////////////////////////////////////////////////////////// LEPTON selection ///////////////////////////////////////////////////////////////////////////////
  const bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  if(channel_==elec )if(event.electrons->at(0).eta()<1.566 && event.electrons->at(0).eta()>1.4442) return false;

  muon_lep1_h->fill(event);
  electron_lep1_h->fill(event);
  eff_lep1_h ->fill(event);
  event_lep1_h->fill(event);
  topjet_lep1_h->fill(event);
  jet_lep1_h->fill(event);
  lumi_lep1_h->fill(event);
  

  if(!event.isRealData){
    if(channel_ == muon) muonscale->process(event);
    if(channel_ == muon) muontracker_sf->process(event);
    if(channel_ == muon) triggerscale->process(event);
    if(channel_ == elec) elecreco->process(event);
    if(channel_ == elec) electrigger->process(event);
  }

  topjetlepton_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  if(berror)  std::cout<<"SelectionModule L:268 vor TopJetSelection"<<std::endl;
  ////////////////////////////////////////////////////////// TopJET selection//////////////////////////////////////////////////////////////////////////////////
  topjet_mass_cleaner->process(event);
  topjet_cleaner->process(event);
  if(b_topjet_cleaner)  topjet_PTcleaner->process(event);
  const bool pass_topjet = topjet_sel->passes(event);
  if(!pass_topjet) return false;

  topjet_topjet2_h->fill(event);
  eff_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);
  muon_topjet2_h->fill(event);
  electron_topjet2_h->fill(event);
  event_topjet2_h->fill(event);
  lumi_topjet2_h->fill(event);
 
  if(berror)   std::cout<<"SelectionModule L:294 vor TwoDCut"<<std::endl;  
  ////////////////////////////////////////////////////////// TwoDCut //////////////////////////////////////////////////////////////////////////////////

  sort_by_pt<Jet>(*event.jets);
  //  cout<<"numer of electrons "<<event.electrons->size()<<endl;
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  jet_cleaner1->process(event);

  // if(channel_==elec){
  //   if(pass_elecTrigger1){
  //     /* 2nd AK4 jet selection */
  //     const bool pass_jet2 = jet2_sel->passes(event);
  //     if(!pass_jet2) return false;
      
  //     /* 1st AK4 jet selection */
  //     const bool pass_jet1 = jet1_sel->passes(event);
  //     if(!pass_jet1) return false;
  //   }
  // } //ZPrime >ttbar
  if(channel_==elec){
    bool pass_met = met_sel->passes(event);
    if(!pass_met) return false;
  }




  topjet_twodcut_h->fill(event);
  jet_twodcut_h->fill(event);
  eff_twodcut_h->fill(event);
  muon_twodcut_h->fill(event);
  electron_twodcut_h->fill(event);
  lumi_twodcut_h->fill(event);
  event_twodcut_h->fill(event);
  if(isMC && channel_==muon) muon_pdf_twodcut_h->fill(event);
  if(isMC && channel_ == elec) electron_pdf_twodcut_h->fill(event);
  if(isMC) jet_pdf_twodcut_h->fill(event);
  if(isMC) topjet_pdf_twodcut_h->fill(event);
  if(isMC) event_pdf_twodcut_h->fill(event);



  
  btag_jet_h->fill(event);
  btag_topjet_h->fill(event);

  lumi_twodcut_h->fill(event);

  if(berror)   std::cout<<"SelectionModule L:294 vor Triangcut"<<std::endl;  
  ///////////////////////////////// /////////////////////////// Triangcut //////////////////////////////////////////////////////////////////////////////////
  // const bool pass_triangc = triangc_sel->passes(event);
  // if(!pass_triangc) return false;
  topjet_triangcut_h->fill(event);
  jet_triangcut_h->fill(event);
  eff_triangcut_h->fill(event);
  muon_triangcut_h->fill(event);
  electron_triangcut_h->fill(event);
  lumi_triangcut_h->fill(event);

  // const bool pass_drele = drele_sel->passes(event);
  // if(!pass_drele) return false;
  event_triangcut_h->fill(event);

  if(berror)  std::cout<<"SelectionModule L:338 vor HiggsTAGGER"<<std::endl;
  /////////////////////////////////////////////////////////// Higgs TAGGER //////////////////////////////////////////////////////////////////////////////////
  if(isMC) btagwAK8->process(event);
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

  bool b_higgs =false;
  bool b_zw =false;
  bool b_top =false;

  for(const auto j: *event.topjets){
    LorentzVector subjet_sum;
    for (const auto s : j.subjets()) {
      subjet_sum += s.v4();
    }
    double  mHiggs_rec=subjet_sum.M();
  
    if(mHiggs_rec < 150 && mHiggs_rec > 100 ) b_higgs = true;
    if(mHiggs_rec < 115 && mHiggs_rec > 60 ) b_zw = true;
    if(mHiggs_rec < 240 && mHiggs_rec > 150 )  b_top = true;
  }
  if(b_higgs)topjet_higgs_masswindow_h->fill(event);
  if(b_zw)topjet_zw_masswindow_h->fill(event);
  if(b_top)topjet_top_masswindow_h->fill(event);


  //Higgstag
  const bool pass_higgstag = higgstag_sel->passes(event);

  input_zwtag_h->fill(event);
  input_higgstag_h->fill(event);
  if(pass_higgstag){
    output_higgstag_h->fill(event);
    topjet_higgstag_h->fill(event);
    eff_higgstag_h->fill(event);
    jet_higgstag_h->fill(event);
    muon_higgstag_h->fill(event);
    electron_higgstag_h->fill(event);
    event_higgstag_h->fill(event);
  }

  if(berror) std::cout<<"vor Z/W tag selection"<<std::endl;
  ///////////////////////////////////////////////////////// HIGGSTAGGER || ZW TAGGER //////////////////////////////////////////////////////////////////////////////////
  bool pass_zwtag=false;
  bool pass_higgs_one_btag=false;
  if(!pass_higgstag){
    pass_higgs_one_btag = higgstag_one_btag_sel->passes(event);
    if(!pass_higgs_one_btag){
      pass_zwtag = ZWtag_sel->passes(event);
      if(!pass_zwtag)return false;
    

      
      topjet_zwtag_h->fill(event);
      eff_zwtag_h->fill(event);
      jet_zwtag_h->fill(event);
      muon_zwtag_h->fill(event);
      electron_zwtag_h->fill(event);
      event_zwtag_h->fill(event);
    }
  }

  if(berror)  std::cout<<"SelectionModule L:338 vor ZWTAGGER"<<std::endl;  
  ///////////////////////////////////////////////////////// ZW TAGGER //////////////////////////////////////////////////////////////////////////////////

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


  ///////////////////////////////////////////////////////// //////////////////////////////////////////////////////////////////////////////////
  if(pass_zwtag){  
    output_zwtag_h->fill(event);
  }


  if(berror)  std::cout<<"SelectionModule L:338 vor Higgstag one btag"<<std::endl;  
  ///////////////////////////////////////////////////////// HIGGS TAGGER with one subjet btag //////////////////////////////////////////////////////////////////////////////////

  std::unique_ptr< std::vector<TopJet> > higgs_one_btag_all(new std::vector<TopJet> (*event.topjets));
  //cleanen der topjets um sie in einen neuen Vektor zu speichern

  higgstag_one_btag_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);
  if(pass_higgs_one_btag){
    ak4_cleaner->process(event);

  }
  std::unique_ptr< std::vector<TopJet>>  topjets_higgstag_one_btag(new std::vector<TopJet> (*event.topjets));
  // zurueck speichern der all topjets
  event.topjets->clear();
  event.topjets->reserve( higgs_one_btag_all->size());
  for(const auto & j : * higgs_one_btag_all) event.topjets->push_back(j); 
  sort_by_pt<TopJet>(*event.topjets);
  //handle auf toptag Topjets
  event.set(h_higgstag_one_btag,*topjets_higgstag_one_btag );
 


  if(berror)  std::cout<<"SelectionModule L:338 vor TopTAGGER"<<std::endl;  
  ///////////////////////////////////////////////////////// TOP TAGGER //////////////////////////////////////////////////////////////////////////////////
  // bool pass_toptag = toptag_sel->passes(event);

  ///////// TEST new toptag
  //clean topjet collection from Higgstagged jet with priority
  std::unique_ptr< std::vector<TopJet> > all (new std::vector<TopJet> (*event.topjets));
 
  TopJet tagged_jet;
  if(pass_higgstag) tagged_jet = topjets_higgstag->at(0);
  else if(pass_zwtag) tagged_jet = topjets_ZWtag->at(0);
  else if(pass_higgs_one_btag) tagged_jet = topjets_higgstag_one_btag->at(0);

  event.topjets->clear();
  for(TopJet topjet:*all){
    if(!(deltaR(topjet,tagged_jet)<0.01)) event.topjets->push_back(topjet);
  }

  
  // require at least one top tag
  bool pass_toptag = toptag_sel->passes(event);


 
 
  ///////// TEST ENDE

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

  event.topjets->clear(); 
  for(const auto & j : *all) event.topjets->push_back(j);
  //Toptag
 
  output_toptag_h->fill(event);
  topjet_toptag_h->fill(event);
  eff_toptag_h->fill(event);
  jet_toptag_h->fill(event);
  muon_toptag_h->fill(event);
  electron_toptag_h->fill(event);
  event_toptag_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_toptag_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_toptag_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_topjet_toptag_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_topjet_toptag_h->fill(event);

  
//////////////////////////////////////////////////////////  dived in tagged and not tagged  ////////////////////////////////////////////////////////
 std::vector<TopJet>* AllTopJets(new std::vector<TopJet> (*event.topjets));
 std::vector<TopJet>* tagged(new std::vector<TopJet> (*event.topjets));
 std::vector<TopJet>* untagged(new std::vector<TopJet> (*event.topjets));
 tagged-> clear();
 tagged->reserve(event.topjets->size());
 untagged-> clear();
 untagged->reserve(event.topjets->size());
 event.topjets->clear();
 event.topjets->reserve(AllTopJets->size());
 for(const auto & topjet: *AllTopJets){
   event.topjets->push_back(topjet);
   bool b_Htag =  higgstag_sel->passes(event);
   bool b_Htag_one_btag =  higgstag_one_btag_sel->passes(event);
   bool b_Ztag =  ZWtag_sel->passes(event);
   
   bool b_tagged = b_Htag || b_Ztag || b_Htag_one_btag;
   if(b_tagged) tagged->push_back(topjet);
   else untagged->push_back(topjet);

   event.topjets->clear();
   event.topjets->reserve(AllTopJets->size());

 }  
 event.set(h_tagged, *tagged);
 event.set(h_other, *untagged);
 for(const auto & j : *AllTopJets) event.topjets->push_back(j); 


  ///////////////////////////////////////////////////////////////// Tagger Hist ///////////////////////////////////////////////////////
  topjet_tagger_comb_h->fill(event);
  topjet_tagger_tagged_h->fill(event);
  topjet_tagger_other_h->fill(event);
  eff_tagger_h->fill(event);
  jet_tagger_h->fill(event);
  muon_tagger_h->fill(event);
  electron_tagger_h->fill(event);
  event_tagger_h->fill(event);
  lumi_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_tagger_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_topjet_tagger_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag)one_higgs_notop_topjet_tagger_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_topjet_tagger_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_topjet_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_jet_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_jet_tagger_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_jet_tagger_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_jet_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_muon_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_muon_tagger_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_muon_tagger_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_muon_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_electron_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_electron_tagger_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_electron_tagger_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_electron_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_event_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_event_tagger_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_event_tagger_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_event_tagger_h->fill(event);

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
  eff_reco_h->fill(event);
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
  sort_by_pt<Jet>(*usedak4);
  event.set(h_used, *usedak4);

 std::vector<Jet>* ak4jets(new std::vector<Jet> (*event.jets));
 std::vector<Jet>* btag_used(new std::vector<Jet> (*event.jets));
 btag_used->clear();
 btag_used ->reserve(ak4jets->size());


  for (auto j : *usedak4){
    event.jets->clear();
    event.jets->reserve(usedak4->size());
    event.jets ->push_back(j);
    bool pass =  btag1_sel->passes(event);
      if(pass) btag_used->push_back(j);
  }
  event.set(h_btag_used, *btag_used);
  for(const auto & j : *ak4jets) event.jets->push_back(j); 

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



  // lumi_h->fill(event);
////////////////////////////////////////////////////////// Btag handle //////////////////////////////////////////////////////////////////////////////////
  //// Set handle on btagged jets
  // save all jets before cleaning
  std::unique_ptr< std::vector<Jet> > jets_all (new std::vector<Jet> (*event.jets));
  //cleaning of jets to get all btagged jets
  btag_cleaner_medium->process(event); 
  sort_by_pt<Jet>(*event.jets);
  //save all btagged jets
  std::unique_ptr< std::vector<Jet>>  jets_btag(new std::vector<Jet> (*event.jets));
  // put original jets back
  event.jets->clear();
  event.jets->reserve(jets_all->size());
  for(const auto & j : *jets_all) event.jets->push_back(j); 
  sort_by_pt<Jet>(*event.jets);
  //handle on btagged jets
  event.set(h_btag_medium,*jets_btag );
  //// END Set handle on btagged jets
 


  // //////////////////////////////////////////////////////////  Eventnumber  ////////////////////////////////////////////////////////
  // // Outputfile of Eventnumber, Runnumber, Lumiblock to compare to Z'>tT'(Wb)
  // if(filename.find("Data_C")!=std::string::npos){
  //   std::fstream g;
  //   g.open("eventnumber_signal_C.txt", ios::out);
  //   g << event.run<<" "<<event.luminosityBlock<<" "<<event.event  << std::endl;
  //   g.close();
  // }

  // if(filename.find("Data_D")!=std::string::npos){
  //   std::fstream f;
  //   f.open("eventnumber_signal_D.txt", ios::out);
  //   f << event.run<<" "<<event.luminosityBlock<<" "<<event.event  << std::endl;
  //   f.close();
  // }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
  //////////////////////////////////////////////////////////  reweighting Eff, mistagrate, sideband  ////////////////////////////////////////////////////////
  
  //SF tagger

  if(isMC){
    /////////////////////////////////   H2b
    if(pass_higgstag){
      bool pass_is_higgs;
      if(filename.Contains("TPrimeToHT")){
	pass_is_higgs = ismistag_higgs ->process(event,"higgs");
	if(pass_is_higgs) higgstagger_mistag_sf->process(event,"higgs");
      }else if(filename.Contains("TTbar")){
	double SF =1.019;
	if(higgs_up_down=="up")SF +=0.1825;
	if(higgs_up_down=="down")SF -=0.1825;
	event.weight *= SF; 
      }else{
	topjet_mistagger_h->fill(event);
	pass_is_higgs=true;
	if(pass_is_higgs) higgstagger_mistag_sf->process(event,"higgs");
      }
      
    }  //higgs mistagrate

    /////////////////////////////////   H1b
    if(pass_higgs_one_btag){
      bool pass_is_higgs_one;
      if(filename.Contains("TPrimeToHT")){	
	pass_is_higgs_one = ismistag_higgs1b ->process(event,"higgs_one");
	if(pass_is_higgs_one) higgs1btagger_mistag_sf->process(event,"higgs_one");
      }else if(filename.Contains("TTbar")){
	double SF =0.9966;
	if(one_higgs_up_down=="up")SF +=0.03;
	if(one_higgs_up_down=="down")SF -=0.03;
	event.weight *= SF; 
      }else{
	topjet_mistagger_h->fill(event);
	pass_is_higgs_one=true;
	if(pass_is_higgs_one) higgs1btagger_mistag_sf->process(event,"higgs_one");
      }
    }
      

   /////////////////////////////////   ZW
    //z/w eff und mistagrate
    if(pass_zwtag){
      bool pass_is_zw;
      pass_is_zw = ismistag_zw ->process(event,"zw");
      if(pass_is_zw) zwtagger_mistag_sf->process(event,"zw");
      else {
	double pT = topjets_ZWtag ->at(0).pt();
	double SF=0.911866;
	if(zw_up_down=="up")SF +=sqrt(pow(0.078,2)+pow(0.01,2) + pow(0.041 * log(pT/200),2));
	if(zw_up_down=="down")SF -=sqrt(pow(0.074,2)+pow(0.01,2) + pow(0.041 * log(pT/200),2));
	event.weight *= SF; 
      }
    }

 /////////////////////////////////   Top
    if(pass_toptag){
      double SF =  1.06;
      if(top_up_down == "up") SF+=0.07;
      if(top_up_down == "down") SF-=0.04;

      double deltaR_Q1=10;
      double deltaR_Q2=10;
      double deltaR_Q3=10;

      if(filename.Contains("TTbar")){
    	auto & zprimegen = event.get(h_ttbargen);
    	if(zprimegen.IsSemiLeptonicDecay()){
    	  deltaR_Q1 = deltaR(topjets_toptag->at(0),zprimegen.Q1().v4());
    	  deltaR_Q2 = deltaR(topjets_toptag->at(0),zprimegen.Q2().v4());
    	  deltaR_Q3 = deltaR(topjets_toptag->at(0),zprimegen.BHad().v4());
    	}
      }else if(filename.Contains("MC_ZPrime")){
    	auto & zprimegen = event.get(h_zprimegen);
    	deltaR_Q1 = deltaR(topjets_toptag->at(0),zprimegen.WHadQ1().v4());
    	deltaR_Q2 = deltaR(topjets_toptag->at(0),zprimegen.WHadQ2().v4());
    	deltaR_Q3 = deltaR(topjets_toptag->at(0),zprimegen.BHad().v4());
      }//filename

      if(deltaR_Q1<=0.8 && deltaR_Q2<=0.8 && deltaR_Q3<=0.8) event.weight = event.weight * SF;
      else{
	SF=0.9445;
	if(top_up_down_miss =="up")SF+=0.01746;
	if(top_up_down_miss =="down")SF-=0.01746;
	event.weight = event.weight * SF;
      }

    }//pass_toptag
  }//isMC
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  chi2min_chi2cut_h->fill(event);
  topjet_chi2cut_comb_h->fill(event);
  topjet_chi2cut_tagged_h->fill(event);
  if(pass_higgstag ) higgs_topjet_chi2cut_tagged_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_topjet_chi2cut_tagged_h->fill(event);
  if(pass_zwtag )zw_topjet_chi2cut_tagged_h->fill(event);
  topjet_chi2cut_other_h->fill(event);
  eff_chi2cut_h->fill(event);
  jet_chi2cut_comb_h->fill(event);
  jet_chi2cut_used_h->fill(event);
  muon_chi2cut_h->fill(event);
  electron_chi2cut_h->fill(event);
  event_chi2cut_h->fill(event);
  lumi_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_topjet_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_topjet_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_topjet_chi2cut_h->fill(event);
  if(pass_zwtag)zw_topjet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_topjet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_topjet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_topjet_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_chi2min_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_chi2min_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_chi2min_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_chi2min_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_chi2min_chi2cut_h->fill(event);


  if(isMC){
  if(pass_higgstag && pass_toptag) higgs_top_pdf_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_pdf_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_pdf_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_pdf_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_pdf_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_pdf_chi2cut_h->fill(event);
  }


  if(pass_higgstag) higgs_chi2min_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_chi2min_chi2cut_h->fill(event);
  if(pass_zwtag )zw_chi2min_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_jet_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_jet_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_jet_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_jet_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_jet_chi2cut_h->fill(event);
  if(pass_zwtag )zw_jet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_jet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_jet_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_jet_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_muon_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_muon_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_muon_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_muon_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_muon_chi2cut_h->fill(event);
  if(pass_zwtag )zw_muon_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_muon_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_muon_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_muon_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_electron_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_electron_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_electron_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_electron_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_electron_chi2cut_h->fill(event);
  if(pass_zwtag )zw_electron_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_electron_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_electron_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_electron_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_event_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_event_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_event_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_event_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_event_chi2cut_h->fill(event);
  if(pass_zwtag )zw_event_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_event_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_event_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_event_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_eff_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_eff_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_eff_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_eff_chi2cut_h->fill(event);
  if(pass_higgstag ) higgs_eff_chi2cut_h->fill(event);
  if(pass_zwtag )zw_eff_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && pass_toptag) one_higgs_top_eff_chi2cut_h->fill(event);
  if(pass_higgs_one_btag && !pass_toptag) one_higgs_notop_eff_chi2cut_h->fill(event);
  if(pass_higgs_one_btag ) one_higgs_eff_chi2cut_h->fill(event);




  if(berror) std::cout<<"Am Ende"<<std::endl;
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeSelectionModule)
