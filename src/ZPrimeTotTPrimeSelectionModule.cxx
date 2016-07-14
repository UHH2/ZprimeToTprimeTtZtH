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
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h>


#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstruction.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHypothesisHists.h" 

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
  std::unique_ptr<ElectronCleaner> ele_cleaner;
  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCleaner>      jet_cleaner2;
  std::unique_ptr<JetCleaner>      jet_cleaner1;
  std::unique_ptr<JetCleaner>      topjet_IDcleaner;
  std::unique_ptr<TopJetCleaner>   topjet_cleaner;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
  std::unique_ptr<JetCleaner>  ak4_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner;

  std::unique_ptr<TopJetCleaner>   toptag_cleaner;
  std::unique_ptr<TopJetCleaner>   higgstag_cleaner;
  std::unique_ptr<TopJetCleaner>   ZWtag_cleaner;

  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;

  //correctors
  std::unique_ptr<SubJetCorrector> subjetcorrector;
 
  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::unique_ptr<uhh2::AnalysisModule> muonscale;
  std::unique_ptr<uhh2::AnalysisModule> btagwAK4;

  //Selections
  std::unique_ptr<AndSelection>  metfilters_selection;
  std::unique_ptr<uhh2::Selection> trigger_sel;
  std::unique_ptr<uhh2::Selection> lumi_sel;
  std::unique_ptr<uhh2::AndSelection> lep1_sel;     //  exactly one lepton(muon) 
  std::unique_ptr<uhh2::Selection> twodcut_sel;     // pt 40 rel 0.4
  std::unique_ptr<uhh2::Selection> topjet_sel;      // one AK8 pT>250
  std::unique_ptr<uhh2::Selection> chi2cut_sel;     // chi2min <50
  std::unique_ptr<uhh2::Selection> btag1_sel;       //for limits
  std::unique_ptr<uhh2::Selection> btag0_sel;

  //TOP TAGGER
  std::unique_ptr<Selection> toptag_sel; 
  uhh2::Event::Handle< std::vector<TopJet> > h_toptag;

  //HIGGS TAGGER
  std::unique_ptr<Selection> higgstag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_higgstag;

  //Z/W TAGGER
  std::unique_ptr<Selection> ZWtag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_ZWtag;

  // Reconstruction ZPrime for Signal
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrimeprod;
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_reco;
  std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_chi;
 

  // Reconstruction TTBar for Background
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

  //Hist Btag1
  std::unique_ptr<Hists> eff_btag1_h;
  std::unique_ptr<Hists> jet_btag1_h;
  std::unique_ptr<Hists> muon_btag1_h;
  std::unique_ptr<Hists> event_btag1_h;
  std::unique_ptr<Hists> topjet_btag1_h;
  std::unique_ptr<Hists> chi2min_btag1_h;

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

  std::unique_ptr<Hists> higgs_top_chi2min_btag0_h;
  std::unique_ptr<Hists> zw_top_chi2min_btag0_h;
  std::unique_ptr<Hists> higgs_notop_chi2min_btag0_h;
  std::unique_ptr<Hists> zw_notop_chi2min_btag0_h;

  //counting signal region events
  std::unique_ptr<Hists> eff_selection_h;
  std::unique_ptr<Hists> chi2min_selection_h;
  std::unique_ptr<Hists> jet_selection_h;
  std::unique_ptr<Hists> muon_selection_h;
  std::unique_ptr<TopJetHists> topjet_selection_h;
  std::unique_ptr<Hists> event_selection_h;

  //general hists
  std::unique_ptr<Hists> lumi_h;
  std::unique_ptr<Hists> btag_jet_h;
  std::unique_ptr<Hists> btag_topjet_h;


  //general
  std::string filename;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ZprimeTotTPrime_hyps;
  uhh2::Event::Handle<double> h_ht;
  uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
  bool berror;
};


ZPrimeTotTPrimeSelectionModule::ZPrimeTotTPrimeSelectionModule(uhh2::Context& ctx){
  //GenParticleprinter
  printer.reset(new GenParticlesPrinter(ctx));

  //choose channel
  const std::string& channel = ctx.get("channel", "");
  if     (channel == "muon") channel_ = muon;
  else if(channel == "elec") channel_ = elec;
  else throw std::runtime_error("ZprimeSelectionModule -- undefined argument for 'channel' key in xml file (must be 'muon' or 'elec'): "+channel);

  const bool isMC = (ctx.get("dataset_type") == "MC");
  //// Data/MC scale
  auto data_dir_path = ctx.get("data_dir_path");

  if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx));
    btagwAK4.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets")); 
    muonscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonID_Z_RunCD_Reco76X_Feb15.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1", 1.));
  }
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

  //JEC
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
  topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));


  //calculator
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));

  //correctors
  if(isMC) subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_MC));
  else subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_DATA));

  
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
 
  // TOP TAGGER
  const TopJetId topjetID = AndId<TopJet>(Type2TopTag(150,240,Type2TopTag::MassType::groomed), Tau32(0.6));
  toptag_sel.reset(new NTopJetSelection(1, -1, topjetID));
  toptag_cleaner.reset(new TopJetCleaner(ctx,topjetID));

  // Higgs TAGGER
  const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150), Tau21(1) );
  higgstag_sel.reset(new NTopJetSelection(1, -1, higgsjetID));
  higgstag_cleaner.reset(new TopJetCleaner(ctx,higgsjetID));

  // W/Z TAGGER
  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
  ZWtag_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
  ZWtag_cleaner.reset(new TopJetCleaner(ctx,ZWjetID));

  /* KINEMATICAL RECO */
  const std::string ZprimeTotTPrime_gen_label ("zprimegen");
  const std::string ttbar_gen_label ("ttbargen");
  const std::string ZprimeTotTPrime_hyps_label("ZPrimeTotTPrimeReconstruction");
  const std::string ZprimeTotTPrime_selection_hyps_label("ZPrimeTotTPrimeReconstruction"); 
  const std::string ZprimeTotTPrime_chi2_label("Chi2");

  //GenParticles
  ZprimeTotTPrimeprod.reset(new ZPrimeGenProducer(ctx, ZprimeTotTPrime_gen_label, false));
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  reco_primlep.reset(new PrimaryLepton(ctx));
  ZprimeTotTPrime_reco.reset(new ZPrimeTotTPrimeReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_hyps_label));
  ZprimeTotTPrime_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ZprimeTotTPrime_hyps_label));
  h_ZprimeTotTPrime_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ZprimeTotTPrime_hyps_label);
 
  chi2cut_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  //BTAG
  btag1_sel.reset(new NBTagSelection(1, -1));
  btag0_sel.reset(new NBTagSelection(0, 0));

  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else                                                    genmttbar_sel.reset(new uhh2::AndSelection(ctx));


  ////////////////////////////////////// Hists  /////////////////////////////////////////////////
  // Hists Input
  input_event_h.reset(new EventHists(ctx, "input"));
  input_lep_h.reset(new MuonHists(ctx, "input_Lep"));
  input_jet_h.reset(new JetHists (ctx, "input_Jet"));
  input_topjet_h.reset(new TopJetHists (ctx, "input_TopJet"));
  input_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
  lumi_h.reset(new LuminosityHists(ctx,"lumi"));
  btag_jet_h.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiencyHistsAK4",CSVBTag::WP_MEDIUM,"jets" ));
  btag_topjet_h.reset(new BTagMCEfficiencyHists(ctx,"BTagMCEfficiencyHistsAK8",CSVBTag::WP_MEDIUM,"topjets" ));


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
  jet_selection_h.reset(new JetHists(ctx, "jet_selection"));
  muon_selection_h.reset(new MuonHists(ctx, "muon_selection"));
  event_selection_h.reset(new EventHists(ctx, "event_selection"));
  chi2min_selection_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_selection",ZprimeTotTPrime_selection_hyps_label,ZprimeTotTPrime_chi2_label ));
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

  //Hists btag0
  topjet_btag0_h.reset(new TopJetHists(ctx, "topjet_btag0"));
  eff_btag0_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag0"));
  jet_btag0_h.reset(new JetHists(ctx, "jet_btag0"));
  muon_btag0_h.reset(new MuonHists(ctx, "muon_btag0"));
  event_btag0_h.reset(new EventHists(ctx, "event_btag0"));
  chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

  higgs_top_chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  higgs_notop_chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_top_chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
  zw_notop_chi2min_btag0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_btag0",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));


  //Higgstag || ZWTag
  // topjet_higgs_masswindow_h.reset(new TopJetHists(ctx, "masswindow_higgs"));
  // topjet_zw_masswindow_h.reset(new TopJetHists(ctx, "masswindow_zw"));
  // topjet_top_masswindow_h.reset(new TopJetHists(ctx, "masswindow_top"));


  //general
  filename =  ctx.get("dataset_version");
  h_toptag = ctx.declare_event_output< std::vector<TopJet> > ("TopTag");
  h_higgstag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag");
  h_ZWtag = ctx.declare_event_output< std::vector<TopJet> > ("ZWTag");
  h_ht = ctx.get_handle<double>("HT");
  h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");

  //Trigger
  const std::string& trigger = ctx.get("trigger", "NULL");
  if(trigger != "NULL") trigger_sel = make_unique<TriggerSelection>(trigger);
  else                  trigger_sel = make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*");

  berror=false;
  
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
  if(filename  == "TTbarAll"){
    ttgenprod->process(event);
    if(!genmttbar_sel->passes(event)) return false;
  }
  

  //ZPrime Genrator Level
  if(filename.find("MC_ZPrime")!=std::string::npos){
    ZprimeTotTPrimeprod->process(event); 
  }

  //printer->process(event);

  if(berror)   std::cout<<"SelectionModule L:232 vor Input Histogrammen"<<std::endl;
  
  /////////////////////////////////////////////////////////// Common Modules   ///////////////////////////////////////////////////////////////////////////////


  //common Modules
  /* luminosity sections from CMS silver-JSON file */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* pileup SF */
  if(!event.isRealData){ pileup_SF->process(event);lumiweight->process(event);btagwAK4->process(event);muonscale->process(event);}
  ////

  //correctors
  subjetcorrector->process(event);

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



  if(berror)  std::cout<<"SelectionModule L:268 vor LeptonSelection"<<std::endl;
  /////////////////////////////////////////////////////////// LEPTON selection ///////////////////////////////////////////////////////////////////////////////
  const bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  muon_lep1_h->fill(event);
  eff_lep1_h ->fill(event);
  event_lep1_h->fill(event);
  topjet_lep1_h->fill(event);
  jet_lep1_h->fill(event);

  


  topjetlepton_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

  if(berror)  std::cout<<"SelectionModule L:268 vor TopJetSelection"<<std::endl;
  ////////////////////////////////////////////////////////// TopJET selection//////////////////////////////////////////////////////////////////////////////////

  const bool pass_topjet = topjet_sel->passes(event);
  if(!pass_topjet) return false;

  topjet_topjet2_h->fill(event);
  eff_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);
  muon_topjet2_h->fill(event);
  event_topjet2_h->fill(event);

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

  btag_jet_h->fill(event);
  btag_topjet_h->fill(event);
 
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
  input_zwtag_h->fill(event);
  input_higgstag_h->fill(event);
  if(pass_higgstag){
    output_higgstag_h->fill(event);
    topjet_higgstag_h->fill(event);
    eff_higgstag_h->fill(event);
    jet_higgstag_h->fill(event);
    muon_higgstag_h->fill(event);
    event_higgstag_h->fill(event);
  }


  ///////////////////////////////////////////////////////// HIGGSTAGGER || ZW TAGGER //////////////////////////////////////////////////////////////////////////////////
  bool pass_zwtag=false;
  if(!pass_higgstag){
    pass_zwtag = ZWtag_sel->passes(event);
    if(!pass_zwtag) return false;
    
    topjet_zwtag_h->fill(event);
    eff_zwtag_h->fill(event);
    jet_zwtag_h->fill(event);
    muon_zwtag_h->fill(event);
    event_zwtag_h->fill(event);
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

  
  //Toptag
 
  output_toptag_h->fill(event);
  topjet_toptag_h->fill(event);
  eff_toptag_h->fill(event);
  jet_toptag_h->fill(event);
  muon_toptag_h->fill(event);
  event_toptag_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_toptag_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_toptag_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_topjet_toptag_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_topjet_toptag_h->fill(event);
 
  ///////////////////////////////////////////////////////////////// Tagger Hist ///////////////////////////////////////////////////////
  topjet_tagger_h->fill(event);
  eff_tagger_h->fill(event);
  jet_tagger_h->fill(event);
  muon_tagger_h->fill(event);
  event_tagger_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_tagger_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_tagger_h->fill(event);
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
  chi2min_chi2cut_h->fill(event);
  topjet_chi2cut_h->fill(event);
  eff_chi2cut_h->fill(event);
  jet_chi2cut_h->fill(event);
  muon_chi2cut_h->fill(event);
  event_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_topjet_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_topjet_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_topjet_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_topjet_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_chi2min_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_chi2min_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_chi2min_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_jet_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_jet_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_jet_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_jet_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_muon_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_muon_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_muon_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_muon_chi2cut_h->fill(event);

  if(pass_higgstag && pass_toptag) higgs_top_event_chi2cut_h->fill(event);
  if(pass_higgstag && !pass_toptag) higgs_notop_event_chi2cut_h->fill(event);
  if(pass_zwtag && pass_toptag)zw_top_event_chi2cut_h->fill(event);
  if(pass_zwtag && !pass_toptag)zw_notop_event_chi2cut_h->fill(event);

  bool pass_btag1 = btag1_sel->passes(event);
  if(pass_btag1){
    topjet_btag1_h->fill(event);
    eff_btag1_h->fill(event);
    jet_btag1_h->fill(event);
    muon_btag1_h->fill(event);
    event_btag1_h->fill(event); 
    chi2min_btag1_h->fill(event);

    if(pass_higgstag && pass_toptag) higgs_top_chi2min_btag1_h->fill(event);
    if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_btag1_h->fill(event);
    if(pass_zwtag && pass_toptag)zw_top_chi2min_btag1_h->fill(event);
    if(pass_zwtag && !pass_toptag)zw_notop_chi2min_btag1_h->fill(event);

  }
 
  bool pass_btag0 = btag0_sel->passes(event);
  if(pass_btag0){
    topjet_btag0_h->fill(event);
    eff_btag0_h->fill(event);
    jet_btag0_h->fill(event);
    muon_btag0_h->fill(event);
    event_btag0_h->fill(event); 
    chi2min_btag0_h->fill(event);

    if(pass_higgstag && pass_toptag) higgs_top_chi2min_btag0_h->fill(event);
    if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_btag0_h->fill(event);
    if(pass_zwtag && pass_toptag)zw_top_chi2min_btag0_h->fill(event);
    if(pass_zwtag && !pass_toptag)zw_notop_chi2min_btag0_h->fill(event);

  }


  lumi_h->fill(event);
  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeSelectionModule)
