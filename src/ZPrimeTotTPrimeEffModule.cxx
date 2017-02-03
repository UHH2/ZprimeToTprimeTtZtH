#include <iostream>
#include <memory>
#include "TH1.h"
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

#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/EffHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h>

#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHypothesisHists.h" 
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstruction.h"
#include <UHH2/ZprimeToTprimeTtZtH/include/TTBarReconstruction.h>
#include <UHH2/common/include/TTbarGen.h>


using namespace uhh2examples;
using namespace uhh2;
using namespace std;

class ZPrimeTotTPrimeEffModule: public uhh2::AnalysisModule {

public:
  explicit ZPrimeTotTPrimeEffModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;

private:
  enum lepton { muon, elec };
  lepton channel_;
  std::unique_ptr<AnalysisModule> printer;

  //cleaner

  std::unique_ptr<MuonCleaner>     muo_cleaner;
  //  std::unique_ptr<ElectronCleaner> ele_cleaner;
  std::unique_ptr<JetCleaner>      jet_IDcleaner;
  std::unique_ptr<JetCleaner>      jet_cleaner2;
  std::unique_ptr<JetCleaner>      jet_cleaner1;
  std::unique_ptr<JetCleaner>      topjet_IDcleaner;
  std::unique_ptr<TopJetCleaner>   topjet_cleaner;
  std::unique_ptr<JetLeptonCleaner> jetlepton_cleaner;
  std::unique_ptr<JetCleaner>  ak4_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner; 
  std::unique_ptr<JetCleaner> jet_btag_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_W_cleaner;
  std::unique_ptr<TopJetCleaner> topjet_PTcleaner;
  std::unique_ptr<JetCleaner>   btag_cleaner_medium;
  std::unique_ptr<JetCleaner>   btag_cleaner_loose;
  std::unique_ptr<JetCleaner>   btag_cleaner_tight;
  std::unique_ptr<TopJetCleaner>      mass_cleaner;

  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;

  //correctors
  //std::unique_ptr<SubJetCorrector> subjetcorrector;

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
  std::unique_ptr<uhh2::Selection> jet_sel;         //at least two ak4 jets
  std::unique_ptr<uhh2::Selection> topjet_sel;      // one AK8 pT>200
  std::unique_ptr<uhh2::Selection> btag2_sel;       // two btags (ak4 jets)
  std::unique_ptr<uhh2::Selection> chi2cut_sel;
  std::unique_ptr<uhh2::Selection> njet3_sel;
  std::unique_ptr<uhh2::Selection> njet2_sel;
  std::unique_ptr<uhh2::Selection> chi2cut20_sel;
  std::unique_ptr<uhh2::Selection> chi2cut10_sel;

  //handle
  uhh2::Event::Handle< std::vector<Jet> > h_btag_medium; // handle of btagged jets (medium) for Reconstruction
  uhh2::Event::Handle< std::vector<Jet> > h_btag_loose; // handle of btagged jets (loose)
  uhh2::Event::Handle< std::vector<Jet> > h_btag_tight; // handle of btagged jets (tight)
  uhh2::Event::Handle<TTbarGen> h_ttbargen;

  uhh2::Event::Handle< std::vector<Jet> > h_btag_jets;
  uhh2::Event::Handle< std::vector<TopJet> > h_WTopjet;
  uhh2::Event::Handle< std::vector<Jet> > h_notused_jets;
  uhh2::Event::Handle< std::vector<TopJet> > h_notused_topjets;

  // Reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod; 
  std::unique_ptr<uhh2::AnalysisModule> reco_primlep;  
  std::unique_ptr<uhh2::AnalysisModule> ttbar_reco;
  std::unique_ptr<uhh2::AnalysisModule> ttbar_chi;
 
  //Tagger
  //Z/W TAGGER
  std::unique_ptr<Selection> ZWtag_sel;
  uhh2::Event::Handle< std::vector<TopJet> > h_ZWtag;


  ////////////////////HISTS///////////////////////////////////////
  //Input Hists
  std::unique_ptr<Hists> input_event_h;
  std::unique_ptr<Hists> input_lep_h;
  std::unique_ptr<Hists> input_eff_h;
  std::unique_ptr<Hists> input_jet_h; 
  std::unique_ptr<Hists> input_topjet_h;
  std::unique_ptr<Hists> input_tagger_h;
  
  //Hists lep1
  std::unique_ptr<Hists> event_lep1_h;
  std::unique_ptr<Hists> topjet_lep1_h;
  std::unique_ptr<Hists> jet_lep1_h;
  std::unique_ptr<Hists> muon_lep1_h;
  std::unique_ptr<Hists> eff_lep1_h;
  std::unique_ptr<Hists> tagger_lep1_h;

  //Hist Topjet2
  std::unique_ptr<Hists> eff_topjet2_h;
  std::unique_ptr<Hists> jet_topjet2_h;
  std::unique_ptr<Hists> muon_topjet2_h;
  std::unique_ptr<Hists> event_topjet2_h;
  std::unique_ptr<Hists> topjet_topjet2_h;
  std::unique_ptr<Hists> tagger_topjet2_h;


  //2DCut Hists
  std::unique_ptr<Hists> eff_twodcut_h;
  std::unique_ptr<Hists> jet_twodcut_h;
  std::unique_ptr<Hists> topjet_twodcut_h;
  std::unique_ptr<Hists> muon_twodcut_h;
  std::unique_ptr<Hists> event_twodcut_h; 
  std::unique_ptr<Hists> tagger_twodcut_h;

  //Jet Hists
  std::unique_ptr<Hists> eff_jet_h;
  std::unique_ptr<Hists> jet_jet_h;
  std::unique_ptr<Hists> topjet_jet_h;
  std::unique_ptr<Hists> muon_jet_h;
  std::unique_ptr<Hists> event_jet_h; 
  std::unique_ptr<Hists> tagger_jet_h;

 //BTag Hists
  std::unique_ptr<Hists> eff_btag_h;
  std::unique_ptr<Hists> jet_btag_h;
  std::unique_ptr<Hists> topjet_btag_h;
  std::unique_ptr<Hists> muon_btag_h;
  std::unique_ptr<Hists> event_btag_h;
  std::unique_ptr<Hists> tagger_btag_h;

 //Reco Hists
  std::unique_ptr<Hists> eff_reco_h;
  std::unique_ptr<Hists> jet_reco_h;
  std::unique_ptr<Hists> topjet_reco_h;
  std::unique_ptr<Hists> muon_reco_h;
  std::unique_ptr<Hists> event_reco_h;
  std::unique_ptr<Hists> chi2min_reco_h;
  std::unique_ptr<Hists> tagger_reco_h;

  //Chi2cut Hists
  std::unique_ptr<Hists> eff_chi2cut_h;
  std::unique_ptr<Hists> jet_chi2cut_comb_h;
  std::unique_ptr<Hists> jet_chi2cut_btag_h;
  std::unique_ptr<Hists> jet_chi2cut_other_h;
  std::unique_ptr<Hists> topjet_chi2cut_comb_h;
  std::unique_ptr<Hists> topjet_chi2cut_w_h;
  std::unique_ptr<Hists> topjet_chi2cut_other_h;
  std::unique_ptr<Hists> muon_chi2cut_h;
  std::unique_ptr<Hists> event_chi2cut_h;
  std::unique_ptr<Hists> chi2min_chi2cut_h;
  std::unique_ptr<Hists> tagger_chi2cut_h;

//ZWTagg Hists
  std::unique_ptr<TopJetHists> input_zwtag_h;
  std::unique_ptr<Hists> eff_zwtag_h;
  std::unique_ptr<Hists> jet_zwtag_comb_h;
  std::unique_ptr<Hists> jet_zwtag_btag_h;
  std::unique_ptr<Hists> jet_zwtag_other_h;
  std::unique_ptr<Hists> topjet_zwtag_comb_h;
  std::unique_ptr<Hists> topjet_zwtag_w_h;
  std::unique_ptr<Hists> topjet_zwtag_other_h;
  std::unique_ptr<Hists> muon_zwtag_h;
  std::unique_ptr<Hists> event_zwtag_h;
  std::unique_ptr<Hists> chi2min_zwtag_h;
  std::unique_ptr<Hists> tagger_zwtag_h;

  //Njet3 Hists
  std::unique_ptr<Hists> eff_njet3_h;
  std::unique_ptr<Hists> jet_njet3_comb_h;
  std::unique_ptr<Hists> jet_njet3_btag_h;
  std::unique_ptr<Hists> jet_njet3_other_h;
  std::unique_ptr<Hists> topjet_njet3_comb_h;
  std::unique_ptr<Hists> topjet_njet3_w_h;
  std::unique_ptr<Hists> topjet_njet3_other_h;
  std::unique_ptr<Hists> muon_njet3_h;
  std::unique_ptr<Hists> event_njet3_h;
  std::unique_ptr<Hists> chi2min_njet3_h;
  std::unique_ptr<Hists> tagger_njet3_h;

//Njet2 Hists
  std::unique_ptr<Hists> eff_njet2_h;
  std::unique_ptr<Hists> jet_njet2_comb_h;
  std::unique_ptr<Hists> jet_njet2_btag_h;
  std::unique_ptr<Hists> jet_njet2_other_h;
  std::unique_ptr<Hists> topjet_njet2_comb_h;
  std::unique_ptr<Hists> topjet_njet2_w_h;
  std::unique_ptr<Hists> topjet_njet2_other_h;
  std::unique_ptr<Hists> muon_njet2_h;
  std::unique_ptr<Hists> event_njet2_h;
  std::unique_ptr<Hists> chi2min_njet2_h;
  std::unique_ptr<Hists> tagger_njet2_h;

// //Chi2cut20 Hists
//   std::unique_ptr<Hists> eff_chi2cut20_h;
//   std::unique_ptr<Hists> jet_chi2cut20_h;
//   std::unique_ptr<Hists> topjet_chi2cut20_h;
//   std::unique_ptr<Hists> muon_chi2cut20_h;
//   std::unique_ptr<Hists> event_chi2cut20_h;
//   std::unique_ptr<Hists> chi2min_chi2cut20_h;
//   std::unique_ptr<Hists> tagger_chi2cut20_h;

// //Chi2cut10 Hists
//   std::unique_ptr<Hists> eff_chi2cut10_h;
//   std::unique_ptr<Hists> jet_chi2cut10_h;
//   std::unique_ptr<Hists> topjet_chi2cut10_h;
//   std::unique_ptr<Hists> muon_chi2cut10_h;
//   std::unique_ptr<Hists> event_chi2cut10_h;
//   std::unique_ptr<Hists> chi2min_chi2cut10_h;
//   std::unique_ptr<Hists> tagger_chi2cut10_h;



// /// Signal Selection /////
//   std::unique_ptr<uhh2::Selection> topjet_sel_signal;      // one AK8 pT>250
//   std::unique_ptr<uhh2::Selection> chi2cut_sel_signal;     // chi2min <50
//   std::unique_ptr<uhh2::Selection> btag1_sel;       //for limits
//   //HIGGS TAGGER
//   std::unique_ptr<Selection> higgstag_sel;
//   uhh2::Event::Handle< std::vector<TopJet> > h_higgstag;
//   std::unique_ptr<TopJetCleaner>   higgstag_cleaner;
//  std::unique_ptr<TopJetCleaner>   ZWtag_cleaner;
//  std::unique_ptr<TopJetCleaner>   toptag_cleaner;
//  //TOP TAGGER
//   std::unique_ptr<Selection> toptag_sel; 
//   uhh2::Event::Handle< std::vector<TopJet> > h_toptag;
//   // Reconstruction ZPrime for Signal
//   std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrimeprod;
//   std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_reco;
//   std::unique_ptr<uhh2::AnalysisModule> ZprimeTotTPrime_chi;
 
//   //Hist Btag1
//   std::unique_ptr<Hists> eff_btag1_h;
//   std::unique_ptr<Hists> jet_btag1_h;
//   std::unique_ptr<Hists> muon_btag1_h;
//   std::unique_ptr<Hists> event_btag1_h;
//   std::unique_ptr<Hists> topjet_btag1_h;
//   std::unique_ptr<Hists> chi2min_btag1_h;

//   std::unique_ptr<Hists> higgs_top_chi2min_btag1_h;
//   std::unique_ptr<Hists> zw_top_chi2min_btag1_h;
//   std::unique_ptr<Hists> higgs_notop_chi2min_btag1_h;
//   std::unique_ptr<Hists> zw_notop_chi2min_btag1_h;

//   uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ZprimeTotTPrime_hyps;

  std::unique_ptr<Hists> lumi_h;

  //systematicen
  std::unique_ptr<AnalysisModule> syst_module;
  bool do_scale_variation;

  //general
  std::string filename;
  double distance;
  double matching;
  uhh2::Event::Handle<double> h_ht;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis> > h_ttbar_hyps;
  bool berror;
};

ZPrimeTotTPrimeEffModule::ZPrimeTotTPrimeEffModule(uhh2::Context& ctx){ 

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
    pileup_SF.reset(new MCPileupReweight(ctx,ctx.get("puVariation"))); 
    lumiweight.reset(new MCLumiWeight(ctx)); 
    btagwAK4.reset(new MCBTagScaleFactor(ctx, CSVBTag::WP_MEDIUM, "jets")); 
    muonscale.reset(new MCMuonScaleFactor(ctx,data_dir_path + "MuonID_Z_RunCD_Reco76X_Feb15.root","MC_NUM_MediumID_DEN_genTracks_PAR_pt_spliteta_bin1", 1.));
  }else     lumi_sel.reset(new LumiSelection(ctx));

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
  muo_cleaner.reset(new MuonCleaner    (AndId<Muon>    (PtEtaCut  (50., 2.1), MuonIDTight())));
  //  ele_cleaner.reset(new ElectronCleaner(AndId<Electron>(PtEtaSCCut(50., 2.4), ElectronID_MVAnotrig_Spring15_25ns_loose)));

  const JetId jetID(JetPFID(JetPFID::WP_LOOSE));
  jet_IDcleaner.reset(new JetCleaner(ctx,jetID));
  jet_cleaner2.reset(new JetCleaner(ctx,15., 2.4));
  jet_cleaner1.reset(new JetCleaner(ctx,30., 2.4));
  topjet_PTcleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(200,2.4,300))));
  jetlepton_cleaner.reset(new JetLeptonCleaner(ctx, JEC_AK4));
  jetlepton_cleaner->set_drmax(.4);
  ak4_cleaner.reset(new JetCleaner(ctx,JetId(ZPrimeTotTPrimeAK4cleaner(1.2))));
  topjetlepton_cleaner.reset(new TopJetLeptonDeltaRCleaner(.8));
 
  const TopJetId massID = Type2TopTag(24,182,Type2TopTag::MassType::groomed);
  mass_cleaner.reset(new TopJetCleaner(ctx,massID));
  //calculator
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));

  // //correctors
  // if(isMC) subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_MC));
  // else subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_DATA));


  ////////////////////////////////////// Selections /////////////////////////////////////////////////
  lep1_sel.reset(new uhh2::AndSelection(ctx));    
 
  if(channel_ == muon){
    lep1_sel->add<NMuonSelection>    ("muoN == 1", 1, 1);
    lep1_sel->add<NElectronSelection>("eleN == 0", 0, 0);
  }
  else if(channel_ == elec){
    lep1_sel->add<NMuonSelection>    ("muoN == 0", 0, 0);
    lep1_sel->add<NElectronSelection>("eleN == 1", 1, 1);
  }


  twodcut_sel.reset(new TwoDCut(.4, 40.));     
  jet_sel.reset(new NJetSelection(2,-1, JetId(PtEtaCut( 50., 2.4))));         
  topjet_sel.reset(new NTopJetSelection(1,-1, TopJetId(PtEtaCut( 200., 2.4)))); 
  njet2_sel.reset(new NJetSelection(2,2, JetId(PtEtaCut( 30., 2.4))));
  njet3_sel.reset(new NJetSelection(2,3, JetId(PtEtaCut( 30., 2.4))));

  const JetId bjetID_medium(CSVBTag(CSVBTag::WP_MEDIUM));
  const JetId bjetID_loose(CSVBTag(CSVBTag::WP_LOOSE));
  const JetId bjetID_tight(CSVBTag(CSVBTag::WP_TIGHT));
  btag2_sel.reset(new NJetSelection(2,-1,bjetID_medium ));
  btag_cleaner_medium.reset(new JetCleaner(ctx,bjetID_medium));
  btag_cleaner_loose.reset(new JetCleaner(ctx,bjetID_loose));
  btag_cleaner_tight.reset(new JetCleaner(ctx,bjetID_tight));

  /////////////////////////////////// Reconstruction /////////////////////////////////////////////////
  const std::string ttbar_gen_label ("ttbargen"); 
  const std::string ttbar_hyps_label("TTbarReconstruction");
  const std::string ttbar_chi2_label("Chi2");


  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  reco_primlep.reset(new PrimaryLepton(ctx));
  ttbar_reco.reset(new LepTopReconstruction(ctx, NeutrinoReconstruction, ttbar_hyps_label));
  ttbar_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ttbar_hyps_label));
  h_ttbar_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ttbar_hyps_label);

  chi2cut_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ttbar_hyps_label,ttbar_chi2_label));
  chi2cut20_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,20,ttbar_hyps_label,ttbar_chi2_label));
  chi2cut10_sel.reset(new ZPrimeTotTPrimeChiCut( ctx,10,ttbar_hyps_label,ttbar_chi2_label));

 // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else genmttbar_sel.reset(new uhh2::AndSelection(ctx));

 //divide event.jets in reco used jets and the rest
  h_btag_jets = ctx.declare_event_output< std::vector<Jet> > ("bjets");
  h_WTopjet = ctx.declare_event_output< std::vector<TopJet> > ("WTopJet");
  h_notused_jets = ctx.declare_event_output< std::vector<Jet> > ("notused_jets");
  h_notused_topjets = ctx.declare_event_output< std::vector<TopJet> > ("notused_topjets");


  //Tagger
  // W/Z TAGGER
  const TopJetId ZWjetID = AndId<TopJet>(Type2TopTag(60,115,Type2TopTag::MassType::groomed), Tau21(0.5));
  ZWtag_sel.reset(new NTopJetSelection(1, -1, ZWjetID));
  

 ////////////////////////////////////// Hists  /////////////////////////////////////////////////
  // Hists Input
  input_event_h.reset(new EventHists(ctx, "input"));
  input_lep_h.reset(new MuonHists(ctx, "input_Lep"));
  input_jet_h.reset(new JetHists (ctx, "input_Jet"));
  input_topjet_h.reset(new TopJetHists (ctx, "input_TopJet"));
  input_eff_h.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
  input_tagger_h.reset(new EffHists(ctx, "input_tagger",ttbar_hyps_label,ttbar_chi2_label));

  lumi_h.reset(new LuminosityHists(ctx,"lumi"));

  // Hists lep1
  topjet_lep1_h.reset(new TopJetHists(ctx, "topjet_lep1"));
  event_lep1_h.reset(new EventHists(ctx, "event_lep1"));
  muon_lep1_h.reset(new MuonHists(ctx, "muon_lep1"));
  jet_lep1_h.reset(new JetHists(ctx, "jet_lep1"));
  eff_lep1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_lep1"));
  tagger_lep1_h.reset(new EffHists(ctx, "tagger_lep1",ttbar_hyps_label,ttbar_chi2_label));

  //Hists topjet2
  topjet_topjet2_h.reset(new TopJetHists(ctx, "topjet_topjet2"));
  eff_topjet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_topjet2"));
  jet_topjet2_h.reset(new JetHists(ctx, "jet_topjet2"));
  muon_topjet2_h.reset(new MuonHists(ctx, "muon_topjet2"));
  event_topjet2_h.reset(new EventHists(ctx, "event_topjet2"));
  tagger_topjet2_h.reset(new EffHists(ctx, "tagger_topjet2",ttbar_hyps_label,ttbar_chi2_label));

  //Hist twodcut
  topjet_twodcut_h.reset(new TopJetHists(ctx, "topjet_twodcut"));
  eff_twodcut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_twodcut"));
  jet_twodcut_h.reset(new JetHists(ctx, "jet_twodcut"));
  muon_twodcut_h.reset(new MuonHists(ctx, "muon_twodcut"));
  event_twodcut_h.reset(new EventHists(ctx, "event_twodcut"));
  tagger_twodcut_h.reset(new EffHists(ctx, "tagger_twodcut",ttbar_hyps_label,ttbar_chi2_label));

  //Hist jet
  topjet_jet_h.reset(new TopJetHists(ctx, "topjet_jet"));
  eff_jet_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_jet"));
  jet_jet_h.reset(new JetHists(ctx, "jet_jet"));
  muon_jet_h.reset(new MuonHists(ctx, "muon_jet"));
  event_jet_h.reset(new EventHists(ctx, "event_jet"));
  tagger_jet_h.reset(new EffHists(ctx, "tagger_jet",ttbar_hyps_label,ttbar_chi2_label));

 //Hist btag
  topjet_btag_h.reset(new TopJetHists(ctx, "topjet_btag"));
  eff_btag_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag"));
  jet_btag_h.reset(new JetHists(ctx, "jet_btag"));
  muon_btag_h.reset(new MuonHists(ctx, "muon_btag"));
  event_btag_h.reset(new EventHists(ctx, "event_btag"));
  tagger_btag_h.reset(new EffHists(ctx, "tagger_btag",ttbar_hyps_label,ttbar_chi2_label));

 //Hist reco
  topjet_reco_h.reset(new TopJetHists(ctx, "topjet_reco"));
  eff_reco_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_reco"));
  jet_reco_h.reset(new JetHists(ctx, "jet_reco"));
  muon_reco_h.reset(new MuonHists(ctx, "muon_reco"));
  event_reco_h.reset(new EventHists(ctx, "event_reco"));
  chi2min_reco_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_reco",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_reco_h.reset(new EffHists(ctx, "tagger_reco",ttbar_hyps_label,ttbar_chi2_label));

  //Hist chi2cut
  topjet_chi2cut_comb_h.reset(new TopJetHists(ctx, "topjet_chi2cut_comb"));
  topjet_chi2cut_w_h.reset(new TopJetHists(ctx, "topjet_chi2cut_W",1,"WTopJet"));
  topjet_chi2cut_other_h.reset(new TopJetHists(ctx, "topjet_chi2cut_other",4,"notused_topjets"));
  eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut"));
  jet_chi2cut_comb_h.reset(new JetHists(ctx, "jet_chi2cut_comb"));
  jet_chi2cut_btag_h.reset(new JetHists(ctx, "jet_chi2cut_btag",4,"bjets"));
  jet_chi2cut_other_h.reset(new JetHists(ctx, "jet_chi2cut_other",4,"notused_jets"));
  muon_chi2cut_h.reset(new MuonHists(ctx, "muon_chi2cut"));
  event_chi2cut_h.reset(new EventHists(ctx, "event_chi2cut"));
  chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_chi2cut_h.reset(new EffHists(ctx, "tagger_chi2cut",ttbar_hyps_label,ttbar_chi2_label));

  //Hist zwtag
  input_zwtag_h.reset(new TopJetHists(ctx, "input_zwtag"));
  input_zwtag_h->set_TopJetId(ZWjetID);
  topjet_zwtag_comb_h.reset(new TopJetHists(ctx, "topjet_zwtag_comb"));
  topjet_zwtag_w_h.reset(new TopJetHists(ctx, "topjet_zwtag_W",4,"WTopJet"));
  topjet_zwtag_other_h.reset(new TopJetHists(ctx, "topjet_zwtag_other",4,"notused_topjets"));
  eff_zwtag_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_zwtag"));
  jet_zwtag_comb_h.reset(new JetHists(ctx, "jet_zwtag_comb"));
  jet_zwtag_btag_h.reset(new JetHists(ctx, "jet_zwtag_btag",4,"bjets"));
  jet_zwtag_other_h.reset(new JetHists(ctx, "jet_zwtag_other",4,"notused_jets"));
  muon_zwtag_h.reset(new MuonHists(ctx, "muon_zwtag"));
  event_zwtag_h.reset(new EventHists(ctx, "event_zwtag"));
  chi2min_zwtag_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_zwtag",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_zwtag_h.reset(new EffHists(ctx, "tagger_zwtag",ttbar_hyps_label,ttbar_chi2_label));


  //Hist njet2
  topjet_njet2_comb_h.reset(new TopJetHists(ctx, "topjet_njet2_comb"));
  topjet_njet2_w_h.reset(new TopJetHists(ctx, "topjet_njet2_W",4,"WTopJet"));
  topjet_njet2_other_h.reset(new TopJetHists(ctx, "topjet_njet2_other",4,"notused_topjets"));
  eff_njet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_njet2"));
  jet_njet2_comb_h.reset(new JetHists(ctx, "jet_njet2_comb"));
  jet_njet2_btag_h.reset(new JetHists(ctx, "jet_njet2_btag",4,"bjets"));
  jet_njet2_other_h.reset(new JetHists(ctx, "jet_njet2_other",4,"notused_jets"));
  muon_njet2_h.reset(new MuonHists(ctx, "muon_njet2"));
  event_njet2_h.reset(new EventHists(ctx, "event_njet2"));
  chi2min_njet2_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_njet2",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_njet2_h.reset(new EffHists(ctx, "tagger_njet2",ttbar_hyps_label,ttbar_chi2_label));

  //Hist njet3
  topjet_njet3_comb_h.reset(new TopJetHists(ctx, "topjet_njet3_comb"));
  topjet_njet3_w_h.reset(new TopJetHists(ctx, "topjet_njet3_W",4,"WTopJet"));
  topjet_njet3_other_h.reset(new TopJetHists(ctx, "topjet_njet3_other",4,"notused_topjets"));
  eff_njet3_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_njet3"));
  jet_njet3_comb_h.reset(new JetHists(ctx, "jet_njet3_comb"));
  jet_njet3_btag_h.reset(new JetHists(ctx, "jet_njet3_btag",4,"bjets"));
  jet_njet3_other_h.reset(new JetHists(ctx, "jet_njet3_other",4,"notused_jets"));
  muon_njet3_h.reset(new MuonHists(ctx, "muon_njet3"));
  event_njet3_h.reset(new EventHists(ctx, "event_njet3"));
  chi2min_njet3_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_njet3",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_njet3_h.reset(new EffHists(ctx, "tagger_njet3",ttbar_hyps_label,ttbar_chi2_label));

 //  //Hist chi2cut20
//   topjet_chi2cut20_h.reset(new TopJetHists(ctx, "topjet_chi2cut20"));
//   eff_chi2cut20_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut20"));
//   jet_chi2cut20_h.reset(new JetHists(ctx, "jet_chi2cut20"));
//   muon_chi2cut20_h.reset(new MuonHists(ctx, "muon_chi2cut20"));
//   event_chi2cut20_h.reset(new EventHists(ctx, "event_chi2cut20"));
//   chi2min_chi2cut20_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut20",ttbar_hyps_label,ttbar_chi2_label ));
//   tagger_chi2cut20_h.reset(new EffHists(ctx, "tagger_chi2cut20",ttbar_hyps_label,ttbar_chi2_label));
 

// //Hist chi2cut10
//   topjet_chi2cut10_h.reset(new TopJetHists(ctx, "topjet_chi2cut10"));
//   eff_chi2cut10_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut10"));
//   jet_chi2cut10_h.reset(new JetHists(ctx, "jet_chi2cut10"));
//   muon_chi2cut10_h.reset(new MuonHists(ctx, "muon_chi2cut10"));
//   event_chi2cut10_h.reset(new EventHists(ctx, "event_chi2cut10"));
//   chi2min_chi2cut10_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut10",ttbar_hyps_label,ttbar_chi2_label ));
//   tagger_chi2cut10_h.reset(new EffHists(ctx, "tagger_chi2cut10",ttbar_hyps_label,ttbar_chi2_label));


// /// Signal selection ///


//   // TOPJET 
//   topjet_sel_signal.reset(new NTopJetSelection(1, -1, TopJetId(PtEtaCut( 250., 2.4))));
 


//   // Higgs TAGGER
//   const TopJetId higgsjetID = AndId<TopJet>(HiggsTag(100,150), Tau21(1) );
//   higgstag_sel.reset(new NTopJetSelection(1, -1, higgsjetID));
//   higgstag_cleaner.reset(new TopJetCleaner(ctx,higgsjetID));
//   h_higgstag = ctx.declare_event_output< std::vector<TopJet> > ("HiggsTag");
//   h_ZWtag = ctx.declare_event_output< std::vector<TopJet> > ("ZWTag");

//   ZWtag_cleaner.reset(new TopJetCleaner(ctx,ZWjetID));
 
// // TOP TAGGER
//   const TopJetId topjetID = AndId<TopJet>(Type2TopTag(150,240,Type2TopTag::MassType::groomed), Tau32(0.6));
//   toptag_sel.reset(new NTopJetSelection(1, -1, topjetID));
//   toptag_cleaner.reset(new TopJetCleaner(ctx,topjetID));
//  h_toptag = ctx.declare_event_output< std::vector<TopJet> > ("TopTag");

//   /* KINEMATICAL RECO */
//   const std::string ZprimeTotTPrime_gen_label ("zprimegen");
//   const std::string ZprimeTotTPrime_hyps_label("ZPrimeTotTPrimeReconstruction");
//   const std::string ZprimeTotTPrime_selection_hyps_label("ZPrimeTotTPrimeReconstruction"); 
//   const std::string ZprimeTotTPrime_chi2_label("Chi2");
// //GenParticles
//   ZprimeTotTPrimeprod.reset(new ZPrimeGenProducer(ctx, ZprimeTotTPrime_gen_label, false));
//   ZprimeTotTPrime_reco.reset(new ZPrimeTotTPrimeReconstruction(ctx, NeutrinoReconstruction, ZprimeTotTPrime_hyps_label));
//   ZprimeTotTPrime_chi.reset(new ZPrimeTotTPrimeChi2Discriminator(ctx, ZprimeTotTPrime_hyps_label));
//   h_ZprimeTotTPrime_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(ZprimeTotTPrime_hyps_label);
 
//   chi2cut_sel_signal.reset(new ZPrimeTotTPrimeChiCut( ctx,50,ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));


//  btag1_sel.reset(new NBTagSelection(1, -1));
//  //Hists btag1
//   topjet_btag1_h.reset(new TopJetHists(ctx, "topjet_btag1"));
//   eff_btag1_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_btag1"));
//   jet_btag1_h.reset(new JetHists(ctx, "jet_btag1"));
//   muon_btag1_h.reset(new MuonHists(ctx, "muon_btag1"));
//   event_btag1_h.reset(new EventHists(ctx, "event_btag1"));
//   chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label ));

//   higgs_top_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_top_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
//   higgs_notop_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "higgs_notop_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
//   zw_top_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_top_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));
//   zw_notop_chi2min_btag1_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "zw_notop_chi2min_btag1",ZprimeTotTPrime_hyps_label,ZprimeTotTPrime_chi2_label));

  //general
  filename =  ctx.get("dataset_version"); 
  distance = string2double( ctx.get("distance"));
  matching = string2double( ctx.get("matching"));
  h_ht = ctx.get_handle<double>("HT");
  h_btag_medium = ctx.declare_event_output< std::vector<Jet> > ("BTag_medium");
  h_btag_loose = ctx.declare_event_output< std::vector<Jet> > ("BTag_loose");
  h_btag_tight = ctx.declare_event_output< std::vector<Jet> > ("BTag_tight");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

  //Trigger
  //  const std::string& trigger = ctx.get("trigger", "NULL");
  //if(trigger != "NULL") trigger_sel = make_unique<TriggerSelection>(trigger);
  //else                  trigger_sel = make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*");


  //systematicen
  syst_module.reset(new MCScaleVariation(ctx));

 do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");

  berror=false;
}

bool ZPrimeTotTPrimeEffModule::process(uhh2::Event& event){
  
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
  if(filename  == "TTbar_Mtt1000toINFT" || filename  == "TTbar_Mtt0700to1000") ttgenprod->process(event);
 
  //common Modules
  /* luminosity sections from CMS silver-JSON file */
  if(event.isRealData && !lumi_sel->passes(event)) return false;
  /* pileup SF */
  if(!event.isRealData){ 
    pileup_SF->process(event);
    lumiweight->process(event);
    btagwAK4->process(event);
    muonscale->process(event);
  }
  ////

 //systematicen
 if(do_scale_variation) syst_module->process(event); 



  //correctors
  // subjetcorrector->process(event);

  // OBJ CLEANING
  muo_cleaner->process(event);
  sort_by_pt<Muon>(*event.muons);

  //  ele_cleaner->process(event);
  // sort_by_pt<Electron>(*event.electrons);

  jet_IDcleaner->process(event);
  jetlepton_cleaner->process(event);
  jet_cleaner2->process(event); 
  
  sort_by_pt<Jet>(*event.jets);
  sort_by_pt<TopJet>(*event.topjets);

  // ///////Trigger///////
  // const bool pass_trigger = trigger_sel->passes(event);
  // if(!pass_trigger) return false;

  if(berror)   std::cout<<"SelectionModule L:232 vor Input Histogrammen"<<std::endl;
  /////////////////////////////////////////////////////////// Input Histogramme ///////////////////////////////////////////////////////////////////////////////
 
  input_eff_h ->fill(event);
  input_event_h->fill(event);
  input_topjet_h->fill(event);
  input_jet_h->fill(event);
  input_lep_h->fill(event);
  input_tagger_h->fill(event);

  if(berror)  std::cout<<"SelectionModule L:268 vor LeptonSelection"<<std::endl;
  /////////////////////////////////////////////////////////// LEPTON selection ///////////////////////////////////////////////////////////////////////////////
  const bool pass_lep1 = lep1_sel->passes(event);
  if(!pass_lep1) return false;
  muon_lep1_h->fill(event);
  eff_lep1_h ->fill(event);
  event_lep1_h->fill(event);
  topjet_lep1_h->fill(event);
  jet_lep1_h->fill(event);
  tagger_lep1_h->fill(event);

  
  topjetlepton_cleaner->process(event);
  sort_by_pt<TopJet>(*event.topjets);

 
  if(berror)  std::cout<<"SelectionModule L:268 vor TopJetSelection"<<std::endl; 
  ////////////////////////////////////////////////////////// TopJET selection//////////////////////////////////////////////////////////////////////////////////
  // topjet_PTcleaner->process(event);
  const bool pass_topjet = topjet_sel->passes(event);
  if(!pass_topjet) return false;

  topjet_topjet2_h->fill(event);
  eff_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);
  muon_topjet2_h->fill(event);
  event_topjet2_h->fill(event);
  tagger_topjet2_h->fill(event);

  if(berror)   std::cout<<"SelectionModule L:294 vor TwoDCut"<<std::endl;  
  ////////////////////////////////////////////////////////// TwoDCut //////////////////////////////////////////////////////////////////////////////////

  sort_by_pt<Jet>(*event.jets);
  const bool pass_twodcut = twodcut_sel->passes(event);
  if(!pass_twodcut) return false;
  topjet_twodcut_h->fill(event);
  jet_twodcut_h->fill(event);
  eff_twodcut_h->fill(event);
  muon_twodcut_h->fill(event);
  tagger_twodcut_h->fill(event);

  jet_cleaner1->process(event);
  event_twodcut_h->fill(event);


  if(berror)   std::cout<<"SelectionModule L:294 vor btag handle"<<std::endl;  
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
  
  //// Set handle on btagged jets
  // save all jets before cleaning
  jets_all->clear();
  jets_all->reserve(event.jets->size());
  for(const auto & j : *event.jets) jets_all->push_back(j); 
  //cleaning of jets to get all btagged jets
  btag_cleaner_loose->process(event); 
  sort_by_pt<Jet>(*event.jets);
  //save all btagged jets
  jets_btag->clear();
  jets_btag->reserve(event.jets->size());
  for(const auto & j : *event.jets) jets_btag->push_back(j); 
  // put original jets back
  event.jets->clear();
  event.jets->reserve(jets_all->size());
  for(const auto & j : *jets_all) event.jets->push_back(j); 
  sort_by_pt<Jet>(*event.jets);
  //handle on btagged jets
  event.set(h_btag_loose,*jets_btag );
  //// END Set handle on btagged jets

  //// Set handle on btagged jets
  // save all jets before cleaning
  jets_all->clear();
  jets_all->reserve(event.jets->size());
  for(const auto & j : *event.jets) jets_all->push_back(j); 
  //cleaning of jets to get all btagged jets
  btag_cleaner_tight->process(event); 
  sort_by_pt<Jet>(*event.jets);
  //save all btagged jets
  jets_btag->clear();
  jets_btag->reserve(event.jets->size());
  for(const auto & j : *event.jets) jets_btag->push_back(j); 
  // put original jets back
  event.jets->clear();
  event.jets->reserve(jets_all->size());
  for(const auto & j : *jets_all) event.jets->push_back(j); 
  sort_by_pt<Jet>(*event.jets);
  //handle on btagged jets
  event.set(h_btag_tight,*jets_btag );
  //// END Set handle on btagged jets



  if(berror)   std::cout<<"SelectionModule L:294 vor Jet Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// Jet Selction //////////////////////////////////////////////////////////////////////////////////
  const bool pass_jet = jet_sel->passes(event);
  if(!pass_jet) return false;
  topjet_jet_h->fill(event);
  jet_jet_h->fill(event);
  eff_jet_h->fill(event);
  muon_jet_h->fill(event);
  event_jet_h->fill(event);
  tagger_jet_h->fill(event);

  if(berror)   std::cout<<"SelectionModule L:294 vor btag Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// btag Selction //////////////////////////////////////////////////////////////////////////////////
  const bool pass_btag = btag2_sel->passes(event);
  if(!pass_btag)return false;
  topjet_btag_h->fill(event);
  jet_btag_h->fill(event);
  eff_btag_h->fill(event);
  muon_btag_h->fill(event);
  event_btag_h->fill(event);
  tagger_btag_h->fill(event);

 
  if(berror) std::cout<<"Number topjet events "<<event.topjets->size() << std::endl;

  if(berror)   std::cout<<"SelectionModule L:303 vor KinReco"<<std::endl;
  //////////////////////////////////////////////////////////// KIN RECO///////////////////////////////////////////////////////////////////////////////////

  reco_primlep->process(event);
  if(!(ttbar_reco->process(event)))return false;
  ttbar_reco->process(event);
  if(berror)std::cout<<"Nach reco"<<std::endl;
  ttbar_chi->process(event); 
  if(berror)std::cout<<"Nach ttbar-Chi"<<std::endl;
  // save only the chi2-best ttbar hypothesis in output sub-ntuple
  std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps = event.get(h_ttbar_hyps);
  if(berror)std::cout <<"Nach best Hyp"<<std::endl;
  const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis(hyps, "Chi2");
  if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
  const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj(*hyp);

  if(berror) std::cout<<"anzahl jets und Top Jets vor cleaning "<<event.jets->size()<<" "<<event.topjets->size()<<std::endl;
  
  if(berror)std::cout<<"In Selection Module: v4 blep " <<hyp->blep_v4() << " v4 bhad "<<hyp->bhad_v4() << "v4 W "<< hyp->W_v4()<<std::endl;

  ///////////////////////////////////// AK8 cleaning from used AK4 //////////////////////////////
  std::vector<Jet>* usedak4(new std::vector<Jet> (*event.jets));
  usedak4->clear();
  usedak4->reserve(event.jets->size());
  for(const Jet jet: hyp->toplep_jets()) usedak4->push_back(jet);
  for(const Jet jet: hyp->tophad_jets()) usedak4->push_back(jet);
  event.set(h_btag_jets,*usedak4);
  if(berror) std::cout<<"SelectionModule L:858 Size usedak4 "<<usedak4->size() <<std::endl;
  
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
      if(deltar < distance) bdeltaR=false;
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
  if(berror) std::cout<<"anzahl jets und Top Jets nach AK8 cleaning "<<event.jets->size()<<" "<<event.topjets->size()<<std::endl;
  ///////////////////////////////////// AK4 cleaning from used AK8 //////////////////////////////
  std::vector<Jet>* AK4Jets(new std::vector<Jet> (*event.jets));
  AK4Jets->clear();
  AK4Jets->reserve(event.jets->size());
  for(const Jet ak4:*event.jets){
    double deltar_ak4hyp = deltaR(ak4,hyp->W_v4());
    if(deltar_ak4hyp > distance) AK4Jets ->push_back(ak4);
  }
  sort_by_pt<Jet>(*AK4Jets);
  event.jets->clear();
  event.jets->reserve(AK4Jets->size());
  for(const auto & j : *AK4Jets) event.jets->push_back(j); 
  sort_by_pt<Jet>(*event.jets);

  std::vector<TopJet>* WBoson(new std::vector<TopJet> (*event.topjets));
  WBoson->clear();
  WBoson->reserve(event.jets->size());
  WBoson->push_back(hyp->HZW());
  event.set(h_WTopjet,*WBoson);

  /////////////////AK4 cleaning end ////////////////
  if(berror)   std::cout<<"anzahl jets und Top Jets nach ganzem cleaning "<<event.jets->size()<<" "<<event.topjets->size()<<std::endl;

  ///////////////////////////////////// divide event.jets in reco used and rest (also for topjets) //////////////////////////////
  std::vector<Jet>* AK4Jets_btagcleaned(new std::vector<Jet> (*event.jets));
  AK4Jets_btagcleaned->clear();
  AK4Jets_btagcleaned->reserve(event.jets->size()-2);
  for(auto j:*event.jets){
    if(!(deltaR(j,usedak4->at(0))==0 || deltaR(j,usedak4->at(1))==0)) AK4Jets_btagcleaned->push_back(j);
  }
  event.set(h_notused_jets,*AK4Jets_btagcleaned );
  
 

  //Topjets
  std::vector<TopJet>* notused_AK8Jets(new std::vector<TopJet> (*event.topjets));
  notused_AK8Jets->clear();
  notused_AK8Jets->reserve(event.topjets->size()-1);
  for(auto j:*event.topjets){
  
    if(!(deltaR(j,hyp->W_v4())==0 ))notused_AK8Jets ->push_back(j);
  }
  event.set(h_notused_topjets, *notused_AK8Jets);

  hyps.clear();
  hyps.push_back(hyp_obj);

  if(berror)std::cout << "ZPrimeTotTPrimeEffModule: nach kin Reco"<<std::endl;
  topjet_reco_h->fill(event);
  jet_reco_h->fill(event);
  eff_reco_h->fill(event);
  muon_reco_h->fill(event);
  event_reco_h->fill(event);
  chi2min_reco_h->fill(event);
  tagger_reco_h->fill(event);

 mass_cleaner->process(event);

  if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// chi2cut Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_chi2cut = chi2cut_sel->passes(event);
  if(!pass_chi2cut) return false;

  lumi_h->fill(event); 

  // topjet_chi2cut_comb_h->fill(event);
  // topjet_chi2cut_w_h->fill(event);
  // topjet_chi2cut_other_h->fill(event);
  // jet_chi2cut_comb_h->fill(event);
  // jet_chi2cut_btag_h->fill(event);
  // jet_chi2cut_other_h->fill(event);
  // eff_chi2cut_h->fill(event);
  // muon_chi2cut_h->fill(event);
  // event_chi2cut_h->fill(event);
  // chi2min_chi2cut_h->fill(event);
  // tagger_chi2cut_h->fill(event);

  // if(berror)   std::cout<<"SelectionModule L:294 vor Z/W tag"<<std::endl;  
  // ////////////////////////////////////////////////////////// zw tag Selction //////////////////////////////////////////////////////////////////////////////////
  //  bool pass_zwtag = ZWtag_sel->passes(event);
  // if(!pass_zwtag) return false;
  // input_zwtag_h->fill(event);
  // topjet_zwtag_comb_h->fill(event);
  // topjet_zwtag_w_h->fill(event);
  // topjet_zwtag_other_h->fill(event);
  // jet_zwtag_comb_h->fill(event);
  // jet_zwtag_btag_h->fill(event);
  // jet_zwtag_other_h->fill(event);
  // eff_zwtag_h->fill(event);
  // muon_zwtag_h->fill(event);
  // event_zwtag_h->fill(event);
  // chi2min_zwtag_h->fill(event);
  // tagger_zwtag_h->fill(event);


 ///////////////////////////////////// hists for wrong and right W //////////////////////////////
 
 hyps = event.get(h_ttbar_hyps);
  hyp = get_best_hypothesis(hyps, "Chi2");
  if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
  const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj3(*hyp);

 
 
  if(event.is_valid(h_ttbargen)){
   

    const auto & ttbargen = event.get(h_ttbargen);
    if(ttbargen.IsSemiLeptonicDecay() ){


      double delta_q1 = deltaR(hyp->W_v4(), ttbargen.Q1());
      double delta_q2 = deltaR(hyp->W_v4(), ttbargen.Q2());
  
      if(delta_q1 < matching && delta_q2 < matching){
	// // 	right
	// topjet_chi2cut_comb_h->fill(event);
	// topjet_chi2cut_w_h->fill(event);
	// topjet_chi2cut_other_h->fill(event);
	// jet_chi2cut_comb_h->fill(event);
	// jet_chi2cut_btag_h->fill(event);
	// jet_chi2cut_other_h->fill(event);
	// eff_chi2cut_h->fill(event);
	// muon_chi2cut_h->fill(event);
	// event_chi2cut_h->fill(event);
	// chi2min_chi2cut_h->fill(event);
	// tagger_chi2cut_h->fill(event);
 	// bool pass_zwtag = ZWtag_sel->passes(event);
 	// if(!pass_zwtag) return false;
 	// input_zwtag_h->fill(event);
 	// topjet_zwtag_comb_h->fill(event);
 	// topjet_zwtag_w_h->fill(event);
 	// topjet_zwtag_other_h->fill(event);
 	// jet_zwtag_comb_h->fill(event);
 	// jet_zwtag_btag_h->fill(event);
 	// jet_zwtag_other_h->fill(event);
 	// eff_zwtag_h->fill(event);
 	// muon_zwtag_h->fill(event);
 	// event_zwtag_h->fill(event);
 	// chi2min_zwtag_h->fill(event);
 	// tagger_zwtag_h->fill(event);
 
      }else{
	// 	wrong
	topjet_chi2cut_comb_h->fill(event);
	topjet_chi2cut_w_h->fill(event);
	topjet_chi2cut_other_h->fill(event);
	jet_chi2cut_comb_h->fill(event);
	jet_chi2cut_btag_h->fill(event);
	jet_chi2cut_other_h->fill(event);
	eff_chi2cut_h->fill(event);
	muon_chi2cut_h->fill(event);
	event_chi2cut_h->fill(event);
	chi2min_chi2cut_h->fill(event);
	tagger_chi2cut_h->fill(event);

	if(berror)   std::cout<<"SelectionModule L:294 vor Z/W tag"<<std::endl;  
	////////////////////////////////////////////////////////////////////////////////zw Selction //////////////////////////////////////////////////////////////////////////////////
	bool pass_zwtag = ZWtag_sel->passes(event);
	if(!pass_zwtag) return false;
	input_zwtag_h->fill(event);
	topjet_zwtag_comb_h->fill(event);
	topjet_zwtag_w_h->fill(event);
	topjet_zwtag_other_h->fill(event);
	jet_zwtag_comb_h->fill(event);
	jet_zwtag_btag_h->fill(event);
	jet_zwtag_other_h->fill(event);
	eff_zwtag_h->fill(event);
	muon_zwtag_h->fill(event);
	event_zwtag_h->fill(event);
	chi2min_zwtag_h->fill(event);
	tagger_zwtag_h->fill(event);

      }
    }
  }
  hyps.clear();
  hyps.push_back(hyp_obj3);


 //  if(berror)   std::cout<<"SelectionModule L:294 vor Signal Selection"<<std::endl;  
//   ////////////////////////////////////////////////////////// Signal Selction //////////////////////////////////////////////////////////////////////////////////
//   std::cout<<"topjet"<<std::endl;  
//   const bool pass_topjet_signal = topjet_sel_signal->passes(event);
//   if(!pass_topjet_signal) return false;

//   std::unique_ptr< std::vector<TopJet> > higgsjets_all (new std::vector<TopJet> (*event.topjets));

  
//   //cleanen der topjets um sie in einen neuen Vektor zu speichern
//   higgstag_cleaner->process(event);
//   sort_by_pt<TopJet>(*event.topjets);
//   ak4_cleaner->process(event);
//   std::unique_ptr< std::vector<TopJet>>  topjets_higgstag(new std::vector<TopJet> (*event.topjets));

//   // zurueck speichern der all topjets
//   event.topjets->clear();
  
//   event.topjets->reserve(higgsjets_all->size());

//   for(const auto & j : *higgsjets_all) event.topjets->push_back(j); 
//   sort_by_pt<TopJet>(*event.topjets);
//   //handle auf toptag Topjets
//   event.set(h_higgstag,*topjets_higgstag );

 
//   //Higgstag
//   const bool pass_higgstag = higgstag_sel->passes(event);
//    pass_zwtag=false;
//   if(!pass_higgstag){
//     pass_zwtag = ZWtag_sel->passes(event);
//     if(!pass_zwtag) return false;
//   }

//   std::cout<<"ZWtag"<<std::endl; 
//   std::unique_ptr< std::vector<TopJet> > ZWjets_all(new std::vector<TopJet> (*event.topjets));
//   //cleanen der topjets um sie in einen neuen Vektor zu speichern
//   ZWtag_cleaner->process(event);
//   sort_by_pt<TopJet>(*event.topjets);
//   if(pass_zwtag){
//     ak4_cleaner->process(event);

//   }
//   std::cout<<"ZWtag abspeichern"<<std::endl; 
//   std::unique_ptr< std::vector<TopJet>>  topjets_ZWtag(new std::vector<TopJet> (*event.topjets));
//   // zurueck speichern der all topjets
//   event.topjets->clear();
//   event.topjets->reserve(ZWjets_all->size());
//   for(const auto & j : *ZWjets_all) event.topjets->push_back(j); 
//   sort_by_pt<TopJet>(*event.topjets);
//   //handle auf toptag Topjets
//   event.set(h_ZWtag,*topjets_ZWtag );


//  if(berror)  std::cout<<"SelectionModule L:338 vor TopTAGGER"<<std::endl;  
//   ///////////////////////////////////////////////////////// TOP TAGGER //////////////////////////////////////////////////////////////////////////////////
//   const bool pass_toptag = toptag_sel->passes(event);

//   //alle topjets
//   std::unique_ptr< std::vector<TopJet> > topjets_all (new std::vector<TopJet> (*event.topjets));
//   //cleanen der topjets um sie in einen neuen Vektor zu speichern
//   toptag_cleaner->process(event);
//   sort_by_pt<TopJet>(*event.topjets);
//   if(pass_toptag){
//     ak4_cleaner->process(event);
//   }
//   std::unique_ptr< std::vector<TopJet>>  topjets_toptag(new std::vector<TopJet> (*event.topjets));
//   // zurueck speichern der all topjets
//   event.topjets->clear();
//   event.topjets->reserve(topjets_all->size());
//   for(const auto & j : *topjets_all) event.topjets->push_back(j); 
//   sort_by_pt<TopJet>(*event.topjets);
//   //handle auf toptag Topjets
//   event.set(h_toptag,*topjets_toptag );



//   std::cout<<"Reco"<<std::endl; 
//   reco_primlep->process(event);
//   std::cout<<"ZPrime reco"<<std::endl; 
//   if(!(ZprimeTotTPrime_reco->process(event)))return false;
//   ZprimeTotTPrime_reco->process(event);
//  std::cout<<"Chi2"<<std::endl; 
//   ZprimeTotTPrime_chi->process(event); 
//   // save only the chi2-best ttbar hypothesis in output sub-ntuple
// std::cout<<"best chi2"<<std::endl; 
//   std::vector<ZPrimeTotTPrimeReconstructionHypothesis>& hyps_signal = event.get(h_ZprimeTotTPrime_hyps);
//   const ZPrimeTotTPrimeReconstructionHypothesis* hyp_signal = get_best_hypothesis(hyps_signal, "Chi2");
//   if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
//   const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj_signal(*hyp_signal);
//   std::cout<<"cleaning"<<std::endl; 

//   ///////////////////////////////////// AK8 cleaning from used AK4 //////////////////////////////
//   std::vector<Jet>* usedak4_signal(new std::vector<Jet> (*event.jets));
//   usedak4_signal->clear();
//   usedak4_signal->reserve(event.jets->size());
//   for(const Jet jet: hyp->toplep_jets()) usedak4_signal->push_back(jet);
//   for(const Jet jet: hyp->tophad_jets()) usedak4_signal->push_back(jet);

//   if(berror) std::cout<<"SelectionModule L:858 Size used ak4 "<<usedak4->size() <<std::endl;

//   //Cleaning AK8 by overlap of AK4
//   std::vector<TopJet>* AK8Jets_signal(new std::vector<TopJet> (*event.topjets));
//   AK8Jets_signal->clear();
//   AK8Jets_signal->reserve(event.topjets->size());

//   if(berror) std::cout<<"SelectionModule L:858 Size AK8 before cleaning "<<event.topjets->size() <<std::endl;
//   for(const TopJet ak8:*event.topjets){
//     bool bdeltaR=true;
//     for(const Jet ak4:*usedak4_signal){
//       double deltar = deltaR(ak4,ak8);
//       if(berror) std::cout<<"SelectionModule L:858 DeltaR(ak4, ak8)<1.2 "<<deltar <<std::endl;
//       if(deltar < 1.2) bdeltaR=false;
//       if(berror) std::cout<<"SelectionModule L:858 bdeltaR  "<<bdeltaR <<std::endl;
//     }
//     if(bdeltaR)AK8Jets_signal ->push_back(ak8);
//   }
//   sort_by_pt<TopJet>(*AK8Jets_signal);
//   ////put cleaned AK8 jets in event.topjet
//   event.topjets->clear();
//   event.topjets->reserve(AK8Jets_signal->size());
//   for(const auto & j : *AK8Jets_signal) event.topjets->push_back(j); 
//   sort_by_pt<TopJet>(*event.topjets);
//   if(berror) std::cout<<"SelectionModule L:858 Size topjets Collection "<<event.topjets->size() <<std::endl;

//   /////////////////AK8 cleaning end ////////////////


//   hyps_signal.clear();
//   hyps_signal.push_back(hyp_obj_signal);
//   if(berror)   std::cout << "SelectionModule L:715 vor Chi2Cut"<<std::endl;
//   //////////////////////////////////////////////////////////  CHi2 Cut  ////////////////////////////////////////////////////////
//   std::cout<<"chi2cut"<<std::endl; 

//   bool pass_chi2cut_signal = chi2cut_sel_signal->passes(event);
//   if(!pass_chi2cut_signal) return false;


//   bool pass_btag1 = btag1_sel->passes(event);
//   if(pass_btag1){
//     topjet_btag1_h->fill(event);
//     eff_btag1_h->fill(event);
//     jet_btag1_h->fill(event);
//     muon_btag1_h->fill(event);
//     event_btag1_h->fill(event); 
//     chi2min_btag1_h->fill(event);
//     bool pass_toptag =false;
//     if(pass_higgstag && pass_toptag) higgs_top_chi2min_btag1_h->fill(event);
//     if(pass_higgstag && !pass_toptag) higgs_notop_chi2min_btag1_h->fill(event);
//     if(pass_zwtag && pass_toptag)zw_top_chi2min_btag1_h->fill(event);
//     if(pass_zwtag && !pass_toptag)zw_notop_chi2min_btag1_h->fill(event);

//   }





















  if(berror)   std::cout<<"SelectionModule L:294 vor njet3 Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// njet3 Selction //////////////////////////////////////////////////////////////////////////////////
  // bool pass_njet3 = njet3_sel->passes(event);
  // if(!pass_njet3) return false;
  
  topjet_njet3_comb_h->fill(event);
  topjet_njet3_w_h->fill(event);
  topjet_njet3_other_h->fill(event);
  jet_njet3_comb_h->fill(event);
  jet_njet3_btag_h->fill(event);
  jet_njet3_other_h->fill(event);
  eff_njet3_h->fill(event);
  muon_njet3_h->fill(event);
  event_njet3_h->fill(event);
  chi2min_njet3_h->fill(event);
  tagger_njet3_h->fill(event);

  if(berror)   std::cout<<"SelectionModule L:294 vor njet2 Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// njet2 Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_njet2 = njet2_sel->passes(event);
  if(!pass_njet2) return false;

  topjet_njet2_comb_h->fill(event);
  topjet_njet2_w_h->fill(event);
  topjet_njet2_other_h->fill(event);
  jet_njet2_comb_h->fill(event);
  jet_njet2_btag_h->fill(event);
  jet_njet2_other_h->fill(event);
  eff_njet2_h->fill(event);
  muon_njet2_h->fill(event);
  event_njet2_h->fill(event);
  chi2min_njet2_h->fill(event);
  tagger_njet2_h->fill(event);
 

//    ///////////////////////////////////////////////////////////////////////////////////////
//  if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut 20 Selection"<<std::endl;  
//   ////////////////////////////////////////////////////////// chi2cut 20 Selction //////////////////////////////////////////////////////////////////////////////////
//   bool pass_chi2cut20 = chi2cut20_sel->passes(event);
//   if(!pass_chi2cut20) return false;

//   topjet_chi2cut20_h->fill(event);
//   jet_chi2cut20_h->fill(event);
//    eff_chi2cut20_h->fill(event);
//   muon_chi2cut20_h->fill(event);
//   event_chi2cut20_h->fill(event);
//   chi2min_chi2cut20_h->fill(event);
//   tagger_chi2cut20_h->fill(event);

// if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut 10 Selection"<<std::endl;  
//   ////////////////////////////////////////////////////////// chi2cut 10 Selction //////////////////////////////////////////////////////////////////////////////////
//   bool pass_chi2cut10 = chi2cut10_sel->passes(event);
//   if(!pass_chi2cut10) return false;

//   topjet_chi2cut10_h->fill(event);
//   jet_chi2cut10_h->fill(event);
//    eff_chi2cut10_h->fill(event);
//   muon_chi2cut10_h->fill(event);
//   event_chi2cut10_h->fill(event);
//   chi2min_chi2cut10_h->fill(event);
//   tagger_chi2cut10_h->fill(event);

 ////////////////////////////////////////////////////////// Lumi Hists //////////////////////////////////////////////////////////////////////////////////


  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeEffModule)
