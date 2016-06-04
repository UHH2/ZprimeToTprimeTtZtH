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
#include <UHH2/ZPrimeTotTPrime/include/EffHists.h>
#include <UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h>



// #include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHypothesisHists.h" 

//#include <UHH2/common/include/HypothesisHists.h>
#include <UHH2/ZPrimeTotTPrime/include/TTBarReconstruction.h>
//#include <UHH2/common/include/ReconstructionHypothesis.h>
//#include <UHH2/common/include/ReconstructionHypothesisDiscriminators.h>
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
  std::unique_ptr<ElectronCleaner> ele_cleaner;
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

  std::unique_ptr<JetCleaner>   btag_cleaner_medium;
  std::unique_ptr<JetCleaner>   btag_cleaner_loose;
  std::unique_ptr<JetCleaner>   btag_cleaner_tight;

  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;

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

//MW0 Hists
  std::unique_ptr<Hists> eff_MW0_h;
  std::unique_ptr<Hists> jet_MW0_h;
  std::unique_ptr<Hists> topjet_MW0_h;
  std::unique_ptr<Hists> muon_MW0_h;
  std::unique_ptr<Hists> event_MW0_h;
  std::unique_ptr<Hists> chi2min_MW0_h;
  std::unique_ptr<Hists> tagger_MW0_h;

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

//wrong W and right W
  std::unique_ptr<Hists> eff_wrong_h;
  std::unique_ptr<Hists> jet_wrong_comb_h;
  std::unique_ptr<Hists> jet_wrong_btag_h;
  std::unique_ptr<Hists> jet_wrong_other_h;
  std::unique_ptr<Hists> topjet_wrong_comb_h;
  std::unique_ptr<Hists> topjet_wrong_w_h;
  std::unique_ptr<Hists> topjet_wrong_other_h;
  std::unique_ptr<Hists> muon_wrong_h;
  std::unique_ptr<Hists> event_wrong_h;
  std::unique_ptr<Hists> tagger_wrong_h;
  std::unique_ptr<Hists> chi2min_wrong_h;

  std::unique_ptr<Hists> eff_right_h;
  std::unique_ptr<Hists> jet_right_comb_h;
  std::unique_ptr<Hists> jet_right_btag_h;
  std::unique_ptr<Hists> jet_right_other_h;
  std::unique_ptr<Hists> topjet_right_comb_h;
  std::unique_ptr<Hists> topjet_right_w_h;
  std::unique_ptr<Hists> topjet_right_other_h;
  std::unique_ptr<Hists> muon_right_h;
  std::unique_ptr<Hists> event_right_h;
  std::unique_ptr<Hists> tagger_right_h;
  std::unique_ptr<Hists> chi2min_right_h;

//Chi2cut20 Hists
  std::unique_ptr<Hists> eff_chi2cut20_h;
  std::unique_ptr<Hists> jet_chi2cut20_h;
  std::unique_ptr<Hists> topjet_chi2cut20_h;
  std::unique_ptr<Hists> muon_chi2cut20_h;
  std::unique_ptr<Hists> event_chi2cut20_h;
  std::unique_ptr<Hists> chi2min_chi2cut20_h;
  std::unique_ptr<Hists> tagger_chi2cut20_h;

//Chi2cut10 Hists
  std::unique_ptr<Hists> eff_chi2cut10_h;
  std::unique_ptr<Hists> jet_chi2cut10_h;
  std::unique_ptr<Hists> topjet_chi2cut10_h;
  std::unique_ptr<Hists> muon_chi2cut10_h;
  std::unique_ptr<Hists> event_chi2cut10_h;
  std::unique_ptr<Hists> chi2min_chi2cut10_h;
  std::unique_ptr<Hists> tagger_chi2cut10_h;



  std::unique_ptr<Hists> lumi_h;

  //general
  std::string filename;
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
  if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx)); 
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

//Hist MW0
  topjet_MW0_h.reset(new TopJetHists(ctx, "topjet_MW0"));
  eff_MW0_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_MW0"));
  jet_MW0_h.reset(new JetHists(ctx, "jet_MW0"));
  muon_MW0_h.reset(new MuonHists(ctx, "muon_MW0"));
  event_MW0_h.reset(new EventHists(ctx, "event_MW0"));
  chi2min_MW0_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_MW0",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_MW0_h.reset(new EffHists(ctx, "tagger_MW0",ttbar_hyps_label,ttbar_chi2_label));

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
  topjet_chi2cut_w_h.reset(new TopJetHists(ctx, "topjet_chi2cut_W",4,"WTopJet"));
  topjet_chi2cut_other_h.reset(new TopJetHists(ctx, "topjet_chi2cut_other",4,"notused_topjets"));
  eff_chi2cut_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut"));
  jet_chi2cut_comb_h.reset(new JetHists(ctx, "jet_chi2cut_comb"));
  jet_chi2cut_btag_h.reset(new JetHists(ctx, "jet_chi2cut_btag",4,"bjets"));
  jet_chi2cut_other_h.reset(new JetHists(ctx, "jet_chi2cut_other",4,"notused_jets"));
  muon_chi2cut_h.reset(new MuonHists(ctx, "muon_chi2cut"));
  event_chi2cut_h.reset(new EventHists(ctx, "event_chi2cut"));
  chi2min_chi2cut_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_chi2cut_h.reset(new EffHists(ctx, "tagger_chi2cut",ttbar_hyps_label,ttbar_chi2_label));

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

  //wrong W and right W
  topjet_wrong_comb_h.reset(new TopJetHists(ctx, "topjet_wrong_comb"));
  topjet_wrong_w_h.reset(new TopJetHists(ctx, "topjet_wrong_W",4,"WTopJet"));
  topjet_wrong_other_h.reset(new TopJetHists(ctx, "topjet_wrong_other",4,"notused_topjets"));
  jet_wrong_comb_h.reset(new JetHists(ctx, "jet_wrong_comb"));
  jet_wrong_btag_h.reset(new JetHists(ctx, "jet_wrong_btag",4,"bjets"));
  jet_wrong_other_h.reset(new JetHists(ctx, "jet_wrong_other",4,"notused_jets"));
  eff_wrong_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_wrong"));
  muon_wrong_h.reset(new MuonHists(ctx, "muon_wrong"));
  event_wrong_h.reset(new EventHists(ctx, "event_wrong"));
  tagger_wrong_h.reset(new EffHists(ctx, "tagger_wrong",ttbar_hyps_label,ttbar_chi2_label));
  chi2min_wrong_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_wrong",ttbar_hyps_label,ttbar_chi2_label ));

  topjet_right_comb_h.reset(new TopJetHists(ctx, "topjet_right_comb"));
  topjet_right_w_h.reset(new TopJetHists(ctx, "topjet_right_W",4,"WTopJet"));
  topjet_right_other_h.reset(new TopJetHists(ctx, "topjet_right_other",4,"notused_topjets"));
  jet_right_comb_h.reset(new JetHists(ctx, "jet_right_comb"));
  jet_right_btag_h.reset(new JetHists(ctx, "jet_right_btag",4,"bjets"));
  jet_right_other_h.reset(new JetHists(ctx, "jet_right_other",4,"notused_jets"));
  eff_right_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_right"));
  muon_right_h.reset(new MuonHists(ctx, "muon_right"));
  event_right_h.reset(new EventHists(ctx, "event_right"));
  tagger_right_h.reset(new EffHists(ctx, "tagger_right",ttbar_hyps_label,ttbar_chi2_label));
  chi2min_right_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_right",ttbar_hyps_label,ttbar_chi2_label ));

  //Hist chi2cut20
  topjet_chi2cut20_h.reset(new TopJetHists(ctx, "topjet_chi2cut20"));
  eff_chi2cut20_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut20"));
  jet_chi2cut20_h.reset(new JetHists(ctx, "jet_chi2cut20"));
  muon_chi2cut20_h.reset(new MuonHists(ctx, "muon_chi2cut20"));
  event_chi2cut20_h.reset(new EventHists(ctx, "event_chi2cut20"));
  chi2min_chi2cut20_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut20",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_chi2cut20_h.reset(new EffHists(ctx, "tagger_chi2cut20",ttbar_hyps_label,ttbar_chi2_label));
 

//Hist chi2cut10
  topjet_chi2cut10_h.reset(new TopJetHists(ctx, "topjet_chi2cut10"));
  eff_chi2cut10_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_chi2cut10"));
  jet_chi2cut10_h.reset(new JetHists(ctx, "jet_chi2cut10"));
  muon_chi2cut10_h.reset(new MuonHists(ctx, "muon_chi2cut10"));
  event_chi2cut10_h.reset(new EventHists(ctx, "event_chi2cut10"));
  chi2min_chi2cut10_h.reset(new ZPrimeTotTPrimeHypothesisHists(ctx, "chi2min_chi2cut10",ttbar_hyps_label,ttbar_chi2_label ));
  tagger_chi2cut10_h.reset(new EffHists(ctx, "tagger_chi2cut10",ttbar_hyps_label,ttbar_chi2_label));

  //general
  filename =  ctx.get("dataset_version"); 
  h_ht = ctx.get_handle<double>("HT");
  h_btag_medium = ctx.declare_event_output< std::vector<Jet> > ("BTag_medium");
  h_btag_loose = ctx.declare_event_output< std::vector<Jet> > ("BTag_loose");
  h_btag_tight = ctx.declare_event_output< std::vector<Jet> > ("BTag_tight");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");

 

  // jet_btag_cleaner.reset(new JetCleaner(ctx,JetId(ZPrimeTotTPrimeBJets(ctx,ttbar_hyps_label))));
  // topjet_W_cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut( 0, 2.4)),"WTopJet"));

  //Trigger
  const std::string& trigger = ctx.get("trigger", "NULL");
  if(trigger != "NULL") trigger_sel = make_unique<TriggerSelection>(trigger);
  else                  trigger_sel = make_unique<TriggerSelection>("HLT_Mu45_eta2p1_v*");

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
   
 
  //common Modules
  /* luminosity sections from CMS silver-JSON file */
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
      if(deltar < 0.8) bdeltaR=false;
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
    if(deltar_ak4hyp > 0.8) AK4Jets ->push_back(ak4);
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

  // ///////////////////////////////////// hists for wrong and right W //////////////////////////////
  // if(event.is_valid(h_ttbargen)){
  //   const auto & ttbargen = event.get(h_ttbargen);
  //   if(ttbargen.IsSemiLeptonicDecay() ){
  //     double deltar_W = deltaR(ttbargen.WHad(),hyp->W_v4());
  //     if(deltar_W > 0.8){
  // 	topjet_wrong_h->fill(event);
  // 	tagger_wrong_h->fill(event);
  // 	jet_wrong_h->fill(event);
  // 	eff_wrong_h->fill(event);
  // 	muon_wrong_h->fill(event);
  // 	event_wrong_h->fill(event);
  //     }else{
  // 	topjet_right_h->fill(event);
  // 	tagger_right_h->fill(event);
  // 	jet_right_h->fill(event);
  // 	eff_right_h->fill(event);
  // 	muon_right_h->fill(event);
  // 	event_right_h->fill(event);
  //     }
  //   }
  // }
 
  //  double mass_W=0;
  //  LorentzVector subjet_sum;
  //  for (const auto s : hyp->W_subjets()) {
  //    subjet_sum += s.v4();
  //  }
  //  mass_W=subjet_sum.M();
  //  if(mass_W<40){
  // 	topjet_MW0_h->fill(event);
  // 	tagger_MW0_h->fill(event);
  // 	jet_MW0_h->fill(event);
  // 	eff_MW0_h->fill(event);
  // 	muon_MW0_h->fill(event);
  // 	event_MW0_h->fill(event);
  //  }

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

  if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// chi2cut Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_chi2cut = chi2cut_sel->passes(event);
  if(!pass_chi2cut) return false;

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

  if(berror)   std::cout<<"SelectionModule L:294 vor njet3 Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// njet3 Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_njet3 = njet3_sel->passes(event);
  if(!pass_njet3) return false;
  
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
 
  if(berror)   std::cout<<"SelectionModule L:294 vor hists wrong anf right"<<std::endl;  
  ///////////////////////////////////// hists for wrong and right W //////////////////////////////
  hyps = event.get(h_ttbar_hyps);
  hyp = get_best_hypothesis(hyps, "Chi2");
  if(!hyp) std::runtime_error("ZprimeTotTPrimeSelectionModule::process -- best hypothesis for reconstruction not found");
  const ZPrimeTotTPrimeReconstructionHypothesis hyp_obj2(*hyp);

 
  if(event.is_valid(h_ttbargen)){
    const auto & ttbargen = event.get(h_ttbargen);
    if(ttbargen.IsSemiLeptonicDecay() ){

      double delta_q1 = deltaR(hyp->W_v4(), ttbargen.Q1());
      double delta_q2 = deltaR(hyp->W_v4(), ttbargen.Q2());
      if(delta_q1 > 0.8 && delta_q2 >0.8){
	topjet_wrong_comb_h->fill(event);
	topjet_wrong_w_h->fill(event);
	topjet_wrong_other_h->fill(event);
	jet_wrong_comb_h->fill(event);
	jet_wrong_btag_h->fill(event);
	jet_wrong_other_h->fill(event);
	eff_wrong_h->fill(event);
	muon_wrong_h->fill(event);
	event_wrong_h->fill(event);
	chi2min_wrong_h->fill(event);
	tagger_wrong_h->fill(event);
 
      }else{
	topjet_right_comb_h->fill(event);
	topjet_right_w_h->fill(event);
	topjet_right_other_h->fill(event);
	jet_right_comb_h->fill(event);
	jet_right_btag_h->fill(event);
	jet_right_other_h->fill(event);
	eff_right_h->fill(event);
	muon_right_h->fill(event);
	event_right_h->fill(event);
	chi2min_right_h->fill(event);
	tagger_right_h->fill(event);
      }
    }
  }

 if(berror)   std::cout<<"SelectionModule L:294 vor M<40"<<std::endl;  
 
   double mass_W=0;
   LorentzVector subjet_sum;
   for (const auto s : hyp->W_subjets()) {
     subjet_sum += s.v4();
   }
   mass_W=subjet_sum.M();
   if(mass_W<40){
	topjet_MW0_h->fill(event);
	tagger_MW0_h->fill(event);
	jet_MW0_h->fill(event);
	eff_MW0_h->fill(event);
	muon_MW0_h->fill(event);
	event_MW0_h->fill(event);
   }
   hyps.clear();
   hyps.push_back(hyp_obj2);

   ///////////////////////////////////////////////////////////////////////////////////////
 if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut 20 Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// chi2cut 20 Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_chi2cut20 = chi2cut20_sel->passes(event);
  if(!pass_chi2cut20) return false;

  topjet_chi2cut20_h->fill(event);
  jet_chi2cut20_h->fill(event);
  eff_chi2cut20_h->fill(event);
  muon_chi2cut20_h->fill(event);
  event_chi2cut20_h->fill(event);
  chi2min_chi2cut20_h->fill(event);
  tagger_chi2cut20_h->fill(event);

if(berror)   std::cout<<"SelectionModule L:294 vor chi2cut 10 Selection"<<std::endl;  
  ////////////////////////////////////////////////////////// chi2cut 10 Selction //////////////////////////////////////////////////////////////////////////////////
  bool pass_chi2cut10 = chi2cut10_sel->passes(event);
  if(!pass_chi2cut10) return false;

  topjet_chi2cut10_h->fill(event);
  jet_chi2cut10_h->fill(event);
  eff_chi2cut10_h->fill(event);
  muon_chi2cut10_h->fill(event);
  event_chi2cut10_h->fill(event);
  chi2min_chi2cut10_h->fill(event);
  tagger_chi2cut10_h->fill(event);

 ////////////////////////////////////////////////////////// Lumi Hists //////////////////////////////////////////////////////////////////////////////////
  lumi_h->fill(event);  

  return true;
}


UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeEffModule)
