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

#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHists.h>
//#include <UHH2/ZprimeToTprimeTtZtH/include/MistagHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h>


using namespace uhh2examples;
using namespace uhh2;
using namespace std;

class ZPrimeTotTPrimePreMisstagModule: public uhh2::AnalysisModule {
  
public:
  explicit ZPrimeTotTPrimePreMisstagModule(uhh2::Context&);
  virtual bool process(uhh2::Event&) override;
  
private:
  enum lepton { muon, elec };
  lepton channel_;
  std::unique_ptr<AnalysisModule> printer;
  
  //cleaner
  std::unique_ptr<TopJetCleaner>      topjet_cleaner;
  std::unique_ptr<TopJetLeptonDeltaRCleaner> topjetlepton_cleaner; 
  
  //calculators
  std::vector<std::unique_ptr<AnalysisModule>> htcalc;
  std::vector<std::unique_ptr<AnalysisModule>> metfilters;
  
  //correctors
  std::unique_ptr<SubJetCorrector> subjetcorrector;
  std::unique_ptr<JetCorrector>     jet_corrector;
  std::unique_ptr<JetResolutionSmearer> jetER_smearer;
  std::unique_ptr<TopJetCorrector>           topjet_corrector;
 
  // Data/MC scale factors
  std::unique_ptr<uhh2::AnalysisModule> pileup_SF;
  std::unique_ptr<uhh2::AnalysisModule> lumiweight;
  std::unique_ptr<uhh2::Selection> lumi_sel;

  //Selections
  std::unique_ptr<uhh2::Selection> muo_sel; //veto on leptons
  std::unique_ptr<uhh2::Selection> ele_sel; //veto on leptons
  std::unique_ptr<uhh2::Selection> topjet2_sel; //two AK8 pT>400GeV eta<2.4 deltaPhi>2.1
  std::unique_ptr<uhh2::Selection> ht_sel; //Ht,jets >1000

  //Hists
  //input
  std::unique_ptr<uhh2::Hists> input_h_event;
  std::unique_ptr<uhh2::Hists> input_h_muo;
  std::unique_ptr<uhh2::Hists> input_h_ele;
  std::unique_ptr<uhh2::Hists> input_h_jet;
  std::unique_ptr<Hists> input_h_eff;
  std::unique_ptr<Hists>input_h_topjet;
  std::unique_ptr<Hists> input_h_mistag;

  //lepton sel
  std::unique_ptr<Hists> eff_lepton_h;
  std::unique_ptr<Hists> jet_lepton_h;
  std::unique_ptr<Hists> muon_lepton_h;
  std::unique_ptr<Hists> event_lepton_h;
  std::unique_ptr<Hists> topjet_lepton_h;
  std::unique_ptr<Hists> mistag_lepton_h;
 
//Hist Topjet2
  std::unique_ptr<Hists> eff_topjet2_h;
  std::unique_ptr<Hists> jet_topjet2_h;
  std::unique_ptr<Hists> muon_topjet2_h;
  std::unique_ptr<Hists> event_topjet2_h;
  std::unique_ptr<Hists> topjet_topjet2_h;
  std::unique_ptr<Hists> tagger_topjet2_h;
 

  //output
  std::unique_ptr<uhh2::Hists> output_h_event;
  std::unique_ptr<uhh2::Hists> output_h_muo;
  std::unique_ptr<uhh2::Hists> output_h_ele;
  std::unique_ptr<uhh2::Hists> output_h_jet;
  std::unique_ptr<Hists> output_h_eff;
  std::unique_ptr<Hists>output_h_topjet;
  std::unique_ptr<Hists> output_h_mistag;
  
  //General
  std::string filename;
  // Reconstruction TTBar for Background
  std::unique_ptr<uhh2::Selection> genmttbar_sel;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;

};

ZPrimeTotTPrimePreMisstagModule::ZPrimeTotTPrimePreMisstagModule(uhh2::Context& ctx){ 

  //GenParticleprinter
  printer.reset(new GenParticlesPrinter(ctx));

  //corrector
  std::vector<std::string> JEC_AK4, JEC_AK8;
  const bool isMC = (ctx.get("dataset_type") == "MC");
  if(isMC){

    JEC_AK4 = JERFiles::Summer15_25ns_L123_AK4PFchs_MC;
    JEC_AK8 = JERFiles::Summer15_25ns_L123_AK8PFchs_MC;
  }
  else {

    JEC_AK4 = JERFiles::Summer15_25ns_L123_AK4PFchs_DATA;
    JEC_AK8 = JERFiles::Summer15_25ns_L123_AK8PFchs_DATA;
  }
  jet_corrector.reset(new JetCorrector(ctx, JEC_AK4));
  jetER_smearer.reset(new JetResolutionSmearer(ctx));
  topjet_corrector.reset(new TopJetCorrector(ctx, JEC_AK8));
  if(isMC) subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_MC));
  else subjetcorrector.reset(new SubJetCorrector(ctx,JERFiles::Fall15_25ns_L123_AK4PFchs_DATA));
 
 
  //// Data/MC scale
   if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx));
   } else     lumi_sel.reset(new LumiSelection(ctx));

   filename =  ctx.get("dataset_version");

  //calculator
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new PrimaryLepton(ctx)));
  htcalc.push_back(std::unique_ptr<AnalysisModule>(new HTlepCalculator(ctx)));

  //cleaner
  topjet_cleaner.reset(new TopJetCleaner(ctx,TopJetId(PtEtaCut(200., 2.5))));

  //Selections
  muo_sel.reset(new NMuonSelection(0,0)); 
  ele_sel.reset(new NElectronSelection(0,0));
  topjet2_sel.reset(new NTopJetSelection(2,-1,TopJetId(PtEtaCut( 250., 2.4)))); 
  ht_sel.reset(new HtSelection(1000,-1));
  
  //Genral
  const std::string ttbar_gen_label ("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  // Select of the inclusiv ttbar sample only events from 0 to 700 GeV
  if(ctx.get("dataset_version") == "TTbarAll"){ genmttbar_sel.reset(new GenMttbarCut(ctx, 0., 700., ttbar_gen_label));}
  else   genmttbar_sel.reset(new uhh2::AndSelection(ctx));

  //Hist
  input_h_event.reset(new EventHists   (ctx, "input_Event"));
  input_h_muo.reset(new MuonHists (ctx, "input_Muons"));
  input_h_ele.reset(new ElectronHists(ctx, "input_Electrons"));
  input_h_jet.reset(new JetHists     (ctx, "input_Jets"));
  input_h_eff.reset(new ZPrimeTotTPrimeHists(ctx, "input_eff"));
  input_h_topjet.reset(new TopJetHists     (ctx, "input_TopJets"));

  //Hists lepton
  topjet_lepton_h.reset(new TopJetHists(ctx, "topjet_lepton"));
  eff_lepton_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_lepton"));
  jet_lepton_h.reset(new JetHists(ctx, "jet_lepton"));
  muon_lepton_h.reset(new MuonHists(ctx, "muon_lepton"));
  event_lepton_h.reset(new EventHists(ctx, "event_lepton"));
 

  //Hists topjet2
  topjet_topjet2_h.reset(new TopJetHists(ctx, "topjet_topjet2"));
  eff_topjet2_h.reset(new ZPrimeTotTPrimeHists(ctx, "eff_topjet2"));
  jet_topjet2_h.reset(new JetHists(ctx, "jet_topjet2"));
  muon_topjet2_h.reset(new MuonHists(ctx, "muon_topjet2"));
  event_topjet2_h.reset(new EventHists(ctx, "event_topjet2"));
 


  output_h_event.reset(new EventHists   (ctx, "output_Event"));
  output_h_muo.reset(new MuonHists (ctx, "output_Muons"));
  output_h_ele.reset(new ElectronHists(ctx, "output_Electrons"));
  output_h_jet.reset(new JetHists     (ctx, "output_Jets"));
  output_h_eff.reset(new ZPrimeTotTPrimeHists(ctx, "output_eff"));
  output_h_topjet.reset(new TopJetHists     (ctx, "output_TopJets"));

}

bool ZPrimeTotTPrimePreMisstagModule::process(Event & event) {
  //  std::cout<<"In PreMisstagModule: Am Anfang"<<std::endl;
  for (auto & mod : htcalc) {
    mod->process(event);
  }

  //correctors
  topjet_corrector->process(event);
  jet_corrector->process(event);
  subjetcorrector->process(event);

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
  //  std::cout<<"In PreMisstagModule: vor Input"<<std::endl;
  input_h_event ->fill(event);
  input_h_muo ->fill(event);
  input_h_ele   ->fill(event);
  input_h_jet   ->fill(event);
  input_h_eff ->fill(event);
  input_h_topjet   ->fill(event);

  //std::cout<<"In PreMisstagModule: vor Input"<<std::endl;
  bool pass_lepton = (ele_sel->passes(event))&&(muo_sel->passes(event));
  if(!pass_lepton)return false;
  muon_lepton_h->fill(event);
  eff_lepton_h ->fill(event);
  event_lepton_h->fill(event);
  topjet_lepton_h->fill(event);
  jet_lepton_h->fill(event);
 
  //std::cout<<"In PreMisstagModule: vor topjet2"<<std::endl;
  bool pass_topjet = topjet2_sel->passes(event);
  if(!pass_topjet)return false;
  muon_topjet2_h->fill(event);
  eff_topjet2_h ->fill(event);
  event_topjet2_h->fill(event);
  topjet_topjet2_h->fill(event);
  jet_topjet2_h->fill(event);

  // std::cout<<"In PreMisstagModule: vor HT"<<std::endl;
  bool pass_ht = ht_sel->passes(event);
  if(!pass_ht)return false;

  output_h_event ->fill(event);
  output_h_muo ->fill(event);
  output_h_ele   ->fill(event);
  output_h_jet   ->fill(event);
  output_h_eff ->fill(event);
  output_h_topjet   ->fill(event);

  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimePreMisstagModule)


