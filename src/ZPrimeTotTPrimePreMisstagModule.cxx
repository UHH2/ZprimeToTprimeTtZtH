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
  bool isMC;
  const int runnr_BCD = 276811;
  const int runnr_EF = 278802;
  const int runnr_G = 280385;

};

ZPrimeTotTPrimePreMisstagModule::ZPrimeTotTPrimePreMisstagModule(uhh2::Context& ctx){ 

  //GenParticleprinter
  printer.reset(new GenParticlesPrinter(ctx));
 
  isMC = (ctx.get("dataset_type") == "MC");
  //corrector
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
 
  //// Data/MC scale
   if(isMC){ 
    pileup_SF.reset(new MCPileupReweight(ctx)); 
    lumiweight.reset(new MCLumiWeight(ctx));
    jet_resolution_smearer.reset(new JetResolutionSmearer(ctx));
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
  ht_sel.reset(new HtJetsSelection(1000,-1));
  
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

  ///correctors
  if(isMC){
    jet_corrector->process(event);
    topjet_corrector->process(event);
    subjet_corrector->process(event);
    jetlepton_cleaner->process(event);

    jet_resolution_smearer->process(event);
  }else{
    if(event.run <= runnr_BCD)  {       
      jet_corrector_BCD->process(event);
      topjet_corrector_BCD->process(event);
      subjet_corrector_BCD->process(event);
      jetlepton_cleaner_BCD->process(event);
    }
    else if(event.run < runnr_EF){       
      jet_corrector_EF->process(event);
      topjet_corrector_EF->process(event);
      subjet_corrector_EF->process(event);
      jetlepton_cleaner_EF->process(event);
    } 
    else if(event.run <= runnr_G) {       
      jet_corrector_G->process(event);
      topjet_corrector_G->process(event);
      subjet_corrector_G->process(event);
      jetlepton_cleaner_G->process(event);
    } 
    else if(event.run > runnr_G) {       
      jet_corrector_H->process(event);
      topjet_corrector_H->process(event);
      subjet_corrector_H->process(event);
      jetlepton_cleaner_H->process(event);
    } 
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
  ti_HT=event.get_trigger_index("HLT_PFHT900_v*");
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


