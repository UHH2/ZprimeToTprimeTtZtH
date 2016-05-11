#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/ElectronHists.h"
#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/MuonHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeSelections.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHists.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

/** \brief Basic analysis example of an AnalysisModule (formerly 'cycle') in UHH2
 * 
 * This is the central class which calls other AnalysisModules, Hists or Selection classes.
 * This AnalysisModule, in turn, is called (via AnalysisModuleRunner) by SFrame.
 */
class ZPrimeTotTPrimeModule: public AnalysisModule {
public:
    
    explicit ZPrimeTotTPrimeModule(Context & ctx);
    virtual bool process(Event & event) override;

private:
    
    std::unique_ptr<CommonModules> common;
    
    std::unique_ptr<JetCleaner> jetcleaner;
   
    // declare the Selections to use. Use unique_ptr to ensure automatic call of delete in the destructor,
    // to avoid memory leaks.
  std::unique_ptr<Selection> njet2_sel, njet3_sel, njet4_sel , dijet_sel, MET_sel, PtJet_cut, btag_sel, Ht_cut;
    
    // store the Hists collection as member variables. Again, use unique_ptr to avoid memory leaks.
  // input hists
  std::unique_ptr<Hists> h_nocuts, input_h_ele, input_h_muon, input_h_event, input_h_jet;

  // selection hists
  std::unique_ptr<Hists> h_njet2, h_dijet, h_ele, h_muon, h_event_met, h_ptjet0,h_btag,h_event_ht;

  std::vector<std::unique_ptr<AnalysisModule>> pre_modules;

  //output hists
  // std::unique_ptr<Hists> output_h_muon;
};


ZPrimeTotTPrimeModule::ZPrimeTotTPrimeModule(Context & ctx){
    // In the constructor, the typical tasks are to initialize the
    // member variables, in particular the AnalysisModules such as
    // CommonModules or some cleaner module, Selections and Hists.
    // But you can do more and e.g. access the configuration, as shown below.
    
    cout << "Hello World from ZPrimeTotTPrimeModule!" << endl;
    
    // If needed, access the configuration of the module here, e.g.:
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    // If running in SFrame, the keys "dataset_version", "dataset_type", "dataset_lumi",
    // and "target_lumi" are set to the according values in the xml file. For CMSSW, these are
    // not set automatically, but can be set in the python config file.
    for(auto & kv : ctx.get_all()){
        cout << " " << kv.first << " = " << kv.second << endl;
    }
    
    // 1. setup other modules. CommonModules and the JetCleaner:
    common.reset(new CommonModules());
    // TODO: configure common here, e.g. by 
    // calling common->set_*_id or common->disable_*
    common->init(ctx);
    //jetcleaner.reset(new JetCleaner(ctx,30.0, 2.4)); 


    pre_modules.push_back(std::unique_ptr<AnalysisModule>(new HTCalculator(ctx)));
    
    // note that the JetCleaner is only kept for the sake of example;
    // instead of constructing a jetcleaner explicitly,
    // the cleaning can also be achieved with less code via CommonModules with:
    // common->set_jet_id(PtEtaCut(30.0, 2.4));
    // before the 'common->init(ctx)' line.
    
    // 2. set up selections
    njet2_sel.reset(new NJetSelection(2)); // see common/include/NSelections.h
    njet3_sel.reset(new NJetSelection(3)); // see common/include/NSelections.h
    njet4_sel.reset(new NJetSelection(4)); // see common/include/NSelections.h
    dijet_sel.reset(new ZPrimeTotTPrimeDijetSelection()); // see ZPrimeTotTPrimeSelections
    MET_sel.reset(new ZPrimeTotTPrimeMETCut(100)); // see ZPrimeTotTPrimeSelections
    Ht_cut.reset(new HtSelection(100,1000)); // see ZPrimeTotTPrimeSelections
    PtJet_cut.reset(new ZPrimeTotTPrimePtJetCut(0,400,1)); // see ZPrimeTotTPrimeSelections
    btag_sel.reset(new NBTagSelection(2,100, btagging::csv_wp::medium )); // see ZPrimeTotTPrimeSelections

    // 3. Set up Hists classes:
    h_nocuts.reset(new ZPrimeTotTPrimeHists(ctx, "NoCuts"));
    input_h_ele.reset(new ElectronHists(ctx, "Input_Ele"));
    input_h_muon.reset(new MuonHists(ctx, "Input_moun"));
    input_h_event.reset(new EventHists(ctx, "Input_event"));
    input_h_jet.reset(new JetHists(ctx, "Input_jet"));

    h_njet2.reset(new JetHists(ctx, "Njet"));
    h_dijet.reset(new ZPrimeTotTPrimeHists(ctx, "Dijet"));
    h_ele.reset(new ElectronHists(ctx, "ele_nocuts"));
    h_muon.reset(new MuonHists(ctx, "moun_cut"));
    h_event_met.reset(new EventHists(ctx, "MET_cut"));
    h_event_ht.reset(new ZPrimeTotTPrimeHists(ctx, "HT_cut"));
    h_ptjet0.reset(new JetHists(ctx, "PT_Cut_Jet0"));
    h_btag.reset(new JetHists(ctx, "Btag"));

    //output_h_muon.reset(new JetHists(ctx, "Output_moun"));
}


bool ZPrimeTotTPrimeModule::process(Event & event) {
    // This is the main procedure, called for each event. Typically,
    // do some pre-processing by calling the modules' process method
    // of the modules constructed in the constructor (1).
    // Then, test whether the event passes some selection and -- if yes --
    // use it to fill the histograms (2).
    // Finally, decide whether or not to keep the event in the output (3);
    // this is controlled by the return value of this method: If it
    // returns true, the event is kept; if it returns false, the event
    // is thrown away.
    
    //cout << "ZPrimeTotTPrimeModule: Starting to process event (runid, eventid) = (" << event.run << ", " << event.event << "); weight = " << event.weight << endl;
    
    // 1. run all modules other modules.
    common->process(event);
    jetcleaner->process(event);


    for (auto & mod : pre_modules) {
      mod->process(event);
    }


    //    double weight = event.weight;
    
    // 2. test selections and fill histograms
    input_h_ele->fill(event);
    h_nocuts->fill(event);
    input_h_muon->fill(event);
    input_h_event->fill(event);
    input_h_jet->fill(event);

    bool njet_selection = njet2_sel->passes(event);
    if(njet_selection){
        h_njet2->fill(event);
    }
    bool dijet_selection = dijet_sel->passes(event);
    if(dijet_selection){
        h_dijet->fill(event);
    }

    bool MET_selection = MET_sel -> passes(event);
    if(MET_selection){
      h_event_met -> fill(event);
    }

    bool PtJet_Selection = PtJet_cut->passes(event);
    if (PtJet_Selection){
      h_ptjet0->fill(event);
    }

    bool btag_selection = btag_sel ->passes(event);
    if(btag_selection){
      h_btag -> fill(event);
    }

    bool Ht_selection = Ht_cut -> passes(event);
    if(Ht_selection){
      h_event_ht->fill(event);
    }

    //  output_h_muon->fill(event);
    // 3. decide whether or not to keep the current event in the output:
    //return njet_selection && dijet_selection;
    return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ZPrimeTotTPrimeModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeModule)

}
