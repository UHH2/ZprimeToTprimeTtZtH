#include <iostream>
#include <memory>
#include <UHH2/core/include/AnalysisModule.h>
#include <UHH2/core/include/Event.h>
#include <UHH2/core/include/Selection.h>

#include "UHH2/common/include/EventVariables.h"
#include "UHH2/common/include/PrimaryLepton.h"
#include <UHH2/common/include/LumiSelection.h>
#include <UHH2/common/include/MCWeight.h>
#include <UHH2/common/include/CleaningModules.h>
#include <UHH2/common/include/NSelections.h>
#include <UHH2/common/include/JetCorrections.h>
#include <UHH2/common/include/ObjectIdUtils.h>
#include <UHH2/common/include/MuonIds.h>
#include <UHH2/common/include/ElectronIds.h>
#include <UHH2/common/include/JetIds.h>
#include <UHH2/common/include/EventHists.h>
#include <UHH2/common/include/MuonHists.h>
#include <UHH2/common/include/ElectronHists.h>
#include <UHH2/common/include/JetHists.h>
#include <UHH2/common/include/Utils.h>

#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimePreSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeSelections.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHists.h>
#include <UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h>


#include "UHH2/common/include/TopPtReweight.h"

using namespace uhh2examples;
using namespace uhh2;

class ZPrimeTotTPrimePreSelectionModule_toppt: public uhh2::AnalysisModule {

public:
  explicit ZPrimeTotTPrimePreSelectionModule_toppt(uhh2::Context& ctx);
  virtual bool process(uhh2::Event&) override;

private:
 
  std::unique_ptr<uhh2::AnalysisModule> topptreweighting_all;
  std::unique_ptr<Hists> topptreweighting_h;
  std::unique_ptr<uhh2::AnalysisModule> ttgenprod;


};

ZPrimeTotTPrimePreSelectionModule_toppt::ZPrimeTotTPrimePreSelectionModule_toppt(uhh2::Context& ctx){

  //TTbar MassSelection 
  const std::string ttbar_gen_label ("ttbargen");
  ttgenprod.reset(new TTbarGenProducer(ctx, ttbar_gen_label, false));
  topptreweighting_h.reset(new TopPtReweightHist(ctx,"topptreweighting","weight_ttbar"));
  topptreweighting_all.reset(new TopPtReweight(ctx, 0.0615 , -0.0005 ,"ttbargen","weight_ttbar", false,1));





}

bool ZPrimeTotTPrimePreSelectionModule_toppt::process(Event & event) {
 
  

  
  ttgenprod->process(event);
  topptreweighting_all->process(event);
  topptreweighting_h->fill(event);


  return true;
}

UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimePreSelectionModule_toppt)
