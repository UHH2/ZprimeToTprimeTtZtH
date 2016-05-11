#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenHists.h"
int ZPrimeMass = 1500;
bool  throw_failure = false;

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  /*
   *ZPrime Generator-Level-Module 
   */
  class ZPrimeTotTPrimeGenModule: public AnalysisModule {
  public:
    
    explicit ZPrimeTotTPrimeGenModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> printer;
    std::unique_ptr<AnalysisModule> zprimegenprod;
    std::unique_ptr<Hists> h_zprimegenhists;
    std::unique_ptr<Hists> h_zprimegenhists_test;
    Event::Handle<ZPrimeGen> h_zprimegen;
  };


  ZPrimeTotTPrimeGenModule::ZPrimeTotTPrimeGenModule(Context & ctx){
    printer.reset(new GenParticlesPrinter(ctx));
    zprimegenprod.reset(new ZPrimeGenProducer(ctx, "zprimegen", false));
    h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
    h_zprimegenhists.reset(new ZPrimeGenHists(ctx, "zprimegenhists"));
    h_zprimegenhists_test.reset(new ZPrimeGenHists(ctx, "test"));
  }


  bool ZPrimeTotTPrimeGenModule::process(Event & event) {
    //Prints a table of all particles
    printer->process(event);
    
//printer->process(event);
    //collects all particles of an event
    zprimegenprod->process(event);

    //const auto & zprimegen = event.get(h_zprimegen);
    // cout << zprimegen.GetDecayChannel(); << endl;

    // if(zprimegen.TPrime_Mass()<1000){
    //   if(zprimegen.invMass2()>= 1000){
    // 	// printer->process(event);
    //   h_zprimegenhists_test->fill(event);
    //   }
    // }

    // if(zprimegen.ZPrime_invMass()<2350 && zprimegen.ZPrime_invMass() ){
    //   if(throw_failure){
    // 	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$  neues Event  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
    //   }

    //   printer->process(event);
    //   h_zprimegenhists_test->fill(event);
    // }

    //Fills all hists of ZPrimeTotTPrimeGenHists.cxx
    h_zprimegenhists->fill(event);

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeGenModule)

}
