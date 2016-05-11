#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/PrintingModules.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenHists.h"
bool bZPrime, bTPrime;
using namespace std;
using namespace uhh2; 
namespace uhh2examples {

  /*
   *ZPrime Generator-Level-Module 
   */
  class ZPrimeTotTPrimeGenModulePre: public AnalysisModule {
  public:
    
    explicit ZPrimeTotTPrimeGenModulePre(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<AnalysisModule> printer;
    std::unique_ptr<AnalysisModule> zprimegenprod;
    std::unique_ptr<Hists> h_zprimegenhists;
    std::unique_ptr<Hists> h_zprimegenhists_ZMass;
    std::unique_ptr<Hists> h_zprimegenhists_TMass;
    std::unique_ptr<Hists> h_zprimegenhists_ZTMass;
    Event::Handle<ZPrimeGen> h_zprimegen;
  };


  ZPrimeTotTPrimeGenModulePre::ZPrimeTotTPrimeGenModulePre(Context & ctx){
    printer.reset(new GenParticlesPrinter(ctx));
    zprimegenprod.reset(new ZPrimeGenProducer(ctx, "zprimegen", false));
    h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
    h_zprimegenhists.reset(new ZPrimeGenHists(ctx, "zprimegenhists"));
    h_zprimegenhists_TMass.reset(new ZPrimeGenHists(ctx, "zprimegenhists_TMass"));
    h_zprimegenhists_ZMass.reset(new ZPrimeGenHists(ctx, "zprimegenhists_ZMass"));
    h_zprimegenhists_ZTMass.reset(new ZPrimeGenHists(ctx, "zprimegenhists_ZTMass"));
  }

 

  bool ZPrimeTotTPrimeGenModulePre::process(Event & event) {
    // double TPrimemass;
    bZPrime = true;
    bTPrime = true;
    //collects all particles of an event
    zprimegenprod->process(event);
    h_zprimegenhists->fill(event);
    //const  auto & zprimegen = event.get(h_zprimegen);

    //double ZPrimemass = zprimegen.ZPrime_Mass();

    // // ZPrime mass 0 GeV -> No ZPrime in event
    // if(ZPrimemass < 100){
    //   //Prints a table of all particles
    //   // printer->process(event);
    //   bZPrime = false;
    //   zprimegenprod->process(event);
    //   ZPrimemass = zprimegen.ZPrime_invMass();
    //   if(ZPrimemass < 300){
    //   	// bZPrime = true;
    //   	// bTPrime = false;
    //   	// zprimegenprod->process(event);
    //   	printer->process(event);
    //   	cout << "Daughter1 " <<zprimegen.pdgDaughter1() <<endl;
    //   	cout << "Daughter2 " <<zprimegen.pdgDaughter2() <<endl;
    //   	cout << "Daughter3 " <<zprimegen.pdgDaughter3() <<endl;
    //   	cout << " daughter123  " <<zprimegen.invMass123() <<endl;
     
    //   		h_zprimegenhists_TMass->fill(event);
    //   	}
    // }


    //   // Fills all hists of ZPrimeTotTPrimeGenHists.cxx
    //   h_zprimegenhists_ZMass->fill(event);
    // }

    //No ZPrime in event and TPrime mass 0 GeV: No event in the sample
    // if(ZPrimemass < 100){
      // bZPrime = false;
      // zprimegenprod->process(event);
      //TPrimemass = zprimegen.TPrime_Mass();
      //   if(TPrimemass < 100){
    	//Prints a table of all particles
	//	printer->process(event);
	// }
      // }
     
    // // TPrime mass 0 GeV because of no TPrime in event
    // bZPrime = true;
    // bTPrime = true;
    // if(ZPrimemass > 100){
    //   TPrimemass = zprimegen.TPrime_Mass();
    //   if(TPrimemass<100){
    // 	//	printer->process(event);
    // 	bTPrime = false;
    // 	zprimegenprod->process(event);
    // 	if(zprimegen.invMass524g() < 1250){
    //   	printer->process(event);
    // 	h_zprimegenhists_ZMass->fill(event);
    // 	}

    	// int pdg = zprimegen.pdgDaughter3();
    	// if(pdg<5){
    	//   printer->process(event);
    	// }

    	//int bzpdgId = zprimegen.pdgDaughter3();
    	// if(bzpdgId < 23){
    	//   printer->process(event);
    	//   // h_zprimegenhists_TMass->fill(event);
    	// }
    	// float mass3 = zprimegen.invMass2();
    	// if(mass3 < 1206 && mass3 > 1194){
    	// h_zprimegenhists_TMass->fill(event);
    	// } 
    	// double invmass123 = zprimegen.invMass123();
    	// if(invmass123 < 1700){
    	// h_zprimegenhists_TMass->fill(event);
    	// printer->process(event);
     	// }
    // 	h_zprimegenhists_TMass->fill(event);
    //   }
    // }
    
     // bZPrime = true;
     // bTPrime = true;
     // //ZPrime and TPrime mass >100 
     // if(ZPrimemass > 100){
     //   TPrimemass = zprimegen.TPrime_Mass();
     //   if(TPrimemass > 100){
     //   	printer->process(event);
     // 	bTPrime = true;
     // 	zprimegenprod->process(event);
     // 	h_zprimegenhists_TMass->fill(event);
     //   }
     // }
  
    return true;
  }



  UHH2_REGISTER_ANALYSIS_MODULE(ZPrimeTotTPrimeGenModulePre)

}
