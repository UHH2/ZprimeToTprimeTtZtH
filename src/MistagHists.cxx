#include "UHH2/core/include/LorentzVector.h"
#include "UHH2/core/include/Utils.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/MistagHists.h"
#include "UHH2/common/include/TTbarGen.h"


using namespace std;
using namespace uhh2;

//Hists for studies to the eff. of the W/Z tagger in MC and Data


MistagHists::MistagHists(Context & ctx, const string & dirname ): Hists(ctx, dirname){

 
  
 
  book<TH1F>("match_genp_Topjet","Match of genparticles to Topjets ",3,0,3);
  hist("match_genp_Topjet")->Fill("matched",0);
  hist("match_genp_Topjet")->Fill("missed",0);
  hist("match_genp_Topjet")->Fill("#topjets",0);


}//constructor definition


void MistagHists::fill(const Event & event){

  assert(event.topjets);
  double weight = event.weight;
    
  
 for(topjet:*event.topjets){
   hist("match_genp_Topjet")->Fill("#topjets",weight);
   for(genp:*event.genparticles){
     double deltar = deltaR(topjet,genp);
     if(deltar<0.8) {
       hist("match_genp_Topjet")->Fill("matched",weight);
     }else{
       hist("match_genp_Topjet")->Fill("missed",weight);
     }
   }
 }
 
  
}//fill function
