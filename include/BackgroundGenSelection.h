#pragma once 


#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

extern bool throw_failure;

class BackgroundGen{
 public:
  explicit BackgroundGen(const std::vector<GenParticle> & genparts);


  LorentzVector Quark1_v4()const;
  LorentzVector Quark2_v4()const;

 private:
  GenParticle Quark1,Quark2;


};



class BackgroundGenProducer: public uhh2::AnalysisModule {
public:
    explicit BackgroundGenProducer(uhh2::Context & ctx, const std::string & name = "backgroundgen", bool throw_on_failure = true);
     virtual bool process(uhh2::Event & event) override;
   
 
private:
    uhh2::Event::Handle<BackgroundGen> h_background;
    //  uhh2::Event::Handle<std::vector<int>> h_n;
    bool throw_on_failure;
    

};
