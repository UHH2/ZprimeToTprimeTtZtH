#pragma once 


#include "UHH2/core/include/GenParticle.h"
#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include <vector>
extern bool bZPrime, bTPrime;
extern int ZPrimeMass;
extern bool throw_failure;
extern   std::vector<int> n_vector;
extern bool b_semileptonic;

class ZPrimeGen{
 public:
  explicit ZPrimeGen(const std::vector<GenParticle> & genparts);
  void assign(GenParticle  daughter1, GenParticle  daughter2,const std::vector<GenParticle> & genparts);
  void recoTPrime(const std::vector<GenParticle> & genparts);
  void recoZPrime(GenParticle TPrime ,GenParticle  main_Top);
  //std::vector<int> DecayChannel( const std::vector<GenParticle> & genparticles);
  void  DecayChannel( const std::vector<GenParticle> & genparticles);
  std::vector<int> GetDecayChannel( );


  GenParticle Top() const;
  GenParticle ZBoson() const;
  GenParticle WBoson() const;
  GenParticle Higgs() const;
  GenParticle ZP() const;
  GenParticle TP() const;
  GenParticle D6()const;
  GenParticle D5()const;
  GenParticle D3()const;
  GenParticle D2()const;
  GenParticle D1()const;
  GenParticle G()const;
  GenParticle ATop()const;
  GenParticle BLep()const;
  GenParticle Lepton()const;
  GenParticle Neutrino()const;
  GenParticle WLep()const;
  GenParticle WHad()const;
  GenParticle WHiggsTop()const;
  GenParticle BHiggsTop()const;
  GenParticle hadTop()const;
  GenParticle lepTop()const;
 GenParticle TopW()const;
 GenParticle muon1()const;
 GenParticle muon2()const;
 GenParticle Quark1()const;
 GenParticle Quark2()const;
 GenParticle BHad()const;
 GenParticle WHadQ1()const;
 GenParticle WHadQ2()const;
 GenParticle ZQ1()const;
 GenParticle ZQ2()const;
 GenParticle HQ1()const;
 GenParticle HQ2()const;
 GenParticle WQ1()const;
 GenParticle WQ2()const;

GenParticle WHiggsTopD1()const;

GenParticle WHiggsTopD2()const;

 GenParticle WTopD1()const;
 GenParticle WTopD2()const;
 GenParticle Btop()const;

  /* GenParticle WHiggsTopD(int i)const; */
  /* GenParticle HiggsD(int i) const; */
  /* GenParticle ZBosonD(int i) const; */
  /* GenParticle WBosonD(int i) const; */

  int pdgDaughter1()const;
  int pdgDaughter2()const;
  int pdgDaughter3()const;
  int pdgDaughter5()const;
  int pdgDaughter6()const;
  
  float Main_Top_PT()const;
  float TPrime_PT()const;
  float ZPrime_PT()const;
  float Daughter3_PT()const;

  float Higgs_PT()const;
  float TopHiggs_PT()const;
  float TopZ_PT()const;
  float WTPrime_PT()const;
  float WTop_PT()const;
  float BTop_PT()const;
  float BTPrime_PT()const;
  float Z_PT()const;


  float ZPrime_Mass()const;
  float TPrime_Mass()const;
  float ZPrime_invMass()const;
  float invMass1()const;
  float invMass2()const;
  float invMass3()const;
  /* float invMass123()const; */
  /* float invMass524()const; */
  /* float invMass621()const; */
  /* float invMass524g()const; */
  /* float invMass5624g()const; */

  float TPrimeTop_Phi()const;
  float HT_Phi()const;
  float ZT_Phi()const;
  float TPrimeWB_Phi()const;
  float TopWB_Phi()const;

  float TPrime_Eta()const;
 float Top_Eta()const;
 float HT_Eta()const;
 float ZT_Eta()const;
 float TPrimeWB_Eta()const;
 float TopWB_Eta()const;

  float TPrimeTop_R()const;
  float HT_R()const;
  float ZT_R()const;
  float TPrimeWB_R()const;
  float TopWB_R()const;
  float TPrimeWG_R()const;
  float TPrimeBG_R()const;
  float TPrimeHT_R()const;
  float TPrimeHAntiT_R()const;

  float Higgs_R()const;
  float ZBoson_R()const;
  float WBoson_R()const;
  bool IsSemiLeptonicDecay()const;
  

 private:

  GenParticle bottom, topHiggs, topZ, Z, W, H, top_W, top_b, main_top, TPrime, ZPrime, daughter1, daughter2, daughter3,up, down, daughter5, daughter6, gluon, Bottom,TPrimedaughter3,Top1,Gluon,Top2, ZPrimeGluon, TPrimereco1,TPrimereco2,TPrimereco3, ZPrimereco, decaydaughterH1,decaydaughterH2,decaydaughterZ1,decaydaughterZ2 , decaydaughterW1, decaydaughterW2,WHiggsTop1,BHiggsTop1,WZTop1, lepton,neutrino, leptop,hadtop, firstH, firstW, firstZ, Wleptop,Whadtop,muon_top, muon_tprime, blep, quark1, quark2, bhad, BZTop1, Whadtopq2,Whadtopq1,WHiggsTop1daughter1,WHiggsTop1daughter2,WTop1daughter1,WTop1daughter2;
 
  int n_TllTll=0,n_TllTh=0,n_ThTll=0,n_ThTh=0;
  int n_TllTllgg=0,n_TllTllZZ=0,n_TllTllWW=0,n_TllTllbb=0, n_TllThgg=0, n_TllThZZ=0, n_TllThWW=0,n_TllThbb=0,n_ThTllgg=0,n_ThTllZZ=0,n_ThTllWW=0,n_ThTllbb=0, n_ThThgg=0, n_ThThZZ=0, n_ThThWW=0,n_ThThbb=0;
  int n_TllWll=0,n_TllWh=0,n_ThWll=0,n_ThWh=0;
  int n_TllTllZll=0,n_TllThZll=0,n_ThTllZll=0,n_ThThZll=0,n_TllTllZh=0,n_TllThZh=0,n_ThTllZh=0,n_ThThZh=0;  
  int n_H=0, n_Z=0, n_W=0;
  

};




class ZPrimeGenProducer: public uhh2::AnalysisModule {
public:
    explicit ZPrimeGenProducer(uhh2::Context & ctx, const std::string & name = "zprimegen", bool throw_on_failure = true);
     virtual bool process(uhh2::Event & event) override;
   
 
private:
    uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
    uhh2::Event::Handle<std::vector<int>> h_n;
    bool throw_on_failure;
    

};
