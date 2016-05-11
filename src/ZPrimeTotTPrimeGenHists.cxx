#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"

#include "TH1F.h"
#include <iostream>

using namespace uhh2;


ZPrimeGenHists::ZPrimeGenHists(Context & ctx, const std::string & dirname): Hists(ctx, dirname){
  
  //Massen von Z' und T'
  book<TH1F>("Masse_ZPrime","Masse_{ZPrime}",1000,0,4000);
  book<TH1F>("Masse_TPrime","Masse_{TPrime}",1000,0,2500);
  book<TH1F>("invMass_ZPrime", "Mass_{ZPrime}",1000,0,4000);
  book<TH1F>("invMass_TPrime1", "Mass_{TPrime1}",1000,0,4000);
  book<TH1F>("invMass_TPrime2", "Mass_{TPrime2}",1000,0,4000);
  book<TH1F>("invMass_TPrime3", "Mass_{TPrime3}",1000,0,4000);
  // book<TH1F>("invMass_TPrime123", "Mass_{TPrime123}",1000,0,4000);
  // book<TH1F>("invMass_TPrime524", "Mass_{TPrime524}",1000,0,4000);
  // book<TH1F>("invMass_TPrime524g", "Mass_{TPrime524g}",1000,0,4000);
  // book<TH1F>("invMass_TPrime621", "Mass_{TPrime621}",1000,0,4000);
  // book<TH1F>("invMass_TPrime5624g", "Mass_{TPrime5624g}",1000,0,4000);

  // //pdgId()of the daughters
  // book<TH1F>("pdgId_Daughter1","pdgId_Daughter1",34,-7,27);
  // book<TH1F>("pdgId_Daughter2","pdgId_Daughter2",34,-7,27);
  // book<TH1F>("pdgId_Daughter3","pdgId_Daughter3",34,-7,27);
  // book<TH1F>("pdgId_Daughter5","pdgId_Daughter5",34,-7,27);
  // book<TH1F>("pdgId_Daughter6","pdgId_Daughter6",34,-7,27);

  //PT von Z', T', T   
  book<TH1F>("PT_ZPrime","P_{T, ZPrime}",1500,0,1500);
  book<TH1F>("PT_TPrime","P_{T, TPrime}",1500,0,1500);
  book<TH1F>("PT_Top","P_{T, T}",1500,0,1500);

  //PT von H,T,W,B,Z,T,W,B
  book<TH1F>("PT_Higgs","P_{T,Higgs}",1500,0,1500);
  book<TH1F>("PT_TopHiggs","P_{T,TopHiggs}",1500,0,1500);
  book<TH1F>("PT_TopZ","P_{T,TopZ}",1500,0,1500);
  book<TH1F>("PT_WTPrime","P_{T,WTPrime}",1500,0,1500);
  book<TH1F>("PT_WTop","P_{T,WTop}",1500,0,1500);
  book<TH1F>("PT_BTop","P_{T,BTop}",1500,0,1500);
  book<TH1F>("PT_BTPrime","P_{T,BTPrime}",1500,0,1500);
  book<TH1F>("PT_Z","P_{T,Z}",1500,0,1500);
  
  //Delta Phi von T' &  T
  book<TH1F>("Phi_TPrime_Top","Delta Phi_{TPrime,Top}",100,0,10);
  
  //Eta von T', T
  book<TH1F>("Eta_TPrime","Eta_{TPrime}",100,-4,4);
  book<TH1F>("Eta_Top","Eta_{Top}",100,-4,4);

  //Delta Phi von Decayparticles
  book<TH1F>("Phi_TPrimeToHT","Phi_{TPrime to HT}",100,0,10);
  book<TH1F>("Phi_TPrimeToZT","Phi_{TPrime to ZT}",100,0,10);
  book<TH1F>("Phi_TPrimeToWB","Phi_{TPrime to WB}",100,0,10);

  book<TH1F>("Phi_TopToWB","Phi_{Top to WB}",100,0,10);

  //Eta von Decay particles
  book<TH1F>("Eta_TPrimeToHT","Eta_{TPrime to HT}",100,-4,4);
  book<TH1F>("Eta_TPrimeToZT","Eta_{TPrime to ZT}",100,-4,4);
  book<TH1F>("Eta_TPrimeToWB","Eta_{TPrime to WB}",100,-4,4);

  book<TH1F>("Eta_TopToWB","Eta_{Top to WB}",100,-4,4);

  //Delta R von Decayparticles
  book<TH1F>("R_TPrimeToHT","R_{TPrime to HT}",100,0,10);
  book<TH1F>("R_TPrimeToZT","R_{TPrime to ZT}",100,0,10);
  book<TH1F>("R_TPrimeToWB","R_{TPrime to WB}",100,0,10);
 book<TH1F>("R_TPrimeToWG","R_{TPrime to WB}",100,0,10);
 book<TH1F>("R_TPrimeToBG","R_{TPrime to BG}",100,0,10);
 book<TH1F>("R_TPrimeToHTreco","R_{TPrime to HTreco}",100,0,10);
book<TH1F>("R_TPrimeToHAntiT","R_{TPrime to HAntiT}",100,0,10);

  book<TH1F>("R_TopToWB","R_{Top to WB}",100,0,10);
  book<TH1F>("R_TPrime_Top","Delta R_{TPrime,Top}",100,0,10);

  book<TH1F>("R_HToxx","R_{H to XX}",100,0,10);
  book<TH1F>("R_ZToxx","R_{Z to XX}",100,0,10);
  book<TH1F>("R_WToxx","R_{W to XX}",100,0,10);

  book<TH1F>("DecayH","final states of H",9, 0, 9);
 book<TH1F>("DecayW","final states of W",3, 0, 3);
 book<TH1F>("DecayZ","final states of Z",4, 0, 4);
 

  //Generator-Teilchen
  h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
}


void ZPrimeGenHists::fill(const uhh2::Event & event){
  //Generator-Teilchen
    if(!event.is_valid(h_zprimegen)){
      return;
    }
    const auto & genparticle = event.get(h_zprimegen);

    //Weightning
  double weight = event.weight;

  //Masses of Z' and T'
  float MZPrime = genparticle.ZPrime_Mass();
  hist("Masse_ZPrime")->Fill(MZPrime, weight);
  float MTPrime = genparticle.TPrime_Mass();
  hist("Masse_TPrime")->Fill(MTPrime, weight);
  float invMZPrime = genparticle.ZPrime_invMass();
  if(invMZPrime){
  hist("invMass_ZPrime")->Fill(invMZPrime, weight);
  }
  float invMTPrime1 = genparticle.invMass1();
  if(invMTPrime1){
  hist("invMass_TPrime1")->Fill(invMTPrime1, weight);
  }
  float invMTPrime2 = genparticle.invMass2();
  if(invMTPrime2){
  hist("invMass_TPrime2")->Fill(invMTPrime2, weight);
  }
  float invMTPrime3 = genparticle.invMass3();
  if(invMTPrime3){
   hist("invMass_TPrime3")->Fill(invMTPrime3, weight);
  }
 //  float invMTPrime524 = genparticle.invMass524();
 //  hist("invMass_TPrime524")->Fill(invMTPrime524, weight);
 //  float invMTPrime524g = genparticle.invMass524g();
 //  hist("invMass_TPrime524g")->Fill(invMTPrime524g, weight);
 //  float invMTPrime123 = genparticle.invMass123();
 //  hist("invMass_TPrime123")->Fill(invMTPrime123, weight);
 // float invMTPrime621 = genparticle.invMass621();
 //  hist("invMass_TPrime621")->Fill(invMTPrime621, weight);
 // float invMTPrime5624g = genparticle.invMass5624g();
 //  hist("invMass_TPrime5624g")->Fill(invMTPrime5624g, weight);

  // //pdgId() of the daughters
  // int pdgIddaughter1 = genparticle.pdgDaughter1();
  // hist("pdgId_Daughter1")->Fill(pdgIddaughter1,weight);
  // int pdgIddaughter2 = genparticle.pdgDaughter2();
  // hist("pdgId_Daughter2")->Fill(pdgIddaughter2,weight);
  // int pdgIddaughter3 = genparticle.pdgDaughter3();
  // hist("pdgId_Daughter3")->Fill(pdgIddaughter3,weight);
  // int pdgIddaughter5 = genparticle.pdgDaughter5();
  // hist("pdgId_Daughter5")->Fill(pdgIddaughter5,weight);
  // int pdgIddaughter6 = genparticle.pdgDaughter6();
  // hist("pdgId_Daughter6")->Fill(pdgIddaughter6,weight);


  //Fill pt Histograms
  float ptTPrime = genparticle.TPrime_PT();
  hist("PT_TPrime")-> Fill(ptTPrime,weight);
  float ptZPrime = genparticle.ZPrime_PT();
  hist("PT_ZPrime")-> Fill(ptZPrime,weight);
  float ptTop = genparticle.Main_Top_PT();
  hist("PT_Top")-> Fill(ptTop,weight);

  //Fill PT of Decay products
  float ptHiggs = genparticle.Higgs_PT();
  hist("PT_Higgs")->Fill(ptHiggs,weight);
  float ptTopHiggs = genparticle.TopHiggs_PT();
  hist("PT_TopHiggs")->Fill(ptTopHiggs,weight);
  float ptTopZ = genparticle.TopZ_PT();
  hist("PT_TopZ")->Fill(ptTopZ,weight);
  float ptWTPrime = genparticle.WTPrime_PT();
  hist("PT_WTPrime")->Fill(ptWTPrime,weight);
  float ptWTop = genparticle.WTop_PT();
  hist("PT_WTop")->Fill(ptWTop,weight);
  float ptBTop = genparticle.BTop_PT();
  hist("PT_BTop")->Fill(ptBTop,weight);
  float ptBTPrime =genparticle.BTPrime_PT();
  hist("PT_BTPrime")->Fill(ptBTPrime,weight);
  float ptZ = genparticle.Z_PT();
  hist("PT_Z")->Fill(ptZ,weight);


  //Fill Phi Histograms 
  float deltaphiTPrimeTop=genparticle.TPrimeTop_Phi() ;
  hist("Phi_TPrime_Top")->Fill(deltaphiTPrimeTop, weight);
  float deltaphiHT=genparticle.HT_Phi() ;
  hist("Phi_TPrimeToHT")->Fill(deltaphiHT, weight);
  float deltaphiZT=genparticle.ZT_Phi() ;
  hist("Phi_TPrimeToZT")->Fill(deltaphiZT, weight);
  float deltaphiWB=genparticle.TPrimeWB_Phi() ;
  hist("Phi_TPrimeToWB")->Fill(deltaphiWB, weight);
  float deltaphiTopWB=genparticle.TopWB_Phi() ;
  hist("Phi_TopToWB")->Fill(deltaphiTopWB, weight);

  //Fill Eta Histograms 
  float etaTPrime = genparticle.TPrime_Eta();
  hist("Eta_TPrime")-> Fill(etaTPrime, weight);
 float etaTop = genparticle.Top_Eta();
  hist("Eta_Top")-> Fill(etaTop, weight);
 float etaHT = genparticle.HT_Eta();
  hist("Eta_TPrimeToHT")-> Fill(etaHT, weight);
 float etaZT = genparticle.ZT_Eta();
  hist("Eta_TPrimeToZT")-> Fill(etaZT, weight);
 float etaTPrimeWB = genparticle.TPrimeWB_Eta();
  hist("Eta_TPrimeToWB")-> Fill(etaTPrimeWB, weight);
 float etaTopWB = genparticle.TopWB_Eta();
  hist("Eta_TopToWB")-> Fill(etaTopWB, weight);

  //Fill R Histograms 
  float deltarTPrimeTop=genparticle.TPrimeTop_R() ;
  hist("R_TPrime_Top")->Fill(deltarTPrimeTop, weight);
  float deltarHT=genparticle.HT_R() ;
  hist("R_TPrimeToHT")->Fill(deltarHT, weight);
  float deltarZT=genparticle.ZT_R() ;
  hist("R_TPrimeToZT")->Fill(deltarZT, weight);
  float deltarWB=genparticle.TPrimeWB_R() ;
  hist("R_TPrimeToWB")->Fill(deltarWB, weight);
 float deltarWG=genparticle.TPrimeWG_R() ;
  hist("R_TPrimeToWG")->Fill(deltarWG, weight);
 float deltarBG=genparticle.TPrimeBG_R() ;
 if(deltarBG){
  hist("R_TPrimeToBG")->Fill(deltarBG, weight);
 }
  float deltarTopWB=genparticle.TopWB_R() ;
  hist("R_TopToWB")->Fill(deltarTopWB, weight);


 float deltarHTreco=genparticle.TPrimeHT_R() ;
float deltarHAntiT=genparticle.TPrimeHAntiT_R() ;
 
    // float deltarHTreco=genparticle.TPrimeHT_R() ;
 if(deltarHTreco){
  hist("R_TPrimeToHTreco")->Fill(deltarHTreco, weight);
 }
 //float deltarHAntiT=genparticle.TPrimeHAntiT_R() ;
 if(deltarHAntiT){
  hist("R_TPrimeToHAntiT")->Fill(deltarHAntiT, weight);
 }

 float deltardecayH  = genparticle.Higgs_R();
 hist("R_HToxx")->Fill(deltardecayH, weight);
 float deltardecayZ  = genparticle.ZBoson_R();
 hist("R_ZToxx")->Fill(deltardecayZ, weight);
 float deltardecayW  = genparticle.WBoson_R();
 hist("R_WToxx")->Fill(deltardecayW, weight); 



 if(n_vector[13]){
   hist("DecayH")->Fill (8, weight);
 }

 if(n_vector[0]){
   hist("DecayH")->Fill(0.,weight);
 }else if(n_vector[1]){
   hist("DecayH")->Fill(1,weight);
 }else if(n_vector[2]){
   hist("DecayH")->Fill(2,weight);
 }else if(n_vector[3]){
   hist("DecayH")->Fill(3,weight);
 }else if(n_vector[4]){
   hist("DecayH")->Fill(4,weight);
 }else if(n_vector[5]){
   hist("DecayH")->Fill(5,weight);
 }else if(n_vector[6]){
   hist("DecayH")->Fill(6,weight);
 }else if(n_vector[7]){
   hist("DecayH")->Fill(7,weight);
 }

 if(n_vector[14]){
   hist("DecayW")-> Fill(2,weight);
 }

 if(n_vector[8]){
   hist("DecayW")->Fill(0.,weight);
 }else if(n_vector[9]){
   hist("DecayW")->Fill(1,weight);
 }


 if(n_vector[15]){
   hist("DecayZ")-> Fill(3, weight);
 }
 if(n_vector[10]){
   hist("DecayZ")->Fill(0.,weight);
 }else if(n_vector[11]){
   hist("DecayZ")->Fill(1,weight);
 }else if(n_vector[12]){
   hist("DecayZ")->Fill(2,weight);
 }


}
