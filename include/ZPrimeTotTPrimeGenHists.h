#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "ZPrimeTotTPrimeGenSelections.h"

/**
 *
 *Histograms for Z' to TT' on generator level
 *
 **/

class ZPrimeGenHists: public uhh2::Hists{

  public:
    ZPrimeGenHists(uhh2::Context & ctx, const std::string & dirname);
    virtual void fill(const uhh2::Event & event) override;


  protected:
    TH1F *M_ZPrime, *M_TPrime,*invMass_ZPrime,*invMass_TPrime1,*invMass_TPrime2,*invMass_TPrime3,*invMass_TPrime524,*invMass_TPrime524g ,*invMass_TPrime5624g, *invMass_TPrime123,*invMass_TPrime621, *PT_ZPrime, *PT_TPrime, *PT_t, *Phi_TPrime_Top,*Phi_TPrimeToHT,*Phi_TPrimeToZT,*Phi_TPrimeToWB, *Phi_TopToWB,  *Eta_TPrime, *Eta_Top, *Eta_TPrimeToHT,*Eta_TPrimeToZT, *Eta_TPrimeToWB, *Eta_TopToWB,*PT_Higgs, *PT_TopHiggs, * PT_TopZ, * PT_WTPrime, *PT_WTop, *PT_BTop, *PT_BTPrime, * PT_Z,*R_TPrimeToHT,*R_TPrimeToZT,*R_TPrimeToWB,*R_TPrimeToWG,*R_TPrimeToBG, *R_TopToWB,*R_TPrime_Top,*R_TPrimeToHTreco,*R_TPrimeToHAntiT,*pdgId_Daughter1,*pdgId_Daughter2,*pdgId_Daughter3,*pdgId_Daughter5,*pdgId_Daughter6, *R_HiggsToxx, *R_ZToxx, *R_WToxx, *DecayH , *DecayW, *DecayZ ; 
    uhh2::Event::Handle<ZPrimeGen> h_zprimegen;

};




