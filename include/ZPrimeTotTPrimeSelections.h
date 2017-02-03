#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"
#include <UHH2/common/include/TTbarGen.h>
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/Utils.h"

namespace uhh2examples {
    
  /* Select events with at least two jets in which the leading two jets have deltaphi > 2.7 and the third jet pt is
   * below 20% of the average of the leading two jets, where the minimum deltaphi and
   * maximum third jet pt fraction can be changed in the constructor.
   * The jets are assumed to be sorted in pt.
   */
  class ZPrimeTotTPrimeDijetSelection: public uhh2::Selection {
  public:
    ZPrimeTotTPrimeDijetSelection(float dphi_min = 2.7f, float third_frac_max = 0.2f);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float dphi_min, third_frac_max;
  };

  class ZPrimeTotTPrimeNJetCut: public uhh2::Selection {
  public:
    ZPrimeTotTPrimeNJetCut(int nmin_, int nmax_, float ptmin1_,float ptmin2_, float etamax_);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    int nmin, nmax;
    float ptmin1, ptmin2, etamax;
  };

 class ZPrimeTotTPrimeNTopJetCut: public uhh2::Selection {
  public:
    ZPrimeTotTPrimeNTopJetCut(int nmin_, int nmax_, float ptmin1_,float ptmin2_, float etamax_);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    int nmin, nmax;
    float ptmin1, ptmin2, etamax;
  };

  class ZPrimeTotTPrimeMassCut: public uhh2::Selection{
  public:
    ZPrimeTotTPrimeMassCut(float mmin_, float mmax_);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float mmin, mmax;
  };

  class ZPrimeTotTPrimeTopMassCut: public uhh2::Selection{
  public:
    ZPrimeTotTPrimeTopMassCut(float mmin_, float mmax_);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float mmin, mmax;
  };

  class ZPrimeTotTPrimeMETCut: public uhh2::Selection{
  public:
    ZPrimeTotTPrimeMETCut(float min_MET=0, float max_MET=1000);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float min_MET, max_MET;
  };


  class ZPrimeTotTPrimePtJetCut: public uhh2::Selection{
  public:
    ZPrimeTotTPrimePtJetCut(float min_pt=0, float max_pt=1000,unsigned int jetnumber=0);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float min_pt, max_pt;
    unsigned int jetnumber;
  };

  class ZPrimeTotTPrimePtTopJetCut: public uhh2::Selection{
  public:
    ZPrimeTotTPrimePtTopJetCut(float min_pt=0, float max_pt=1000,unsigned int jetnumber=0);
    virtual bool passes(const uhh2::Event & event) override;
  private:
    float min_pt, max_pt;
    unsigned int jetnumber;
  };

}

/*
class MuonptSlection:public uhh2::Selection{
 public:
  MuonptSelection(float ptmin = 2);
  virtual bool passes(const uhh2::Event & event) override;
private:
  float ptmin;
}
*/



namespace uhh2{

namespace btagging {
    
enum class csv_wp {
    loose, medium, tight
};

/// convert a CSV working point to a numerical threshold of the discriminator.
float csv_threshold(const csv_wp & wp);

}

/// Select events with certain minimum / maximum number of b-tagged jets using the CSV tagger
class NBTagSelection: public Selection {
public:
    /// In case nmax=-1, no cut on the maximum is applied.
    explicit NBTagSelection(int nmin, int nmax = -1, btagging::csv_wp wp = btagging::csv_wp::medium);
    virtual bool passes(const Event & event) override;
    
private:
    int nmin, nmax;
    float min_csv;
};

 /////
 class HtSelection: public uhh2::Selection {
  public:
    explicit HtSelection(double ht_min=0., double ht_max=-1);
    virtual bool passes(const uhh2::Event & event);
 private:
    double ht_min, ht_max;
 };


/////

 class GenMttbarCut: public Selection {
 public:
  explicit GenMttbarCut(Context&, const float, const float, const std::string&);
  virtual bool passes(const Event&) override;

 private:
  float mttbar_min_, mttbar_max_;
  Event::Handle<TTbarGen> h_ttbargen_;
};
/////
// ////////////////////////////////////////////////////////
 class ZPrimeTotTPrimeDeltaRCut:public Selection{
 public:
  explicit ZPrimeTotTPrimeDeltaRCut(uhh2::Context& ctx ,float deltar_,const std::string & discriminator_name );
  virtual bool passes(const Event& event) override;

 private:
  float deltar;
  std::string name;
  uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>  h_hyps;
  uhh2::Event::Handle<ZPrimeGen> h_zprimegen;

};


// ////////
// ////////////////////////////////////////////////////////
 class ZPrimeTotTPrimeAK4cleaner{
 public:
 ZPrimeTotTPrimeAK4cleaner(float deltar_ , int num_=1):deltar(deltar_),num(num_){};

   bool operator()(const Particle & p, const Event & event)const{
     assert(event.topjets);
     std::vector<TopJet> *Topjets = event.topjets;    
     if(Topjets->size()<1)return true;  
     TopJet jet= Topjets->at(0);
     if(num>=2){
       jet= Topjets->at(1);
     }

    float delta1 = deltaR(jet, p);
    return delta1 > deltar;
  }

 private:
  float deltar;
  int num;
};

// ////////


// ////////////////////////////////////////////////////////
 class ZPrimeTotTPrimeUnHiggstag{
 public:
 ZPrimeTotTPrimeUnHiggstag(float m_min_,float m_max_):m_min(m_min_),m_max(m_max_){};

   bool operator()(const TopJet & topjet , const Event & event)const{
     assert(event.topjets);
  LorentzVector subjet_sum;
  for (const auto s : topjet.subjets()) {
    subjet_sum += s.v4();
  }

  return (subjet_sum.M()<m_min ||subjet_sum.M()>m_max );
  }

 private:
   float m_min,m_max;
};




// ////////////////////////////////////////////////////////
 /* class ZPrimeTotTPrimeAK8cleaner{ */
/*  public: */
/*  ZPrimeTotTPrimeAK8cleaner(float deltar_ ,ReconstructionHypothesis hyp_):deltar(deltar_),hyp(hyp_){}; */

/*    bool operator()(const Particle & p, const Event & event)const{ */
/*      assert(event.topjets); */
/*      std::vector<TopJet> *Topjets = event.topjets;     */
/*      if(Topjets->size()<1)return true;   */
/*      TopJet jet= Topjets->at(0); */
/*      if(num>=2){ */
/*        jet= Topjets->at(1); */
/*      } */

/*     float delta1 = deltaR(jet, p); */
/*     return delta1 > deltar; */
/*   } */

/*  private: */
/*   float deltar; */
/*   int num; */
/* }; */

/* // //////// */
// ////////////////////////////////////////////////////////
 class ZPrimeTotTPrimenumbersub:public Selection{
 public:
   explicit ZPrimeTotTPrimenumbersub(unsigned int num_);
   virtual bool passes(const Event& event) override;
 private:
   unsigned int num;
 };

// ////////////////////////////////////////////////////////
 class ZPrimeTotTPrimedrmin:public Selection{
 public:
   explicit ZPrimeTotTPrimedrmin(float dr_);
   virtual bool passes(const Event& event) override;
 private:
   float dr, drmin;
 };

// ////////////////////////////////////////////////////////
 class TwoDCut : public Selection {
   public:
    explicit TwoDCut(float min_deltaR, float min_pTrel): min_deltaR_(min_deltaR), min_pTrel_(min_pTrel) {}
    virtual bool passes(const Event&) override;

   private:
    float min_deltaR_, min_pTrel_;
  };

// ////////////////////////////////////////////////////////
class ZPrimeTotTPrimeChiCut:public Selection{
 public:
   explicit ZPrimeTotTPrimeChiCut(uhh2::Context& ctx,float chi_,const std::string & hyps_name,const std::string & discriminator_name);
   virtual bool passes(const Event& event) override;
 private:
   float chi;
   std::string m_discriminator_name;
   uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>> h_hyps;
 };


// ////////////////////////////////////////////////////////
class ZPrimeTotTPrimeRelIsoCut:public Selection{
 public:
   explicit ZPrimeTotTPrimeRelIsoCut(double relmin_);
   virtual bool passes(const Event& event) override;
 private:
   double relisomin;

 };
}
// ////////////////////////////////////////////////////////
class TopJetLeptonDeltaRCleaner : public uhh2::AnalysisModule {
 public:
  explicit TopJetLeptonDeltaRCleaner(float mindr=0.8): minDR_(mindr) {}
  virtual bool process(uhh2::Event&) override;

 private:
  float minDR_;
};

namespace uhh2{
// ////////////////////////////////////////////////////////
class ZPrimeTotTPrimeMuonPT:public Selection{
 public:
   explicit ZPrimeTotTPrimeMuonPT(double ptmax_);
   virtual bool passes(const Event& event) override;
 private:
   double ptmax;

 };
}
// ////////////////////////////////////////////////////////
class Tau32_inverted {
public:
    explicit Tau32_inverted(double threshold_ = 0.7): threshold(threshold_){}
    
    bool operator()(const TopJet & topjet, const uhh2::Event & event) const;
    
private:
    double threshold;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////
class UnHiggsTag {
public:
    explicit UnHiggsTag(float minmass = 60.f, float maxmass = std::numeric_limits<float>::infinity(), JetId const & id = CSVBTag(CSVBTag::WP_MEDIUM)) :
        minmass_(minmass), maxmass_(maxmass), btagid_(id) {}

    bool operator()(TopJet const & topjet, uhh2::Event const & event) const;

private:
    float minmass_, maxmass_;
    JetId btagid_;

};

//////////////////////////////////////////////////////////////////////////////////////////////////////////
class UnType2TopTag {
public:
  
  enum class MassType {groomed, ungroomed};
  explicit UnType2TopTag(double mjetLower=60., double mjetUpper=100., MassType typeOfMass = MassType::ungroomed, boost::optional<JetId> SubjetId=boost::none);
  explicit UnType2TopTag(MassType typeOfMass);
  
  bool operator()(const TopJet & topjet, const uhh2::Event & event) const;

 private:
   double m_mjetLower;
   double m_mjetUpper;
   MassType m_typeOfMass;
   boost::optional<JetId> m_SubjetId;
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////

/* class ZPrimeTotTPrimeEff: public uhh2::AnalysisModule { */
/* public: */
/*   explicit ZPrimeTotTPrimeEff(uhh2::Context & ctx);  */

/*   virtual bool process(uhh2::Event & event) override; */

/* private: */
/*   double SF; */
/* }; */
//////////////////////////////////////////////////////////////////////////////////////////////////////////
