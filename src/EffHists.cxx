#include "UHH2/core/include/LorentzVector.h"

#include "UHH2/ZPrimeTotTPrime/include/EffHists.h"
#include "UHH2/common/include/TTbarGen.h"
#include "TH1F.h"


using namespace std;
using namespace uhh2;

//Hists for studies to the eff. of the W/Z tagger in MC and Data


EffHists::EffHists(Context & ctx, const string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  ///////////////////////////      GenLevel Studies for ttbar  ///////////////////////////////////
  
  //Analysis idea: W has to be a topjet (also tagger requierment)
  book<TH1F>("match_gen_WTopjet","Match of W(GenLevel) and Topjets (recolevel)",3,0,3);
  hist("match_gen_WTopjet")->Fill("matched",0);
  hist("match_gen_WTopjet")->Fill("missed",0);
  hist("match_gen_WTopjet")->Fill("top matched",0);

  //Analysis idea: Are b (GenLevel) fullfilling b jet criteria (loose, medium, tight)
  book<TH1F>("match_gen_blepbjet","Match of b(GenLevel) and b jet criteria",6,0,6);
  hist("match_gen_blepbjet")->Fill("matched loose",0);
  hist("match_gen_blepbjet")->Fill("missed loose",0);
  hist("match_gen_blepbjet")->Fill("matched medium",0);
  hist("match_gen_blepbjet")->Fill("missed medium",0);
  hist("match_gen_blepbjet")->Fill("matched tight",0);
  hist("match_gen_blepbjet")->Fill("missed tight",0);

  book<TH1F>("match_gen_bhadbjet","Match of b(GenLevel) and b jet criteria",6,0,6);
  hist("match_gen_bhadbjet")->Fill("matched loose",0);
  hist("match_gen_bhadbjet")->Fill("missed loose",0);
  hist("match_gen_bhadbjet")->Fill("matched medium",0);
  hist("match_gen_bhadbjet")->Fill("missed medium",0);
  hist("match_gen_bhadbjet")->Fill("matched tight",0);
  hist("match_gen_bhadbjet")->Fill("missed tight",0);

  //Analysis idea: the distance between blep und mu is smaller than bhad and mu
  book<TH1F>("deltar_gen_blepmu","DeltaR between b (lep. top) and mu",50,0,5);
  book<TH1F>("deltaphi_gen_blepmu","DeltaPhi between b (lep. top) and mu",50,0,M_PI);
  book<TH1F>("deltar_gen_bhadmu","DeltaR between b (had. top) and mu",50,0,5);
  book<TH1F>("deltaphi_gen_bhadmu","DeltaPhi between b (had. top) and mu",50,0,M_PI);

  book<TH1F>("deltar_gen_blepW","DeltaR between b (lep. top) and W",50,0,5);
  book<TH1F>("deltaphi_gen_blepW","DeltaPhi between b (lep. top) and W",50,0,M_PI);
  book<TH1F>("deltar_gen_bhadW","DeltaR between b (had. top) and W",50,0,5);
  book<TH1F>("deltaphi_gen_bhadW","DeltaPhi between b (had. top) and W",50,0,M_PI);

  //difference between deltaR and deltaPhi between lep top and had top
  book<TH1F>("deltar_gen_lephadtop","DeltaR between lep. top and had. top (GenLevel)",50,0,5);
  book<TH1F>("deltaphi_gen_lephadtop","DeltaPhi between lep. top and had. top (GenLevel)",50,0,M_PI);


   ///////////////////////////      Reco Studies for ttbar  ///////////////////////////////////

  // Matching test

  //b jets matching to GenLevel
  book<TH1F>("match_reco_blep","Matching reco b jet (lep top) match to GenLevel b",2,0,2);
  hist("match_reco_blep")->Fill("matched",0);
  hist("match_reco_blep")->Fill("missed",0);

  book<TH1F>("match_reco_bhad","Matching reco b jet (had top) match to GenLevel b",2,0,2);
  hist("match_reco_bhad")->Fill("matched",0);
  hist("match_reco_bhad")->Fill("missed",0);

  //matching v4 lep top (Analysis idea> back to back)
  book<TH1F>("match_reco_leptop","Matching v4 lep top (reco vs GenLevel)",2,0,2);
  hist("match_reco_leptop")->Fill("matched",0);
  hist("match_reco_leptop")->Fill("missed",0);

  //bjets on lep hemisphere (loose, medium, tight)
  book<TH1F>("reco_binleptop_loose","bjets (loose) in lep. hemisphere",3,0,3);
  hist("reco_binleptop_loose")->Fill("otherbjet",0);
  hist("reco_binleptop_loose")->Fill("no_otherbjet",0);
  hist("reco_binleptop_loose")->Fill("otherbjet_bhad",0);
  
  book<TH1F>("reco_binleptop_medium","bjets (medium) in lep. hemisphere",3,0,3);
  hist("reco_binleptop_medium")->Fill("otherbjet",0);
  hist("reco_binleptop_medium")->Fill("no_otherbjet",0);
  hist("reco_binleptop_medium")->Fill("otherbjet_bhad",0);

  
  book<TH1F>("reco_binleptop_tight","bjets (tight) in lep. hemisphere",3,0,3);
  hist("reco_binleptop_tight")->Fill("otherbjet",0);
  hist("reco_binleptop_tight")->Fill("no_otherbjet",0);
  hist("reco_binleptop_tight")->Fill("otherbjet_bhad",0);

  //matched W
  book<TH1F>("matched_reco_W","Matching reco W to GenLevel W",2,0,2);
  hist("matched_reco_W")->Fill("W matched",0);
  hist("matched_reco_W")->Fill("W missed",0);

  //controll variables

  //mass reco W
  book<TH1F>("reco_mass_W", "mass of reco W Boson",100,0,200);
  
  //mass reco lep top
  book<TH1F>("reco_mass_leptop", "mass of reco lep top",100,100,300);
  
  //mass reco had top
  book<TH1F>("reco_mass_hadtop", "mass of reco had top",100,100,300);
  
  //deltaR and deltaPhi between lep top and had top
  book<TH1F>("reco_deltaphi_top","DeltaPhi between reco lep top and had top",50,0,M_PI);
  book<TH1F>("reco_deltar_top","DeltaR between reco lep top and had top",50,0,5);
    
  //deltaR and deltaPhi between decay products
  book<TH1F>("reco_deltaphi_blepmu","DeltaPhi between blep and mu (reco)",50,0,M_PI);
  book<TH1F>("reco_deltar_blepmu","DeltaR between blep and mu (reco)",50,0,5);
  book<TH1F>("reco_deltaphi_bhadmu","DeltaPhi between bhad and mu (reco)",50,0,M_PI);
  book<TH1F>("reco_deltar_bhadmu","DeltaR between bhad and mu (reco)",50,0,5);

  book<TH1F>("reco_deltaphi_blepW","DeltaPhi between blep and W (reco)",50,0,M_PI);
  book<TH1F>("reco_deltar_blepW","DeltaR between blep and W (reco)",50,0,5);
  book<TH1F>("reco_deltaphi_bhadW","DeltaPhi between bhad and W (reco)",50,0,M_PI);
  book<TH1F>("reco_deltar_bhadW","DeltaR between bhad and W (reco)",50,0,5);

  //////////////////////////////////    Handles   //////////////////////////////////////////////////////

  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
  h_btag_medium = ctx.get_handle<std::vector<Jet> >("BTag_medium");
  h_btag_loose = ctx.get_handle<std::vector<Jet> >("BTag_loose");
  h_btag_tight = ctx.get_handle<std::vector<Jet> >("BTag_tight");
  h_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(hyps_name);

  m_discriminator_name = discriminator_name;
}//constructor definition


  void EffHists::fill(const Event & event){

    assert(event.met);
    assert(event.jets);
    assert(event.topjets);
    double weight = event.weight;
    
    
    ///////////////////////////      GenLevel Studies for ttbar  ///////////////////////////////////
  
    //Analysis idea: W has to be a topjet (also tagger requierment)
    if(event.is_valid(h_ttbargen)){  
      const auto & ttbargen = event.get(h_ttbargen);
 if(ttbargen.IsSemiLeptonicDecay()){
      for(const auto & s : *event.topjets){
	double delta_Wtopjet = deltaR(s, ttbargen.WHad());
	if(delta_Wtopjet <= 0.8)hist("match_gen_WTopjet")->Fill("matched",weight);
	else {
	  hist("match_gen_WTopjet")->Fill("missed",weight);
	  double delta_TopTopjet = deltaR(s, ttbargen.TopHad());
	  if(delta_TopTopjet<=0.8) hist("match_gen_WTopjet")->Fill("top matched",weight);
	}
      }// over all topjets
   
      //Analysis idea: Are b (GenLevel) fullfilling b jet criteria (loose, medium, tight)   BLep
      if(event.is_valid(h_btag_medium)){
	const auto & btag_medium = event.get(h_btag_medium);
	for(Jet s:btag_medium){
	  double delta_blepbjet = deltaR(s,ttbargen.BLep());
	  if(delta_blepbjet <= 0.4) hist("match_gen_blepbjet")->Fill("matched medium",weight);
	else  hist("match_gen_blepbjet")->Fill("missed medium",weight);
      }//over all medium btags
    }//valid h_medium_btag

    if(event.is_valid(h_btag_loose)){
      const auto & btag_loose = event.get(h_btag_loose);
      for(Jet s:btag_loose){
	double delta_blepbjet = deltaR(s,ttbargen.BLep());
	if(delta_blepbjet <= 0.4) hist("match_gen_blepbjet")->Fill("matched loose",weight);
	else  hist("match_gen_blepbjet")->Fill("missed loose",weight);
      }//over all loose btags
    }//valid h_loose_btag


    if(event.is_valid(h_btag_tight)){
      const auto & btag_tight = event.get(h_btag_tight);
      for(Jet s:btag_tight){
	double delta_blepbjet = deltaR(s,ttbargen.BLep());
	if(delta_blepbjet <= 0.4) hist("match_gen_blepbjet")->Fill("matched tight",weight);
	else  hist("match_gen_blepbjet")->Fill("missed tight",weight);
      }//over all tight btags
    }//valid h_tight_btag


 //Analysis idea: Are b (GenLevel) fullfilling b jet criteria (loose, medium, tight)   BHad
if(event.is_valid(h_btag_medium)){
      const auto & btag_medium = event.get(h_btag_medium);
      for(Jet s:btag_medium){
	double delta_bhadbjet = deltaR(s,ttbargen.BHad());
	if(delta_bhadbjet <= 0.4) hist("match_gen_bhadbjet")->Fill("matched medium",weight);
	else  hist("match_gen_bhadbjet")->Fill("missed medium",weight);
      }//over all medium btags
    }//valid h_medium_btag

    if(event.is_valid(h_btag_loose)){
      const auto & btag_loose = event.get(h_btag_loose);
      for(Jet s:btag_loose){
	double delta_bhadbjet = deltaR(s,ttbargen.BHad());
	if(delta_bhadbjet <= 0.4) hist("match_gen_bhadbjet")->Fill("matched loose",weight);
	else  hist("match_gen_bhadbjet")->Fill("missed loose",weight);
      }//over all loose btags
    }//valid h_loose_btag


    if(event.is_valid(h_btag_tight)){
      const auto & btag_tight = event.get(h_btag_tight);
      for(Jet s:btag_tight){
	double delta_bhadbjet = deltaR(s,ttbargen.BHad());
	if(delta_bhadbjet <= 0.4) hist("match_gen_bhadbjet")->Fill("matched tight",weight);
	else  hist("match_gen_bhadbjet")->Fill("missed tight",weight);
      }//over all tight btags
    }//valid h_tight_btag

   
    //Analysis idea: the distance between blep und mu is smaller than bhad and mu
    double deltaR_gen_blepmu = deltaR(ttbargen.ChargedLepton(),ttbargen.BLep());
    double deltaPhi_gen_blepmu = deltaPhi(ttbargen.ChargedLepton(),ttbargen.BLep());

    double deltaR_gen_bhadmu = deltaR(ttbargen.ChargedLepton(),ttbargen.BHad());
    double deltaPhi_gen_bhadmu = deltaPhi(ttbargen.ChargedLepton(),ttbargen.BHad());

    hist("deltar_gen_blepmu")->Fill(deltaR_gen_blepmu,weight);
    hist("deltaphi_gen_blepmu")->Fill(deltaPhi_gen_blepmu,weight);
    hist("deltar_gen_bhadmu")->Fill(deltaR_gen_bhadmu,weight);
    hist("deltaphi_gen_bhadmu")->Fill(deltaPhi_gen_bhadmu,weight);


    double deltaR_gen_blepW = deltaR(ttbargen.WHad(),ttbargen.BLep());
    double deltaPhi_gen_blepW = deltaPhi(ttbargen.WHad(),ttbargen.BLep());

    double deltaR_gen_bhadW = deltaR(ttbargen.WHad(),ttbargen.BHad());
    double deltaPhi_gen_bhadW = deltaPhi(ttbargen.WHad(),ttbargen.BHad());
 
    hist("deltar_gen_blepW")->Fill(deltaR_gen_blepW,weight);
    hist("deltaphi_gen_blepW")->Fill(deltaPhi_gen_blepW,weight);
    hist("deltar_gen_bhadW")->Fill(deltaR_gen_bhadW,weight);
    hist("deltaphi_gen_bhadW")->Fill(deltaPhi_gen_bhadW,weight);
    

    //difference between deltaR and deltaPhi between lep top and had top
    double deltaR_lephadtop = deltaR(ttbargen.TopLep(),ttbargen.TopHad());
    double deltaPhi_lephadtop = deltaPhi(ttbargen.TopLep(),ttbargen.TopHad());

    hist("deltar_gen_lephadtop")->Fill(deltaR_lephadtop,weight);
    hist("deltaphi_gen_lephadtop")->Fill(deltaPhi_lephadtop,weight);

   
    ///////////////////////////      Reco Studies for ttbar  ///////////////////////////////////

    // Matching test
    //b jets matching to GenLevel
    if(event.is_valid(h_hyps)){

      std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = event.get(h_hyps);
      const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

      double deltaR_blep = deltaR(hyp->blep_v4(),ttbargen.BLep());
      if(deltaR_blep <= 0.4) hist("match_reco_blep")->Fill("matched",weight);
      else   hist("match_reco_blep")->Fill("missed",weight);


      double deltaR_bhad = deltaR(hyp->bhad_v4(),ttbargen.BLep());
      if(deltaR_bhad <= 0.4) hist("match_reco_bhad")->Fill("matched",weight);
      else   hist("match_reco_bhad")->Fill("missed",weight);


      double deltaR_leptop = deltaR(hyp->toplep_v4(), ttbargen.TopLep());
      if(deltaR_leptop <= 0.4)  hist("match_reco_leptop")->Fill("matched",weight);
      else  hist("match_reco_leptop")->Fill("missed",weight);

      //bjets on lep hemisphere (medium, medium, tight)
      if(event.is_valid(h_btag_medium)){
	const auto & btag_medium = event.get(h_btag_medium);
	for(Jet s:btag_medium){
	  double delta_bleptop = deltaPhi(s,hyp->toplep_v4());
	  if(delta_bleptop <= M_PI/2) hist("reco_binleptop_medium")->Fill("no_otherbjet",weight);
	  else  {
	    hist("reco_binleptop_medium")->Fill("otherbjet",weight);
	    double deltaR_bbhad = deltaR(s,hyp->bhad_v4());
	    if(deltaR_bbhad <= 0.4 ) hist("reco_binleptop_medium")->Fill("otherbjet_bhad",0);
	  }
	}//over all medium btags
      }//valid h_medium_btag

      //matched W
      double delta_W = deltaR(hyp->W_v4(), ttbargen.WHad());      
      if(delta_W < 0.8) hist("matched_reco_W")->Fill("W matched",weight);
      else hist("matched_reco_W")->Fill("W missed",weight);

    }//valid h_hyps
    }//isSemiLeptonicDecay
    }// valid h_ttbargen



    ///////////////////////////      Reco Studies for ttbar  ///////////////////////////////////

    //controll variables

    //mass reco W
    if(event.is_valid(h_hyps)){

      std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = event.get(h_hyps);
      const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
      
      double mass_W=0;
      if( (hyp->W_v4()).isTimelike() ){
	LorentzVector subjet_sum;
	for (const auto s : hyp->W_subjets()) {
	  subjet_sum += s.v4();
	}
	mass_W=subjet_sum.M();
      }else{
	LorentzVector subjet_sum;
	for (const auto s : hyp->W_subjets()) {
	  subjet_sum += s.v4();
	}
	mass_W=sqrt( -subjet_sum.mass2());
      }

      hist("reco_mass_W")->Fill(mass_W,weight);

      double mass_leptop=0;
      if(hyp->toplep_v4().isTimelike()) mass_leptop = hyp->toplep_v4().M();
      else mass_leptop = sqrt(hyp->toplep_v4().mass2());

      hist("reco_mass_leptop")->Fill(mass_leptop,weight);

      double mass_hadtop=0;
      if(hyp->tophad_v4().isTimelike()) mass_hadtop = hyp->tophad_v4().M();
      else mass_hadtop = sqrt(hyp->tophad_v4().mass2());

      hist("reco_mass_hadtop")->Fill(mass_hadtop,weight);

      double deltaphi_top = deltaPhi(hyp->toplep_v4(), hyp->tophad_v4());
      hist("reco_deltaphi_top")->Fill(deltaphi_top,weight);

      double deltar_top = deltaR(hyp->toplep_v4(), hyp->tophad_v4());
      hist("reco_deltar_top")->Fill(deltar_top,weight);

      //deltaR and deltaPhi between decay products
      double deltaphi_blepmu = deltaPhi(hyp->blep_v4(), hyp->lepton());
      hist("reco_deltaphi_blepmu")->Fill(deltaphi_blepmu,weight);

      double deltar_blepmu = deltaR(hyp->blep_v4(), hyp->lepton());
      hist("reco_deltar_blepmu")->Fill(deltar_blepmu,weight);

      double deltaphi_bhadmu = deltaPhi(hyp->bhad_v4(), hyp->lepton());
      hist("reco_deltaphi_bhadmu")->Fill(deltaphi_bhadmu,weight);

      double deltar_bhadmu = deltaR(hyp->bhad_v4(), hyp->lepton());
      hist("reco_deltar_bhadmu")->Fill(deltar_bhadmu,weight);

      double deltaphi_blepW = deltaPhi(hyp->blep_v4(), hyp->W_v4());
      hist("reco_deltaphi_blepW")->Fill(deltaphi_blepW,weight);

      double deltar_blepW = deltaR(hyp->blep_v4(), hyp->W_v4());
      hist("reco_deltar_blepW")->Fill(deltar_blepW,weight);

      double deltaphi_bhadW = deltaPhi(hyp->bhad_v4(), hyp->W_v4());
      hist("reco_deltaphi_bhadW")->Fill(deltaphi_bhadW,weight);

      double deltar_bhadW = deltaR(hyp->bhad_v4(), hyp->W_v4());
      hist("reco_deltar_bhadW")->Fill(deltar_bhadW,weight);


    } //valid h_hyps

 
  
}//fill function
