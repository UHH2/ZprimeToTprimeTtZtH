#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/TTbarGen.h"

#include "TH1F.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;

ZPrimeTotTPrimeHists::ZPrimeTotTPrimeHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  // book all histograms here
  // jets
  book<TH1F>("N_jets", "N_{jets}", 20, 0, 20); 
  book<TH1F>("N_jets_pt20", "N_{jets}, p_T > 20", 20, 0, 20); 
  book<TH1F>("N_jets_pt50", "N_{jets}, p_T > 50", 20, 0, 20); 
  book<TH1F>("N_jets_pt70", "N_{jets}, p_T > 70", 20, 0, 20); 
  book<TH1F>("N_jets_pt100", "N_{jets}, p_T > 100", 20, 0, 20); 
 
  book<TH1F>("eta_jet1", "#eta^{jet 1}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet2", "#eta^{jet 2}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet3", "#eta^{jet 3}", 40, -2.5, 2.5);
  book<TH1F>("eta_jet4", "#eta^{jet 4}", 40, -2.5, 2.5);

  // leptons
  book<TH1F>("N_mu", "N^{#mu}", 10, 0, 10);
  book<TH1F>("pt_mu", "p_{T}^{#mu} [GeV/c]", 40, 0, 200);
  book<TH1F>("eta_mu", "#eta^{#mu}", 40, -2.1, 2.1);
  book<TH1F>("reliso_mu", "#mu rel. Iso", 40, 0, 0.5);

  // primary vertices
  book<TH1F>("N_pv", "N^{PV}", 50, 0, 50);

  //general
  book<TH1F>("H_T", "H_{T}", 200, 0, 1000);
  book<TH1F>("E_Tmiss", "missingET", 200, 0, 1000);
  book<TH1F>("H_T1", "H_{T}", 200, 0, 1000);
 // book<TH1F>("H_T", "H_{T}", 200, 0, 1000);

  //Matched
  book<TH1F>("Higgs_matched","Higgs_matched",100,0,5 );
  book<TH1F>("Wlep_matched","Wlep_matched",100,0,5 );
  book<TH1F>("T_matched","T_matched",100,0,5 );
  book<TH1F>("Blep_matched","Blep_matched",100,0,5 );
  book<TH1F>("Top_matched","Top_matched",100,0,5 );
  book<TH1F>("Bh_matched","Bh_matched",100,0,5 );
  book<TH1F>("Wh_matched","Wh_matched",100,0,5 );
  book<TH1F>("MissMatched","MissMatched",22,0,22 );

  //TOPTAG
  book<TH1F>("MatchedTopTag","Matched Toptag Jets",5,0,5);
  hist("MatchedTopTag")->Fill("Toptag matched",0);
  hist("MatchedTopTag")->Fill("Toptag miss",0);
  hist("MatchedTopTag")->Fill("Higgs matched",0);
  hist("MatchedTopTag")->Fill("Z matched",0);
  hist("MatchedTopTag")->Fill("Wlep matched",0);

 book<TH1F>("pdgIdTopTag", "pdgId Toptag",26, 0,26);


  //HIGGSTAG
  book<TH1F>("MatchedHiggsTag","Matched Higgstag Jets",5,0,5);

  hist("MatchedHiggsTag")->Fill("Higgstag matched",0);
  hist("MatchedHiggsTag")->Fill("W matched",0);
  hist("MatchedHiggsTag")->Fill("Top matched",0);
  hist("MatchedHiggsTag")->Fill("Z matched",0);
  hist("MatchedHiggsTag")->Fill("Higgstag miss",0);

  book<TH1F>("pdgIdHiggsTag", "pdgId Higgstag",26, 0,26);


  //ZWTAG
  book<TH1F>("MatchedZWTag","Matched ZWtag Jets",14,0,14);
  hist("MatchedZWTag")->Fill("Higgs matched",0);
  hist("MatchedZWTag")->Fill("W matched",0);
  hist("MatchedZWTag")->Fill("Top matched",0);
  hist("MatchedZWTag")->Fill("Ztag matched",0);
  hist("MatchedZWTag")->Fill("Ztag miss",0);
 hist("MatchedZWTag")->Fill("Wlep matched",0);

 book<TH1F>("pdgIdZWTag", "pdgId zwTag",26, 0,26);

 //matched W in ttbar
 book<TH1F>("matchedWTTbar","W machted ttbar",2,0,2);
 hist("matchedWTTbar")->Fill("matched W",0);
 hist("matchedWTTbar")->Fill("not matched",0);


  //TAU
 book<TH1F>("Tau32H","Matched Higgs tau32 ",50,0,1);
 book<TH1F>("Tau32T","Matched Top tau32 ",50,0,1);
 book<TH1F>("Tau21H","Matched Higgs tau21 ",50,0,1);
 book<TH1F>("Tau21T","Matched Top tau21 ",50,0,1);

  //efficiency
  book<TH1F>("eff", "Eff",10,0,9);

  h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
  h_toptag = ctx.get_handle<std::vector<TopJet>>("TopTag");
  h_higgstag = ctx.get_handle<std::vector<TopJet>>("HiggsTag");
  h_ZWtag = ctx.get_handle<std::vector<TopJet>>("ZWTag");
  h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");


  //to get the ordering correct
 
  hist("MissMatched")->Fill("H matched",0);
  hist("MissMatched")->Fill("H miss",0);
  hist("MissMatched")->Fill("Z matched",0);
  hist("MissMatched")->Fill("Z miss",0);
  hist("MissMatched")->Fill("Wl matched",0);
  hist("MissMatched")->Fill("Wl miss",0);
  hist("MissMatched")->Fill("T matched",0);
  hist("MissMatched")->Fill("T miss",0);
  hist("MissMatched")->Fill("Bl matched",0);
  hist("MissMatched")->Fill("Bl miss",0);
  hist("MissMatched")->Fill("Bh matched",0);
  hist("MissMatched")->Fill("Bh miss",0);
  hist("MissMatched")->Fill("Wh matched",0);
  hist("MissMatched")->Fill("Wh miss",0);
  hist("MissMatched")->Fill("Top matched",0);
  hist("MissMatched")->Fill("Top miss",0);
  hist("MissMatched")->Fill("H2 matched",0);
  hist("MissMatched")->Fill("H2 miss",0);
  hist("MissMatched")->Fill("Top1 matched",0);
  hist("MissMatched")->Fill("Top1 miss",0);
  hist("MissMatched")->Fill("MainTop1 matched",0);
  hist("MissMatched")->Fill("MainTop1 miss",0);

  // Verzweigungsverhaeltnis Muon
 book<TH1F>("MuonMatch","Matched Muon",3,0,3);
 hist("MuonMatch")->Fill("Muon aus top",0);
 hist("MuonMatch")->Fill("Muon aus TPrime",0);
 hist("MuonMatch")->Fill("Muon aus Higgs",0);

 //Abstandsbestimmung GenMuon - RecMuon
 book<TH1F>("MuonDistance1","Distance GenMuon-RecMuon aus top",40,0,0.8);
 book<TH1F>("MuonDistance2","Distance GenMuon-RecMuon aus TPrime",40,0,0.8);


 //Anzahl der gebtaggeten jets
 book<TH1D>("number_btag","Number of btagged Jets",10,0,10);
 


}


void ZPrimeTotTPrimeHists::fill(const Event & event){
  // fill the histograms. Please note the comments in the header file:
  // 'hist' is used here a lot for simplicity, but it will be rather
  // slow when you have many histograms; therefore, better
  // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
  assert(event.met);
  assert(event.pvs);
  assert(event.jets);
  // Don't forget to always use the weight when filling.
  double weight = event.weight;
 
  int Njets20 =0;
  int Njets50 =0;
  int Njets70 =0;
  int Njets100 =0;


  std::vector<TopJet>* jets = event.topjets;
  int Njets = event.jets->size();
  hist("N_jets")->Fill(Njets, weight);
  
  for(const Jet &jet: *event.jets){
    if(jet.pt()>20) Njets20++;
    if(jet.pt()>50) Njets50++;
    if(jet.pt()>70) Njets70++;
    if(jet.pt()>100) Njets100++;
  }

 hist("N_jets_pt20")->Fill(Njets20, weight);
 hist("N_jets_pt50")->Fill(Njets50, weight);
 hist("N_jets_pt70")->Fill(Njets70, weight);
 hist("N_jets_pt100")->Fill(Njets100, weight);

  Njets = event.topjets->size();

  if(Njets>=1){
    hist("eta_jet1")->Fill(jets->at(0).eta(), weight);
  }
  if(Njets>=2){
    hist("eta_jet2")->Fill(jets->at(1).eta(), weight);
  }
  if(Njets>=3){
    hist("eta_jet3")->Fill(jets->at(2).eta(), weight);
  }
  if(Njets>=4){
    hist("eta_jet4")->Fill(jets->at(3).eta(), weight);
  }
 
  int Nmuons = event.muons->size();
  hist("N_mu")->Fill(Nmuons, weight);
  for (const Muon & thismu : *event.muons){
    hist("pt_mu")->Fill(thismu.pt(), weight);
    hist("eta_mu")->Fill(thismu.eta(), weight);
    hist("reliso_mu")->Fill(thismu.relIso(), weight);
  }
  
  int Npvs = event.pvs->size();
  hist("N_pv")->Fill(Npvs, weight);

  //HT

  auto met = event.met->pt();
  hist("E_Tmiss")->Fill(met, weight);
  double ht = 0.0;
  double ht_jets = 0.0;
  double ht_lep = 0.0;
  for(const auto & jet : *jets){
    //ht += jet.pt();
    ht_jets += jet.pt();
  }

  for(const auto & electron : *event.electrons){
    ht_lep += electron.pt();
  }
  for(const auto & muon : *event.muons){
    ht_lep += muon.pt();
  }
  /*for(const auto & tau : *event.taus){
    ht_lep += tau.pt();
    }*/
  ht = ht_lep + ht_jets + met;
  //ht = met;
  hist("H_T")->Fill(ht, weight);
  hist("H_T1")->Fill(ht, weight);

  hist("eff")-> Fill(1, weight);
 
  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);

    if(Njets>=1){ 
      double deltaRH = deltaR(jets->at(0),  zprimegen.Higgs());
      hist ("Higgs_matched") ->Fill(deltaRH, weight);

      double deltaRTop1 = deltaR(jets->at(0),  zprimegen.ATop());
      double deltaRMainTop1 = deltaR(jets->at(0),  zprimegen.Top());

      if(deltaRH <= 0.8){
	hist("MissMatched")->Fill("H matched",weight);
	hist("Tau32H")->Fill(jets->at(0).tau3() /jets->at(0).tau2() ,weight);
	hist("Tau21H")->Fill(jets->at(0).tau2() /jets->at(0).tau1() ,weight);
	//	std::cout << "(Hists) L:163 Masse Higgs und DeltaRH "<<jets->at(0).v4().M() <<"    "<< deltaRH << std::endl;
	//	std::cout << "(Hists) L:163 Und zum ATOp "<< deltaRTop1 << std::endl;
      }else{
	hist("MissMatched")->Fill("H miss",weight);
      }
      if(deltaRTop1 <= 0.8){
	hist("MissMatched")->Fill("Top1 matched",weight);
      }else{
	hist("MissMatched")->Fill("Top1 miss",weight);
      }
      if(deltaRMainTop1 <= 0.8){
	hist("MissMatched")->Fill("MainTop1 matched",weight);
      }else{
	hist("MissMatched")->Fill("MainTop1 miss",weight);
      }
    }

    if(Njets>=2){
      double deltaRWlep = deltaR(jets->at(1),  zprimegen.WLep());
      hist ("Wlep_matched") ->Fill(deltaRWlep, weight);

      double deltaRWh = deltaR(jets->at(1),  zprimegen.WHiggsTop());
      hist ("Wh_matched") ->Fill(deltaRWh, weight);
 
      double deltaRT = deltaR(jets->at(1),  zprimegen.Top());
      hist ("T_matched") ->Fill(deltaRT, weight);
 
      double deltaRBlep = deltaR(jets->at(1),  zprimegen.BLep());
      hist ("Blep_matched") ->Fill(deltaRBlep, weight);
 
      double deltaRBh = deltaR(jets->at(1),  zprimegen.BHiggsTop());
      hist ("Bh_matched") ->Fill(deltaRBh, weight);
 
      double deltaRTop = deltaR(jets->at(1),  zprimegen.hadTop());
      hist ("Top_matched") ->Fill(deltaRTop, weight);

      double deltaRH2 = deltaR(jets->at(1),  zprimegen.Higgs());
  
     

    
      if(deltaRWlep <= 0.8){
	hist("MissMatched")->Fill("Wl matched",weight);
      }else{
	hist("MissMatched")->Fill("Wl miss",weight);
      }
      if(deltaRT <= 0.8){
	hist("MissMatched")->Fill("T matched",weight);
      }else{
	hist("MissMatched")->Fill("T miss",weight);
      }
      if(deltaRBlep <= 0.8){
	hist("MissMatched")->Fill("Bl matched",weight);
      }else{
	hist("MissMatched")->Fill("Bl miss",weight);
      }
      if(deltaRBh <= 0.8){
	hist("MissMatched")->Fill("Bh matched",weight);
      }else{
	hist("MissMatched")->Fill("Bh miss",weight);
      }
      if(deltaRWh <= 0.8){
	hist("MissMatched")->Fill("Wh matched",weight);
      }else{
	hist("MissMatched")->Fill("Wh miss",weight);
      }
      if(deltaRTop <= 0.8){
	hist("MissMatched")->Fill("Top matched",weight);
	hist("Tau32T")->Fill(jets->at(1).tau3() /jets->at(1).tau2() ,weight);
	hist("Tau21T")->Fill(jets->at(1).tau2() /jets->at(1).tau1() ,weight);
      }else{
	hist("MissMatched")->Fill("Top miss",weight);
      }
      if(deltaRH2 <= 0.6){
	hist("MissMatched")->Fill("H2 matched",weight);
      }else{
	hist("MissMatched")->Fill("H2 miss",weight);
      }
    }


  }//h_zprimegen


  //TOPTAGEFF
  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);
    if(event.is_valid(h_toptag)){
      const auto & toptagevt = event.get(h_toptag);
      if( toptagevt.size() >=1){

	const GenParticle & genp = event.genparticles->at(2);
	float drmin = deltaR(genp.v4(),toptagevt.at(0));
	int jmin =2; 
	for(unsigned int i = 2; i<event.genparticles->size();i++){
	  const auto & genp = event.genparticles->at(i);
	  float dr = deltaR(genp.v4(),toptagevt.at(0));
	  if(dr < drmin){drmin = dr;jmin = i;}
	}

	hist("pdgIdTopTag")->Fill(abs(event.genparticles->at(jmin).pdgId()),1);


	if(deltaR( toptagevt.at(0).v4(), zprimegen.hadTop())<= 0.8){
	  hist("MatchedTopTag")->Fill("Toptag matched",1);
	}else hist("MatchedTopTag")->Fill("Toptag miss",1);
	if(toptagevt.size()>1){
	if(deltaR( toptagevt.at(1).v4(), zprimegen.hadTop())<= 0.8){
	  hist("MatchedTopTag")->Fill("Toptag1 matched",1);
	}else hist("MatchedTopTag")->Fill("Toptag1 miss",1);
	}
	if( deltaR(toptagevt.at(0).v4(), zprimegen.Higgs())<=0.8)  hist("MatchedTopTag")->Fill("Higgs matched",1);
	if( deltaR(toptagevt.at(0).v4(), zprimegen.ZBoson())<=0.8)  hist("MatchedTopTag")->Fill("Z matched",1);
	if(deltaR(toptagevt.at(0).v4(), zprimegen.WBoson())<= 0.8) hist("MatchedTopTag")->Fill("W matched",1);
      }
    }//h_toptag
  }
  
  //HIGGSTAG
  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);
   
    if(event.is_valid(h_higgstag)){
      const auto & higgstagevt = event.get(h_higgstag);
     
      if(higgstagevt.size()>=1){
	const GenParticle & genp = event.genparticles->at(2);
	float drmin = deltaR(genp.v4(),higgstagevt.at(0));
	int jmin =2; 
	for(unsigned int i = 2; i<event.genparticles->size();i++){
	  const auto & genp = event.genparticles->at(i);
	  float dr = deltaR(genp.v4(),higgstagevt.at(0));
	  if(dr < drmin){drmin = dr;jmin = i;}
	}

	hist("pdgIdHiggsTag")->Fill(abs(event.genparticles->at(jmin).pdgId()),1);

	if(deltaR( higgstagevt.at(0).v4(), zprimegen.Higgs())<= 0.8){
	  hist("MatchedHiggsTag")->Fill("Higgstag matched",1);
	}else hist("MatchedHiggsTag")->Fill("Higgstag miss",1);

	if( deltaR(higgstagevt.at(0).v4(), zprimegen.hadTop())<=0.8)  hist("MatchedHiggsTag")->Fill("Top matched",1);
	if( deltaR(higgstagevt.at(0).v4(), zprimegen.WBoson())<=0.8)  hist("MatchedHiggsTag")->Fill("WBoson (T') matched",1);
	if( deltaR(higgstagevt.at(0).v4(), zprimegen.ZBoson())<=0.8)  hist("MatchedHiggsTag")->Fill("Z matched",1);

      }
    }//h_higgstag
  }

  
  //ZWTAG
  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);
 
    if(event.is_valid(h_ZWtag)){
      const auto & ZWtagevt = event.get(h_ZWtag);
     
      if( ZWtagevt.size() >=1){

	
	const GenParticle & genp = event.genparticles->at(2);
	float drmin = deltaR(genp.v4(),ZWtagevt.at(0));
	int jmin =2; 
	for(unsigned int i = 2; i<event.genparticles->size();i++){
	  const auto & genp = event.genparticles->at(i);
	  float dr = deltaR(genp.v4(),ZWtagevt.at(0));
	  if(dr < drmin){drmin = dr;jmin = i;}
	}

	hist("pdgIdZWTag")->Fill(abs(event.genparticles->at(jmin).pdgId()),1);



	if(deltaR( ZWtagevt.at(0).v4(), zprimegen.Higgs())<= 0.8)	hist("MatchedZWTag")->Fill("Higgs matched",1);
	if(deltaR( ZWtagevt.at(0).v4(), zprimegen.WBoson())<= 0.8)	hist("MatchedZWTag")->Fill("WBoson (TPrime) matched",1);
	if(deltaR(ZWtagevt.at(0).v4(), zprimegen.hadTop())<=0.8)  hist("MatchedZWTag")->Fill("had. Top matched",1);
	if(deltaR( ZWtagevt.at(0).v4(), zprimegen.TopW())<= 0.8)	hist("MatchedZWTag")->Fill("W aus Top matched",1);
	if(deltaR( ZWtagevt.at(0).v4(), zprimegen.WHiggsTop())<= 0.8)	hist("MatchedZWTag")->Fill("W aus top vom Higgs matched",1);
	if(ZWtagevt.size()>1){	
	  if(deltaR(ZWtagevt.at(1).v4(), zprimegen.WLep())<=0.8)  hist("MatchedZWTag")->Fill("Wlep2 matched",1);
	if(deltaR(ZWtagevt.at(1).v4(), zprimegen.WHad())<=0.8)  hist("MatchedZWTag")->Fill("Whad2 matched",1);
	  hist("MatchedZWTag")->Fill("2nd zwtag",1);
	  if(deltaR( ZWtagevt.at(1).v4(), zprimegen.TopW())<= 0.8)	hist("MatchedZWTag")->Fill("W aus Top(2 ztag) matched",1);
	  if(deltaR( ZWtagevt.at(1).v4(), zprimegen.WHiggsTop())<= 0.8)	hist("MatchedZWTag")->Fill("W aus top vom Higgs(2 ztag) matched",1);
	  if(deltaR(ZWtagevt.at(1).v4(), zprimegen.ZBoson())<=0.8) {
	    hist("MatchedZWTag")->Fill("Ztag (2tag)matched",1);
	  }else hist("MatchedZWTag")->Fill("Ztag (2tag)miss",1);
	}	
	if(deltaR(ZWtagevt.at(0).v4(), zprimegen.WLep())<=0.8)  hist("MatchedZWTag")->Fill("Wlep matched",1);
	if(deltaR(ZWtagevt.at(0).v4(), zprimegen.WHad())<=0.8)  hist("MatchedZWTag")->Fill("Whad matched",1);
	if(deltaR(ZWtagevt.at(0).v4(), zprimegen.ZBoson())<=0.8) {
	  hist("MatchedZWTag")->Fill("Ztag matched",1);
	}else hist("MatchedZWTag")->Fill("Ztag miss",1);

      }
    }//h_ZWtag
  }


  if(event.is_valid(h_zprimegen)){
    const auto & zprimegen = event.get(h_zprimegen);
    if( event.muons->size()>0){
    const auto & muon = event.muons->at(0);
    double distance1 = deltaR(muon, zprimegen.muon1());
    double distance2 = deltaR(muon, zprimegen.muon2());
    hist("MuonDistance1")->Fill(distance1,weight);
    hist("MuonDistance2")->Fill(distance2,weight);
    if(deltaR(muon, zprimegen.muon1())<0.1) hist("MuonMatch")->Fill("Muon aus top",weight);
    else if(deltaR(muon, zprimegen.muon2())<0.1) hist("MuonMatch")->Fill("Muon aus TPrime",weight);
    else  hist("MuonMatch")->Fill("Muon aus Higgs",weight);
    }
  }

  int number_of_btag=0;
  for(auto const & jet: *event.jets){
    if( jet.btag_combinedSecondaryVertex()>0.800f) number_of_btag++;
  }

  hist("number_btag")->Fill(number_of_btag, weight);


  //matched W in ttbar
  if(event.is_valid(h_ttbargen)&&event.is_valid(h_ZWtag) ){
    const auto & ttbargen = event.get(h_ttbargen);
    if( ttbargen.IsSemiLeptonicDecay()){
      const auto & zwtag = event.get(h_ZWtag);
      if(zwtag.size()){
	double distance_WLep = deltaR(zwtag.at(0), ttbargen.WLep());
	double distance_WHad = deltaR(zwtag.at(0), ttbargen.WHad());
	if(distance_WLep <= 0.8 || distance_WHad <= 0.8 ) hist("matchedWTTbar")->Fill("matched W",1);
	else hist("matchedWTTbar")->Fill("not matched",1);
      }
    }
  }
}

ZPrimeTotTPrimeHists::~ZPrimeTotTPrimeHists(){}
