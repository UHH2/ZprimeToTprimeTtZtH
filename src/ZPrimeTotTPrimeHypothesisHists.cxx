#include "TH1F.h"
#include "TH2F.h"
 #include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeHypothesisHists.h" 
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"

using namespace uhh2;

ZPrimeTotTPrimeHypothesisHists::ZPrimeTotTPrimeHypothesisHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name ): Hists(ctx, dirname){

  TString name = discriminator_name;
    double min=0;
    double max=200;
    if(discriminator_name=="Chi2"){
      name = "#Chi^{2}";
    }
    else{
      name += " discriminator";
    }

    if( discriminator_name=="CorrectMatch"){
      min=0;
      max=2;
    }
    if( discriminator_name=="TopDRMC"){
      min=0;
      max=6;
    }

    Discriminator = book<TH1F>("Discriminator",name,100,min,max);
    Discriminator_2 = book<TH1F>("Discriminator_2",name,50,0,50);
    Discriminator_3 = book<TH1F>("Discriminator_3",name,300,0,30); 
    DiscriminatorH = book<TH1F>("DiscriminatorH",name+"H",100,min,max);
    Discriminator_mismatchedHadTop = book<TH1F>("Discriminator_mismatchedHadTops",name,100,min,max);
    Discriminator_mismatchedHiggs = book<TH1F>("Discriminator_mismatchedHiggs",name,100,min,max);

    M_ttbar_rec = book<TH1F>( "M_ttbar_rec", "M_{t#bar{t}}^{rec} [GeV/c^{2}]", 100, 0, 5000 ) ;
    M_Higgs_rec = book<TH1F>( "M_Higgs_rec", "M_{Higgs}^{rec} [GeV/c^{2}]", 60, 50, 170 ) ;
    M_ZPrime_rec = book<TH1F>( "M_ZPrime_rec", "M_{ZPrime}^{rec} [GeV/c^{2}]", 50, 600, 5000 ) ;
    M_ZPrime_scaled = book<TH1F>( "M_ZPrime_scaled", "M_{ZPrime}^{rec} [GeV/c^{2}]", 40, 600, 3000 ) ;
    M_TPrime_rec = book<TH1F>( "M_TPrime_rec", "M_{TPrime}^{rec} [GeV/c^{2}]", 50, 500, 1600 ) ;
    M_TPrime_rec2 = book<TH1F>( "M_TPrime_rec2", "M_{TPrime}^{rec} [GeV/c^{2}]", 50, 500, 1600 ) ;
    M_TPrime_rec3 = book<TH1F>( "M_TPrime_rec3", "M_{TPrime}^{rec} [GeV/c^{2}]", 50, 500, 1600 ) ;
    M_TPrime_rec4 = book<TH1F>( "M_TPrime_rec4", "M_{TPrime}^{rec} [GeV/c^{2}]", 50, 500, 1600 ) ;

 
    M_ttbar_gen = book<TH1F>( "M_ttbar_gen", "M_{t#bar{t}}^{gen} [GeV/c^{2}]", 100, 0, 5000 ) ;
    M_Higgs_gen = book<TH1F>( "M_Higgs_gen", "M_{ H}^{gen} [GeV/c^{2}]", 300, 0, 300 ) ;
    M_Z_gen = book<TH1F>( "M_Z_gen", "M_{Z}^{gen} [GeV/c^{2}]", 300, 0, 300 ) ;
    M_W_gen = book<TH1F>( "M_W_gen", "M_{W}^{gen} [GeV/c^{2}]", 300, 0, 300 ) ;
    M_ZPrime_gen = book<TH1F>( "M_ZPrime_gen", "M_{ZPrime}^{gen} [GeV/c^{2}]", 200, 500, 3000 ) ;
    M_TPrime_gen = book<TH1F>( "M_TPrime_gen", "M_{TPrime}^{gen} [GeV/c^{2}]", 500, 0, 3000 ) ;


    M_toplep_rec = book<TH1F>( "M_toplep_rec", "M^{top,lep} [GeV/c^{2}]", 100, 0, 300 ) ;
    M_tophad_rec = book<TH1F>( "M_tophad_rec", "M^{top,had} [GeV/c^{2}]", 100, 0, 500 ) ;
    M_toplep_gen = book<TH1F>( "M_toplep_gen", "M^{top,lep}_{gen} [GeV/c^{2}]", 100, 0, 300 ) ;
    M_tophad_gen = book<TH1F>( "M_tophad_gen", "M^{top,had}_{gen} [GeV/c^{2}]", 100, 0, 500 ) ;
  
    M_tophad_rec_1jet = book<TH1F>( "M_tophad_rec_1jet", "M^{top,had}_{1jet} [GeV/c^{2}]", 100, 0, 500 ) ;
    M_tophad_rec_2jet = book<TH1F>( "M_tophad_rec_2jet", "M^{top,had}_{2jet} [GeV/c^{2}]", 100, 0, 500 ) ;
    M_tophad_rec_3jet = book<TH1F>( "M_tophad_rec_3jet", "M^{top,had}_{3jet} [GeV/c^{2}]", 100, 0, 500 ) ;

    M_tophad_rec_1jet_pdgid = book<TH1F>( "M_tophad_rec_1jet_pdgid", "M^{top,had}_{1jet} [GeV/c^{2}]", 22, 0, 22 ) ;
    M_tophad_rec_2jet_pdgid = book<TH1F>( "M_tophad_rec_2jet_pdgid", "M^{top,had}_{2jet} [GeV/c^{2}]", 22, 0, 22 ) ;
    M_tophad_rec_3jet_pdgid = book<TH1F>( "M_tophad_rec_3jet_pdgid", "M^{top,had}_{3jet} [GeV/c^{2}]", 22, 0, 22 ) ;
    M_tophad_rec_1jet_pdgid->Fill("Gluon",0);
    M_tophad_rec_1jet_pdgid->Fill("Quark",0);
    M_tophad_rec_2jet_pdgid->Fill("Gluon",0);
    M_tophad_rec_2jet_pdgid->Fill("Quark",0);
    M_tophad_rec_3jet_pdgid->Fill("Gluon",0);
    M_tophad_rec_3jet_pdgid->Fill("Quark",0);
    
    Pt_toplep_rec = book<TH1F>( "Pt_toplep_rec", "P_{T}^{top,lep} [GeV/c]", 60, 0, 1200 ) ;
    Pt_tophad_rec = book<TH1F>( "Pt_tophad_rec", "P_{T}^{top,had} [GeV/c]", 60, 0, 1200 ) ;
    //ttbar PT
    Pt_ttbar_rec = book<TH1F>( "Pt_ttbar_rec", "P_{T,t#bar{t}}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_ttbar_gen = book<TH1F>( "Pt_ttbar_gen", "P_{T,t#bar{t}}^{gen} [GeV/c]", 60, 0, 600 ) ;
   //Higgs PT
    Pt_Higgs_rec = book<TH1F>( "Pt_Higgs_rec", "P_{T,Higgs}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_Higgs_gen = book<TH1F>( "Pt_Higgs_gen", "P_{T,Higgs}^{gen} [GeV/c]", 60, 0, 600 ) ;
   //ZPrime PT
    Pt_ZPrime_rec = book<TH1F>( "Pt_ZPrime_rec", "P_{T,ZPrime}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_ZPrime_gen = book<TH1F>( "Pt_ZPrime_gen", "P_{T,ZPrime}^{gen} [GeV/c]", 60, 0, 600 ) ;
   //TPrime PT
    Pt_TPrime_rec = book<TH1F>( "Pt_TPrime_rec", "P_{T,TPrime}^{rec} [GeV/c]", 60, 0, 600 ) ;
    Pt_TPrime_gen = book<TH1F>( "Pt_TPrime_gen", "P_{T,TPrime}^{gen} [GeV/c]", 60, 0, 600 ) ;


    ///////////////////////////////////ETA//////////////////////////////////////////////////////////////////////////////////
    eta_toplep_rec = book<TH1F>( "eta_toplep_rec", "#eta^{top,lep}_rec [GeV/c]", 100, -5, 5 ) ;
    eta_tophad_rec = book<TH1F>( "eta_tophad_rec", "#eta^{top,had}_rec [GeV/c]", 100, -5, 5 ) ;
    eta_toplep_gen = book<TH1F>( "eta_toplep_gen", "#eta^{top,lep}_gen [GeV/c]", 100, -5, 5 ) ;
    eta_tophad_gen = book<TH1F>( "eta_tophad_gen", "#eta^{top,had}_gen [GeV/c]", 100, -5, 5 ) ;
   
   
   //Higgs PT
    eta_Higgs_rec = book<TH1F>( "eta_Higgs_rec", "#eta_{Higgs}^{rec} [GeV/c]", 100, -5, 5 ) ;
    eta_Higgs_gen = book<TH1F>( "eta_Higgs_gen", "#eta_{Higgs}^{gen} [GeV/c]", 100, -5, 5 ) ;
   //ZPrime PT
    eta_ZPrime_rec = book<TH1F>( "eta_ZPrime_rec", "#eta_{ZPrime}^{rec} [GeV/c]", 100, -5, 5 ) ;
    eta_ZPrime_gen = book<TH1F>( "eta_ZPrime_gen", "#eta_{ZPrime}^{gen} [GeV/c]", 100, -5, 5 ) ;
   //TPrime PT
    eta_TPrime_rec = book<TH1F>( "eta_TPrime_rec", "#eta_{TPrime}^{rec} [GeV/c]", 100, -5, 5 ) ;
    eta_TPrime_gen = book<TH1F>( "eta_TPrime_gen", "#eta_{TPrime}^{gen} [GeV/c]", 100, -5, 5 ) ;

    //////////////////////////////////////////////////////////Rec vs. Gen/////////////////////////////////////////////////////
    M_Higgs_vs = book<TH1F>( "M_Higgs_vs", "(M_{H}^{rec} - M_{H}^{gen})/ M_{H}^{gen}  [GeV/c^{2}]", 100, -1, 1 ) ;
    M_Z_vs = book<TH1F>( "M_Z_vs", "(M_{Z}^{rec} - M_{Z}^{gen})/ M_{Z}^{gen}  [GeV/c^{2}]", 100, -1, 1 ) ;
    M_W_vs = book<TH1F>( "M_W_vs", "(M_{W}^{rec} - M_{W}^{gen})/ M_{W}^{gen}  [GeV/c^{2}]", 100, -1, 1 ) ;
    M_ZPrime_vs = book<TH1F>( "M_ZPrime_vs", "(M_{Z'}^{rec} - M_{Z'}^{gen})/ M_{Z'}^{gen} [GeV/c^{2}]", 100,-1, 1 ) ;
    M_TPrime_vs = book<TH1F>( "M_TPrime_vs", "(M_{T'}^{rec} - M_{T'}^{gen})/ M_{T'}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    M_toplep_vs = book<TH1F>( "M_toplep_vs", "(M_{top,lep}^{rec} - M_{top,lep}^{gen})/ M_{top,lep}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    M_tophad_vs = book<TH1F>( "M_tophad_vs", "(M_{top,had}^{rec} - M_{top,had}^{gen})/ M_{top,had}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;

    Pt_Higgs_vs = book<TH1F>( "Pt_Higgs_vs", "(P_{T,H}^{rec} - P_{T,H}^{gen})/ P_{T,H}^{gen}  [GeV/c^{2}]", 100, -1, 1 ) ;
    Pt_ZPrime_vs = book<TH1F>( "Pt_ZPrime_vs", "(P_{T,Z'}^{rec} - P_{T,Z'}^{gen})/ P_{T,Z'}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    Pt_TPrime_vs = book<TH1F>( "Pt_TPrime_vs", "(P_{T,T'}^{rec} - P_{T,T'}^{gen})/ P_{T,T'}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    Pt_toplep_vs = book<TH1F>( "Pt_toplep_vs", "(P_{T,top,lep}^{rec} - P_{T,top,lep}^{gen})/ P_{T,top,lep}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    Pt_tophad_vs = book<TH1F>( "Pt_tophad_vs", "(P_{T,top,had}^{rec} - P_{T,top,had}^{gen})/ P_{T,top,had}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;

    eta_Higgs_vs = book<TH1F>( "eta_Higgs_vs", "(#eta_{H}^{rec} - #eta_{H}^{gen})/ #eta_{H}^{gen}  [GeV/c^{2}]", 100, -1, 1 ) ;
    eta_ZPrime_vs = book<TH1F>( "eta_ZPrime_vs", "(#eta_{Z'}^{rec} - #eta_{Z'}^{gen})/ #eta_{Z'}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    eta_TPrime_vs = book<TH1F>( "eta_TPrime_vs", "(#eta_{T'}^{rec} - #eta_{T'}^{gen})/ #eta_{T'}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    eta_toplep_vs = book<TH1F>( "eta_toplep_vs", "(#eta_{top,lep}^{rec} - #eta_{top,lep}^{gen})/ #eta_{top,lep}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;
    eta_tophad_vs = book<TH1F>( "eta_tophad_vs", "(#eta_{top,had}^{rec} - #eta_{top,had}^{gen})/ #eta_{top,had}^{gen} [GeV/c^{2}]", 100, -1, 1 ) ;


    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    Pt_ttbar_rec_vs_Pt_ttbar_gen = book<TH2F>( "Pt_ttbar_rec_vs_Pt_ttbar_gen", "P_{T,t#bar{t}}^{rec} [GeV/c] vs P_{T,t#bar{t}}^{gen} [GeV/c]", 60, 0, 600 ,60, 0, 600);
    Pt_Higgs_rec_vs_Pt_Higgs_gen = book<TH2F>( "Pt_Higgs_rec_vs_Pt_Higgs_gen", "P_{T,Higgs}^{rec} [GeV/c] vs P_{T,Higgs}^{gen} [GeV/c]", 60, 0, 600 ,60, 0, 600);
    Pt_ZPrime_rec_vs_Pt_ZPrime_gen = book<TH2F>( "Pt_ZPrime_rec_vs_Pt_ZPrime_gen", "P_{T,ZPrime}^{rec} [GeV/c] vs P_{T,ZPrime}^{gen} [GeV/c]", 60, 0, 600 ,60, 0, 600);
    Pt_TPrime_rec_vs_Pt_TPrime_gen = book<TH2F>( "Pt_TPrime_rec_vs_Pt_ZPrime_gen", "P_{T,TPrime}^{rec} [GeV/c] vs P_{T,TPrime}^{gen} [GeV/c]", 60, 0, 600 ,60, 0, 600);


    Discriminator_vs_M_ttbar = book<TH2F>("Discriminator_vs_M_ttbar", name+" vs M_{t#bar{t}}^{rec}" , 50, min,max, 100,0,5000) ;

    M_ttbar_rec_vs_M_ttbar_gen = book<TH2F>("M_ttbar_rec_vs_M_ttbar_gen","M_{t#bar{t}}^{rec} [GeV/c^{2}] vs M_{t#bar{t}}^{gen} [GeV/c^{2}]",100,0,5000,100,0,5000);
    M_ttbar_gen_vs_M_ttbar_rec = book<TH2F>("M_ttbar_gen_vs_M_ttbar_rec","M_{t#bar{t}}^{gen} [GeV/c^{2}] vs M_{t#bar{t}}^{rec} [GeV/c^{2}]",100,0,5000,100,0,5000);

    M_Higgs_rec_vs_M_Higgs_gen = book<TH2F>("M_Higgs_rec_vs_M_Higgs_gen","M_{Higgs}^{rec} [GeV/c^{2}] vs M_{Higgs}^{gen} [GeV/c^{2}]",100,0,5000,100,0,5000);
    M_Higgs_gen_vs_M_Higgs_rec = book<TH2F>("M_Higgs_gen_vs_M_Higgs_rec","M_{Higgs}^{gen} [GeV/c^{2}] vs M_{Higgs}^{rec} [GeV/c^{2}]",100,0,5000,100,0,5000);

    M_ZPrime_rec_vs_M_ZPrime_gen = book<TH2F>("M_ZPrime_rec_vs_M_ZPrime_gen","M_{ZPrime}^{rec} [GeV/c^{2}] vs M_{ZPrime}^{gen} [GeV/c^{2}]",100,0,5000,100,0,5000);
    M_ZPrime_gen_vs_M_ZPrime_rec = book<TH2F>("M_ZPrime_gen_vs_M_ZPrime_rec","M_{ZPrime}^{gen} [GeV/c^{2}] vs M_{ZPrime}^{rec} [GeV/c^{2}]",100,0,5000,100,0,5000);

    M_TPrime_rec_vs_M_TPrime_gen = book<TH2F>("M_TPrime_rec_vs_M_TPrime_gen","M_{TPrime}^{rec} [GeV/c^{2}] vs M_{TPrime}^{gen} [GeV/c^{2}]",100,0,5000,100,0,5000);
    M_TPrime_gen_vs_M_TPrime_rec = book<TH2F>("M_TPrime_gen_vs_M_TPrime_rec","M_{TPrime}^{gen} [GeV/c^{2}] vs M_{TPrime}^{rec} [GeV/c^{2}]",100,0,5000,100,0,5000);
    
    
    M_ttbar_resolution = book<TH1F>("M_ttbar_resolution", "(M_{t#bar{t}}^{gen} - M_{t#bar{t}}^{rec})/M_{t#bar{t}}^{rec}", 100, -5,5) ;
  
    Pt_toplep_rec_vs_Pt_toplep_gen = book<TH2F>("Pt_toplep_rec_vs_Pt_toplep_gen","P_{T}^{top,lep,rec} [GeV/c] vs P_{T}^{top,lep,gen} [GeV/c]",60, 0, 1200,60, 0, 1200);
    Pt_tophad_rec_vs_Pt_tophad_gen = book<TH2F>("Pt_tophad_rec_vs_Pt_tophad_gen","P_{T}^{top,had,rec} [GeV/c] vs P_{T}^{top,had,gen} [GeV/c]",60, 0, 1200,60, 0, 1200);
    Pt_lepton_rec_vs_Pt_lepton_gen = book<TH2F>("Pt_lepton_rec_vs_Pt_lepton_gen","P_{T}^{lepton,rec} [GeV/c] vs P_{T}^{lepton,gen} [GeV/c]",60, 0, 800,60, 0, 800);
    Pt_blep_rec_vs_Pt_blep_gen = book<TH2F>("Pt_blep_rec_vs_Pt_blep_gen","P_{T}^{blep,rec} [GeV/c] vs P_{T}^{blep,gen} [GeV/c]",60, 0, 800,60, 0, 800);
    Pt_neutrino_rec_vs_Pt_neutrino_gen = book<TH2F>("Pt_neutrino_rec_vs_Pt_neutrino_gen","P_{T}^{#nu,rec} [GeV/c] vs P_{T}^{#nu,gen} [GeV/c]",60, 0, 800,60, 0, 800);
    
    eta_toplep_rec_vs_eta_toplep_gen = book<TH2F>("eta_toplep_rec_vs_eta_toplep_gen","#eta^{top,lep,rec} vs #eta^{top,lep,gen}",100,-5,5,100,-5,5);
    eta_tophad_rec_vs_eta_tophad_gen = book<TH2F>("eta_tophad_rec_vs_eta_tophad_gen","#eta^{top,had,rec} vs #eta^{top,had,gen}",100,-5,5,100,-5,5);
    eta_lepton_rec_vs_eta_lepton_gen = book<TH2F>("eta_lepton_rec_vs_eta_lepton_gen","#eta^{lepton,rec} vs #eta^{lepton,gen}",100,-5,5,100,-5,5);
    eta_blep_rec_vs_eta_blep_gen = book<TH2F>("eta_blep_rec_vs_eta_blep_gen","#eta^{blep,rec} vs #eta^{blep,gen}",100,-5,5,100,-5,5 );
    eta_neutrino_rec_vs_eta_neutrino_gen = book<TH2F>("eta_neutrino_rec_vs_eta_neutrino_gen","#eta^{#nu,rec} vs #eta^{#nu,gen}",100,-5,5,100,-5,5 );


    ///////////////////////////////////////////////////         Matching tests //////////////////////////////////////////////////////////////////////////////////////////////
    massfitHiggs = book<TH1F>("massfitHiggs","Mass reco Higgs",200, 50, 400);
    massfitZ = book <TH1F>("massfitZ","Mass reco Z",200, 50, 400);
    massfithadTop = book<TH1F>("massfithadTop","Mass reco hadTop",200, 50, 400);
    massfitlepTop = book<TH1F>("massfitlepTop","Mass reco lepTop",200, 50, 400);
 
    missmatch = book<TH1F>("missmatch","Match Rate",15,0,15);

    //distance between lepton and bjet of the lep top
    deltar_mu_blep = book<TH1F>("deltar_mu_blep","Distance between the lepton and the jet used for top lep",50,0,5);
    number_closet_jet = book<TH1F>("number_closet_jet","Number of jets which are closest, number of all jets",3,0,3);
    number_closet_jet->Fill("#closest",0);
    number_closet_jet->Fill("#reconstruction",0);

    number_closet_jet_with_matching = book<TH1F>("number_closet_jet_with_matching","Number of jets which are closest, number of all jets",3,0,3);
    number_closet_jet_with_matching->Fill("#closest",0);
    number_closet_jet_with_matching->Fill("#reconstruction",0);

       
    h_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(hyps_name);
    h_zprimegen = ctx.get_handle<ZPrimeGen>("zprimegen");
    h_ttbargen = ctx.get_handle<TTbarGen>("ttbargen");
    m_discriminator_name = discriminator_name;
}


void ZPrimeTotTPrimeHypothesisHists::fill(const uhh2::Event & e){

  std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = e.get(h_hyps);
  const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );
  //  const ZPrimeTotTPrimeReconstructionHypothesis* hypHiggs = get_best_hypothesis( hyps, m_discriminator_name+"H" );
  double weight = e.weight;

  double mttbar_rec = 0;
  double mHiggs_rec = 0;
  double mZPrime_rec = 0;
  double mTPrime_rec = 0;
  double mTPrime_rec2 = 0;
  double mTPrime_rec3 = 0;
 double mTPrime_rec4 = 0;

  //ttbar Mass
  if( (hyp->top_v4()+hyp->antitop_v4()).isTimelike() )
    mttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).M();
  else{
    mttbar_rec = sqrt( -(hyp->top_v4()+hyp->antitop_v4()).mass2());
  }
  //Higgs Mass
  if( (hyp->HZW_v4()).isTimelike() ){
    LorentzVector subjet_sum;
    for (const auto s : hyp->HZW_subjets()) {
      subjet_sum += s.v4();
    }
    mHiggs_rec=subjet_sum.M();
  }else{
    LorentzVector subjet_sum;
    for (const auto s : hyp->HZW_subjets()) {
      subjet_sum += s.v4();
    }
    mHiggs_rec=sqrt( -subjet_sum.mass2());
  }
  //ZPrime Mass
  if( (hyp->HZW_v4() + hyp->top_v4()+ hyp->antitop_v4()).isTimelike() )
    mZPrime_rec = (hyp->HZW_v4()+ hyp->top_v4()+ hyp->antitop_v4()).M();
  else{
    mZPrime_rec = sqrt( -(hyp->HZW_v4()+ hyp->top_v4()+ hyp->antitop_v4()).mass2());
  }
  //TPrime Mass
  if( (hyp->HZW_v4() + hyp->tophad_v4() ).isTimelike() ){
    mTPrime_rec = (hyp->HZW_v4()+ hyp->tophad_v4()).M();

    if((hyp->HZW_v4()+ hyp->tophad_v4()).M()>=(hyp->HZW_v4()+ hyp->toplep_v4()).M()) mTPrime_rec2 = (hyp->HZW_v4()+ hyp->tophad_v4()).M();
    else (hyp->HZW_v4()+ hyp->toplep_v4()).M();
    if(hyp->tophad_v4().Pt()>= hyp->toplep_v4().Pt()) mTPrime_rec3 = (hyp->HZW_v4()+ hyp->tophad_v4()).M();
    else  mTPrime_rec3 = (hyp->HZW_v4()+ hyp->toplep_v4()).M();
    if(deltaR((hyp->tophad_v4()+hyp->HZW_v4()),hyp->toplep_v4()) >= deltaR((hyp->toplep_v4()+hyp->HZW_v4()),hyp->tophad_v4()) ) mTPrime_rec4 = (hyp->HZW_v4()+ hyp->tophad_v4()).M();
    else  mTPrime_rec4 = (hyp->HZW_v4()+ hyp->toplep_v4()).M();
  } else{
    mTPrime_rec = sqrt( -(hyp->HZW_v4()+ hyp->tophad_v4()).mass2());
    if((hyp->HZW_v4()+ hyp->tophad_v4()).M()>=(hyp->HZW_v4()+ hyp->toplep_v4()).M())  mTPrime_rec2 = sqrt( -(hyp->HZW_v4()+ hyp->tophad_v4()).mass2());
    else  mTPrime_rec2 = sqrt( -(hyp->HZW_v4()+ hyp->toplep_v4()).mass2());
    if(hyp->tophad_v4().Pt()>= hyp->toplep_v4().Pt()) mTPrime_rec3 = sqrt( -(hyp->HZW_v4()+ hyp->tophad_v4()).mass2());
    else mTPrime_rec3 = sqrt( -(hyp->HZW_v4()+ hyp->toplep_v4()).mass2());
    if(deltaR((hyp->tophad_v4()+hyp->HZW_v4()),hyp->toplep_v4()) >= deltaR((hyp->toplep_v4()+hyp->HZW_v4()),hyp->tophad_v4()) )mTPrime_rec4 = sqrt( -(hyp->HZW_v4()+ hyp->tophad_v4()).mass2());
    else mTPrime_rec4 = sqrt( -(hyp->HZW_v4()+ hyp->toplep_v4()).mass2());
  }

  M_ttbar_rec->Fill(mttbar_rec, weight);
  M_Higgs_rec->Fill(mHiggs_rec, weight);
  M_ZPrime_rec->Fill(mZPrime_rec, weight);
  //hier
  // if(mZPrime_rec>4000){
  //   std::fstream g;
  //   g.open("masses.txt", ios::out);
  //   std::cout << e.run<<" "<<e.luminosityBlock<<" "<<e.event  << std::endl;
  //   g.close();
  // }

  M_ZPrime_scaled->Fill(mZPrime_rec, (0.048 - 0.000038 * hyp->HZW_v4().Pt() ) * weight);
  M_TPrime_rec->Fill(mTPrime_rec, weight);
 M_TPrime_rec2->Fill(mTPrime_rec2, weight);
 M_TPrime_rec3->Fill(mTPrime_rec3, weight);
M_TPrime_rec4->Fill(mTPrime_rec4, weight);

  double mtoplep_rec = (hyp->toplep_v4()).M();
  double mtophad_rec = (hyp->tophad_v4()).M();
  double mtoplep_gen = 0;
  double mtophad_gen = 0;

  double mttbar_gen = 0;
  double mHiggs_gen = 0;
  double mZ_gen = 0;
  double mW_gen = 0;
  double mZPrime_gen = 0;
  double mTPrime_gen = 0;


  double ptttbar_rec = (hyp->top_v4()+hyp->antitop_v4()).Pt();
  double ptttbar_gen = 0;

  double pttoplep_rec=(hyp->toplep_v4()).Pt();
  double pttophad_rec=(hyp->tophad_v4()).Pt();
  double pttoplep_gen=0;
  double pttophad_gen=0;
  double ptHiggs_rec = (hyp->HZW_v4()).Pt();
  double ptHiggs_gen = 0;
  double ptZPrime_rec = (hyp->HZW_v4() + hyp->top_v4()+ hyp->antitop_v4()).Pt();
  double ptZPrime_gen = 0;
  double ptTPrime_rec = (hyp->HZW_v4()+ hyp->tophad_v4()).Pt();
  double ptTPrime_gen = 0;


  //////ETA Variables/////
  double etatoplep_rec=(hyp->toplep_v4()).Eta();
  double etatoplep_gen=0;
  double etatophad_rec=(hyp->tophad_v4()).Eta();
  double etatophad_gen=0;
  double etaHiggs_rec = (hyp->HZW_v4()).Eta();
  double etaHiggs_gen = 0;
  double etaZPrime_rec = (hyp->HZW_v4() + hyp->top_v4()+ hyp->antitop_v4()).Eta();
  double etaZPrime_gen = 0;
  double etaTPrime_rec = (hyp->HZW_v4()+ hyp->tophad_v4()).Eta();
  double etaTPrime_gen = 0;

  //////

  Pt_ttbar_rec->Fill ( ptttbar_rec, weight);
  Pt_Higgs_rec->Fill ( ptHiggs_rec, weight);
  Pt_ZPrime_rec->Fill ( ptZPrime_rec, weight);
  Pt_TPrime_rec->Fill ( ptTPrime_rec, weight);
  
  if(e.is_valid(h_zprimegen)){
    const auto & zprimegen = e.get(h_zprimegen);
    mttbar_gen = ( zprimegen.Top().v4() + zprimegen.ATop().v4()).M();
    mHiggs_gen = ( zprimegen.Higgs().v4()).M();
    mZ_gen = ( zprimegen.ZBoson().v4()).M();
    mW_gen = ( zprimegen.WBoson().v4()).M();
    mZPrime_gen = ( zprimegen.ZP().v4()).M();
    mTPrime_gen = ( zprimegen.TP().v4()).M();
    mtoplep_gen = (zprimegen.lepTop().v4()).M();
    mtophad_gen =(zprimegen.hadTop().v4()).M();

    ptttbar_gen = (zprimegen.Top().v4() + zprimegen.ATop().v4()).Pt();

    pttoplep_gen=(zprimegen.lepTop().v4()).Pt();
    pttophad_gen=(zprimegen.hadTop().v4()).Pt();
    ptHiggs_gen = (zprimegen.Higgs().v4()).Pt();
    ptZPrime_gen = (zprimegen.ZP().v4()).Pt();
    ptTPrime_gen = (zprimegen.TP().v4()).Pt();

    //Eta
    etatoplep_gen=(zprimegen.lepTop().v4()).Eta();
    etatophad_gen=(zprimegen.hadTop().v4()).Eta();
    etaHiggs_gen = (zprimegen.Higgs().v4()).Eta();
    etaZPrime_gen = (zprimegen.ZP().v4()).Eta();
    etaTPrime_gen = (zprimegen.TP().v4()).Eta();


    M_toplep_gen->Fill(mtoplep_gen,weight);
    M_tophad_gen->Fill(mtophad_gen,weight);

    //ttbar Mass comparison
    M_ttbar_gen->Fill(mttbar_gen, weight);
    M_ttbar_rec_vs_M_ttbar_gen->Fill(mttbar_rec, mttbar_gen,weight);
    M_ttbar_gen_vs_M_ttbar_rec->Fill(mttbar_gen, mttbar_rec,weight);
    //Higgs Mass comparison
    M_Higgs_gen->Fill(mHiggs_gen, weight);
    M_Z_gen->Fill(mZ_gen, weight);
    M_W_gen->Fill(mW_gen, weight);
    M_Higgs_rec_vs_M_Higgs_gen->Fill(mHiggs_rec, mHiggs_gen,weight);
    M_Higgs_gen_vs_M_Higgs_rec->Fill(mHiggs_gen, mHiggs_rec,weight);
    //ZPrime Mass comparison
    M_ZPrime_gen->Fill(mZPrime_gen, weight);
    M_ZPrime_rec_vs_M_ZPrime_gen->Fill(mZPrime_rec, mZPrime_gen,weight);
    M_ZPrime_gen_vs_M_ZPrime_rec->Fill(mZPrime_gen, mZPrime_rec,weight);
    //TPrime Mass comparison
    M_TPrime_gen->Fill(mTPrime_gen, weight);
    M_TPrime_rec_vs_M_TPrime_gen->Fill(mTPrime_rec, mTPrime_gen,weight);
    M_TPrime_gen_vs_M_TPrime_rec->Fill(mTPrime_gen, mTPrime_rec,weight);

    M_ttbar_resolution->Fill( (mttbar_gen-mttbar_rec)/mttbar_gen, weight);
    //ttbar PT comparison
    Pt_ttbar_gen->Fill ( ptttbar_gen, weight);
    Pt_ttbar_rec_vs_Pt_ttbar_gen->Fill(ptttbar_rec, ptttbar_gen,weight);
    //Higgs PT comparison
    Pt_Higgs_gen->Fill ( ptHiggs_gen, weight);
    Pt_Higgs_rec_vs_Pt_Higgs_gen->Fill(ptHiggs_rec, ptHiggs_gen,weight);
    //ZPrime PT comparison
    Pt_ZPrime_gen->Fill ( ptZPrime_gen, weight);
    Pt_ZPrime_rec_vs_Pt_ZPrime_gen->Fill(ptZPrime_rec, ptZPrime_gen,weight);
    //TPrime PT comparison
    Pt_TPrime_gen->Fill ( ptTPrime_gen, weight);
    Pt_TPrime_rec_vs_Pt_TPrime_gen->Fill(ptTPrime_rec, ptTPrime_gen,weight);

    //Eta
    eta_toplep_rec  ->Fill(etatoplep_rec,weight);
    eta_tophad_rec ->Fill(etatophad_rec,weight);
    eta_toplep_gen ->Fill(etatoplep_gen,weight);
    eta_tophad_gen ->Fill(etatoplep_gen,weight);
    eta_Higgs_rec ->Fill(etaHiggs_rec,weight);
    eta_Higgs_gen ->Fill(etaHiggs_gen,weight);
    eta_ZPrime_rec ->Fill(etaZPrime_rec,weight);
    eta_ZPrime_gen ->Fill(etaZPrime_gen,weight);
    eta_TPrime_rec ->Fill(etaTPrime_rec,weight);
    eta_TPrime_gen ->Fill(etaTPrime_gen,weight);

    ////////////////////////////////////////////////// Rec vs. Gen ///////////////////////////////////////////////////////////////////
    M_Higgs_vs ->Fill((mHiggs_rec-mHiggs_gen)/mHiggs_gen,weight);
    M_Z_vs ->Fill((mHiggs_rec-mZ_gen)/mZ_gen,weight);
    M_W_vs ->Fill((mHiggs_rec-mW_gen)/mW_gen,weight);
    M_ZPrime_vs ->Fill((mZPrime_rec-mZPrime_gen)/mZPrime_gen,weight);
    M_TPrime_vs ->Fill((mTPrime_rec-mTPrime_gen)/mTPrime_gen,weight);
    M_toplep_vs ->Fill((mtoplep_rec-mtoplep_gen)/mtoplep_gen,weight);
    M_tophad_vs ->Fill((mtophad_rec-mtophad_gen)/mtophad_gen,weight);


    Pt_Higgs_vs ->Fill((ptHiggs_rec-ptHiggs_gen)/ptHiggs_gen,weight);
    Pt_ZPrime_vs ->Fill((ptZPrime_rec-ptZPrime_gen)/ptZPrime_gen,weight);
    Pt_TPrime_vs ->Fill((ptTPrime_rec-ptTPrime_gen)/ptTPrime_gen,weight);
    Pt_toplep_vs ->Fill((pttoplep_rec-pttoplep_gen)/pttoplep_gen,weight);
    Pt_tophad_vs ->Fill((pttophad_rec-pttophad_gen)/pttophad_gen,weight);

    eta_Higgs_vs ->Fill((etaHiggs_rec-etaHiggs_gen)/etaHiggs_gen,weight);
    eta_ZPrime_vs ->Fill((etaZPrime_rec-etaZPrime_gen)/etaZPrime_gen,weight);
    eta_TPrime_vs ->Fill((etaTPrime_rec-etaTPrime_gen)/etaTPrime_gen,weight);
    eta_toplep_vs ->Fill((etatoplep_rec-etaHiggs_gen)/etatoplep_gen,weight);
    eta_tophad_vs ->Fill((etatophad_rec-etaHiggs_gen)/etatophad_gen,weight);

    


    ////////////////////////////////////////////////////////

    if(n_vector[16]){
      double lep_pt_gen = zprimegen.Lepton().pt();
      double lep_eta_gen = zprimegen.Lepton().eta();
      double nu_pt_gen = zprimegen.Neutrino().pt();
      double nu_eta_gen = zprimegen.Neutrino().eta();
      double blep_pt_gen = zprimegen.BLep().pt();
      double blep_eta_gen = zprimegen.BLep().eta();
      double toplep_pt_gen = zprimegen.Top().pt();
      double toplep_eta_gen = zprimegen.Top().eta();
      double tophad_pt_gen = zprimegen.ATop().pt();
      double tophad_eta_gen = zprimegen.ATop().eta();

      Pt_tophad_rec_vs_Pt_tophad_gen->Fill(hyp->tophad_v4().Pt(), tophad_pt_gen, weight);
      Pt_toplep_rec_vs_Pt_toplep_gen->Fill(hyp->toplep_v4().Pt(), toplep_pt_gen, weight);
      Pt_lepton_rec_vs_Pt_lepton_gen->Fill(hyp->lepton().pt(), lep_pt_gen, weight);
      Pt_neutrino_rec_vs_Pt_neutrino_gen->Fill(hyp->neutrino_v4().Pt(), nu_pt_gen, weight);
      Pt_blep_rec_vs_Pt_blep_gen->Fill(hyp->blep_v4().Pt(), blep_pt_gen, weight);
      eta_tophad_rec_vs_eta_tophad_gen->Fill(hyp->tophad_v4().eta(), tophad_eta_gen, weight);
      eta_toplep_rec_vs_eta_toplep_gen->Fill(hyp->toplep_v4().eta(), toplep_eta_gen, weight);
      eta_lepton_rec_vs_eta_lepton_gen->Fill(hyp->lepton().eta(), lep_eta_gen, weight);
      eta_neutrino_rec_vs_eta_neutrino_gen->Fill(hyp->neutrino_v4().eta(), nu_eta_gen, weight);
      eta_blep_rec_vs_eta_blep_gen->Fill(hyp->blep_v4().eta(), blep_eta_gen, weight);
    }
    //hier
    if(deltaR(zprimegen.Higgs(), hyp->HZW_v4() )<=0.8){
      massfitHiggs -> Fill(hyp->HZW_v4().M(), weight);
      missmatch->Fill("H matched",1);
    }else missmatch->Fill("H miss",1);
    if(zprimegen.IsSemiLeptonicDecay()){
	if(deltaR(zprimegen.hadTop(), hyp->HZW_v4() )<=0.8){
	  missmatch->Fill("had Top matched",weight);
	}else missmatch->Fill("had Top miss",weight);
      if(deltaR(zprimegen.lepTop(), hyp->HZW_v4() )<=0.8){
	missmatch->Fill("lep Top matched",weight);
      }else missmatch->Fill("lep Top miss",weight);
    }

    if(deltaR(zprimegen.ZBoson(), hyp->HZW_v4() )<=0.8){
      massfitZ -> Fill(hyp->HZW_v4().M(), weight);
      missmatch->Fill("Z matched",1);
    }else{missmatch->Fill("Z miss",1);}
    if(deltaR(zprimegen.WBoson(), hyp->HZW_v4() )<=0.8){
      massfitZ -> Fill(hyp->HZW_v4().M(), weight);
      missmatch->Fill("W matched",1);
    }else{missmatch->Fill("W miss",1);}


    bool b_dquark1 = false;
    bool b_dquark2 = false;
    bool b_dbhad = false;
    //massfit had top  
    //hier
    for(const auto hadtop : hyp->tophad_jets()){
      if(deltaR(zprimegen.Quark1(),hadtop)<=0.8) b_dquark1=true;
      if(deltaR(zprimegen.Quark2(),hadtop)<=0.8) b_dquark2=true;
      if(deltaR(zprimegen.BHad(),hadtop)<=0.8) b_dbhad=true;
    }

    bool all_with_tophad = (deltaR(zprimegen.Quark1(),hyp->tophad_v4())<0.8)&&(deltaR(zprimegen.Quark2(),hyp->tophad_v4())<0.8)&&(deltaR(zprimegen.BHad(),hyp->tophad_v4())<0.8);

    if(zprimegen.IsSemiLeptonicDecay()){
      if(b_dquark1 && b_dquark2 && b_dbhad){
	massfithadTop ->Fill(hyp->tophad_v4().M(), weight);
	missmatch->Fill("had TopR matched",weight);
      }else  missmatch->Fill("had TopR miss",weight);
    }

    //Massfit lep. top

    double distance_blep = deltaR(zprimegen.BLep(),hyp->blep_v4()) ;
    double distance_lepton = deltaR(hyp->lepton(),zprimegen.Lepton());

    if(zprimegen.IsSemiLeptonicDecay()){
      if(distance_blep<=0.8&& distance_lepton<= 0.8){
	massfitlepTop ->Fill(hyp->toplep_v4().M(), weight);
	missmatch->Fill("lep TopR matched",weight);
      }else{missmatch->Fill("lep TopR miss",weight);}

      if(distance_blep<=0.8){
	missmatch->Fill("lep Top b  matched",weight);
      }else{missmatch->Fill("lep Top b miss",weight);}
      
      if(distance_lepton<=0.8){
	missmatch->Fill("lep Top lepton  matched",weight);
      }else{missmatch->Fill("lep Top lepton miss",weight);}


    }
    // go over each jet to define the matching by the closest jet to the generator b quark and not over deltaR (might be different due to hadronisation)
    double deltaR_bgen_breco_min=1000;
    Jet closets_reco_jet;
    for(auto jet:*e.jets){
      double tmp_deltar = deltaR(zprimegen.BLep(), jet);
      if(tmp_deltar < deltaR_bgen_breco_min){
	deltaR_bgen_breco_min = tmp_deltar;
	closets_reco_jet = jet;
      }
    }
    if(deltaR(closets_reco_jet,hyp->blep_v4()))  missmatch->Fill("lep Top matched with closet jet",weight);
    else missmatch->Fill("lep Top missmatched with closest jet",weight);


    double deltar_WQ1= deltaR(hyp->HZW_v4(), zprimegen.WHadQ1());
    double deltar_WQ2= deltaR(hyp->HZW_v4(), zprimegen.WHadQ2());
    double deltar_bhad= deltaR(hyp->HZW_v4(), zprimegen.BHad());

    if(deltar_WQ1<=0.8 && deltar_WQ2<=0.8  && deltar_bhad<=0.8) Discriminator_mismatchedHadTop->Fill(hyp->discriminator(m_discriminator_name) ,weight);

    if(deltaR(hyp->HZW_v4(), zprimegen.Higgs())>0.8) Discriminator_mismatchedHiggs->Fill(hyp->discriminator(m_discriminator_name) ,weight);

    // with bjet matching of the closest jet
    //for each jet in event check if it is the closest, if it is check deltar with hyp->blep_v4() if match count number_closest up and for each event fill number_closet_jet->Fill("#reconstruction",weight);
    if(deltaR(hyp->blep_v4(),zprimegen.BLep())<0.4){
      number_closet_jet_with_matching->Fill("#reconstruction",weight);
      double deltar_min = 1000;
      Jet closest_jet;
      for(jet:*e.jets){
	double deltar = deltaR(hyp->lepton().v4(),jet);
	if(deltar < deltar_min){
	  deltar_min = deltar;
	  closest_jet = jet;
	}
      }
      if(deltaR(closest_jet,hyp->blep_v4())<0.01) number_closet_jet_with_matching->Fill("#closest",weight);
    }



  }

 
  
 
  
  Discriminator->Fill(hyp->discriminator(m_discriminator_name) ,weight);
  Discriminator_2->Fill(hyp->discriminator(m_discriminator_name) ,weight);
  Discriminator_3->Fill(hyp->discriminator(m_discriminator_name) ,weight);
  DiscriminatorH->Fill(hyp->discriminator(m_discriminator_name+"H") ,weight);
  Discriminator_vs_M_ttbar->Fill(hyp->discriminator(m_discriminator_name), mttbar_rec, weight);

  double mtoplep=0;
  double mtophad=0;
  if(hyp->toplep_v4().isTimelike()) mtoplep = hyp->toplep_v4().M();
  if(hyp->tophad_v4().isTimelike()) mtophad = hyp->tophad_v4().M();
  M_toplep_rec->Fill(mtoplep,weight);
  M_tophad_rec->Fill(mtophad,weight);

  if(hyp->tophad_jets().size()==1) M_tophad_rec_1jet->Fill(mtophad,weight);
  if(hyp->tophad_jets().size()==2) M_tophad_rec_2jet->Fill(mtophad,weight);
  if(hyp->tophad_jets().size()>=3) M_tophad_rec_3jet->Fill(mtophad,weight);

  if(e.is_valid(h_ttbargen)){
    const auto & ttbargen = e.get(h_ttbargen);  

    if(hyp->tophad_jets().size()==1){
      double deltar_gluon = deltaR(ttbargen.Gluon(),hyp->tophad_jets()[0]);
      double deltar_quark = deltaR(ttbargen.Quark(),hyp->tophad_jets()[0]);
      if(deltar_gluon<=0.8) M_tophad_rec_1jet_pdgid->Fill("Gluon",weight);
      else if(deltar_quark<=0.8) M_tophad_rec_1jet_pdgid->Fill("Quark",weight);
    }

    if(hyp->tophad_jets().size()==2){
      double deltar_gluon = deltaR(ttbargen.Gluon(),hyp->tophad_jets()[1]);
      double deltar_quark = deltaR(ttbargen.Quark(),hyp->tophad_jets()[1]);
      if(deltar_gluon<=0.8) M_tophad_rec_2jet_pdgid->Fill("Gluon",weight);
      else if(deltar_quark<=0.8) M_tophad_rec_2jet_pdgid->Fill("Quark",weight);
    }

    if(hyp->tophad_jets().size()==3){
      double deltar_gluon = deltaR(ttbargen.Gluon(),hyp->tophad_jets()[2]);
      double deltar_quark = deltaR(ttbargen.Quark(),hyp->tophad_jets()[2]);
      if(deltar_gluon<=0.8) M_tophad_rec_3jet_pdgid->Fill("Gluon",weight);
      else if(deltar_quark<=0.8) M_tophad_rec_3jet_pdgid->Fill("Quark",weight);
    }
  }



  Pt_toplep_rec->Fill( hyp->toplep_v4().Pt(),weight );
  Pt_tophad_rec->Fill( hyp->tophad_v4().Pt(),weight );

//////////// distance between lepton and bjet of the lep top
    double deltaR_lepton_bjet = deltaR(hyp->lepton().v4(),hyp->blep_v4());
    deltar_mu_blep->Fill(deltaR_lepton_bjet,weight);
    
    //for each jet in event check if it is the closest, if it is check deltar with hyp->blep_v4() if match count number_closest up and for each event fill number_closet_jet->Fill("#reconstruction",weight);
    number_closet_jet->Fill("#reconstruction",weight);
    double deltar_min = 1000;
    Jet closest_jet;
    for(jet:*e.jets){
      double deltar = deltaR(hyp->lepton().v4(),jet);
      if(deltar < deltar_min){
	deltar_min = deltar;
	closest_jet = jet;
      }
    }
    if(deltaR(closest_jet,hyp->blep_v4())<0.01) number_closet_jet->Fill("#closest",weight);





}
