#include "UHH2/ZprimeToTprimeTtZtH/include/SR_PDFHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetIds.h"
#include <math.h>
#include <sstream>

#include "TH1F.h"
#include "TH2D.h"
#include <iostream>

using namespace std;
using namespace uhh2;
using namespace uhh2examples;


SR_PDFHists::SR_PDFHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name, bool use_pdf_weights_): Hists(ctx, dirname), use_pdf_weights(use_pdf_weights_),dirnames(dirname){   

  TString name = discriminator_name;

  is_mc = ctx.get("dataset_type") == "MC";
 
  m_discriminator_name ="Chi2"; 
  TString m_oname = ctx.get("dataset_version");

  TString m_pdfname = "NNPDF30_lo_as_0130";
  if(m_oname.Contains("TTbar")) m_pdfname = "PDF4LHC15_nlo_mc";
  if(m_oname.Contains("WJets")) m_pdfname = "PDF4LHC15_nlo_mc";
  if(m_oname.Contains("ST")) m_pdfname = "PDF4LHC15_nlo_mc";
  TString weightpath = ctx.get("PDFWeightPath");  cout << "File: " << weightpath+m_oname << endl;
  if(m_oname.Contains("MC_ZPrime")){
  m_pdfweights.reset(new PDFWeights(m_pdfname,weightpath+m_oname)); 
  }else{
    m_pdfweights.reset(new PDFWeights(m_pdfname));
  } 


  if(dirnames.find("chi2cut")!= std::string::npos|| dirnames.find("btag")!= std::string::npos){
    std::vector<TString> hists = {"M_ZPrime_rec", "Pt_toplep_rec"};

    for(TString obs:hists){

    for(int i=0; i<100; i++){
      stringstream ss_name;
      ss_name << obs<<"_PDF_"  << i+1 ;
    
      stringstream ss_title;
      ss_title << obs<<" [GeV/c^{2}] for PDF No. "  << i+1 << " out of 100" ;
   

      string s_name = ss_name.str();
      string s_title = ss_title.str();
   
      const char* char_name = s_name.c_str();
      const char* char_title = s_title.c_str();
   
      if(obs=="M_ZPrime_rec")histo_names[i] = s_name;
      else if(obs=="Pt_toplep_rec")histo_names_toppt[i] = s_name;

      if(obs=="M_ZPrime_rec")book<TH1F>(char_name, char_title, 50,600,5000);
      else if(obs=="Pt_toplep_rec")book<TH1F>(char_name, char_title,60, 0, 1200);
    }//for loop
    }//observables

    h_hyps = ctx.get_handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>>(hyps_name);
  }else if(dirnames.find("twodcut")!= std::string::npos){
    std::vector<TString> hists = {"pt", "eta","mass"};

    for(TString obs:hists){

      for(int i=0; i<100; i++){
	stringstream ss_name;
	ss_name <<obs<< "_PDF_"  << i+1 ;
	
	stringstream ss_title;
	ss_title << obs <<" [GeV] for PDF No. "  << i+1 << " out of 100" ;


	string s_name = ss_name.str();
	string s_title = ss_title.str();

	const char* char_name = s_name.c_str();
	const char* char_title = s_title.c_str();

	if(obs=="pt")histo_names_pt[i] = s_name;
	if(obs=="eta")histo_names_eta[i] = s_name;
	if(obs=="mass")histo_names_mass[i]=s_name;

	if(obs=="pt")book<TH1F>(char_name, char_title, 200,0,1500);
	else if(obs=="eta")book<TH1F>(char_name,char_title,100,-3,3);
	else if(obs=="mass")book<TH1F>(char_name,char_title,100,0,300);
	else book<TH1F>(char_name, char_title, 200,0,1500);
      }//for loop
    }//observables
  }else if(dirnames.find("event")!= std::string::npos){
    std::vector<TString> hists = {"MET"};

    for(TString obs:hists){
      for(int i=0; i<100; i++){
	stringstream ss_name;
	ss_name <<obs<< "_PDF_"  << i+1 ;
	
	stringstream ss_title;
	ss_title << obs <<" [GeV] for PDF No. "  << i+1 << " out of 100" ;

	string s_name = ss_name.str();
	string s_title = ss_title.str();

	const char* char_name = s_name.c_str();
	const char* char_title = s_title.c_str();

	histo_names_MET[i] = s_name;

	if(obs=="MET")book<TH1F>(char_name,char_title,200,0,1000);
	else book<TH1F>(char_name, char_title, 200,0,1500);
      }//for loop
    }//observables
  }//dirname event
}

void SR_PDFHists::fill(const Event & event){
  double weight = event.weight;

  if(dirnames.find("chi2cut")!= std::string::npos|| dirnames.find("btag")!= std::string::npos){
    std::vector<ZPrimeTotTPrimeReconstructionHypothesis> hyps = event.get(h_hyps);
  
    const ZPrimeTotTPrimeReconstructionHypothesis* hyp = get_best_hypothesis( hyps, m_discriminator_name );

    double mZPrime_rec = 0;
    if(is_mc){
    
      if( (hyp->HZW_v4() + hyp->top_v4()+ hyp->antitop_v4()).isTimelike() )
	mZPrime_rec = (hyp->HZW_v4()+ hyp->top_v4()+ hyp->antitop_v4()).M();
      else{
	mZPrime_rec = sqrt( -(hyp->HZW_v4()+ hyp->top_v4()+ hyp->antitop_v4()).mass2());
      }

      
      double pt_toplep_rec=(hyp->toplep_v4()).Pt();
      

      std::vector<double> weights = m_pdfweights->GetWeightList(event);
      for(int i=0; i<100; i++){
	if(use_pdf_weights){
	  double fillweight = weight*weights[i]; // LQ PDF
	  const char* name = histo_names[i].c_str();
	  const char* name_toppt = histo_names_toppt[i].c_str();
	  hist(name)->Fill(mZPrime_rec, fillweight);
	  hist(name_toppt)->Fill(pt_toplep_rec, fillweight);
	}//if use pdf weights
      }//over all 100 possibilities
    } //is_mc
  }else if(dirnames.find("twodcut")!= std::string::npos){
    double pt=0;
    double eta=0;
    if(dirnames.find("muon")!= std::string::npos) pt= event.muons->at(0).pt();
    if(dirnames.find("muon")!= std::string::npos) eta=event.muons->at(0).eta();

    if(dirnames.find("elec")!= std::string::npos) pt= event.electrons->at(0).pt();
    if(dirnames.find("elec")!= std::string::npos) eta=event.electrons->at(0).eta();

    if(dirnames.find("topjet")!= std::string::npos){

      pt= event.topjets->at(0).pt();
      eta=event.topjets->at(0).eta();


    }else if(dirnames.find("jet")!= std::string::npos){
      //Fuer alle jets nicht nur fuer den ersten
      //      pt= event.jets->at(0).pt();
      //      eta=event.jets->at(0).eta();
    }


    std::vector<double> weights = m_pdfweights->GetWeightList(event);
    for(int i=0; i<100; i++){
      if(use_pdf_weights){
	double fillweight = weight*weights[i]; 
	const char* pt_name = histo_names_pt[i].c_str();
	const char* eta_name = histo_names_eta[i].c_str();


	if(dirnames.find("topjet")!=std::string::npos){
	  const char* mass_name = histo_names_mass[i].c_str();
	  //mass_subjet_sum
	  for(auto topjet:*event.topjets){
	    LorentzVector subjet_sum;
	    for (const auto s : topjet.subjets()) {
	      subjet_sum += s.v4();
	    }//over all subjets
	    hist(mass_name)->Fill(subjet_sum.M(),fillweight);
	  }//over all topjets
	  //pt
	  hist(pt_name)->Fill(pt, fillweight);	
	  hist(eta_name)->Fill(eta, fillweight);
	}else if((dirnames.find("jet")!= std::string::npos) && (dirnames.find("topjet")== std::string::npos)){
	  for(auto ak4:*event.jets){
	    hist(pt_name)->Fill(ak4.pt(), fillweight);	
	    hist(eta_name)->Fill(ak4.eta(), fillweight);
	  }
	}else{
	  hist(pt_name)->Fill(pt, fillweight);	
	  hist(eta_name)->Fill(eta, fillweight);
	}//ak4 jets or electron, muons, ak8 jets
      }//if use pdf weights                                                                                                                                                        
    }//over all 100 possibilities    
  }else if(dirnames.find("event")!= std::string::npos){
    double met = event.met->pt();

    std::vector<double> weights = m_pdfweights->GetWeightList(event);
    for(int i=0; i<100; i++){
      if(use_pdf_weights){
	double fillweight = weight*weights[i]; 

	const char* met_name = histo_names_MET[i].c_str();

	hist(met_name)->Fill(met, fillweight);

      }//if use pdf weights                                                                                                                                                        
    }//over all 100 possibilities    
  }// dirname event

}

SR_PDFHists::~SR_PDFHists(){}














