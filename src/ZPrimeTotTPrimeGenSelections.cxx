#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"

#include <stdexcept>
#include <cmath>

// extern variables
std::vector<int> n_vector= {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

using namespace uhh2;
using namespace std;

ZPrimeGen::ZPrimeGen(const vector<GenParticle> & genparticles){
  
  //allparticles.reseve(genparticles->size());

  //find the first gluon and top and not one later on
   bool oncetop = true;  
   bool oncegluon = true;
 
  //Go over all Particles of the event to find the first Z', T', Top and Gluon(which should belong to a ZPrime if the ZPrime is not listet)
  for(unsigned int i=0; i<genparticles.size(); i++){
    const GenParticle & genp = genparticles[i];
    // allparticles.push_back(*genp);
    if(abs(genp.pdgId()) ==9900113){ZPrime = genp;}
    if(abs(genp.pdgId())==8000001){TPrime = genp;}
    if(abs(genp.pdgId())==6 && oncetop){main_top = genp;oncetop=false;}
    if(abs(genp.pdgId())==21 && i>1 && oncegluon && genp.mother(&genparticles,1)-> index() == 0){ZPrimeGluon = genp;oncegluon= false;}
  }

  //Check if there is a ZPrime in the event
  if(abs(ZPrime.pdgId()) ==9900113){

    //Check if one of the daughters of the Z' is the T'
    if(abs(ZPrime.daughter(&genparticles,1)->pdgId())==8000001 || abs(ZPrime.daughter(&genparticles,2)->pdgId())==8000001){

      //Assign all Particles
      ZPrimeGen::assign(*ZPrime.daughter(&genparticles,1), *ZPrime.daughter(&genparticles,2), genparticles);

      //if there is no T' in event reconstrutct it out of the decay products and then assign all particles
    }else{
    
      ZPrimeGen::recoTPrime(genparticles);

      if(throw_failure){
	cout << "**********************  T' found **********************" << endl;
	cout << " T' (daughter W,Z,H || b,t,t )"<< TPrime.daughter(&genparticles,1)->pdgId() << endl;
	cout << " T' (v4())"<< TPrime.v4() << endl;
	cout << " T' (Mass)"<< TPrime.v4().mass() << endl;
      }

      ZPrimeGen::assign(*TPrime.daughter(&genparticles,1), *TPrime.daughter(&genparticles,2), genparticles);
    }
    //if no Z' was found check if a T' was found
  }else if(abs(TPrime.pdgId()) ==8000001 ){
 
    //assign all particles and then recostruct the Z'
    ZPrimeGen::assign(TPrime,main_top , genparticles);
    ZPrimeGen::recoZPrime(TPrime , main_top);

    if(throw_failure){
      cout << "**********************  Z' found ********************** " << ZPrime.v4() <<endl;
      cout << "Z' (v4()) " << ZPrime.v4() <<endl;
      cout << "Z' (Mass) " << ZPrime.v4().mass() <<endl;
    }
 

    //event in which is no T' and no Z'. So first recostruct the TPrime and then assign all particles, so one can reconstruct the ZPrime
  }else{
 
    ZPrimeGen::recoTPrime(genparticles);

    if(throw_failure){
      cout << "**********************  T' found with no Z'  **********************" << endl;
      cout << " T' (daughter W,Z,H || b,t,t )"<< TPrime.daughter(&genparticles,1)->pdgId() << endl;
      cout << " T' (v4())"<< TPrime.v4() << endl;
      cout << " T' (Mass)"<< TPrime.v4().mass() << endl;
    }

    ZPrimeGen::assign(*TPrime.daughter(&genparticles,1), *TPrime.daughter(&genparticles,2), genparticles);
    ZPrimeGen::recoZPrime(TPrime , main_top);

    if(throw_failure){
      cout << "**********************  Z' found with no T' ********************** " << ZPrime.v4() <<endl;
      cout << "Z' (v4()) " << ZPrime.v4() <<endl;
      cout << "Z' (Mass) " << ZPrime.v4().mass() <<endl;
    }
  }

  ZPrimeGen::DecayChannel(genparticles);
}

      
// assign all decay particles to T' and T
void ZPrimeGen::assign(GenParticle  daughter1, GenParticle  daughter2,const vector<GenParticle> & genparticles){
  //Initiate all countervariables
  int  n_TPrime=0, n_t =0;

  if (abs(daughter1.pdgId()) ==8000001 ){
    TPrime = daughter1;
    n_TPrime++;
    //Check the Decay chanel of the TPrime. If it is none of these there is an error message
    switch(abs(TPrime.daughter(&genparticles,1)->pdgId())){
    case 5: bottom = *TPrime.daughter(&genparticles,1);break;
    case 6: if(abs(TPrime.daughter(&genparticles,2)->pdgId())== 23){ topZ =  *TPrime.daughter(&genparticles,1);break;}else {topHiggs =  *TPrime.daughter(&genparticles,1);break;};
    case 23:Z   =  *TPrime.daughter(&genparticles,1);decaydaughterZ1 = *Z.daughter(&genparticles,1) ; decaydaughterZ2 = *Z.daughter(&genparticles,2);break; 
    case 24:W   =  *TPrime.daughter(&genparticles,1);decaydaughterW1 = *W.daughter(&genparticles,1) ; decaydaughterW2 = *W.daughter(&genparticles,2);break;
    case 25:H   =  *TPrime.daughter(&genparticles,1);if(H.daughter(&genparticles,1)&&H.daughter(&genparticles,2) ){decaydaughterH1 = *H.daughter(&genparticles,1) ; decaydaughterH2 = *H.daughter(&genparticles,2);};break;
    default :;// cout << "L106: Zerfall der nicht bedacht wurde" <<endl;
    }
    switch(abs(TPrime.daughter(&genparticles,2)->pdgId())){
    case 5: bottom = *TPrime.daughter(&genparticles,2);break;
    case 6: if(abs(TPrime.daughter(&genparticles,1)->pdgId())== 23){ topZ =  *TPrime.daughter(&genparticles,1);break;}else {topHiggs =  *TPrime.daughter(&genparticles,1);break;};
    case 23:Z   =  *TPrime.daughter(&genparticles,2);decaydaughterZ1 = *Z.daughter(&genparticles,1) ; decaydaughterZ2 = *Z.daughter(&genparticles,2);break; 
    case 24:W   =  *TPrime.daughter(&genparticles,2);decaydaughterW1 = *W.daughter(&genparticles,1) ; decaydaughterW2 = *W.daughter(&genparticles,2);break;
    case 25:H   =  *TPrime.daughter(&genparticles,2);if(H.daughter(&genparticles,1)&&H.daughter(&genparticles,2)){decaydaughterH1 = *H.daughter(&genparticles,1) ; decaydaughterH2 = *H.daughter(&genparticles,2);};break;
    default :;// cout << "L114: Zerfall der nicht bedacht wurde" <<endl;
    }
  }

  if(abs(daughter2.pdgId()) == 8000001){ 
    TPrime = daughter2;
    n_TPrime++;
    switch(abs(TPrime.daughter(&genparticles,1)->pdgId())){
    case 5: bottom = *TPrime.daughter(&genparticles,1);break;
    case 6: if(abs(TPrime.daughter(&genparticles,2)->pdgId())== 23){ topZ =  *TPrime.daughter(&genparticles,1);break;}else {topHiggs =  *TPrime.daughter(&genparticles,1);break;};
    case 23:Z   =  *TPrime.daughter(&genparticles,1);decaydaughterZ1 = *Z.daughter(&genparticles,1) ; decaydaughterZ2 = *Z.daughter(&genparticles,2);break; 
    case 24:W   =  *TPrime.daughter(&genparticles,1);decaydaughterW1 = *W.daughter(&genparticles,1) ; decaydaughterW2 = *W.daughter(&genparticles,2);break;
    case 25:H   =  *TPrime.daughter(&genparticles,1);if(H.daughter(&genparticles,1)&&H.daughter(&genparticles,2)){decaydaughterH1 = *H.daughter(&genparticles,1) ; decaydaughterH2 = *H.daughter(&genparticles,2);};break;
    default : ;//cout << "L127: Zerfall der nicht bedacht wurde" <<endl;
    }
    switch(abs(TPrime.daughter(&genparticles,2)->pdgId())){
    case 5: bottom = *TPrime.daughter(&genparticles,2);break;
    case 6: if(abs(TPrime.daughter(&genparticles,1)->pdgId())== 23){ topZ =  *TPrime.daughter(&genparticles,1);break;}else {topHiggs =  *TPrime.daughter(&genparticles,1);break;};
    case 23:Z   =  *TPrime.daughter(&genparticles,2);decaydaughterZ1 = *Z.daughter(&genparticles,1) ; decaydaughterZ2 = *Z.daughter(&genparticles,2);break; 
    case 24:W   =  *TPrime.daughter(&genparticles,2);decaydaughterW1 = *W.daughter(&genparticles,1) ; decaydaughterW2 = *W.daughter(&genparticles,2);break;
    case 25:H   =  *TPrime.daughter(&genparticles,2);if(H.daughter(&genparticles,1)&&H.daughter(&genparticles,2)){decaydaughterH1 = *H.daughter(&genparticles,1) ; decaydaughterH2 = *H.daughter(&genparticles,2);};break;
    default :;// cout << "L135: Zerfall der nicht bedacht wurde" <<endl;
    }}

  //Check if the decay product is a top and collect the decay products of it
  if(abs(daughter1.pdgId())==6){
    main_top = daughter1;
    n_t++;
    //Check the decay chanel of the top.If it is none of these there is an error message 
    if(abs(main_top.daughter(&genparticles,1)->pdgId())==24){
      top_W = *main_top.daughter(&genparticles,1);
      top_b = *main_top.daughter(&genparticles,2);
    }

    if(abs(main_top.daughter(&genparticles,2)->pdgId())==24){
      top_W = *main_top.daughter(&genparticles,2);
      top_b = *main_top.daughter(&genparticles,1);
    }
  }

  if(abs(daughter2.pdgId())==6){
    main_top = daughter2;
    n_t++;
    if(abs(main_top.daughter(&genparticles,1)->pdgId())==24){
      top_W = *main_top.daughter(&genparticles,1);
      top_b = *main_top.daughter(&genparticles,2);
    }

    if(abs(main_top.daughter(&genparticles,2)->pdgId())==24){
      top_W = *main_top.daughter(&genparticles,2);
      top_b = *main_top.daughter(&genparticles,1);
    }
  }
}


//reconstruct the T' by finding all decay products and then calculate the invariant mass
void ZPrimeGen::recoTPrime( const vector<GenParticle> & genparticles){
  if(throw_failure){
    cout << "----------------in TPrimeReco--------------" <<endl;
  }
  //To Reconstruct the T' we need first and second top, gluon if there, Bottom if there and a Higgs, W or Z
  bool firstTPD3 = true;
  bool firstTop = true;
  bool secondTop = true;
  bool firstGluon = true;
  bool firstBottom = true;
 
  for(unsigned int k=0; k<genparticles.size(); k++){
    const GenParticle & genpa = genparticles[k];
    if((abs(genpa.pdgId())==23 || abs(genpa.pdgId())==24 || abs(genpa.pdgId())==25) && firstTPD3 ){
      TPrimedaughter3 = genpa;
      firstTPD3=false;
    }
    if(abs(genpa.pdgId())== 6 && secondTop && !firstTop){
      Top2 = genpa;
      secondTop = false;
    }
    if(abs(genpa.pdgId())== 6 && firstTop){
      Top1 = genpa;
      firstTop = false;
    }
    if(abs(genpa.pdgId())== 21 && k>1 &&( abs(genpa.mother(&genparticles,1)->pdgId())==9900113 || abs(genpa.mother(&genparticles,1)->index()) == 0 ) && firstGluon){
      Gluon = genpa;
      firstGluon = false;
    }
    if(abs(genpa.pdgId())== 5 && (abs(genpa.mother(&genparticles,1)->pdgId())==9900113 || abs(genpa.mother(&genparticles,1)->index()) == 0 ) && firstBottom){
      Bottom = genpa;
      firstBottom = false;
    }
  }

  //if H or Z Decay we need to fin the right top by the distance between the particles
  if(abs(TPrimedaughter3.pdgId())== 23 || abs(TPrimedaughter3.pdgId())== 25 ){
    float deltaR1 = deltaR(TPrimedaughter3,Top1);
    float deltaR2 = deltaR(TPrimedaughter3,Top2) ;

    if(deltaR1 < deltaR2){
      
      if(throw_failure){
	cout << "TPrimereco v4 " << TPrimereco1.v4() << endl;
      }
      
      TPrimereco1.set_v4(TPrimedaughter3.v4()+ Top1.v4());
      TPrime.set_v4(TPrimedaughter3.v4()+ Top1.v4()); 
      TPrime.set_daughter1(TPrimedaughter3.index()); 
      TPrime.set_daughter2(Top1.index());
    }else{
      TPrimereco2.set_v4(TPrimedaughter3.v4()+ Top2.v4());
      TPrime.set_v4(TPrimedaughter3.v4()+ Top2.v4());
      TPrime.set_daughter1(TPrimedaughter3.index()); 
      TPrime.set_daughter2(Top2.index());
    }
    // Else we need the W, B and sometimes an G
  }else{
    float deltaRWG = deltaR(TPrimedaughter3, Gluon);
    float deltaRBG = deltaR(Bottom, Gluon);
   
    if(throw_failure){
      cout << "Delta R WG "<<deltaRWG << endl;
      cout << "pdgId Bottom und Gluon " << Bottom.pdgId() << " " << Gluon.pdgId()<<endl; 
      cout << "Delta R BG "<<deltaRBG << endl;
    }
 
    if(deltaRWG > 2){
      TPrimereco1.set_v4(TPrimedaughter3.v4()+ Bottom.v4());
      TPrime.set_v4(TPrimedaughter3.v4() + Bottom.v4());
      TPrime.set_daughter1(TPrimedaughter3.index()); 
      TPrime.set_daughter2(Bottom.index());
    }else{

      TPrimereco3.set_v4(TPrimedaughter3.v4()+ Bottom.v4()+ Gluon.v4());
      TPrime.set_v4(TPrimedaughter3.v4() + Bottom.v4() + Gluon.v4());
      TPrime.set_daughter1(TPrimedaughter3.index()); 
      TPrime.set_daughter2(Bottom.index());
    } 
  }
}

//reconstruct the Z' by calculating the invariant mass of the T' Top
void ZPrimeGen::recoZPrime(GenParticle TPrime ,GenParticle  main_top){
 
  if(throw_failure){
    cout << "----------------------------- reco ZPrime ----------------------"<<endl;
  }
  
  ZPrime.set_v4( TPrime.v4() + main_top.v4());
  ZPrimereco.set_v4(TPrime.v4()+ main_top.v4());
}



void ZPrimeGen::DecayChannel( const vector<GenParticle> & genparticles ){
  
  if(topHiggs.daughter(&genparticles,1) && topHiggs.daughter(&genparticles,2) ){
    if(abs(topHiggs.daughter(&genparticles,1)-> pdgId()) == 24){
      WHiggsTop1 = *topHiggs.daughter(&genparticles,1);
      BHiggsTop1 = *topHiggs.daughter(&genparticles,2);
    }else{
      WHiggsTop1 = *topHiggs.daughter(&genparticles,2);
      BHiggsTop1 = *topHiggs.daughter(&genparticles,1);
    }
  }



  if(topZ.daughter(&genparticles,1) && topZ.daughter(&genparticles,2) ){
    if(abs(topZ.daughter(&genparticles,1)-> pdgId()) == 24){
      WZTop1 = *topZ.daughter(&genparticles,1);
      BZTop1 = *topZ.daughter(&genparticles,2);
    }else{
      WZTop1 = *topZ.daughter(&genparticles,2);
      BZTop1 = *topZ.daughter(&genparticles,1);

    }
  }

  if(throw_failure){
    cout << "----------------------------- in DecayChannel --------------------------------------------" << endl;
  }
  //  cout << "L:290 (GenSelections) "<< endl;
  //top decays leptonic
  if(top_W.daughter(&genparticles,1) && top_W.daughter(&genparticles,2)){
    if(throw_failure){
      cout << "Top W daughter pdgId (<=6 || <= 18 ) " << abs(top_W.daughter(&genparticles,1)->pdgId())<< endl;
    }
    if ((abs(top_W.daughter(&genparticles,1)->pdgId())<=18 && abs(top_W.daughter(&genparticles,1)->pdgId()) >6)  || (abs(top_W.daughter(&genparticles,2)->pdgId())<=18&& abs(top_W.daughter(&genparticles,2)->pdgId()) >6)){
      leptop = main_top;

      //muon zum Vergleich des Verzeigungsverhaeltnisses
      if(abs(top_W.daughter(&genparticles,1)->pdgId())==13)muon_top = *top_W.daughter(&genparticles,1);
	else muon_top = *top_W.daughter(&genparticles,2);
      blep = top_b;

      if(throw_failure){
	cout << "leptonischer Zerfall von Top W" << endl;
      }
      // cout << "L:301 (GenSelections) "<< endl;  
      //top (TPrime)  decays  leptonic
      if( H.daughter(&genparticles,1)&&WHiggsTop1.daughter(&genparticles,1) && WHiggsTop1.daughter(&genparticles,2)){
	if((abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())<=18 && abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())>6) || (abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())<=18 && abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())>6)){

	  //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *WHiggsTop1.daughter(&genparticles,1);
	  else muon_tprime =*WHiggsTop1.daughter(&genparticles,2); 

	  if(throw_failure){
	    cout << "top leptonisch, Top Higgs leptonisch, higgs daughter pdgId()" << abs(H.daughter(&genparticles,1)->pdgId()) << endl;
	  } 
	  n_H++;
	  switch(abs(H.daughter(&genparticles,1)->pdgId())){
	  case 21: n_TllTllgg++;n_TllTll++; break;
	  case 23: n_TllTllZZ++;n_TllTll++; break;
	  case 24: n_TllTllWW++;n_TllTll++; break;
	  case 5: n_TllTllbb++;n_TllTll++; break;
	  default: ;//cout<<"L461: Zerfall nicht bedacht, pdgID "<< abs(H.daughter(&genparticles,1)->pdgId()) <<endl;break;
	  }
	}
      }
      // cout << "L:318 (GenSelections) "<< endl;
      //top (TPrime)  decays  leptonic
      if(Z.daughter(&genparticles,1) &&WZTop1.daughter(&genparticles,1) && WZTop1.daughter(&genparticles,2) ){
	// cout << "L:321 (GenSelections) "<< endl;
	if((abs(WZTop1.daughter(&genparticles,1)->pdgId())<=18 && abs(WZTop1.daughter(&genparticles,1)->pdgId())>6) || (abs(WZTop1.daughter(&genparticles,2)->pdgId())<=18&& abs(WZTop1.daughter(&genparticles,2)->pdgId())>6)){
	  n_Z++;
	  //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(WZTop1.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *WZTop1.daughter(&genparticles,1);
	  else muon_tprime =*WZTop1.daughter(&genparticles,2); 

	  if( abs(Z.daughter(&genparticles,1)->pdgId())<=6){
	    // cout << "L:325 (GenSelections) "<< endl;
 	    n_TllTllZh++;
	    n_TllTll++;  
	  }else if(abs(Z.daughter(&genparticles,1)->pdgId())<=18){
	    // cout << "L:329 (GenSelections) "<< endl;
 	    n_TllTllZll++;
	    n_TllTll++;

	  }else{
	    // cout<<"L477: Z Zerfall nicht bedacht"<<endl;
	  }
	}
      }
      //top (TPrime)  decays  hadronic
      if( H.daughter(&genparticles,1) && WHiggsTop1.daughter(&genparticles,1) && WHiggsTop1.daughter(&genparticles,2)  ){
	if(abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())<=6 || abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())<=6){
	 
	  quark1 = *WHiggsTop1.daughter(&genparticles,1);
	  quark2 = *WHiggsTop1.daughter(&genparticles,2);
	  bhad = BHiggsTop1;

	  if(throw_failure){
	    cout << "Top Higgs hadronisch" << endl;
	  } 
	  n_H++;
	  hadtop = topHiggs;
	  switch(abs(H.daughter(&genparticles,1)->pdgId())){
	  case 21: n_TllThgg++;	n_TllTh++; break;
	  case 23: n_TllThZZ++;	n_TllTh++; break;
	  case 24: n_TllThWW++;	n_TllTh++;break;
	  case 5: n_TllThbb++;	n_TllTh++;break;
	  default:;// cout<<"L493: Zerfall nicht bedacht, pdgID "<< abs(H.daughter(&genparticles,1)->pdgId())<<endl;break;
	  }
	}
      }
      //cout << "L:350 (GenSelections) "<< endl;  
      //top (TPrime)  decays  hadronic
      if(Z.daughter(&genparticles,1)&&WZTop1.daughter(&genparticles,1) && WZTop1.daughter(&genparticles,2) ){
	if(abs(WZTop1.daughter(&genparticles,1)->pdgId())<=6 || abs(WZTop1.daughter(&genparticles,2)->pdgId())<=6){

	  quark1 = *WZTop1.daughter(&genparticles,1);
	  quark2 = *WZTop1.daughter(&genparticles,2);
	  bhad = BZTop1;

	  n_Z++;
	  hadtop = topZ;
	  if(abs(Z.daughter(&genparticles,1)->pdgId())<=6){
	    n_TllThZh++;
	    n_TllTh++;
	  }else if(abs(Z.daughter(&genparticles,1)->pdgId())<=18){
	    n_TllThZll++;
	    n_TllTh++;
	  }else{
	    // cout<<"L508: Z Zerfall nicht bedacht"<<endl;
	  }
	}
      }
      //top decays leptonic
      if(W.daughter(&genparticles,1)){
	if(throw_failure){
	  cout << "(Top leptonisch) WTochter pdgId()" << abs(W.daughter(&genparticles,1)->pdgId()) << endl;
	}
	  n_W++;
	if(abs(W.daughter(&genparticles,1)->pdgId())<=6){
	  n_TllWh++;
	  
	  n_TllTh++;n_TllTll++;
	}else if(abs(W.daughter(&genparticles,1)->pdgId())<=18){

 //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(W.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *W.daughter(&genparticles,1);
	  else muon_tprime =*W.daughter(&genparticles,2); 


	  if(abs(W.daughter(&genparticles,1)->pdgId()) == 12 || abs(W.daughter(&genparticles,1)->pdgId()) == 14  ){
	    neutrino = *W.daughter(&genparticles,1);
	    lepton = *W.daughter(&genparticles,2);
	  } else{
	    neutrino = *W.daughter(&genparticles,2);
	    lepton = *W.daughter(&genparticles,1);
	  }
 	  n_TllWll++;
	  n_TllTh++;n_TllTll++;
	  
	}else{
	  // cout<<"L534: W Zerfall nicht bedacht"<<endl;
	}
      }

      //top decays hadronic
    }else if(abs(top_W.daughter(&genparticles,1)->pdgId())<=6 || abs(top_W.daughter(&genparticles,2)->pdgId())<=6){
      if(throw_failure){
	cout << "Top W hadronisch" << endl;
      }
      hadtop = main_top;
      quark1 = *top_W.daughter(&genparticles,1);
      quark2 = *top_W.daughter(&genparticles,2);
      bhad = top_b;

      //top (TPrime)  decays  leptonic
      if( H.daughter(&genparticles,1) &&WHiggsTop1.daughter(&genparticles,1) && WHiggsTop1.daughter(&genparticles,2) ){
	if((abs((WHiggsTop1.daughter(&genparticles,1))->pdgId())<=18 &&abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())>6)  || (abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())<=18&&abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())>6)){

 //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *WHiggsTop1.daughter(&genparticles,1);
	  else muon_tprime =*WHiggsTop1.daughter(&genparticles,2); 

	  blep = BHiggsTop1;

	  n_H++;
	  leptop = topHiggs;
 	  switch(abs(H.daughter(&genparticles,1)->pdgId())){
	  case 21: n_ThTllgg++;  n_ThTll++; break;
	  case 23: n_ThTllZZ++;  n_ThTll++; break;
	  case 24: n_ThTllWW++;  n_ThTll++;break;
	  case 5: n_ThTllbb++;  n_ThTll++;break;
	  default:;// cout<<"L554: Zerfall nicht bedacht, pdgID "<< abs(H.daughter(&genparticles,1)->pdgId())<<endl;break;
	  }
	}
      }
      //  cout << "L:410 (GenSelections) "<< endl;
      //top (TPrime)  decays  leptonic
      if(Z.daughter(&genparticles,1)&&WZTop1.daughter(&genparticles,1) && WZTop1.daughter(&genparticles,2) ){
  
	if((abs(WZTop1.daughter(&genparticles,1)->pdgId())<=18 && abs(WZTop1.daughter(&genparticles,1)->pdgId())>6)|| (abs(WZTop1.daughter(&genparticles,2)->pdgId())<=18&& abs(WZTop1.daughter(&genparticles,2)->pdgId())>6)){
	  n_Z++;
	  leptop = topZ;
	  blep = BZTop1;

	  //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(WZTop1.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *WZTop1.daughter(&genparticles,1);
	  else muon_tprime =*WZTop1.daughter(&genparticles,2); 


	  if( abs(Z.daughter(&genparticles,1)->pdgId())<=6){
  
	    n_ThTllZh++;
	    n_ThTll++;
	  }else if(abs(Z.daughter(&genparticles,1)->pdgId())<=18){
  
	    n_ThTllZll++;
	    n_ThTll++;

	  }else{
	    // cout<<"L570: Z Zerfall nicht bedacht"<<endl;
	  }
	}
      }
       //top (TPrime)  decays  hadronic
      if(H.daughter(&genparticles,1)&&WHiggsTop1.daughter(&genparticles,1) && WHiggsTop1.daughter(&genparticles,2)  ){
	if(abs(WHiggsTop1.daughter(&genparticles,1)->pdgId())<=6 || abs(WHiggsTop1.daughter(&genparticles,2)->pdgId())<=6){
	  n_H++;
	  switch(abs(H.daughter(&genparticles,1)->pdgId())){
	  case 21: n_ThThgg++;	  n_ThTh++; break;
	  case 23: n_ThThZZ++;	  n_ThTh++; break;
	  case 24: n_ThThWW++;	  n_ThTh++;break;
	  case 5: n_ThThbb++;	  n_ThTh++;break;
	  default: ;//cout<<"L583: Zerfall nicht bedacht, pdgID "<< abs(H.daughter(&genparticles,1)->pdgId())<<endl;break;
	  }
	}
      }
      //top (TPrime)  decays  hadronic
      if(Z.daughter(&genparticles,1)&&WZTop1.daughter(&genparticles,1)&&WZTop1.daughter(&genparticles,2) ){
	if(abs(WZTop1.daughter(&genparticles,1)->pdgId()<=6) || abs(WZTop1.daughter(&genparticles,2)->pdgId())<=6){
	  n_Z++;
	  if( abs(Z.daughter(&genparticles,1)->pdgId())<=6){
	    n_ThThZh++;
	    n_ThTh++;
	  }else if(abs(Z.daughter(&genparticles,1)->pdgId())<=18){
	    n_ThThZll++;
	    n_ThTh++;

	  }else{
	    // cout<<"L599: Z Zerfall nicht bedacht"<<endl;
	  }
	}
      }
      //top  decays  hadronic
      if(W.daughter(&genparticles,1)){
	if(throw_failure){
	  cout << "(Top W hadronsich) WTochter pdgId()" << abs(W.daughter(&genparticles,1)->pdgId()) << endl;
	}
	  n_W++;
	if(abs(W.daughter(&genparticles,1)->pdgId())<=6){
	  n_ThWh++;
	  n_ThTh++;n_ThTll++;

	}else if(abs(W.daughter(&genparticles,1)->pdgId())<=18){
	  //muon zum Vergleich des Verzeigungsverhaeltnisses
	  if(abs(W.daughter(&genparticles,1)->pdgId())==13)muon_tprime = *W.daughter(&genparticles,1);
	  else muon_tprime =*W.daughter(&genparticles,2); 

	  n_ThWll++;
	  n_ThTh++;n_ThTll++;
	}else{
	  // cout<<"L617: W Zerfall nicht bedacht"<<endl;
	}
      }
    }
  }
  //  cout << "L:479 (GenSelections) "<< endl; 
  //vector<int> 
  n_vector = { n_TllThgg, n_TllThZZ,n_TllThWW, n_TllThbb,n_ThTllgg,n_ThTllZZ, n_ThTllWW,n_ThTllbb,n_TllWh, n_ThWll, n_ThThZll, n_TllThZh,n_ThTllZh,n_H, n_W, n_Z, n_TllTh  } ;
  if(leptop.daughter(&genparticles,1)&& leptop.daughter(&genparticles,2)){
  if(abs(leptop.daughter(&genparticles,1)->pdgId()) == 24 ) Wleptop =  *leptop.daughter(&genparticles,1);
  else Wleptop = *leptop.daughter(&genparticles,2);
  }
 if(hadtop.daughter(&genparticles,1)&& hadtop.daughter(&genparticles,2)){
  if(abs(hadtop.daughter(&genparticles,1)->pdgId()) == 24 ) Whadtop =  *hadtop.daughter(&genparticles,1);
  else Whadtop = *hadtop.daughter(&genparticles,2);
  }

  //return n_vector;
}
vector<int> ZPrimeGen::GetDecayChannel(){return n_vector;}



//Get the particle
GenParticle ZPrimeGen::Top() const{
  return main_top;
}
GenParticle ZPrimeGen::TopW() const{
  return top_W;
}
GenParticle ZPrimeGen::WBoson() const{
  return W;
}
GenParticle ZPrimeGen::ZBoson() const{
  return Z;
}
GenParticle ZPrimeGen::Higgs() const{
  return H;
}
GenParticle ZPrimeGen::TP() const{
  return TPrime;
}
GenParticle ZPrimeGen::ZP() const{
  return ZPrime;
}
GenParticle ZPrimeGen::D1() const{
  return daughter1;
}
GenParticle ZPrimeGen::D2() const{
  return daughter2;
}
GenParticle ZPrimeGen::D3() const{
  return daughter3;
}
GenParticle ZPrimeGen::D5() const{
  return daughter5;
}
GenParticle ZPrimeGen::D6() const{
  return daughter6;
}

GenParticle ZPrimeGen::G() const{
  return gluon;
}
GenParticle ZPrimeGen::ATop()const{
  if(topHiggs.pdgId()){return topHiggs;}else{return topZ;}
}
GenParticle ZPrimeGen::BLep()const{
  return blep;
}
GenParticle ZPrimeGen::WLep()const{
  return Wleptop;
}
GenParticle ZPrimeGen::WHad()const{
  return Whadtop;
}
GenParticle ZPrimeGen::Neutrino()const{
  return neutrino;
}
GenParticle ZPrimeGen::Lepton()const{
  return lepton;
}
GenParticle ZPrimeGen::WHiggsTop()const{
  return WHiggsTop1;
}
GenParticle ZPrimeGen::BHiggsTop()const{
  return BHiggsTop1;
}
GenParticle ZPrimeGen::lepTop()const{
  return leptop;
}
GenParticle ZPrimeGen::hadTop()const{
  return hadtop;
}

GenParticle ZPrimeGen::muon1()const{
  return muon_top;
}

GenParticle ZPrimeGen::muon2()const{
  return muon_tprime;
}
GenParticle ZPrimeGen::Quark1()const{
  return quark1;
}
GenParticle ZPrimeGen::Quark2()const{
  return quark2;
}
GenParticle ZPrimeGen::BHad()const{
  return bhad;
}


//Get the PT of the leading particles
float ZPrimeGen::Main_Top_PT()const{return main_top.pt();}
float ZPrimeGen::TPrime_PT()const{return TPrime.pt();}
float ZPrimeGen::ZPrime_PT()const{return ZPrime.pt();}
float ZPrimeGen::Daughter3_PT()const{return daughter3.pt();}

//Get the PT of the decay particle
float ZPrimeGen::Higgs_PT()const{return H.pt();}
float ZPrimeGen::TopHiggs_PT()const{return topHiggs.pt();}
float ZPrimeGen::TopZ_PT()const{return topZ.pt();}
float ZPrimeGen::WTPrime_PT()const{return W.pt();}
float ZPrimeGen::WTop_PT()const{return top_W.pt();}
float ZPrimeGen::BTop_PT()const{return top_b.pt();}
float ZPrimeGen::BTPrime_PT()const{return bottom.pt();}
float ZPrimeGen::Z_PT()const{return Z.pt();}

//Get the mass of the particles
float ZPrimeGen::ZPrime_Mass()const{return ZPrime.v4().mass();}
float ZPrimeGen::TPrime_Mass()const{return TPrime.v4().mass();}
float ZPrimeGen::ZPrime_invMass()const{return ZPrimereco.v4().mass();}
float ZPrimeGen::invMass1()const{return TPrimereco3.v4().mass();}
float ZPrimeGen::invMass2()const{return TPrimereco1.v4().mass();}
float ZPrimeGen::invMass3()const{return TPrimereco2.v4().mass();}

//Get delta phi
float ZPrimeGen::TPrimeTop_Phi()const{return abs(TPrime.phi()- main_top.phi());}
float ZPrimeGen::HT_Phi()const{return abs(H.phi()- topHiggs.phi());}
float ZPrimeGen::ZT_Phi()const{return abs(Z.phi()- topZ.phi());}
float ZPrimeGen::TopWB_Phi()const{return abs(W.phi()- bottom.phi());}
float ZPrimeGen::TPrimeWB_Phi()const{return abs(top_W.phi()- top_b.phi());}

//Get eta
float ZPrimeGen::TPrime_Eta()const{return TPrime.eta();}
float ZPrimeGen::Top_Eta()const{return main_top.eta();}
float ZPrimeGen::HT_Eta()const{return H.eta();}
float ZPrimeGen::ZT_Eta()const{return Z.eta();}
float ZPrimeGen::TPrimeWB_Eta()const{return W.eta();}
float ZPrimeGen::TopWB_Eta()const{return top_W.eta();}

//Get delta R
float ZPrimeGen::TPrimeTop_R()const{return abs(deltaR(main_top,TPrime));}
float ZPrimeGen::HT_R()const{return abs(deltaR(topHiggs,H));}
float ZPrimeGen::ZT_R()const{return abs(deltaR(Z,topZ));}
float ZPrimeGen::TopWB_R()const{return abs(deltaR(top_W, top_b));}
float ZPrimeGen::TPrimeWB_R()const{return abs(deltaR(W,bottom));}
float ZPrimeGen::TPrimeWG_R()const{return abs(deltaR(W,Gluon));}
float ZPrimeGen::TPrimeBG_R()const{return abs(deltaR(Bottom, Gluon));}
float ZPrimeGen::TPrimeHT_R()const{return deltaR(TPrimedaughter3,Top2);}
float ZPrimeGen::TPrimeHAntiT_R()const{return deltaR(TPrimedaughter3,Top1);}

float ZPrimeGen::Higgs_R()const{return deltaR(decaydaughterH1,decaydaughterH2);}
float ZPrimeGen::ZBoson_R()const{return deltaR(decaydaughterZ1,decaydaughterZ2);}
float ZPrimeGen::WBoson_R()const{return deltaR(decaydaughterW1,decaydaughterW2);}





ZPrimeGenProducer::ZPrimeGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
   h_zprimegen = ctx.get_handle<ZPrimeGen>(name);
}


bool ZPrimeGenProducer::process(Event & event){
    event.set(h_zprimegen, ZPrimeGen(*event.genparticles));
    return true;
}
