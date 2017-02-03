#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/TTbarGen.h"

/** \brief Common histograms for reconstruction hypotheses
 *
 * hyps_name is the name of the reconstruction hypothesis collection, for instance "HighMassReconstruction"
 * discriminator_name is the name of the discriminator used to choose the best reconstruction hypothesis, for instance "Chi2"
 */
class ZPrimeTotTPrimeHypothesisHists: public uhh2::Hists {
public:
    ZPrimeTotTPrimeHypothesisHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F *Discriminator, *Discriminator_2, *Discriminator_3,*DiscriminatorH;
    TH1F *M_ttbar_rec, *M_ttbar_gen,*M_Higgs_rec, *M_Higgs_gen, *M_Z_gen, *M_W_gen,*M_ZPrime_rec, *M_ZPrime_gen,*M_TPrime_rec,*M_ZPrime_scaled, *M_TPrime_gen, *M_toplep_rec, *M_tophad_rec,  *M_toplep_gen, *M_tophad_gen,*M_tophad_rec_1jet, *M_tophad_rec_2jet, *M_tophad_rec_3jet,*M_TPrime_rec2,*M_TPrime_rec3,*M_TPrime_rec4;
    TH1F *Pt_toplep_rec, *Pt_tophad_rec, *Pt_ttbar_rec, *Pt_ttbar_gen, *Pt_Higgs_rec, *Pt_Higgs_gen, *Pt_ZPrime_rec, *Pt_ZPrime_gen, *Pt_TPrime_rec, *Pt_TPrime_gen;
    TH2F *Pt_ttbar_rec_vs_Pt_ttbar_gen,*Pt_Higgs_rec_vs_Pt_Higgs_gen,*Pt_ZPrime_rec_vs_Pt_ZPrime_gen,*Pt_TPrime_rec_vs_Pt_TPrime_gen, *Discriminator_vs_M_ttbar, *M_ttbar_rec_vs_M_ttbar_gen,*M_Higgs_rec_vs_M_Higgs_gen,*M_TPrime_rec_vs_M_TPrime_gen,*M_ZPrime_rec_vs_M_ZPrime_gen, *M_ttbar_gen_vs_M_ttbar_rec, *M_Higgs_gen_vs_M_Higgs_rec, *M_ZPrime_gen_vs_M_ZPrime_rec, *M_TPrime_gen_vs_M_TPrime_rec;
    TH1F *M_ttbar_resolution;
    TH2F *Pt_toplep_rec_vs_Pt_toplep_gen, *Pt_tophad_rec_vs_Pt_tophad_gen, *Pt_lepton_rec_vs_Pt_lepton_gen, *Pt_blep_rec_vs_Pt_blep_gen, *Pt_neutrino_rec_vs_Pt_neutrino_gen;
    TH2F *eta_toplep_rec_vs_eta_toplep_gen, *eta_tophad_rec_vs_eta_tophad_gen, *eta_lepton_rec_vs_eta_lepton_gen, *eta_blep_rec_vs_eta_blep_gen, *eta_neutrino_rec_vs_eta_neutrino_gen;
    TH1F *massfitHiggs, *massfitZ, *missmatch,*massfitlepTop,*massfithadTop;


    /////////////////////////Eta/////////////////////////////
    TH1F *eta_toplep_rec, *eta_tophad_rec,*eta_toplep_gen, *eta_tophad_gen, *eta_Higgs_rec, *eta_Higgs_gen, *eta_ZPrime_rec, *eta_ZPrime_gen, *eta_TPrime_rec, *eta_TPrime_gen;
    ////////////////////////Rec vs. Gen///////////////////////////////////
    TH1F *M_Higgs_vs,*M_Z_vs,*M_W_vs, *M_ZPrime_vs, *M_TPrime_vs, *M_toplep_vs, *M_tophad_vs, *Pt_Higgs_vs, *Pt_ZPrime_vs, *Pt_TPrime_vs, *Pt_toplep_vs, *Pt_tophad_vs, *eta_Higgs_vs, *eta_ZPrime_vs, *eta_TPrime_vs, *eta_toplep_vs, *eta_tophad_vs;



    uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::string m_discriminator_name;
    
};
