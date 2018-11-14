#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesis.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeReconstructionHypothesisDiscriminators.h"
#include "UHH2/ZprimeToTprimeTtZtH/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/PDFWeights.h" 
#include "UHH2/common/include/TTbarGen.h"


namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class SR_PDFHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    SR_PDFHists(uhh2::Context & ctx, const std::string & dirname, const std::string & hyps_name, const std::string & discriminator_name, bool use_pdf_weights_ = false);

    virtual void fill(const uhh2::Event & ev) override;

    std::string histo_names[100];
    std::string histo_names_toppt[100];
    std::string histo_names_pt[100];
    std::string histo_names_mass[100];
    std::string histo_names_eta[100];
    std::string histo_names_MET[100];


  protected: 

    uhh2::Event::Handle<std::vector<ZPrimeTotTPrimeReconstructionHypothesis>> h_hyps;
    uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::string m_discriminator_name;
    bool use_pdf_weights, is_mc;
    std::string dirnames;
    
    std::unique_ptr<PDFWeights> m_pdfweights;


    virtual ~SR_PDFHists();
};

}
