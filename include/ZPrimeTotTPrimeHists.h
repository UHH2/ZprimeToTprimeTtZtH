#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/ZPrimeTotTPrime/include/ZPrimeTotTPrimeGenSelections.h"
#include "UHH2/common/include/TTbarGen.h"

namespace uhh2examples {

/**  \brief Example class for booking and filling histograms
 * 
 * NOTE: This class uses the 'hist' method to retrieve histograms.
 * This requires a string lookup and is therefore slow if you have
 * many histograms. Therefore, it is recommended to use histogram
 * pointers as member data instead, like in 'common/include/ElectronHists.h'.
 */
class ZPrimeTotTPrimeHists: public uhh2::Hists {
public:
    // use the same constructor arguments as Hists for forwarding:
    ZPrimeTotTPrimeHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~ZPrimeTotTPrimeHists();
  private:
    uhh2::Event::Handle<ZPrimeGen> h_zprimegen;
    uhh2::Event::Handle<std::vector<TopJet> > h_toptag;
    uhh2::Event::Handle<std::vector<TopJet> > h_higgstag;
    uhh2::Event::Handle<std::vector<TopJet> > h_ZWtag;
    uhh2::Event::Handle<TTbarGen> h_ttbargen;
    std::string btopjets;
    bool isMC=false;;
};

}
