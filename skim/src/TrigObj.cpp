#include "analysis_suite/skim/interface/TrigObj.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

#include "analysis_suite/skim/interface/Muon.h"
#include "analysis_suite/skim/interface/Electron.h"

#include <bitset>

void TrigObj::setup(TTreeReader& fReader)
{
    id.setup(fReader, "TrigObj_id");
    filterBits.setup(fReader, "TrigObj_filterBits");

    m_pt.setup(fReader, "TrigObj_pt");
    m_eta.setup(fReader, "TrigObj_eta");
    m_phi.setup(fReader, "TrigObj_phi");

    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) {
        hlt_isomu.setup(fReader, "HLT_IsoMu24");
        isomu_cut = 27.;
        hlt_mu50.setup(fReader, "HLT_Mu50");
        hlt_el.setup(fReader, "HLT_Ele27_WPTight_Gsf");
        el_cut = 30.;
    } else if (year_ == Year::yr2017) {
        hlt_isomu.setup(fReader, "HLT_IsoMu27");
        isomu_cut = 30.;
        hlt_mu50.setup(fReader, "HLT_Mu50");
        hlt_el.setup(fReader, "HLT_Ele35_WPTight_Gsf");
        el_cut = 38.;
    } else if (year_ == Year::yr2018) {
        hlt_isomu.setup(fReader, "HLT_IsoMu24");
        isomu_cut = 27.;
        hlt_mu50.setup(fReader, "HLT_Mu50");
        hlt_el.setup(fReader, "HLT_Ele32_WPTight_Gsf");
        el_cut = 35.;
    }

    // setup_map(Level::Loose);
    // setup_map(Level::Fake);
    // setup_map(Level::Tight);
}

void TrigObj::match_muon(Muon& muon)
{
    for (size_t i = 0; i < size(); i++) {
        std::bitset<16> bits(filterBits.at(i));
        if (id.at(i) != static_cast<int>(PID::Muon)) continue;
        if (!(*hlt_isomu && (bits & std::bitset<16>(10)).any() && pt(i) > isomu_cut)
            && !(*hlt_mu50 && bits[10]==1 && pt(i) > 50)) continue;

        size_t minIdx = 0;
        float minDr = 100;
        float trig_phi = phi(i);
        float trig_eta = eta(i);
        for (auto mu : muon.list(Level::Loose)) {
            float dr2 = deltaR2(muon.eta(mu), trig_eta, muon.phi(mu), trig_phi);
            if (minDr > dr2) {
                minIdx = mu;
                minDr = dr2;
            }
        }
        if (minDr < 0.01
            && muon.tightCharge.at(minIdx) == 2
            && muon.sip3d.at(minIdx) < 4) {
            muon_match.insert(minIdx);
            // muon_match[i] = minIdx;
                // std::cout << "Muon: " << muon.pt(minIdx) << " " << bits << " " << filterBits.at(i) << " " << *hlt_isomu << " " << *hlt_mu50 << " " << *hlt_dimu << std::endl;
        }
    }
}

/*
 * 1 = TrkIsoVVL,
 * 2 = Iso,
 * 3 = OverlapFilter PFTau,
 * 4 = 1mu,
 * 5 = 2mu,
 * 6 = 1mu-1e,
 * 7 = 1mu-1tau,
 * 8 = 3mu,
 * 9 = 2mu-1e,
 * 10 = 1mu-2e,
 * 11 = 1mu (Mu50),
 * 12 = 1mu (Mu100) for Muon;
 */

void TrigObj::match_electron(Electron& elec)
{
    for (size_t i = 0; i < size(); i++) {
        std::bitset<16> bits(filterBits.at(i));
        if (id.at(i) != static_cast<int>(PID::Electron)) continue;
        if (!(*hlt_el && bits[1] == 1 && pt(i) > el_cut)) continue; // Doesn't match trigger
        size_t minIdx = 0;
        float minDr = 100;
        float trig_phi = phi(i);
        float trig_eta = eta(i);
        for (auto el : elec.list(Level::Loose)) {
            float dr2 = deltaR2(elec.eta(el), trig_eta, elec.phi(el), trig_phi);
            if (minDr > dr2) {
                minIdx = el;
                minDr = dr2;
            }
        }
        if (minDr < 0.01
            && elec.lostHits.at(minIdx) == 0
            && elec.tightCharge.at(minIdx) == 2
            && elec.sip3d.at(minIdx) < 4
            && elec.passTriggerRequirements(minIdx)) {
            elec_match.insert(minIdx);
            // elec_match.push_back(std::make_pair(i, minIdx));
                // std::cout << "Electron: " << elec.pt(minIdx) << " " <<  bits << " " << filterBits.at(i) << " " << *hlt_el << std::endl;
        }
    }
}

/*
 * 1 = CaloIdL_TrackIdL_IsoVL,
 * 2 = 1e (WPTight),
 * 3 = 1e (WPLoose),
 * 4 = OverlapFilter PFTau,
 * 5 = 2e,
 * 6 = 1e-1mu,
 * 7 = 1e-1tau,
 * 8= 3e,
 * 9 = 2e-1mu,
 * 10= 1e-2mu,
 * 11 = 1e (32_L1DoubleEG_AND_L1SingleEGOr),
 * 12 = 1e (CaloIdVT_GsfTrkIdT),
 * 13 = 1e (PFJet),
 * 14 = 1e (Photon175_OR_Photon200) for Electron (PixelMatched e/gamma);
 */
