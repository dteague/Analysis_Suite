#include "analysis_suite/skim/interface/Met.h"
#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

#include "analysis_suite/skim/interface/XYMETCorrection.h"

void Met::setup(TTreeReader& fReader, MET_Type type)
{
    // name = met_name[type];
    name = "PuppiMET";

    ispuppi = (type == MET_Type::PUPPI);
    m_pt.setup(fReader, (name+"_pt").c_str());
    m_phi.setup(fReader, (name+"_phi").c_str());
    // m_raw_pt.setup(fReader, "RawMET_pt");
    // m_raw_phi.setup(fReader, "RawMET_phi");

    if (ispuppi && isMC) {
        m_jer_pt_up.setup(fReader, "PuppiMET_ptJERUp");
        m_jer_pt_down.setup(fReader, "PuppiMET_ptJERDown");
        m_jer_phi_up.setup(fReader, "PuppiMET_phiJERUp");
        m_jer_phi_down.setup(fReader, "PuppiMET_phiJERDown");
    }

    m_corr_pt[Systematic::Nominal] = {{eVar::Nominal, 0}};
    m_corr_phi[Systematic::Nominal] = {{eVar::Nominal, 0}};
    for (auto syst : jec_systs) {
        m_corr_pt[syst] = {{eVar::Up, 0}, {eVar::Down, 0}};
        m_corr_phi[syst] = {{eVar::Up, 0}, {eVar::Down, 0}};
    }
}

void Met::setSyst()
{
    if (currentSyst == Systematic::Nominal || isJECSyst()) {
        corr_pt = &m_corr_pt[currentSyst][currentVar];
        corr_phi = &m_corr_phi[currentSyst][currentVar];
    } else {
        corr_pt = &m_corr_pt[Systematic::Nominal][eVar::Nominal];
        corr_phi = &m_corr_phi[Systematic::Nominal][eVar::Nominal];
    }
}

void Met::setupMet(Jet& jet, UInt_t run, int nVertices)
{
    if (currentSyst == Systematic::Nominal || isJECSyst()) {
        // Get change in Met from change in Jet momentum
        if (ispuppi && currentSyst == Systematic::Jet_JER) {
            (*corr_pt) = (currentVar == eVar::Up) ? *m_jer_pt_up : *m_jer_pt_down;
            (*corr_phi) = (currentVar == eVar::Up) ? *m_jer_phi_up : *m_jer_phi_down;
        }
        else {
            // auto met_vec = std::polar(*m_raw_pt, *m_raw_phi) - jet.get_momentum_change();
            auto met_vec = std::polar(*m_pt, *m_phi) - jet.get_momentum_change();
            (*corr_pt) = std::abs(met_vec);
            (*corr_phi) = std::arg(met_vec);
        }
        fix_xy(run, nVertices);
    }
}

void Met::fix_xy(UInt_t run, int nVertices)
{
    auto met_corr = METXYCorr_Met_MetPhi(pt(), phi(), run, yearMap.at(year_), isMC, nVertices, true, ispuppi);
    (*corr_pt) = met_corr.first;
    (*corr_phi) = met_corr.second;

}

float Met::mt(float l_pt, float l_phi)
{
    return sqrt(2*pt()*l_pt*(1 - cos(phi() - l_phi)));
}
