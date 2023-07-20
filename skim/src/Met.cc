#include "analysis_suite/skim/interface/Met.h"
#include "analysis_suite/skim/interface/Jet.h"

#include "analysis_suite/skim/interface/XYMETCorrection.h"

void Met::setup(TTreeReader& fReader, MET_Type type)
{
    name = met_name[type];
    ispuppi = (type == MET_Type::PUPPI);
    m_pt.setup(fReader, (name+"_pt").c_str());
    m_phi.setup(fReader, (name+"_phi").c_str());

    if (ispuppi && isMC) {
        m_puppi_pt_syst[Systematic::Jet_JES] = Puppi_pair();
        m_puppi_pt_syst[Systematic::Jet_JES].first.setup(fReader, "PuppiMET_ptJESDown");
        m_puppi_pt_syst[Systematic::Jet_JES].second.setup(fReader, "PuppiMET_ptJESUp");
        m_puppi_phi_syst[Systematic::Jet_JES] = Puppi_pair();
        m_puppi_phi_syst[Systematic::Jet_JES].first.setup(fReader, "PuppiMET_phiJESDown");
        m_puppi_phi_syst[Systematic::Jet_JES].second.setup(fReader, "PuppiMET_phiJESUp");

        m_puppi_pt_syst[Systematic::Jet_JER] = Puppi_pair();
        m_puppi_pt_syst[Systematic::Jet_JER].first.setup(fReader, "PuppiMET_ptJERDown");
        m_puppi_pt_syst[Systematic::Jet_JER].second.setup(fReader, "PuppiMET_ptJERUp");
        m_puppi_phi_syst[Systematic::Jet_JER] = Puppi_pair();
        m_puppi_phi_syst[Systematic::Jet_JER].first.setup(fReader, "PuppiMET_phiJERDown");
        m_puppi_phi_syst[Systematic::Jet_JER].second.setup(fReader, "PuppiMET_phiJERUp");
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
        if (!ispuppi) {
            auto met_vec = std::polar(*m_pt, *m_phi) - jet.get_momentum_change();
            (*corr_pt) = std::abs(met_vec);
            (*corr_phi) = std::arg(met_vec);
        } else if (isJECSyst()) {
            if (currentVar == eVar::Down) {
                (*corr_pt) = *m_puppi_pt_syst[currentSyst].first;
                (*corr_phi) = *m_puppi_phi_syst[currentSyst].first;
            } else {
                (*corr_pt) = *m_puppi_pt_syst[currentSyst].second;
                (*corr_phi) = *m_puppi_phi_syst[currentSyst].second;
            }
        } else {
            (*corr_pt) = *m_pt;
            (*corr_phi) = *m_phi;
        }
        pt_unfix = pt();
        phi_unfix = phi();
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
