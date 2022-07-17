#include "analysis_suite/Analyzer/interface/Met.h"
#include "analysis_suite/Analyzer/interface/Jet.h"

void Met::setup(TTreeReader& fReader)
{
    m_pt.setup(fReader, "type1_pfMETEt");
    m_phi.setup(fReader, "type1_pfMETPhi");

    auto corr_set = getScaleFile("USER", "met_xy");

    xcorr = WeightHolder(corr_set->at("x_correction"));
    ycorr = WeightHolder(corr_set->at("y_correction"));

    m_corr_pt[Systematic::Nominal] = {{eVar::Nominal, 0}};
    m_corr_phi[Systematic::Nominal] = {{eVar::Nominal, 0}};
    for (auto syst : jec_systs) {
        m_corr_pt[syst] = {{eVar::Up, 0}, {eVar::Down, 0}};
        m_corr_phi[syst] = {{eVar::Up, 0}, {eVar::Down, 0}};
    }
}

void Met::setupJEC(Jet& jet)
{
    if (currentSyst == Systematic::Nominal || jet.isJECSyst()) {
        corr_pt = &m_corr_pt[currentSyst][currentVar];
        corr_phi = &m_corr_phi[currentSyst][currentVar];
        auto delta_jet = jet.get_momentum_change();
        PolarVector met_vec(*m_pt, *m_phi);
        met_vec -= delta_jet;
        (*corr_pt) = met_vec.R();
        (*corr_phi) = met_vec.Phi();
    } else {
        corr_pt = &m_corr_pt[Systematic::Nominal][eVar::Nominal];
        corr_phi = &m_corr_pt[Systematic::Nominal][eVar::Nominal];
    }
}

void Met::fix_xy(UInt_t run, int nVertices)
{
    float corr_metx = *m_pt*cos(*m_phi)+xcorr.evaluate({"MET", (float)run, (float)nVertices});
    float corr_mety = *m_pt*sin(*m_phi)+ycorr.evaluate({"MET", (float)run, (float)nVertices});

    (*corr_pt) = sqrt(pow(corr_metx, 2) + pow(corr_mety, 2));
    (*corr_phi) = atan2(corr_mety, corr_metx);
}
