#include "analysis_suite/skim/interface/JEC.hpp"
#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void JEC::init(TTreeReader& fReader, std::vector<Systematic> used_jec_systs_)
{
    used_jec_systs = used_jec_systs_;
    rho.setup(fReader, "fixedGridRhoFastjetAll");

    m_jet_scales[Systematic::Nominal] = {{eVar::Nominal, std::vector<float>()}};
    for (auto syst : used_jec_systs) {
        m_jet_scales[syst] = {
            {eVar::Up, std::vector<float>()},
            {eVar::Down, std::vector<float>()}
        };
    }

    // JEC Weights
    auto corr_set = getScaleFile("JME", "jet_jerc");
    if(isMC) {
        jet_resolution = WeightHolder(corr_set->at(jer_source[year_]+"_PtResolution_AK4PFchs"));
        jer_scale = WeightHolder(corr_set->at(jer_source[year_]+"_ScaleFactor_AK4PFchs"),
                                 Systematic::Jet_JER, {"nom","up","down"});

        // JES uncertatiny
        auto unc_set = getScaleFile("JME", "jec_unc");
        for (auto syst: used_jec_systs) {
            if (syst == Systematic::Jet_JER) continue;
            jec_unc_vec[syst] = unc_set->at(jes_source[year_]+unc_by_syst.at(syst)+"_AK4PFchs");
        }

    }
}


void JEC::setupJEC(Jet& jet, GenericParticle& genJet)
{
    for (auto& [syst, var_scales] : m_jet_scales ) {
        for (auto& [var, scales] : var_scales) {
            scales.clear();
        }
    }

    m_met_change[Systematic::Nominal][eVar::Nominal] = std::polar(0.f, 0.f);

    if (!isMC) {
        m_jet_scales[Systematic::Nominal][eVar::Nominal].assign(jet.size(), 1);
        return;
    } else {
        for (auto& [syst, var_scales] : m_jet_scales ) {
            for (auto& [var, scales] : var_scales) {
                scales.resize(jet.size());
                m_met_change[syst][var] = std::polar(0.f, 0.f);
            }
        }
    }

    for(size_t i = 0; i < jet.size(); ++i) {
        float pt = jet.nomPt(i);
        float eta = jet.eta(i);
        float phi = jet.phi(i);

        auto jer = get_jer(pt, eta, phi, jet.genJetIdx.at(i), genJet);
        float central = jer.at(0);

        m_jet_scales[Systematic::Nominal][eVar::Nominal][i] = central;
        m_met_change[Systematic::Nominal][eVar::Nominal] += std::polar((central-1)*pt, phi);
        for (auto syst: used_jec_systs) {
            if (syst == Systematic::Jet_JER) {
                m_jet_scales[syst][eVar::Up][i] = jer.at(1);
                m_met_change[syst][eVar::Up] += std::polar((jer.at(1)-1)*pt, phi);

                m_jet_scales[syst][eVar::Down][i] = jer.at(2);
                m_met_change[syst][eVar::Down] += std::polar((jer.at(2)-1)*pt, phi);
            } else {
                float jec_unc = jec_unc_vec[syst]->evaluate({eta, pt});
                m_jet_scales[syst][eVar::Down][i] = central*(1-jec_unc);
                m_met_change[syst][eVar::Down] += std::polar((central*(1-jec_unc)-1)*pt, phi);

                m_jet_scales[syst][eVar::Up][i] = central*(1+jec_unc);
                m_met_change[syst][eVar::Up] += std::polar((central*(1+jec_unc)-1)*pt, phi);
            }
        }
    }
}


std::vector<float> JEC::get_jer(float pt, float eta, float phi, int gIdx, GenericParticle& genJets)
{
    if (!isMC) {
        return {1.0};
    }

    float resolution = jet_resolution.evaluate({eta, pt, *rho});
    float rand = 1.;
    float pt_ratio = (gIdx != -1) ? 1-genJets.pt(gIdx)/pt : 0.;

    bool matched_gen = (gIdx != -1)
        && (deltaR2(eta, genJets.eta(gIdx), phi, genJets.phi(gIdx)) < 0.04) // (0.4/2)**2
        && (fabs(pt_ratio) < 3*resolution);

    if (!matched_gen) {
        std::normal_distribution<> gaussian{0, resolution};
        rand = gaussian(gen);
    }

    std::vector<float> output;
    for (auto var : { eVar::Nominal, eVar::Up, eVar::Down }) {
        float scale = jer_scale.evaluate({eta, jer_scale.name_by_var.at(var)});
        if (matched_gen) {
            output.push_back(1 + (scale-1)*pt_ratio);
        } else if (scale > 1.) {
            output.push_back(1 + rand*sqrt(pow(scale,2) - 1));
        } else {
            output.push_back(1);
        }
    }
    return output;
}

void JEC::set_syst(bool isJECSyst)
{
    if (currentSyst == Systematic::Nominal || isJECSyst) {
        m_jec = &m_jet_scales[currentSyst][currentVar];
        m_met_diff = &m_met_change[currentSyst][currentVar];
    } else {
        m_jec = &m_jet_scales[Systematic::Nominal][eVar::Nominal];
        m_met_diff = &m_met_change[Systematic::Nominal][eVar::Nominal];
    }
}
