#ifndef JEC_H_
#define JEC_H_

#include "analysis_suite/skim/interface/Systematic.h"
#include "analysis_suite/skim/interface/Particle.h"

#include <complex>
#include <unordered_map>
#include <random>

class Jet;

class JEC : SystematicWeights {
public:
    void init(TTreeReader& fReader, std::vector<Systematic> used_jec_systs_);
    float get_shift(size_t i) const { return m_jec->at(i); }
    float get_shift(size_t i, Systematic syst, eVar var) { return m_jet_scales[syst][var].at(i); }
    void setupJEC(Jet& jet, GenericParticle& genJet);
    void set_syst(bool isJECSyst);
    std::complex<float> get_met_change() { return *m_met_diff; }

private:
    std::vector<float> get_jer(float pt, float eta, float phi, int gIdx, GenericParticle& genJets);

    TRVariable<Float_t> rho;

    std::unordered_map<Systematic, std::unordered_map<eVar, std::vector<float>>> m_jet_scales;
    std::unordered_map<Systematic, std::unordered_map<eVar, std::complex<float>>> m_met_change;
    std::vector<float>* m_jec;
    std::complex<float>* m_met_diff;

    std::vector<Systematic> used_jec_systs;
    WeightHolder jer_scale, jet_resolution, jes_scale;
    std::unordered_map<Systematic, correction::Correction::Ref> jec_unc_vec;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    std::unordered_map<Year, std::string> jes_source = {
        {Year::yr2016pre, "Summer19UL16APV_V7_MC"},
        {Year::yr2016post, "Summer19UL16_V7_MC"},
        {Year::yr2017, "Summer19UL17_V5_MC"},
        {Year::yr2018, "Summer19UL18_V5_MC"},
    };

    std::unordered_map<Year, std::string> jer_source = {
        {Year::yr2016pre, "Summer20UL16APV_JRV3_MC"},
        {Year::yr2016post, "Summer20UL16_JRV3_MC"},
        {Year::yr2017, "Summer19UL17_JRV2_MC"},
        {Year::yr2018, "Summer19UL18_JRV2_MC"},
    };

    const std::unordered_map<Systematic, std::string> unc_by_syst = {
        { Systematic::Jet_JES, "_Total" },
        { Systematic::Jet_JEC_Absolute, "_Absolute_uncorr" },
        { Systematic::Jet_JEC_Absolute_corr, "_Absolute" },
        { Systematic::Jet_JEC_BBEC1, "_BBEC1_uncorr" },
        { Systematic::Jet_JEC_BBEC1_corr, "_BBEC1" },
        { Systematic::Jet_JEC_EC2, "_EC2_uncorr" },
        { Systematic::Jet_JEC_EC2_corr, "_EC2" },
        { Systematic::Jet_JEC_HF, "_HF_uncorr" },
        { Systematic::Jet_JEC_HF_corr, "_HF" },
        { Systematic::Jet_JEC_RelativeBal, "_RelativeBal" },
        { Systematic::Jet_JEC_RelativeSample, "_RelativeSample_uncorr" },
        { Systematic::Jet_JEC_FlavorQCD, "_FlavorQCD" },
    };
};

#endif // JEC_H_
