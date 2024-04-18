#ifndef __JET_H_
#define __JET_H_

#include <unordered_map>
#include <random>

#include "analysis_suite/skim/interface/Particle.h"

#include <complex>

enum PUID { PU_Tight = 7};

class Jet : public Particle {
public:
    void setup(TTreeReader& fReader, std::vector<Systematic> used_jec_systs_);

    virtual float getScaleFactor() override;
    float getTotalBTagWeight(std::string btag_wp = "L");
    float getTotalShapeWeight();
    float getBJetWeight(size_t idx, std::string lvl);

    float pt_(size_t idx) const override { return m_pt.at(idx)*m_jec->at(idx); }
    float nompt(size_t idx) const { return m_pt.at(idx); }
    float rawPt(size_t idx) const { return (1-rawFactor.at(idx))*m_pt.at(idx); }

    float getHT(Level level, size_t syst) { return getHT(list(level, syst)); };
    float getHT(Level level) { return getHT(list(level)); };
    float getMHT();

    float getCentrality(Level level, size_t syst) { return getCentrality(list(level, syst)); };
    float getCentrality(Level level) { return getCentrality(list(level)); };

    std::complex<float> get_momentum_change();

    virtual void setupGoodLists() override
    {
        n_loose_bjet.push_back(0);
        n_medium_bjet.push_back(0);
        n_tight_bjet.push_back(0);

        createLooseList();
        createBJetList();
        createTightList();
    }

    virtual void clear() override
    {
        Particle::clear();
        closeJetDr_by_index.clear();
        n_loose_bjet.clear();
        n_medium_bjet.clear();
        n_tight_bjet.clear();

        for (auto& [syst, var_scales] : m_jet_scales ) {
            for (auto& [var, scales] : var_scales) {
                scales.clear();
            }
        }
    }

    std::unordered_map<size_t, size_t> closeJetDr_by_index;

    std::vector<Int_t> n_loose_bjet, n_medium_bjet, n_tight_bjet;

    TRArray<Int_t> jetId;
    TRArray<Int_t> hadronFlavour;
    TRArray<Float_t> btag;
    TRArray<Int_t> genJetIdx;
    TRArray<Float_t> area;
    TRArray<Float_t> rawFactor;
    TRArray<Int_t> puId;
    TRVariable<Float_t> rho;

    TRArray<Float_t> neHEF, neEmEF, muEF, chHEF, chEmEF;
    TRArray<UChar_t> nConstituents;


    void setSyst(size_t syst) override;
    void setupJEC(GenericParticle& genJet);

    std::pair<Float_t, Float_t> get_JEC_pair(Systematic syst, size_t idx) const
    {
        if (m_jet_scales.find(syst) == m_jet_scales.end() ||
            m_jet_scales.at(syst).size() == 0) {
            return std::make_pair(1., 1.);
        }
        const auto scales = m_jet_scales.at(syst);
        const auto central = m_jet_scales.at(Systematic::Nominal).at(eVar::Nominal);
        return std::make_pair(scales.at(eVar::Down).at(idx)/central.at(idx), scales.at(eVar::Up).at(idx)/central.at(idx));
    }
    float loose_bjet_cut, medium_bjet_cut, tight_bjet_cut;

    void fillJet(JetOut& output, Level level, const Bitmap& event_bitmap);
    void fillJetEff(BEffOut& output, Level level, const Bitmap& event_bitmap);

private:
    int looseId = 0b11;
    int tightId = 0b10;
    float jet_dr = 0.4;
    std::unordered_map<Systematic, std::unordered_map<eVar, std::vector<float>>> m_jet_scales;

    void setup_jes(size_t i);
    void setup_jer(size_t i, GenericParticle& genJets);

    std::pair<float, float> get_jec_unc(size_t i, correction::Correction::Ref& jec_unc);
    std::vector<float> get_jer(size_t i, GenericParticle& genJets);

    std::vector<float>* m_jec;

    void createLooseList();
    void createBJetList();
    void createTightList();

    float getHT(const std::vector<size_t>& jet_list);
    float getCentrality(const std::vector<size_t>& jet_list);

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

    const std::unordered_map<Systematic, std::string> shape_by_syst = {
        { Systematic::BJet_Shape_hf, "hf" },
        { Systematic::BJet_Shape_hfstats1, "hfstats1" },
        { Systematic::BJet_Shape_hfstats2, "hfstats2" },
        { Systematic::BJet_Shape_lf, "lf" },
        { Systematic::BJet_Shape_lfstats1, "lfstats1" },
        { Systematic::BJet_Shape_lfstats2, "lfstats2" },
        { Systematic::BJet_Shape_cferr1, "cferr1" },
        { Systematic::BJet_Shape_cferr2, "cferr2" },
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

    std::unordered_map<Systematic, std::string> bjet_jec_syst = {
        { Systematic::Jet_JEC_Absolute, "jesAbsolute_20" },
        { Systematic::Jet_JEC_Absolute_corr, "jesAbsolute" },
        { Systematic::Jet_JEC_BBEC1, "jesBBEC1_20" },
        { Systematic::Jet_JEC_BBEC1_corr, "jesBBEC1" },
        { Systematic::Jet_JEC_EC2, "jesEC2_20" },
        { Systematic::Jet_JEC_EC2_corr, "jesEC2" },
        { Systematic::Jet_JEC_HF, "jesHF_20" },
        { Systematic::Jet_JEC_HF_corr, "jesHF" },
        { Systematic::Jet_JEC_RelativeBal, "jesRelativeBal" },
        { Systematic::Jet_JEC_RelativeSample, "jesRelativeSample_20" },
        { Systematic::Jet_JEC_FlavorQCD, "jesFlavorQCD" },
    };


    const std::vector<Systematic> charm_systs = {
        Systematic::BJet_Shape_cferr1,
        Systematic::BJet_Shape_cferr2,
    };

    WeightHolder veto_map;
    WeightHolder jer_scale, jet_resolution, jes_scale;
    WeightHolder puid_scale;
    WeightHolder btag_bc_scale, btag_udsg_scale, btag_eff, btag_shape_scale;
    std::unordered_map<Systematic, correction::Correction::Ref> jec_unc_vec;

    bool use_shape_btag = false;
    std::vector<Systematic> used_jec_systs;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::unordered_map<std::string, float> bjet_cuts;

};

#endif // __JET_H_
