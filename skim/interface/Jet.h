#ifndef __JET_H_
#define __JET_H_

#include <unordered_map>
#include <random>
#include <complex>

#include "analysis_suite/skim/interface/Particle.h"
#include "analysis_suite/skim/interface/JEC.hpp"

enum PUID { PU_Tight = 7};

class JEC;

class Jet : public Particle {
public:
    void setup(TTreeReader& fReader, std::string data_run, std::vector<Systematic> used_jec_systs_);

    virtual float getScaleFactor() override;
    float getTotalBTagWeight(std::string btag_wp = "L");
    float getTotalShapeWeight();
    float getBJetWeight(size_t idx, std::string lvl);
    float getCutBasedBTagWeight(std::string btag_wp="M");

    float mass_(size_t idx) const override { return jec.get_shift(idx)*m_mass.at(idx); }
    float pt_(size_t idx) const override { return jec.get_shift(idx)*m_pt.at(idx); }
    float pt_shift(size_t idx, Systematic syst, eVar var) { return jec.get_shift(idx, syst, var)*m_pt.at(idx); }
    float mass_shift(size_t idx, Systematic syst, eVar var) { return jec.get_shift(idx, syst, var)*m_mass.at(idx); }
    float nomPt(size_t idx) const { return m_pt.at(idx); }
    float jec_shift(size_t i) { return jec.get_shift(i, Systematic::Nominal, eVar::Nominal); }
    void setup_shift_output(TTree* tree);

    float getHT(Level level, size_t syst) { return getHT(list(level, syst)); };
    float getHT(Level level) { return getHT(list(level)); };
    float getMHT();
    bool passVeto(float eta, float phi);

    float getCentrality(Level level, size_t syst) { return getCentrality(list(level, syst)); };
    float getCentrality(Level level) { return getCentrality(list(level)); };

    std::complex<float> get_momentum_change();

    void setup_jets(GenericParticle& genJet, size_t run);
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

    }

    std::unordered_map<size_t, size_t> closeJetDr_by_index;
    std::vector<float> closeLep_dr;

    std::vector<Int_t> n_loose_bjet, n_medium_bjet, n_tight_bjet;

    TRArray<Int_t> jetId;
    TRArray<Int_t> hadronFlavour;
    TRArray<Float_t> btag;
    TRArray<Int_t> genJetIdx;
    TRArray<Float_t> area;
    TRArray<Float_t> rawFactor;
    TRArray<Int_t> puId;
    TRArray<Float_t> muonSubtrFactor;
    TRVariable<Float_t> rho;

    TRArray<Float_t> neHEF, neEmEF, muEF, chHEF, chEmEF;
    TRArray<UChar_t> nConstituents;

    // Met stuff
    TRArray<Float_t> t1_rawpt, t1_eta, t1_phi, t1_area, t1_muonSub;
    correction::Correction::Ref jec_l1_total;
    correction::CompoundCorrection::Ref jec_total;

    void setSyst(size_t syst) override;
    void setupJEC(GenericParticle& genJet, size_t run);

    float loose_bjet_cut, medium_bjet_cut, tight_bjet_cut;

    void fillJet(JetOut& output, Level level, const Bitmap& event_bitmap);
    void fillJetEff(BEffOut& output, Level level, const Bitmap& event_bitmap);

private:
    int looseId = 0b11;
    int tightId = 0b10;

    JEC jec;
    std::vector<JECShiftOut*> o_jet_shifts;
    std::vector<JECShiftOut*> o_bjet_shifts;

    void createLooseList();
    void createBJetList();
    void createTightList();

    float getHT(const std::vector<size_t>& jet_list);
    float getCentrality(const std::vector<size_t>& jet_list);

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
    WeightHolder puid_scale;
    WeightHolder btag_bc_scale, btag_udsg_scale, btag_eff, btag_shape_scale;

    std::vector<Systematic> used_jec_systs;
    bool applyHEMVeto = false;
    bool inHEM = false;
    const float pi = 3.141592;

    std::unordered_map<std::string, float> bjet_cuts;
};

#endif // __JET_H_
