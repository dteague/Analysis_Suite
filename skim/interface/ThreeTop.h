#ifndef THREETOP_H
#define THREETOP_H

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"
#include "analysis_suite/skim/interface/ResolvedTop.h"

class ThreeTop : public BaseSelector {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void ApplyChannelScaleFactors() override;
    void clearParticles() override;
    void clearOutputs() override;
    void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(ThreeTop, 0);

private:
    void printStuff();
    float getLeadPt(size_t idx = 0);
    void setSubChannel();
    bool isSameSign(Level level);
    bool getTriggerValue();
    void xor_lists(Particle& part, std::vector<size_t>& outlist);

    bool baseline_cuts();
    bool signal_cuts();
    bool nonprompt_cuts();
    bool charge_misid_cuts();
    bool opposite_sign();

    float muPt(size_t idx) { return muon.size(Level::Fake) > idx ? muon.pt(Level::Fake, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Fake) > idx ? elec.pt(Level::Fake, idx) : -1; }

    ResolvedTop rTop;

    ParticleOut* o_tightMuons;
    ElectronOut* o_tightElectrons;
    JetOut* o_jets;
    JetOut* o_bJets;
    TopOut* o_top;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality, o_ZMass;
    std::vector<Int_t> o_nb_loose, o_nb_medium, o_nb_tight;
    std::vector<Float_t> o_bwgt_loose, o_bwgt_medium, o_bwgt_tight;
    std::vector<Bool_t> o_hlt_dilepton_ht, o_hlt_dilepton, o_hlt_trilepton;
    std::vector<Bool_t> o_hltind_dilepton_ht, o_hltind_dilepton, o_hltind_trilepton;
    std::vector<std::vector<Float_t>> o_nb_weight;
    std::vector<Int_t> o_nloose_mu,o_nloose_el;
    std::vector<Bool_t> o_pass_zveto;

    std::vector<size_t> fake_el, fake_mu;

    bool use_nonprompt = false;
    // TrigEff trigEff_leadPt;
};

#endif
