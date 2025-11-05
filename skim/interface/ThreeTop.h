#ifndef THREETOP_H
#define THREETOP_H

// #include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

class ThreeTop : public DileptonBase {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void ApplyChannelScaleFactors() override;
    void clearOutputs() override;
    bool passPreselection() override;
    ClassDefOverride(ThreeTop, 0);

private:
    void printStuff();
    void xor_lists(Particle& part, std::vector<size_t>& outlist);

    bool baseline_cuts();
    bool signal_cuts();
    bool nonprompt_cuts();
    bool charge_misid_cuts();
    // bool opposite_sign(CutInfo cuts);

    float muPt(size_t idx) { return muon.size(Level::Fake) > idx ? muon.pt(Level::Fake, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Fake) > idx ? elec.pt(Level::Fake, idx) : -1; }

    LeptonOut* o_tightMuons;
    LeptonOut* o_tightElectrons;
    JetOut* o_jets;
    JetOut* o_bJets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality, o_ZMass;
    std::vector<Int_t> o_nb_loose, o_nb_medium, o_nb_tight;
    std::vector<Float_t> o_bwgt_loose, o_bwgt_medium, o_bwgt_tight;
    std::vector<std::vector<Float_t>> o_nb_weight;
    std::vector<Int_t> o_nloose_mu,o_nloose_el;
    std::vector<Bool_t> o_pass_zveto;
    std::vector<Float_t> o_wgt_nobtag;
    std::vector<Float_t> masses, os_masses;
    float o_raw_met, o_raw_metphi;

    bool o_hasVetoJet;

    bool use_nonprompt = false;
    // TrigEff trigEff_leadPt;
};

#endif
