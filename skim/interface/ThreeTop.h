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
    void clearParticles() override;
    void clearOutputs() override;
    void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(ThreeTop, 0);

private:
    void ApplyDataScaleFactors();
    void printStuff();
    float getLeadPt(size_t idx = 0);
    void setSubChannel();
    bool isSameSign(Level level);

    bool baseline_cuts();
    bool signal_cuts();
    bool nonprompt_cuts();
    bool charge_misid_cuts();
    bool ttz_CR_cuts();

    void applyNonprompt(Particle& part, PID pid);

    float muPt(size_t idx) { return muon.size(Level::Tight) > idx ? muon.pt(Level::Tight, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Tight) > idx ? elec.pt(Level::Tight, idx) : -1; }

    ResolvedTop rTop;

    ParticleOut* o_tightMuons;
    ParticleOut* o_tightElectrons;
    JetOut* o_jets;
    JetOut* o_bJets;
    TopOut* o_top;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality;
    std::vector<Int_t> o_nb_loose, o_nb_tight, o_nloose_muon, o_nloose_elec;
    std::vector<Bool_t> o_pass_zveto;
    // TrigEff trigEff_leadPt;
};

#endif
