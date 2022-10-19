#ifndef FAKERATE_H_
#define FAKERATE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

#include <set>

class FakeRate : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void ApplyScaleFactors() override;
    void ApplyDataSpecifics();
    virtual void clearParticles() override;
    virtual void clearOutputs() override;
    virtual void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(FakeRate, 0);

private:
    void setSubChannel();
    void set_leadlep();
    bool measurement_cuts();
    bool sideband_cuts();
    bool single_lep_cuts(CutInfo& cuts);

    LeptonOut_Fake* o_fakeMuons;
    LeptonOut_Fake* o_tightMuons;
    LeptonOut_Fake* o_fakeElectrons;
    LeptonOut_Fake* o_tightElectrons;
    ParticleOut* o_muon;
    JetOut* o_jets;
    JetOut* o_bJets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi;
    std::vector<size_t> o_nb_loose, o_nb_tight;
    std::vector<int> o_lo_hlt, o_hi_hlt;

    LorentzVector lead_lep;
    float lead_rawpt;
};


#endif // FAKERATE_H_
