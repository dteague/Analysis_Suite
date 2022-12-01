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

    LeptonOut_Fake *o_looseMuons, *o_fakeMuons, *o_tightMuons;
    LeptonOut_Fake *o_looseElectrons, *o_fakeElectrons, *o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_met, o_metphi;

    LorentzVector lead_lep;
};


#endif // FAKERATE_H_
