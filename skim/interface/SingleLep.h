#ifndef SINGLELEP_TRIGGER_H_
#define SINGLELEP_TRIGGER_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

#include <set>

class SingleLep : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void ApplyScaleFactors() override;
    virtual void clearParticles() override;
    virtual void clearOutputs() override;
    virtual void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(SingleLep, 0);

 private:
    void setSubChannel();
    bool measurement_cuts();

    LeptonOut_small *o_fakeElectrons, *o_fakeMuons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;
    TRVariable<Bool_t> hlt_lo_mu, hlt_hi_mu, hlt_lo_el, hlt_hi_el;

    std::vector<Float_t> o_ht, o_met, o_metphi;
    Bool_t o_hlt_lo_mu, o_hlt_hi_mu, o_hlt_lo_el, o_hlt_hi_el;

};


#endif // SINGLELEP_TRIGGER_H_
