#ifndef TWOLEPTON_H_
#define TWOLEPTON_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

#include <set>

class TwoLepton : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void ApplyScaleFactors() override;
    virtual void clearParticles() override;
    virtual void clearOutputs() override;
    virtual void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(TwoLepton, 0);

private:
    void setSubChannel();
    bool measurement_cuts();
    float getLeadPt();

    LeptonOut_Fake *o_fakeMuons, *o_tightMuons;
    LeptonOut_Fake *o_fakeElectrons, *o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_met, o_metphi;
    std::vector<size_t> o_nloose_mu, o_nloose_el;

};

#endif // TWOLEPTON_H_
