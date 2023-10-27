#ifndef NONPROMPT_CLOSURE_H_
#define NONPROMPT_CLOSURE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

#include <set>

class Nonprompt_Closure : public BaseSelector {
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
    ClassDefOverride(Nonprompt_Closure, 0);

private:
    void setSubChannel();
    bool closure_cuts();
    bool dy_closure_cuts();
    float getLeadPt();
    bool isSameSign();
    bool getTriggerValue();

    LeptonOut_Fake *o_looseMuons, *o_fakeMuons, *o_tightMuons;
    LeptonOut_Fake* o_looseElectrons, *o_fakeElectrons, *o_tightElectrons;
    JetOut* o_jets;
    JetOut* o_bJets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi;
    std::vector<size_t> o_nb_loose, o_nb_tight;

};


#endif // NONPROMPT_CLOSURE_H_
