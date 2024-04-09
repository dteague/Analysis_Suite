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
    bool measurement_cuts();
    bool sideband_cuts();
    bool single_lep_cuts(CutInfo& cuts);
    size_t get_trigger();
    bool apply_trigger();
    LeptonOut_Fake *o_fakeMuons, *o_tightMuons;
    ElectronOut_Fake *o_fakeElectrons, *o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    Dataset highpt_e_dataset;

    std::vector<Float_t> o_ht, o_met, o_metphi;
    std::vector<Float_t> bjet_scales;
    Int_t hlt_loPt_prescale, hlt_hiPt_prescale;
};


#endif // FAKERATE_H_
