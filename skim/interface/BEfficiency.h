#ifndef BEFFICIENCY_H_
#define BEFFICIENCY_H_

#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

class BEfficiency : public DileptonBase {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void ApplyScaleFactors() override;
    void clearOutputs() override;
    ClassDefOverride(BEfficiency, 0);

private:
    bool signal_cuts();

    TRVariable<Float_t> Pileup_nTrueInt;
    TRVariable<Bool_t> hlt_mu, hlt_el, hlt_mm, hlt_me, hlt_em, hlt_ee, hlt_ee_high;

    std::vector<Bool_t> o_hlt_mu, o_hlt_el, o_hlt_mm, o_hlt_me, o_hlt_em, o_hlt_ee, o_hlt_ee_high;


    ElectronOut_Endcap* o_elec;
    ParticleOut* o_muon;
    JetOut* o_jets;
};

#endif // BEFFICIENCY_H_
