#ifndef BEFFICIENCY_H_
#define BEFFICIENCY_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

class BEfficiency : public BaseSelector {
 public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void ApplyScaleFactors() override;
    void clearOutputs() override;
    ClassDefOverride(BEfficiency, 0);

private:
    bool isSameSign();
    bool signal_cuts();

    float muPt(size_t idx) { return muon.size(Level::Fake) > idx ? muon.pt(Level::Fake, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Fake) > idx ? elec.pt(Level::Fake, idx) : -1; }

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Bool_t> o_hlt_dilepton_ht, o_hlt_dilepton, o_hlt_trilepton;
    std::vector<Bool_t> o_pass_zveto_loose, o_pass_zveto_fake;
    BEffOut* o_beff;
};

#endif // BEFFICIENCY_H_
