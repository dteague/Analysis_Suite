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
    ClassDefOverride(BEfficiency, 0);

 private:
    bool isSameSign();
    bool signal_cuts();

    TRVariable<Float_t> Pileup_nTrueInt;

    BEffOut* o_beff;
};

#endif // BEFFICIENCY_H_
