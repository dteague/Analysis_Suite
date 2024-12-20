#ifndef BSCALE_H_
#define BSCALE_H_

#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

class BScale : public DileptonBase {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void clearOutputs() override;
    ClassDefOverride(BScale, 0);

private:
    bool baseline_cuts();
    bool signal_cuts();

    float muPt(size_t idx) { return muon.size(Level::Fake) > idx ? muon.pt(Level::Fake, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Fake) > idx ? elec.pt(Level::Fake, idx) : -1; }

    TRVariable<Float_t> Pileup_nTrueInt;

    JetOut* o_jets;
    std::vector<Float_t> o_wgt_nobtag;

};


#endif // BSCALE_H_
