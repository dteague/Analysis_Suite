#ifndef LEPTON_MISID_CLOSURE_H_
#define LEPTON_MISID_CLOSURE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

class Closure_MisId : public BaseSelector {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void clearParticles() override;
    void clearOutputs() override;
    ClassDefOverride(Closure_MisId, 0);

private:
    void printStuff();
    bool isSameSign();
    bool measurement_cuts();
    bool closure_cuts();
    void setSubChannel();
    float getLeadPt();
    float get_mass();
    bool getTriggerValue();

    LeptonOut* o_tightMuons;
    LeptonOut* o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality;
};

#endif // LEPTON_MISID_CLOSURE_H_
