#ifndef LEPTON_MISID_CLOSURE_H_
#define LEPTON_MISID_CLOSURE_H_

#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

class Charge_MisId : public DileptonBase {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void clearOutputs() override;
    ClassDefOverride(Charge_MisId, 0);

private:
    void printStuff();
    float get_mass();

    bool measurement_cuts();
    bool closure_cuts();

    LeptonOut* o_tightMuons;
    LeptonOut* o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality;
};

#endif // LEPTON_MISID_CLOSURE_H_
