#ifndef DY_TEST_H_
#define DY_TEST_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"

class DY_test : public BaseSelector {
 public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void clearParticles() override;
    void clearOutputs() override;
    virtual void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(DY_test, 0);

 private:
    void printStuff();
    bool isSameSign();
    bool closure_cuts();
    void setSubChannel();
    float getLeadPt();
    float get_mass();

    LeptonOut_Fake *o_tightMuons, *o_fakeMuons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality;
    std::vector<int> o_nlooseMu, o_nlooseEl;
};


#endif // DY_TEST_H_
