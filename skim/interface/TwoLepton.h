#ifndef TWOLEPTON_H_
#define TWOLEPTON_H_

#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

#include <set>

class TwoLepton : public DileptonBase {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void ApplyScaleFactors() override;
    ClassDefOverride(TwoLepton, 0);


private:
    bool measurement_cuts();

    LeptonOut_Fake *o_fakeMuons, *o_tightMuons;
    LeptonOut_Fake *o_fakeElectrons, *o_tightElectrons;
    JetOut* o_jets;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_met, o_metphi;
    std::vector<size_t> o_nloose_mu, o_nloose_el;

};

#endif // TWOLEPTON_H_
