#ifndef ZZANALYSIS_H
#define ZZANALYSIS_H

#include "analysis_suite/Analyzer/interface/BaseSelector.h"
#include "analysis_suite/Analyzer/interface/Output.h"

class ZZAnalysis : public BaseSelector {
public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const std::vector<bool>& passVec) override;
    void ApplyScaleFactors() override;
    void clearOutputs() override;
    void SetupOutTreeBranches(TTree* tree) override;
    ClassDefOverride(ZZAnalysis, 0);

private:
    void printStuff();
    void setSubChannel();

    bool signal_cuts();

    float getLeadPt();
    float getSubLeadPt();

    ParticleOut* o_looseMuons;
    ParticleOut* o_looseElectrons;
    ParticleOut* o_tightMuons;
    ParticleOut* o_tightElectrons;
    std::vector<float> o_met, o_metphi;

    TRVariable<Float_t> Pileup_nTrueInt;
};

#endif
