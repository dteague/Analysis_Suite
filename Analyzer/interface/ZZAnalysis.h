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
    void setSubChannel() {};

    bool signal_cuts();

    ParticleOut* o_looseMuons;
    ParticleOut* o_looseElectrons;
    LeptonOut* o_tightMuons;
    LeptonOut* o_tightElectrons;
    ParticleOut* o_tightLeptons;
    JetOut* o_jets;
    JetOut* o_bJets;
    TopOut* o_resolvedTop;

    TRVariable<Float_t> Pileup_nTrueInt;
};

#endif
