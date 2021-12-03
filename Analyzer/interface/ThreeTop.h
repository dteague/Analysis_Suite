#ifndef THREETOP_H
#define THREETOP_H

#include "analysis_suite/Analyzer/interface/BaseSelector.h"
#include "analysis_suite/Analyzer/interface/Output.h"
#include "analysis_suite/Analyzer/interface/ResolvedTop.h"

class ThreeTop : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow(cut_info& cuts) override;
    virtual bool getTriggerCut(cut_info& cuts) override;
    virtual void fillTriggerEff(bool passCuts, bool passTrigger) override;
    virtual void FillValues(const std::vector<bool>& passVec) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void setupChannel() override;
    virtual void ApplyScaleFactors() override;
    virtual void clearParticles() override;
    virtual void clearOutputs() override;
    virtual void setOtherGoodParticles(size_t syst) override;
    ClassDefOverride(ThreeTop, 0);

private:
    void printStuff();
    float getLeadPt(size_t idx = 0);
    void setSubChannel();
    bool isSameSign();

    TTree* treeFakeRate_;

    ResolvedTop rTop;

    ParticleOut* o_looseMuons;
    ParticleOut* o_tightMuons;
    ParticleOut* o_looseElectrons;
    ParticleOut* o_tightElectrons;
    ParticleOut* o_tightLeptons;
    JetOut* o_jets;
    BJetOut* o_bJets;
    TopOut* o_resolvedTop;

    TTRValue<ULong64_t>* event;
    TTRValue<Bool_t>* Flag_goodVertices;
    TTRValue<Bool_t>* Flag_globalSuperTightHalo2016Filter;
    TTRValue<Bool_t>* Flag_HBHENoiseFilter;
    TTRValue<Bool_t>* Flag_HBHENoiseIsoFilter;
    TTRValue<Bool_t>* Flag_EcalDeadCellTriggerPrimitiveFilter;
    TTRValue<Bool_t>* Flag_BadPFMuonFilter;
    TTRValue<Bool_t>* Flag_ecalBadCalibFilter;
    TTRValue<Float_t>* Met_pt;
    TTRValue<Float_t>* Met_phi;
    TTRValue<Float_t>* Pileup_nTrueInt;

    std::vector<Float_t> o_ht, o_htb, o_met, o_metphi, o_centrality;

    TH2F *passTrigger_leadPt, *failTrigger_leadPt;
};

#endif
