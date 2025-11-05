#ifndef TRIGGER_EFF_H_
#define TRIGGER_EFF_H_

#include "analysis_suite/skim/interface/DileptonBase.h"
#include "analysis_suite/skim/interface/Output.h"

class Trigger_Eff : public DileptonBase {
 public:
    void Init(TTree* tree) override;
    bool getCutFlow() override;
    void FillValues(const Bitmap& event_bitmap) override;
    void SetupOutTreeBranches(TTree* tree) override;
    void ApplyScaleFactors() override;
    void clearOutputs() override;
    ClassDefOverride(Trigger_Eff, 0);

private:
    void setSubChannel();
    bool passMetTrigger();
    bool passDileptonTrigger();

    bool signal_cuts();

    float muPt(size_t idx) { return muon.size(Level::Fake) > idx ? muon.pt(Level::Fake, idx) : -1; }
    float elPt(size_t idx) { return elec.size(Level::Fake) > idx ? elec.pt(Level::Fake, idx) : -1; }

    ParticleOut* o_tightMuons;
    ParticleOut* o_tightElectrons;

    TRVariable<Float_t> Pileup_nTrueInt;

    std::vector<Float_t>  o_met, o_metphi, o_ht, o_mht;
    std::vector<Bool_t> o_dilepton_trigger, o_met_trigger;

};


#endif // TRIGGER_EFF_H_
