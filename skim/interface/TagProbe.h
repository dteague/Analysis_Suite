#ifndef TAGPROBE_H_
#define TAGPROBE_H_

#include "analysis_suite/skim/interface/BaseSelector.h"
#include "analysis_suite/skim/interface/Output.h"
#include "analysis_suite/skim/interface/TrigObj.h"

#include <set>

class TagProbe : public BaseSelector {
public:
    virtual void Init(TTree* tree) override;
    virtual bool getCutFlow() override;
    virtual void FillValues(const Bitmap& event_bitmap) override;
    virtual void SetupOutTreeBranches(TTree* tree) override;
    virtual void setOtherGoodParticles(size_t syst) override;
    virtual void ApplyScaleFactors() override;
    virtual void clearParticles() override;
    virtual void clearOutputs() override;
    ClassDefOverride(TagProbe, 0);

private:
    void setSubChannel();
    bool closure_cuts();
    bool isOppositeSign(Level level);
    float getLeadPt();
    float get_mass(size_t t, size_t p, bool useRaw);

    std::unordered_map<size_t, size_t> muon_pair, elec_pair;

    TrigObj trigobj;

    TRVariable<Float_t> Pileup_nTrueInt;
    TRVariable<Bool_t> trig_mm, trig_ee;

    Float_t tag_pt, tag_abseta, tag_iso, probe_pt, probe_fakept, probe_abseta, probe_iso, probe_mva;
    Float_t mass_fake, mass, wgt;
    Int_t  njets;
    bool isMuon, tag_tight, probe_tight, pass_mm, pass_ee;

    Subchannel file_chan;
    TTree* final_tree;

};





#endif // TAGPROBE_H_
