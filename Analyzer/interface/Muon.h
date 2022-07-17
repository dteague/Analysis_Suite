#ifndef __MUON_H_
#define __MUON_H_

#include "analysis_suite/Analyzer/interface/Lepton.h"

class Muon : public Lepton {
public:
    void setup(TTreeReader& fReader, bool isMC);
    virtual void createTightList(Particle& jets) override;
    virtual float getScaleFactor() override;

    TRArray<Bool_t> isGlobal;
    TRArray<Bool_t> isTracker;
    TRArray<Bool_t> isPFcand;
    TRArray<Int_t> tightCharge;
    TRArray<Bool_t> mediumId;
    TRArray<Float_t> sip3d;

    WeightHolder muon_scale;

    Float_t ptMax = 119;
    Float_t ptMin = 20;
};

#endif // __MUON_H_
