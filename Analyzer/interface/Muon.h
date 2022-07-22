#ifndef __MUON_H_
#define __MUON_H_

#include "analysis_suite/Analyzer/interface/Lepton.h"

class Muon : public Lepton {
public:
    void setup(TTreeReader& fReader, bool isMC);
    void createLooseList() override;
    void createTightList() override;
    void createIsolatedList() override;
    float getScaleFactor() override;

private:
    NTupleArray<Bool_t> isGlobal, isTracker;
    NTupleArray<UInt_t> nMatches, bestTrackType;
    NTupleArray<Bool_t> isPF, highPtId;

    WeightHolder muon_scale;

    Float_t ptMax = 119;
    Float_t ptMin = 20;
};

#endif // __MUON_H_
