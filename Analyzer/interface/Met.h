#ifndef MET_H_
#define MET_H_

#include "analysis_suite/Analyzer/interface/Systematic.h"
#include "analysis_suite/Analyzer/interface/Variable.h"

#include <TTreeReader.h>

enum class MET_Type {
    PUPPI, PF
};

class Jet;

class Met : SystematicWeights {
public:
    float pt() { return *corr_pt; }
    float phi() { return *corr_phi; }

    void setup(MET_Type type, TTreeReader& fReader);
    void setupJEC(Jet& jet);
    void fix_xy(UInt_t run, int nVertices);

private:
    TRVariable<float> m_pt;
    TRVariable<float> m_phi;
    std::unordered_map<MET_Type, std::string> met_name = {
        {MET_Type::PUPPI, "PuppiMET"},
        {MET_Type::PF, "type1_pfMET"},
    };

    float *corr_pt, *corr_phi;
    std::unordered_map<Systematic, std::unordered_map<eVar, float>> m_corr_pt, m_corr_phi;
    WeightHolder xcorr, ycorr;
    std::string name;
};


#endif // MET_H_
