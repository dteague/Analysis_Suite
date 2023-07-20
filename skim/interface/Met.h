#ifndef MET_H_
#define MET_H_

#include "analysis_suite/skim/interface/Systematic.h"
#include "analysis_suite/skim/interface/Variable.h"

#include <TTreeReader.h>

enum class MET_Type {
    PUPPI, PF
};

class Jet;

typedef std::pair<TRVariable<float>, TRVariable<float>> Puppi_pair;

class Met : SystematicWeights {
public:
    float pt() { return *corr_pt; }
    float phi() { return *corr_phi; }

    void setup(TTreeReader& fReader, MET_Type type=MET_Type::PF);
    void setSyst();
    void setupMet(Jet& jet, UInt_t run, int nVertices);
    void fix_xy(UInt_t run, int nVertices);

    float mt(float pt, float phi);
    float pt_unfix, phi_unfix;
private:
    TRVariable<float> m_pt;
    TRVariable<float> m_phi;
    std::unordered_map<MET_Type, std::string> met_name = {
        {MET_Type::PUPPI, "PuppiMET"},
        {MET_Type::PF, "MET"},
    };

    float *corr_pt, *corr_phi;
    std::unordered_map<Systematic, std::unordered_map<eVar, float>> m_corr_pt, m_corr_phi;
    std::string name;
    bool ispuppi;

    std::unordered_map<Systematic, Puppi_pair> m_puppi_pt_syst, m_puppi_phi_syst;
};


#endif // MET_H_
