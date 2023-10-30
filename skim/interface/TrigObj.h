#ifndef TRIGOBJ_H_
#define TRIGOBJ_H_

#include "analysis_suite/skim/interface/Particle.h"
#include "analysis_suite/skim/interface/GenParticle.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

#include <unordered_map>
#include <set>

class Electron;
class Muon;

class TrigObj : public Particle {
public:
    void setup(TTreeReader& fReader);


    void setupGoodLists(Electron& elec, Muon& muon)
    {
        match_muon(muon);
        match_electron(elec);
    }

    virtual void clear() override
    {
        Particle::clear();
        muon_match.clear();
        elec_match.clear();
    }

    void match_muon(Muon& muons);
    void match_electron(Electron& elecs);

    bool is_match() { return muon_match.size() > 0 || elec_match.size() > 0; }

    TRArray<Int_t> id;
    TRArray<Int_t> filterBits;
    TRVariable<Bool_t> hlt_isomu, hlt_mu50, hlt_dimu, hlt_el;
    float isomu_cut, el_cut;

    std::set<size_t> muon_match, elec_match;
};

#endif // TRIGOBJ_H_
