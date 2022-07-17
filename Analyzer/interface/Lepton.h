#ifndef __LEPTON_H_
#define __LEPTON_H_

#include "analysis_suite/Analyzer/interface/Particle.h"
#include "analysis_suite/Analyzer/interface/CommonFuncs.h"

#include <unordered_map>

class Lepton : public Particle {
public:
    size_t _size() const override { return m_pt.size(); }
    float _pt(size_t idx) const override { return m_pt.at(idx); }
    float _eta(size_t idx) const override { return m_eta.at(idx); }
    float _phi(size_t idx) const override { return m_phi.at(idx); }
    float _mass(size_t idx) const override { return m_mass.at(idx); }

    virtual void createTightList(Particle& jets){};
    bool passZVeto();
    bool passZCut(float low, float high);
    void setup(std::string name, TTreeReader& fReader, bool isMC);
    std::pair<size_t, float> getCloseJet(size_t lidx, const Particle& jet);
    bool passJetIsolation(size_t idx, const Particle& jets);

    Int_t charge(size_t idx) { return m_charge.at(idx); };
    Int_t charge(Level level, size_t i) { return charge(idx(level, i)); };

    float getMT(size_t idx, float met, float met_phi) const
    {
        return sqrt(2*met*(1-cos(phi(idx) - met_phi)));
    }
    float getMT(Level level, size_t i, float met, float met_phi) const { return getMT(idx(level, i), met, met_phi); }

    virtual void setupGoodLists(Particle& jets) override
    {
        createTightList(jets);
    }

    virtual void clear() override
    {
        Particle::clear();
        closeJet_by_lepton.clear();
    }

    std::unordered_map<size_t, size_t> closeJet_by_lepton;

    float isoCut, ptRatioCut, ptRelCut;

    TRArray<Float_t> ptRel;
    TRArray<Float_t> ptRatio;
 
protected:
    NTupleArray<Float_t> m_pt, m_eta, m_phi, m_mass;
    NTupleArray<Int_t> m_charge;
    TRArray<Float_t> iso;
    NTupleArray<Float_t> dz;
    NTupleArray<Float_t> dxy;
    NTupleArray<Int_t> genPartIdx;

    const float ZMASS = 91.188;
    const float ZWINDOW = 15;
    const float LOW_ENERGY_CUT = 12;

    PID id;

};

#endif // __LEPTON_H_
