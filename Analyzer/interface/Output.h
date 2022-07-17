#ifndef __OUTPUT_H_
#define __OUTPUT_H_

#include "analysis_suite/Analyzer/interface/Jet.h"
#include "analysis_suite/Analyzer/interface/Lepton.h"
#include "analysis_suite/Analyzer/interface/Particle.h"

struct ParticleOut {
    std::vector<Float_t> pt;
    std::vector<Float_t> eta;
    std::vector<Float_t> phi;
    std::vector<Float_t> mass;
    std::vector<Int_t> syst_bitMap;
    void clear()
    {
        pt.clear();
        eta.clear();
        phi.clear();
        mass.clear();
        syst_bitMap.clear();
    }
};

struct JetOut {
    std::vector<Float_t> pt;
    std::vector<Float_t> eta;
    std::vector<Float_t> phi;
    std::vector<Float_t> mass;
    std::vector<Int_t> syst_bitMap;
    std::vector<std::pair<Float_t, Float_t>> jer;
    std::vector<std::pair<Float_t, Float_t>> jes;
    std::vector<Float_t> discriminator;
    void clear() {
        pt.clear();
        eta.clear();
        phi.clear();
        mass.clear();
        syst_bitMap.clear();
        jer.clear();
        jes.clear();
        discriminator.clear();
    }
};

struct LeptonOut_Fake {
    std::vector<Float_t> pt;
    std::vector<Float_t> eta;
    std::vector<Float_t> phi;
    std::vector<Float_t> mass;
    std::vector<Float_t> ptRatio;
    std::vector<Float_t> ptRel;
    std::vector<Int_t> syst_bitMap;
    void clear()
    {
        pt.clear();
        eta.clear();
        phi.clear();
        mass.clear();
        ptRatio.clear();
        ptRel.clear();
        syst_bitMap.clear();
    }
};


struct BEffOut {
    std::vector<Float_t> pt;
    std::vector<Float_t> eta;
    std::vector<Float_t> phi;
    std::vector<Float_t> mass;
    std::vector<Int_t> flavor;
    std::vector<Int_t> pass_btag;
    std::vector<Int_t> syst_bitMap;
    void clear()
    {
        pt.clear();
        eta.clear();
        phi.clear();
        mass.clear();
        flavor.clear();
        pass_btag.clear();
        syst_bitMap.clear();
    }
};

template <class T>
size_t fillParticle(const Particle& particle, Level level, T& fillObject, size_t idx, size_t pass_bitmap);

template <class T>
void fillParticle(const Particle& particle, Level level, T& fillObject, size_t pass_bitmap);

void fillJet(const Jet& jet, Level level, JetOut& fillObject, size_t pass_bitmap);

void fillAllLeptons(const Lepton& muon, const Lepton& elec, ParticleOut& fillObject, size_t pass_bitmap);

void fillLepton_Fake(const Lepton& lep, Level level, LeptonOut_Fake& fillObject, size_t pass_bitmap);

void fillBEff(const Jet& jet, Level level, BEffOut& fillObject, size_t pass_bitmap);

#include "analysis_suite/Analyzer/interface/Output.hxx"

#endif // __OUTPUT_H_
