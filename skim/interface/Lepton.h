#ifndef __LEPTON_H_
#define __LEPTON_H_

#include "analysis_suite/skim/interface/Particle.h"
#include "analysis_suite/skim/interface/GenParticle.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

#include <unordered_map>

class Jet;

class Lepton : public Particle {
public:
    virtual void createLooseList(){};
    virtual void createFakeList(Particle& jets){};
    virtual void createTightList(Particle& jets){};
    float massInRange(Level level, float low=ZMASS-ZWINDOW, float high=ZMASS+ZWINDOW);
    bool isInMassRange(Level level, float low=ZMASS-ZWINDOW, float high=ZMASS+ZWINDOW);
    void setup(std::string name, TTreeReader& fReader);
    std::pair<size_t, float> getCloseJet(size_t lidx, const Particle& jet);
    virtual bool passJetIsolation(size_t idx) const;
    void fillFlippedCharge(GenParticle& gen);
    float getFakePtFactor(size_t idx) const;

    float pt_(size_t idx) const override { return m_pt.at(idx)*fakePtFactor.at(idx); }
    float rawpt(size_t i) const { return m_pt.at(i); }
    float rawpt(Level level, size_t i) const { return m_pt.at(idx(level, i)); }

    Int_t charge(size_t idx) { return m_charge.at(idx); };
    Int_t charge(Level level, size_t i) { return charge(idx(level, i)); };

    float getMT(size_t idx, float met, float met_phi) const
    {
        return sqrt(2*pt(idx)*met*(1-cos(phi(idx) - met_phi)));
    }
    float getMT(Level level, size_t i, float met, float met_phi) const { return getMT(idx(level, i), met, met_phi); }

    float ptRatio(size_t i) const { return 1/(1+ptRatio_.at(i)); }
    float ptRatio2(size_t i, Jet& jets);
    float ptRatio3(size_t i, Particle& jets) const {
        if (jetIdx.at(i) < 0) {
            return 1/(1+ptRatio_.at(i));
        } else {
            return m_pt.at(i)/jets.pt(jetIdx.at(i));
        }
    }
    virtual float tthMVA_(size_t i) { return mvaTTH.at(i); }
    virtual void setupGoodLists(Particle& jets, GenParticle& gen) override
    {
        fakePtFactor.assign(size(), 1.);
        createLooseList();
        createFakeList(jets);
        createTightList(jets);
        fillFlippedCharge(gen);
    }

    virtual void clear() override
    {
        Particle::clear();
        closeJet_by_lepton.clear();
        flips.clear();
        fakePtFactor.clear();
    }

    std::unordered_map<size_t, size_t> closeJet_by_lepton;
    std::vector<bool> flips;
    std::vector<float> fakePtFactor;

    float isoCut, ptRatioCut, ptRelCut, mvaCut;
    float cone_correction;

    TRArray<Float_t> mvaTTH;
    TRArray<Float_t> ptRel;
    TRArray<Float_t> ptRatio_;
    TRArray<Float_t> iso;
    TRArray<Bool_t> tid;
    TRArray<Int_t> jetIdx;

    void fillLepton(LeptonOut& output, Level level, const Bitmap& event_bitmap);
    void fillLepton_Iso(LeptonOut_Fake& output, Jet& jet, Level level, const Bitmap& event_bitmap);

    TRArray<Float_t> dz;
    TRArray<Float_t> dxy;
    TRArray<Float_t> sip3d;

protected:
    TRArray<Int_t> m_charge;

    TRArray<Int_t> genPartIdx;
    TRArray<UChar_t> genPartFlav;

    PID id;

};

#endif // __LEPTON_H_
