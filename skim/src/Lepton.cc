#include "analysis_suite/skim/interface/Lepton.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

#include <limits>

// #define NOMVATTH

void Lepton::setup(std::string name, TTreeReader& fReader)
{
    m_charge.setup(fReader, name + "_charge");
    dz.setup(fReader, name + "_dz");
    dxy.setup(fReader, name + "_dxy");
    mvaTTH.setup(fReader, name + "_mvaTTH");
    sip3d.setup(fReader, name+"_sip3d");
    if (isMC) {
        genPartIdx.setup(fReader, name+"_genPartIdx");
    }

    // Isolation variables
    ptRel.setup(fReader, name + "_jetPtRelv2");
    ptRatio_.setup(fReader, name + "_jetRelIso");
    iso.setup(fReader, name + "_miniPFRelIso_all");

    isoCut = 0.1;
    // ptRatioCut = 0.75;
    // ptRelCut = 8.;
    mvaCut = 0.4;

    GenericParticle::setup(name, fReader);
    setup_map(Level::Loose);
    setup_map(Level::Fake);
    setup_map(Level::Tight);
}

void Lepton::fillLepton(LeptonOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            output.flip.push_back(flips.at(idx));
        }
    }

}

void Lepton::fillLepton_Iso(LeptonOut_Fake& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            output.rawPt.push_back(m_pt.at(idx));
            output.ptRatio.push_back(ptRatio(idx));
            output.ptRel.push_back(ptRel.at(idx));
            output.mvaTTH.push_back(mvaTTH.at(idx));
            output.iso.push_back(iso.at(idx));
        }
    }
}

std::pair<size_t, float> Lepton::getCloseJet(size_t lidx, const Particle& jet)
{
    size_t minIdx = SIZE_MAX;
    float minDr = 100;
    float lphi = phi(lidx);
    float leta = eta(lidx);
    for (size_t jidx = 0; jidx < jet.size(); jidx++) {
        float dr2 = deltaR(jet.eta(jidx), leta, jet.phi(jidx), lphi);
        if (minDr > dr2) {
            minIdx = jidx;
            minDr = dr2;
        }
    }

    if (minIdx != SIZE_MAX)
        closeJet_by_lepton[lidx] = minIdx;

    return std::make_pair(minIdx, minDr);
}

bool Lepton::passZVeto()
{
    for (auto tidx : list(Level::Loose)) {
        LorentzVector tlep = p4(tidx);
        for (auto lidx : list(Level::Loose)) {
            if (tidx >= lidx || charge(tidx) * charge(lidx) > 0)
                continue;
            float mass = (p4(lidx) + tlep).M();
            if (mass < LOW_ENERGY_CUT || (fabs(mass - ZMASS) < ZWINDOW))
                return false;
        }
    }
    return true;
}

bool Lepton::passZCut(float low, float high)
{
    for (auto tidx : list(Level::Fake)) { //tightList
        LorentzVector tlep = p4(tidx);
        for (auto lidx : list(Level::Fake)) {
            if (tidx >= lidx || charge(tidx) * charge(lidx) > 0)
                continue;
            float mass = (p4(lidx) + tlep).M();
            if (mass > low && mass < high)
                return true;
        }
    }
    return false;
}

bool Lepton::passJetIsolation(size_t idx) const
{
#ifdef NOMVATTH
    return iso.at(idx) < isoCut && ( ptRatio(idx) > ptRatioCut || ptRel.at(idx) > ptRelCut );
#else
    return mvaTTH.at(idx) > mvaCut;
#endif
}

float Lepton::getFakePtFactor(size_t idx) const
{
#ifdef NOMVATTH
    if (ptRel.at(idx) > ptRelCut) {
        if (iso.at(idx) > isoCut)
            return (1 + iso.at(idx)-isoCut);
    } else if (ptRatio(idx) <= ptRatioCut) {
        return ptRatioCut/ptRatio(idx);
    }
    return 1.;
#else
    if (passJetIsolation(idx)) {
        return 1.;
    } else {
        // return 1.;
        return cone_correction/ptRatio(idx);
        // return 1./ptRatio(idx);
    }
#endif
}

void Lepton::fillFlippedCharge(GenParticle& gen)
{
    for (size_t i = 0; i < size(); i++) {
        int pdg = gen.pdgId.at(genPartIdx.at(i));
        // remember, pdg is postive for electron, negative for positron
        flips.push_back(abs(pdg) == static_cast<int>(id) && pdg*charge(i) > 0);
    }
}
