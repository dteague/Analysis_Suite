#include "analysis_suite/Analyzer/interface/Lepton.h"
#include "analysis_suite/Analyzer/interface/CommonFuncs.h"

#include <limits>

void Lepton::setup(std::string name, TTreeReader& fReader, bool isMC)
{
    m_pt.setup(fReader, name, "Pt");
    m_eta.setup(fReader, name, "Eta");
    m_phi.setup(fReader, name, "Phi");
    m_mass.setup(fReader, name, "Mass");

    m_charge.setup(fReader, name, "Charge");
    dz.setup(fReader, name, "PVDZ");
    dxy.setup(fReader, name, "PVDXY");
    sips.setup(fReader, name, "SIP3D");

    zzTightId.setup(fReader, name, "ZZTightIDNoVtx");
    zzIso.setup(fReader, name, "ZZIsoPass");
    // if (isMC) {
    //     genPartIdx.setup(fReader, name+"_genPartIdx");
    // }

    // // Isolation variables
    // ptRel.setup(fReader, name + "_jetPtRelv2");
    // ptRatio.setup(fReader, name + "_jetRelIso");
    // iso.setup(fReader, name + "_miniPFRelIso_all");


    setup_map(Level::Loose);
    setup_map(Level::Tight);
    setup_map(Level::Iso);
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
    for (auto tidx : list(Level::Loose)) { //tightList
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

// bool Lepton::passZCut(float low, floa,t high)
// {
//     for (auto tidx : list(Level::Fake)) { //tightList
//         LorentzVector tlep = p4(tidx);
//         for (auto lidx : list(Level::Fake)) {
//             if (tidx >= lidx || charge(tidx) * charge(lidx) > 0)
//                 continue;
//             float mass = (p4(lidx) + tlep).M();
//             if (mass > low && mass < high)
//                 return true;
//         }
//     }
//     return false;
// }

bool Lepton::passJetIsolation(size_t idx, const Particle& jets)
{
    // if (closeJet_by_lepton.find(idx) == closeJet_by_lepton.end())
    //     return true; /// no close jet (probably no jets)
    return iso.at(idx) < isoCut && ( 1/(1+ptRatio.at(idx)) > ptRatioCut || ptRel.at(idx) > ptRelCut );
}
