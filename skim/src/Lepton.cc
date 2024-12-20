#include "analysis_suite/skim/interface/Lepton.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"
#include "analysis_suite/skim/interface/Jet.h"

#include <limits>

void Lepton::setup(std::string name, TTreeReader& fReader)
{
    m_charge.setup(fReader, name + "_charge");
    dz.setup(fReader, name + "_dz");
    dxy.setup(fReader, name + "_dxy");
    mvaTTH.setup(fReader, name + "_mvaTTH");
    sip3d.setup(fReader, name+"_sip3d");
    if (isMC) {
        genPartIdx.setup(fReader, name+"_genPartIdx");
        genPartFlav.setup(fReader, name+"_genPartFlav");
    }
    // Isolation variables
    jetIdx.setup(fReader, name + "_jetIdx");
    ptRel.setup(fReader, name + "_jetPtRelv2");
    ptRatio_.setup(fReader, name + "_jetRelIso");
    iso.setup(fReader, name + "_miniPFRelIso_all");

    isoCut = 0.1;
    ptRatioCut = 0.4;
    mvaCut = 0.65;
    // mvaCut = -2;

    GenericParticle::setup(name, fReader);
    setup_map(Level::Loose);
    setup_map(Level::Fake);
    setup_map(Level::Tight);
}

float Lepton::ptRatio2(size_t i, Jet& jets) {
    if (jetIdx.at(i) < 0 && closeJet_by_lepton.find(i) == closeJet_by_lepton.end()) {
        return ptRatio(i);
    } else if (jetIdx.at(i) < 0) {
        return ptRatio(i)/jets.jec_shift(closeJet_by_lepton[i]);
    } else {
        return ptRatio(i)/jets.jec_shift(jetIdx.at(i));
    }
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

void Lepton::fillLepton_Iso(LeptonOut_Fake& output, Jet& jet, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            output.rawPt.push_back(m_pt.at(idx));
            output.ptRatio.push_back(ptRatio(idx));
            output.ptRatio2.push_back(ptRatio2(idx, jet));
            output.ptRatio3.push_back(ptRatio3(idx, jet));
            output.mvaTTH.push_back(tthMVA_(idx));
            output.iso.push_back(iso.at(idx));
            if (jetIdx.at(idx) != -1) {
                output.jet_btag.push_back(jet.btag.at(jetIdx.at(idx)));
            } else {
                output.jet_btag.push_back(-1.);
            }
            if (isMC) {
                output.jet_flav.push_back(static_cast<Int_t>(genPartFlav.at(idx)));
            }
        }
    }
}

std::pair<size_t, float> Lepton::getCloseJet(size_t lidx, const Particle& jet)
{
    size_t minIdx = SIZE_MAX;
    float minDr = 0.16; // 0.4^2
    float lphi = phi(lidx);
    float leta = eta(lidx);
    for (size_t jidx = 0; jidx < jet.size(); jidx++) {
        float dr2 = deltaR2(jet.eta(jidx), leta, jet.phi(jidx), lphi);
        if (minDr > dr2) {
            minIdx = jidx;
            minDr = dr2;
        }
    }

    if (minIdx != SIZE_MAX) {
        closeJet_by_lepton[lidx] = minIdx;
    }
    return std::make_pair(minIdx, minDr);
}

float Lepton::massInRange(Level level, float low, float high)
{
    auto lep_list = list(level);
    for (auto i = lep_list.begin(); i != lep_list.end(); ++i) {
        LorentzVector lep = p4(*i);
        for (auto j = i; j != lep_list.end(); ++j) {
            if (*i == *j || charge(*i) * charge(*j) > 0) {
                continue;
            }
            float mass = (p4(*j) + lep).M();
            if (mass >= low && mass < high)
                return mass;
        }
    }
    return -1.;
}

bool Lepton::isInMassRange(Level level, float low, float high)
{
    auto lep_list = list(level);
    for (auto i = lep_list.begin(); i != lep_list.end(); ++i) {
        LorentzVector lep = p4(*i);
        for (auto j = i; j != lep_list.end(); ++j) {
            if (*i == *j || charge(*i) * charge(*j) > 0) {
                continue;
            }
            float mass = (p4(*j) + lep).M();
            if (mass >= low && mass < high)
                return true;
        }
     }
    return false;
}

bool Lepton::passJetIsolation(size_t idx) const
{
     return mvaTTH.at(idx) > mvaCut;
}

float Lepton::getFakePtFactor(size_t idx) const
{
    if (passJetIsolation(idx)) {
        return 1.;
    } else {
        return cone_correction/ptRatio(idx);
    }
}

void Lepton::fillFlippedCharge(GenParticle& gen)
{
    for (size_t i = 0; i < size(); i++) {
        int pdg = gen.pdgId.at(genPartIdx.at(i));
        // remember, pdg is postive for electron, negative for positron
        flips.push_back(abs(pdg) == static_cast<int>(id) && pdg*charge(i) > 0);
    }
}
