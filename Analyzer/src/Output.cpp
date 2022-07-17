#include "analysis_suite/Analyzer/interface/Output.h"
#include"analysis_suite/Analyzer/interface/logging.h"

void fillJet(const Jet& jet, Level level, JetOut& fillObject, size_t pass_bitmap)
{
    LOG_FUNC << "Start fillJet";
    fillObject.clear();

    for (size_t idx = 0; idx < jet.size(); ++idx) {
        size_t final_bitmap = fillParticle(jet, level, fillObject, idx, pass_bitmap);
        if (final_bitmap != 0) {
            fillObject.discriminator.push_back(jet.btag.at(idx));
            fillObject.jer.push_back(jet.get_JEC_pair(Systematic::Jet_JER, idx));
            fillObject.jes.push_back(jet.get_JEC_pair(Systematic::Jet_JES, idx));
        }
    }
    LOG_FUNC << "End of fillJet";
}

void fillAllLeptons(const Lepton& muon, const Lepton& elec, ParticleOut& fillObject, size_t pass_bitmap)
{
    LOG_FUNC << "Start of fillLeptons";
    fillObject.clear();
    size_t midx = 0;
    size_t eidx = 0;
    while (midx != muon.size() || eidx != elec.size()) {
        if (midx != muon.size() && (eidx == elec.size() || muon.pt(midx) > elec.pt(eidx))) {
            fillParticle(muon, Level::Tight, fillObject, midx, pass_bitmap);
            midx++;
        } else {
            fillParticle(elec, Level::Tight, fillObject, eidx, pass_bitmap);
            eidx++;
        }
    }
    LOG_FUNC << "End of fillLeptons";
}

void fillLepton_Fake(const Lepton& lep, Level level, LeptonOut_Fake& fillObject, size_t pass_bitmap)
{
    LOG_FUNC << "Start fillLepton";
    fillObject.clear();

    for (size_t idx = 0; idx < lep.size(); ++idx) {
        size_t final_bitmap = fillParticle(lep, level, fillObject, idx, pass_bitmap);
        if (final_bitmap != 0) {
            fillObject.ptRatio.push_back(lep.ptRatio.at(idx));
            fillObject.ptRel.push_back(lep.ptRel.at(idx));
        }
    }
    LOG_FUNC << "End of fillLepon";
}

void fillBEff(const Jet& jet, Level level, BEffOut& fillObject, size_t pass_bitmap) {
    LOG_FUNC << "Start fillJet";
    fillObject.clear();

    for (size_t idx = 0; idx < jet.size(); ++idx) {
        size_t final_bitmap = fillParticle(jet, level, fillObject, idx, pass_bitmap);
        if (final_bitmap != 0) {
            float btag = jet.btag.at(idx);
            if (btag > jet.tight_bjet_cut) fillObject.pass_btag.push_back(3);
            else if (btag > jet.medium_bjet_cut) fillObject.pass_btag.push_back(2);
            else if (btag > jet.loose_bjet_cut) fillObject.pass_btag.push_back(1);
            else fillObject.pass_btag.push_back(0);
            fillObject.flavor.push_back(jet.hadronFlavour.at(idx));
        }
    }
    LOG_FUNC << "End of fillJet";

}
