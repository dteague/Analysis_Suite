#include "analysis_suite/Analyzer/interface/Electron.h"
#include "analysis_suite/Analyzer/interface/Jet.h"

void Electron::setup(TTreeReader& fReader, bool isMC)
{
    Lepton::setup("e", fReader, isMC);
    // eCorr.setup(fReader, "Electron_eCorr");
    // lostHits.setup(fReader, "Electron_lostHits");
    // convVeto.setup(fReader, "Electron_convVeto");
    // sieie.setup(fReader, "Electron_sieie");
    // hoe.setup(fReader, "Electron_hoe");
    // eInvMinusPInv.setup(fReader, "Electron_eInvMinusPInv");
    // sip3d.setup(fReader, "Electron_sip3d");
    // tightCharge.setup(fReader, "Electron_tightCharge");
    // ecalSumEt.setup(fReader, "Electron_dr03EcalRecHitSumEt");
    // hcalSumEt.setup(fReader, "Electron_dr03HcalDepth1TowerSumEt");
    // tkSumPt.setup(fReader, "Electron_dr03TkSumPt");
    id = PID::Electron;
    // mva_l.setup(fReader, "Electron_mvaFall17V2noIso_WPL");
    // mva_80.setup(fReader, "Electron_mvaFall17V2noIso_WP80");
    // mva_90.setup(fReader, "Electron_mvaFall17V2noIso_WP90");

    ptRatioCut = 0.8;
    ptRelCut = 5.0;
    // ptRatioCut = 0.78;
    // ptRelCut = 8.0;
    isoCut = 0.1;

    // if(isMC) {
    //     auto corr_set = getScaleFile("EGM", "electron");
    //     electron_scale = WeightHolder(corr_set->at("UL-Electron-ID-SF"), Systematic::Electron_Scale,
    //                                   {"sf", "sfup", "sfdown"});
    // }

}

void Electron::createLooseList()
{
    for (size_t i = 0; i < size(); ++i) {
        if (pt(i) > 7
            && fabs(eta(i) < 2.5)
            && fabs(dxy.at(i)) < 0.5
            && fabs(dz.at(i) < 1)
            && sips.at(i) < 4)
            m_partList[Level::Loose]->push_back(i);
    }
}


void Electron::createIsolatedList()
{
    for (auto i : list(Level::Loose)) {
        if (zzTightId.at(i)
            && zzIso.at(i)) {
            m_partList[Level::Tight]->push_back(i);
            m_partList[Level::Iso]->push_back(i);
        }

    }
}

float Electron::getScaleFactor()
{
    float weight = 1.;
    // std::string syst = systName(electron_scale);
    // for (auto eidx : list(Level::Tight)) {
    //     weight *= electron_scale.evaluate({yearMap.at(year_), syst, "wp90noiso", fabs(eta(eidx)), pt(eidx)});
    // }
    return weight;
}

bool Electron::passTriggerRequirements(size_t i)
{
    return ecalSumEt.at(i) / pt(i) < 0.5
        && hcalSumEt.at(i) / pt(i) < 0.3
        && tkSumPt.at(i) / pt(i) < 0.2
        && ((fabs(eta(i)) < BARREL_ETA && sieie.at(i) < 0.013) || (fabs(eta(i)) >= BARREL_ETA && sieie.at(i) < 0.035))
        && hoe.at(i) < 0.13;

}
