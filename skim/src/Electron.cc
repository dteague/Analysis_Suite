#include "analysis_suite/skim/interface/Electron.h"
#include "analysis_suite/skim/interface/Jet.h"

void Electron::setup(TTreeReader& fReader)
{
    Lepton::setup("Electron", fReader);
    lostHits.setup(fReader, "Electron_lostHits");
    convVeto.setup(fReader, "Electron_convVeto");
    sieie.setup(fReader, "Electron_sieie");
    hoe.setup(fReader, "Electron_hoe");
    eInvMinusPInv.setup(fReader, "Electron_eInvMinusPInv");
    tightCharge.setup(fReader, "Electron_tightCharge");
    sc_delta.setup(fReader, "Electron_deltaEtaSC");

    ecalSumEt.setup(fReader, "Electron_dr03EcalRecHitSumEt");
    hcalSumEt.setup(fReader, "Electron_dr03HcalDepth1TowerSumEt");
    tkSumPt.setup(fReader, "Electron_dr03TkSumPt");
    id = PID::Electron;
    mva_l.setup(fReader, "Electron_mvaFall17V2noIso_WPL");
    mva_80.setup(fReader, "Electron_mvaFall17V2noIso_WP80");
    mva_90.setup(fReader, "Electron_mvaFall17V2noIso_WP90");
    tid.setup(fReader, "Electron_mvaFall17V2noIso_WP80");

    // TTH MVA variables
    iso_chg.setup(fReader, "Electron_miniPFRelIso_chg");
    jetNDauCharged.setup(fReader, "Electron_jetNDauCharged");
    mva_vals.setup(fReader, "Electron_mvaFall17V2noIso");

    if(isMC) {
        auto corr_set = getScaleFile("EGM", "electron");
        electron_scale = WeightHolder(corr_set->at("UL-Electron-ID-SF"), Systematic::Electron_Scale,
                                      {"sf", "sfup", "sfdown"});
        // genPartFlav.setup(fReader, "Electron_genPartFlav");
        auto tth_set = getScaleFile("USER", "tthMVA_scales");
        electron_tthMVA = WeightHolder(tth_set->at("NUM_mvaTTH_DEN_isElectron"), Systematic::Electron_tthMVA,
                                       {"value", "systup", "systdown"});
    }

    std::string mva_xmlfile = scaleDir_+"/POG/USER/"+yearMap.at(year_)+"_UL/electron_tth_weights.xml";
    tth_mva.setup(mva_xmlfile, this);
    // if (year_ == Year::yr2016pre)        cone_correction = 0.825;
    // else if (year_ == Year::yr2016post)  cone_correction = 0.825;
    // else if (year_ == Year::yr2017)      cone_correction = 0.775;
    // else if (year_ == Year::yr2018)      cone_correction = 0.800;
    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) cone_correction = 0.850;
    else if (year_ == Year::yr2017)                            cone_correction = 0.800;
    else if (year_ == Year::yr2018)                            cone_correction = 0.800;


}

void Electron::fillElectron_Endcap(ElectronOut_Endcap& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            output.dxy.push_back(fabs(dxy.at(idx)));
            output.dz.push_back(fabs(dz.at(idx)));
            output.convVeto.push_back(convVeto.at(idx));
            output.lostHits.push_back(lostHits.at(idx));
            output.sip3d.push_back(sip3d.at(idx));
            output.mvaTTH.push_back(tth_mva_vec[idx]);
            output.truth.push_back(genPartFlav.at(idx));
        }
    }
}

void Electron::fillElectron_small(LeptonOut_small& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        setSyst(0); // Fill values with nominal values and corrections done on top
        Bitmap final_bitmap = bitmap(level).at(idx) & event_bitmap;
        if (final_bitmap.any()) {
            output.syst_bitMap.push_back(final_bitmap.to_ulong());
            output.rawPt.push_back(m_pt.at(idx));
            output.ptRatio.push_back(ptRatio(idx));
            output.mvaTTH.push_back(tth_mva_vec[idx]);
        } else {
            continue;
        }
    }
}

bool Electron::inCrack(size_t i)
{
    float abseta = fabs(eta(i)+sc_delta.at(i));
    return abseta > 1.444 && abseta < 1.566;
}

void Electron::createLooseList()
{
    tth_mva_vec.clear();
    tth_mva_vec.resize(size());
    for (size_t i = 0; i < size(); i++) {
        if (pt(i) > 10
            && fabs(eta(i)) < 2.5
            && convVeto.at(i)
            && !inCrack(i)
            && lostHits.at(i) <= 1
            && fabs(dz.at(i)) < 0.1
            && fabs(dxy.at(i)) < 0.05
            && mva_l.at(i)
            && iso.at(i) < 0.4)
            {
                m_partList[Level::Loose]->push_back(i);
            }
    }
}

/// closejet w iso
void Electron::createFakeList(Particle& jets)
{
    for (auto i : list(Level::Loose)) {
        // std::cout << tth_mva.get_mva(i, dynamic_cast<Jet&>(jets)) << " " << mvaTTH.at(i) << std::endl;
        tth_mva_vec[i] = tth_mva.get_mva(i, dynamic_cast<Jet&>(jets));
        if (m_pt.at(i) > 15 // Trigger pt threshold
            && mva_90.at(i)
            && sip3d.at(i) < 4
            && lostHits.at(i) == 0
            && tightCharge.at(i) == 2
            && passTriggerRequirements(i) // Nonmva
            && m_pt.at(i)*getFakePtFactor(i) > 20
            && (ptRatio(i) > ptRatioCut || passJetIsolation(i)) // MVA
            )
            {
                m_partList[Level::Fake]->push_back(i);
                dynamic_cast<Jet&>(jets).closeJetDr_by_index.insert(getCloseJet(i, jets));
            }
    }
}

// need iso
void Electron::createTightList(Particle& jets)
{
    for (auto i : list(Level::Fake)) {
        if (pt(i) > 20 &&
            passJetIsolation(i))
            {
                m_partList[Level::Tight]->push_back(i);
        } else {
            fakePtFactor[i] = getFakePtFactor(i);
        }
    }
}

float Electron::getScaleFactor()
{
    float weight = 1.;
    std::string syst = systName(electron_scale);
    std::string syst_mva = systName(electron_tthMVA);
    for (auto eidx : list(Level::Fake)) {
        float pt = (m_pt.at(eidx) >= 20) ? m_pt.at(eidx) : 20.;
        weight *= electron_scale.evaluate({yearMap.at(year_), syst, "wp90noiso", eta(eidx), m_pt.at(eidx)});
        weight *= electron_tthMVA.evaluate({fabs(eta(eidx)), pt, syst_mva});
    }
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

bool Electron::passJetIsolation(size_t idx) const
{
    return tth_mva_vec.at(idx) > mvaCut;
}
