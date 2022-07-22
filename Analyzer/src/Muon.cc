#include "analysis_suite/Analyzer/interface/Muon.h"
#include "analysis_suite/Analyzer/interface/Jet.h"

void Muon::setup(TTreeReader& fReader, bool isMC)
{
    Lepton::setup("m", fReader, isMC);
    // isGlobal.setup(fReader, "Muon_isGlobal");
    // isTracker.setup(fReader, "Muon_isTracker");
    // isPFcand.setup(fReader, "Muon_isPFcand");
    // iso.setup(fReader, "Muon_miniPFRelIso_all");
    // tightCharge.setup(fReader, "Muon_tightCharge");
    // mediumId.setup(fReader, "Muon_mediumId");
    // sip3d.setup(fReader, "Muon_sip3d");
    isGlobal.setup(fReader, "m", "IsGlobal");
    isTracker.setup(fReader, "m", "IsTracker");
    nMatches.setup(fReader, "m", "NoOfMatches");
    bestTrackType.setup(fReader, "m", "BestTrackType");

    isPF.setup(fReader, "m", "IsPFMuon");
    highPtId.setup(fReader, "m", "HighPtID");


    id = PID::Muon;
    isoCut = 0.1;
    ptRatioCut = 0.8;
    ptRelCut = 5.0;
    // if (year_ == Year::yr2016pre || year_ == Year::yr2016post) {
    //     isoCut = 0.16;
    //     ptRatioCut = 0.76;
    //     ptRelCut = 7.2;
    // } else {
    //     isoCut = 0.11;
    //     ptRatioCut = 0.74;
    //     ptRelCut = 6.8;
    // }

    // if (isMC) {
    //     auto corr_set = getScaleFile("MUO", "muon_Z");
    //     muon_scale = WeightHolder(corr_set->at("NUM_MediumID_DEN_TrackerMuons"), Systematic::Muon_Scale,
    //                               {"sf", "systup", "systdown"});
    // }
}

void Muon::createLooseList()
{
    for (size_t i = 0; i < size(); ++i) {
        if (pt(i) > 5
            && fabs(eta(i) < 2.4)
            && fabs(dxy.at(i)) < 0.5
            && fabs(dz.at(i) < 1)
            && sips.at(i) < 4
            && (isGlobal.at(i) || (isTracker.at(i) && nMatches.at(i) > 0))
            && bestTrackType.at(i) != 2) {
            m_partList[Level::Loose]->push_back(i);
        }

    }
}

void Muon::createTightList()
{
    for (auto i : list(Level::Loose)) {
        if ((year_ == Year::yr2018 && (true)) // 2018 cuts
            ||
            (isPF.at(i) || (pt(i) > 200 && highPtId.at(i)))) // 2016/2017 cuts
            m_partList[Level::Tight]->push_back(i);
    }
}

void Muon::createIsolatedList()
{
    for (auto i : list(Level::Tight)) {
        if (zzTightId.at(i)
            && zzIso.at(i))
            m_partList[Level::Iso]->push_back(i);
    }
}

float Muon::getScaleFactor()
{
    float weight = 1.;
    // std::string syst = systName(muon_scale);
    // for (auto midx : list(Level::Tight)) {
    //     weight *= muon_scale.evaluate({yearMap.at(year_)+"_UL", fabs(eta(midx)), pt(midx), syst});
    // }
    return weight;
}
