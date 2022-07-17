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

void Muon::createTightList(Particle& jets)
{
    for (size_t i = 0; i < size(); ++i) {
        if (pt(i) > 15)
            m_partList[Level::Tight]->push_back(i);
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
