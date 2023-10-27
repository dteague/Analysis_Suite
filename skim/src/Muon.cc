#include "analysis_suite/skim/interface/Muon.h"
#include "analysis_suite/skim/interface/Jet.h"

enum ROC_ERROR {
    eDefault = 0,
    eStat = 1,
    eZpt = 2,
    eEwk = 3,
    eDeltaM = 4,
    eEwk2 = 5,
};

void Muon::setup(TTreeReader& fReader)
{
    Lepton::setup("Muon", fReader);
    isGlobal.setup(fReader, "Muon_isGlobal");
    isTracker.setup(fReader, "Muon_isTracker");
    isPFcand.setup(fReader, "Muon_isPFcand");
    tightCharge.setup(fReader, "Muon_tightCharge");
    mediumId.setup(fReader, "Muon_mediumId");
    nTrackerLayers.setup(fReader, "Muon_nTrackerLayers");
    tid.setup(fReader, "Muon_tightId");

    id = PID::Muon;

    if (isMC) {
        auto corr_set = getScaleFile("MUO", "muon_Z");
        muon_scale = WeightHolder(corr_set->at("NUM_MediumID_DEN_genTracks"), Systematic::Muon_Scale,
                                  {"sf", "systup", "systdown"});
        muon_tthMVA = WeightHolder(corr_set->at("NUM_MediumID_DEN_genTracks"), Systematic::Muon_tthMVA,
                                   {"sf", "systup", "systdown"});
    }
    std::string roccor_file = scaleDir_+"/POG/USER/"+yearMap.at(year_)+"_UL/RoccoR.txt";
    roc_corr.init(roccor_file);

    // if (year_ == Year::yr2016pre)        cone_correction = 0.750;
    // else if (year_ == Year::yr2016post)  cone_correction = 0.750;
    // else if (year_ == Year::yr2017)      cone_correction = 0.700;
    // else if (year_ == Year::yr2018)      cone_correction = 0.725;
    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) cone_correction = 0.775;
    else if (year_ == Year::yr2017)                            cone_correction = 0.725;
    else if (year_ == Year::yr2018)                            cone_correction = 0.75;
}

void Muon::createLooseList()
{
    for (size_t i = 0; i < size(); i++) {
        if (pt(i) > 5
            && fabs(eta(i)) < 2.4
            && (isGlobal.at(i) || isTracker.at(i))
            && isPFcand.at(i)
            && iso.at(i) < 0.4
            && fabs(dz.at(i)) < 0.1
            && fabs(dxy.at(i)) < 0.05)
            m_partList[Level::Loose]->push_back(i);
    }
}

void Muon::createFakeList(Particle& jets)
{
    for (auto i : list(Level::Loose)) {
        if (m_pt.at(i) > 12 // Trigger pt threshold
            && tightCharge.at(i) == 2
            && mediumId.at(i)
            && sip3d.at(i) < 4
            && (ptRatio(i) > ptRatioCut || passJetIsolation(i))
            && m_pt.at(i)*getFakePtFactor(i) > 15 // Total pt threshold
            )
            {
                m_partList[Level::Fake]->push_back(i);
                dynamic_cast<Jet&>(jets).closeJetDr_by_index.insert(getCloseJet(i, jets));
            }
    }
}

void Muon::createTightList(Particle& jets)
{
    for (auto i : list(Level::Fake)) {
        if (pt(i) > 15
            && passJetIsolation(i)
            ) {
            m_partList[Level::Tight]->push_back(i);
        } else {
            fakePtFactor[i] = getFakePtFactor(i);
        }

    } 
}

float Muon::getScaleFactor()
{
    float weight = 1.;
    std::string syst = systName(muon_scale);
    for (auto midx : list(Level::Fake)) {
        float pt = (m_pt.at(midx) > 15) ? m_pt.at(midx) : 15.;
        weight *= muon_scale.evaluate({yearMap.at(year_)+"_UL", fabs(eta(midx)), pt, syst});
    }
    return weight;
}

float Muon::getRocCorrection(GenParticle& gen)
{
    float weight = 1.;
    for (auto i : list(Level::Fake)) {
        if (isMC) {
            if (genPartIdx.at(i) != -1) { //(recommended), MC scale and resolution correction when matched gen muon is available
                float genPt = gen.pt(genPartIdx.at(i));
                weight *= roc_corr.kSpreadMC(charge(i), pt(i), eta(i), phi(i), genPt, eDefault, 0);
            } else { //MC scale and extra smearing when matched gen muon is not available
                float rand_uniform = static_cast<float>(rand())/static_cast <float>(RAND_MAX);
                weight *= roc_corr.kSmearMC(charge(i), pt(i), eta(i), phi(i), nTrackerLayers.at(i), rand_uniform, eDefault, 0);
            }
        } else { //data
            weight *= roc_corr.kScaleDT(charge(i), pt(i), eta(i), phi(i), eDefault, 0);
        }
    }
    return weight;
}
