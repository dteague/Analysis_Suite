#include "analysis_suite/Analyzer/interface/Jet.h"

void Jet::setup(TTreeReader& fReader, bool isMC)
{
    m_pt.setup(fReader, "jetPt");
    m_eta.setup(fReader, "jetEta");
    m_phi.setup(fReader, "jetPhi");


    isMC_ = isMC;
    // jetId.setup(fReader, "Jet_jetId");
    btag.setup(fReader, "jetDeepCSV");
    puId.setup(fReader, "jetPUID");
    if (isMC) {
        hadronFlavour.setup(fReader, "jetHadronFlavor");
        // genJetIdx.setup(fReader, "Jet_genJetIdx");
        // rawFactor.setup(fReader, "Jet_rawFactor");
        rho.setup(fReader, "pvRho");
    }

    // setup_map(Level::Loose);
    // setup_map(Level::Bottom);
    setup_map(Level::Tight);

    // if (year_ == Year::yr2016pre) {
    //     loose_bjet_cut =  0.0508;
    //     medium_bjet_cut = 0.2598;
    //     tight_bjet_cut =  0.6502;
    // } else if (year_ == Year::yr2016post) {
    //     loose_bjet_cut =  0.0480;
    //     medium_bjet_cut = 0.2489;
    //     tight_bjet_cut =  0.6377;
    // } else if (year_ == Year::yr2017) {
    //     loose_bjet_cut =  0.0532;
    //     medium_bjet_cut = 0.3040;
    //     tight_bjet_cut =  0.7476;
    // } else if (year_ == Year::yr2018) {
    //     loose_bjet_cut =  0.0490;
    //     medium_bjet_cut = 0.2783;
    //     tight_bjet_cut =  0.7100;
    // }

    // use_shape_btag = true;
    if (isMC) {
    //     // JEC Weights
    //     auto corr_set = getScaleFile("JME", "jet_jerc");
    //     jec_scale = WeightHolder(corr_set->at(jec_source[year_]+"_Total_AK4PFchs"));
    //     jet_resolution = WeightHolder(corr_set->at(jer_source[year_]+"_PtResolution_AK4PFchs"));
    //     jer_scale = WeightHolder(corr_set->at(jer_source[year_]+"_ScaleFactor_AK4PFchs"),
    //                              Systematic::Jet_JER, {"nom","up","down"});

        // Pileup Weights
        auto jmar_set = getScaleFile("JME", "jmar");
        puid_scale = WeightHolder(jmar_set->at("PUJetID_eff"),
                                  Systematic::Jet_PUID, {"nom","up","down"});

        //     // BTagging Weights
    //     auto btag_set = getScaleFile("BTV", "btagging");
    //     btag_bc_scale = WeightHolder(btag_set->at("deepJet_comb"),
    //                                  Systematic::Jet_PUID, {"central","up","down"});
    //     btag_udsg_scale = WeightHolder(btag_set->at("deepJet_incl"),
    //                                    Systematic::Jet_PUID, {"central","up","down"});

    //     // // BTagging Efficiencies
    //     // auto beff_set = getScaleFile("BTV", "tagging_eff");
    //     // btag_eff = WeightHolder(btag_set->at("SS"),
    //     //                         Systematic::Jet_PUID, {"central","up","down"});

    }

    m_jet_scales[Systematic::Nominal] = {{eVar::Nominal, std::vector<float>()}};
    for (auto syst : jec_systs) {
        m_jet_scales[syst] = {
            {eVar::Up, std::vector<float>()},
            {eVar::Down, std::vector<float>()}
        };
    }
}

void Jet::createTightList()
{
    for (size_t i = 0; i < size(); i++) {
        if (pt(i) > 25)
            m_partList[Level::Tight]->push_back(i);
    }
}

float Jet::getHT(const std::vector<size_t>& jet_list)
{
    float ht = 0;
    for (auto i : jet_list) {
        ht += pt(i);
    }
    return ht;
}

float Jet::getScaleFactor()
{
    float weight = 1.;
    // std::string syst = systName(puid_scale);
    // for (auto idx : list(Level::Tight)) {
    //     if (pt(idx) < 50)
    //         weight *= puid_scale.evaluate({eta(idx), pt(idx), syst, "M"});
    // }
    return weight;
}

float Jet::getCentrality(const std::vector<size_t>& jet_list)
{
    float etot = 0;
    for (auto i : jet_list) {
        etot += p4(i).E();
    }
    return getHT(jet_list) / etot;
}

void Jet::setupJEC(GenericParticle& genJet) {
    if (currentSyst == Systematic::Nominal || isJECSyst()) {
        m_jec = &m_jet_scales[currentSyst][currentVar];
        m_jec->assign(size(), 1);
        if (!isMC_) return;
        // for(size_t i = 0; i < size(); ++i) {
        //     (*m_jec)[i] *= get_jec(i);
        //     (*m_jec)[i] *= get_jer(i, genJet);
        // }
    } else {
        m_jec = &m_jet_scales[Systematic::Nominal][eVar::Nominal];
    }
}

float Jet::getTotalBTagWeight() {
    float weight = 1;
    std::string syst = systName(puid_scale);
    const auto& goodBJets = list(Level::Bottom);
    for (auto bidx : goodBJets) {
        auto scaler = (hadronFlavour.at(bidx) == 0) ? btag_udsg_scale : btag_bc_scale;
        weight *= scaler.evaluate({syst, "M", hadronFlavour.at(bidx), fabs(eta(bidx)), pt(bidx)});
    }
    for (auto jidx : list(Level::Tight)) {
        if (std::find(goodBJets.begin(), goodBJets.end(), jidx) != goodBJets.end()) {
            continue; // is a bjet, weighting already taken care of
        }
        auto scaler = (hadronFlavour.at(jidx) == 0) ? btag_udsg_scale : btag_bc_scale;
        float bSF = scaler.evaluate({syst, "M", hadronFlavour.at(jidx), fabs(eta(jidx)), pt(jidx)});
        float eff = btag_eff.evaluate({syst, "M", hadronFlavour.at(jidx), fabs(eta(jidx)), pt(jidx)});
        weight *= (1 - bSF * eff) / (1 - eff);
    }
    return weight;
}

float Jet::get_jec(size_t i) {
    if (currentSyst != Systematic::Jet_JER || currentVar == eVar::Nominal) {
        return 1.;
    } else {
        float delta = jec_scale.evaluate({eta(i), pt(i)});
        return (currentVar == eVar::Up) ? (1+delta) : (1-delta);
    }
}

float Jet::get_jer(size_t i, GenericParticle& genJets) {
    using namespace ROOT::Math::VectorUtil;
    float resolution = jet_resolution.evaluate({eta(i), pt(i), *rho});
    float scale = jer_scale.evaluate({eta(i), systName(jer_scale)});

    bool hasGenJet = genJetIdx.at(i) != -1;
    auto genJet = (hasGenJet) ? genJets.p4(genJetIdx.at(i)) : LorentzVector();
    float pt_ratio = (hasGenJet) ? 1-pt(i)/genJet.Pt() : 0.;

    if (hasGenJet && DeltaR(p4(i), genJet) < jet_dr/2 && abs(pt_ratio) < 3*resolution) {
        return 1 + (scale-1)*pt_ratio;
    } else if (scale > 1.) {
        std::normal_distribution<> gaussian{0, resolution};
        return 1 + gaussian(gen)*sqrt(pow(scale,2) - 1);
    }

    return 1.;
}

PolarVector Jet::get_momentum_change()
{
    PolarVector change;
    for(size_t i = 0; i < size(); ++i) {
        change += PolarVector(pt(i)-nompt(i), phi(i));
    }


    return change;
}
