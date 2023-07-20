#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void Jet::setup(TTreeReader& fReader)
{
    GenericParticle::setup("Jet", fReader);
    jetId.setup(fReader, "Jet_jetId");
    btag.setup(fReader, "Jet_btagDeepFlavB");
    puId.setup(fReader, "Jet_puId");
    if (isMC) {
        hadronFlavour.setup(fReader, "Jet_hadronFlavour");
        genJetIdx.setup(fReader, "Jet_genJetIdx");
        rawFactor.setup(fReader, "Jet_rawFactor");
        rho.setup(fReader, "fixedGridRhoFastjetAll");
    }
    neHEF.setup(fReader, "Jet_neHEF");
    neEmEF.setup(fReader, "Jet_neEmEF");
    muEF.setup(fReader, "Jet_muEF");
    chHEF.setup(fReader, "Jet_chHEF");
    chEmEF.setup(fReader, "Jet_chEmEF");
    nConstituents.setup(fReader, "Jet_nConstituents");

    setup_map(Level::Loose_NoPUID);
    setup_map(Level::Loose);
    setup_map(Level::Bottom);
    setup_map(Level::Tight);

    if (year_ == Year::yr2016pre) {
        loose_bjet_cut =  0.0508;
        medium_bjet_cut = 0.2598;
        tight_bjet_cut =  0.6502;
    } else if (year_ == Year::yr2016post) {
        loose_bjet_cut =  0.0480;
        medium_bjet_cut = 0.2489;
        tight_bjet_cut =  0.6377;
    } else if (year_ == Year::yr2017) {
        loose_bjet_cut =  0.0532;
        medium_bjet_cut = 0.3040;
        tight_bjet_cut =  0.7476;
    } else if (year_ == Year::yr2018) {
        loose_bjet_cut =  0.0490;
        medium_bjet_cut = 0.2783;
        tight_bjet_cut =  0.7100;
    }
    bjet_cuts["L"] = loose_bjet_cut;
    bjet_cuts["M"] = medium_bjet_cut;
    bjet_cuts["T"] = tight_bjet_cut;


    // use_shape_btag = true;
    if (isMC) {
        // JEC Weights
        auto corr_set = getScaleFile("JME", "jet_jerc");
        jes_scale = WeightHolder(corr_set->at(jes_source[year_]+"_Total_AK4PFchs"));
        jet_resolution = WeightHolder(corr_set->at(jer_source[year_]+"_PtResolution_AK4PFchs"));
        jer_scale = WeightHolder(corr_set->at(jer_source[year_]+"_ScaleFactor_AK4PFchs"),
                                 Systematic::Jet_JER, {"nom","up","down"});

        // Pileup Weights
        auto jmar_set = getScaleFile("JME", "jmar");
        puid_scale = WeightHolder(jmar_set->at("PUJetID_eff"),
                                  Systematic::Jet_PUID, {"nom","up","down"});

        // BTagging Weights
        auto btag_set = getScaleFile("BTV", "btagging");
        btag_bc_scale = WeightHolder(btag_set->at("deepJet_comb"),
                                     Systematic::BJet_BTagging, {"central","up","down"});
        btag_udsg_scale = WeightHolder(btag_set->at("deepJet_incl"),
                                       Systematic::BJet_BTagging, {"central","up","down"});

        // BTagging Efficiencies
        try {
            auto beff_set = getScaleFile("USER", "beff");
            btag_eff = WeightHolder(beff_set->at("SS"),
                                    Systematic::BJet_Eff, {"central","up","down"});
        } catch (...) {
            LOG_WARN << "BTagging Efficiencies not found for this year. May not be necessary, will continue";
        }
    }

    m_jet_scales[Systematic::Nominal] = {{eVar::Nominal, std::vector<float>()}};
    for (auto syst : jec_systs) {
        m_jet_scales[syst] = {
            {eVar::Up, std::vector<float>()},
            {eVar::Down, std::vector<float>()}
        };
    }
    for(auto var: all_vars) {
        jer[var] = std::vector<float>();
        jes[var] = std::vector<float>();
    }
}

void Jet::fillJet(JetOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            if (btag.at(idx) > tight_bjet_cut) output.pass_btag.push_back(3);
            else if (btag.at(idx) > medium_bjet_cut) output.pass_btag.push_back(2);
            else if (btag.at(idx) > loose_bjet_cut) output.pass_btag.push_back(1);
            else output.pass_btag.push_back(0);
            output.discriminator.push_back(btag.at(idx));
            output.jer.push_back(get_JEC_pair(Systematic::Jet_JER, idx));
            output.jes.push_back(get_JEC_pair(Systematic::Jet_JES, idx));
        }
    }
}

void Jet::fillJetEff(BEffOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            if (btag.at(idx) > tight_bjet_cut) output.pass_btag.push_back(3);
            else if (btag.at(idx) > medium_bjet_cut) output.pass_btag.push_back(2);
            else if (btag.at(idx) > loose_bjet_cut) output.pass_btag.push_back(1);
            else output.pass_btag.push_back(0);
            output.flavor.push_back(hadronFlavour.at(idx));
        }
    }
}

void Jet::createLooseList()
{
    for (size_t i = 0; i < size(); i++) {
        if (pt(i) > 25
            && fabs(eta(i)) < 2.4
            && (jetId.at(i) & tightId) != 0
            && (neHEF.at(i) < 0.9 && neEmEF.at(i) < 0.9 && muEF.at(i) < 0.8 && chEmEF.at(i) < 0.8 && chHEF.at(i) > 0. && nConstituents.at(i) > 1)
            && (closeJetDr_by_index.find(i) == closeJetDr_by_index.end() || closeJetDr_by_index.at(i) >= pow(0.4, 2))) {
            if (pt(i) > 50 || (puId.at(i) >> PU_Medium) & 1) {
                m_partList[Level::Loose]->push_back(i);
            } else {
                m_partList[Level::Loose_NoPUID]->push_back(i);
            }
        }
    }
}

void Jet::createBJetList()
{
    for (auto i : list(Level::Loose)) {
        if (btag.at(i) > medium_bjet_cut)
            m_partList[Level::Bottom]->push_back(i);
        n_loose_bjet.back() += (btag.at(i) > loose_bjet_cut) ? 1 : 0;
        n_medium_bjet.back() += (btag.at(i) > medium_bjet_cut) ? 1 : 0;
        n_tight_bjet.back() += (btag.at(i) > tight_bjet_cut) ? 1 : 0;
    }
}

void Jet::createTightList()
{
    for (auto i : list(Level::Loose)) {
        if (pt(i) > 40)
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
    std::string syst = systName(puid_scale);
    for (auto idx : list(Level::Loose)) {
        if (pt(idx) < 50)
            weight *= puid_scale.evaluate({eta(idx), pt(idx), syst, "M"});
    }
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

void Jet::setSyst(size_t syst)
{
    Particle::setSyst(syst);
    if (currentSyst == Systematic::Nominal || isJECSyst()) {
        m_jec = &m_jet_scales[currentSyst][currentVar];
    } else {
        m_jec = &m_jet_scales[Systematic::Nominal][eVar::Nominal];
    }
}

void Jet::setupJEC(GenericParticle& genJet)
{
    for (auto& [syst, var_scales] : m_jet_scales ) {
        for (auto& [var, scales] : var_scales) {
            scales.assign(size(), 1);
        }
    }
    if (!isMC)
        return;

    for(size_t i = 0; i < size(); ++i) {
        setup_jes(i);
        setup_jer(i, genJet);
        m_jet_scales[Systematic::Nominal][eVar::Nominal][i] = jer[eVar::Nominal].at(i)*jes[eVar::Nominal].at(i);
        for (auto var: syst_vars) {
            m_jet_scales[Systematic::Jet_JES][var][i] = jer[eVar::Nominal].at(i)*jes[var].at(i);
            m_jet_scales[Systematic::Jet_JER][var][i] = jer[var].at(i)*jes[eVar::Nominal].at(i);
        }
    }
}

float Jet::getPileupIDWeight()
{
    float weight = 1;
    std::string syst = systName(puid_scale);
    for (auto idx : list(Level::Loose)) {
        if (pt(idx) < 50)
            weight *= puid_scale.evaluate({eta(idx), pt(idx), syst, "M"});
    }
    for (auto idx : list(Level::Loose_NoPUID)) {
        float eff = puid_scale.evaluate({eta(idx), pt(idx), "MCEff", "M"});
        float sf = puid_scale.evaluate({eta(idx), pt(idx), syst, "M"});
        weight *= (1 - sf*eff) / (1 -eff);
    }
    return weight;
}

float Jet::getBJetWeight(size_t idx, std::string lvl)
{
    auto scaler = (hadronFlavour.at(idx) == 0) ? btag_udsg_scale : btag_bc_scale;
    return scaler.evaluate({"central", lvl, hadronFlavour.at(idx), fabs(eta(idx)), pt(idx)});
}

float Jet::getTotalBTagWeight(std::string btag_wp) {
    float weight = 1;
    if (!isMC) {
        return weight;
    }
    std::string tag_syst = systName(btag_bc_scale);
    std::string eff_syst = systName(btag_eff);
    for(auto i : list(Level::Loose)) {
        auto scaler = (hadronFlavour.at(i) == 0) ? btag_udsg_scale : btag_bc_scale;
        float bSF = scaler.evaluate({tag_syst, btag_wp, hadronFlavour.at(i), fabs(eta(i)), pt(i)});
        if (btag.at(i) > bjet_cuts[btag_wp]) {
            weight *= bSF;
        } else {
            float eff = btag_eff.evaluate({eff_syst, btag_wp, hadronFlavour.at(i), fabs(eta(i)), pt(i)});
            weight *= (1 - bSF * eff) / (1 - eff);
        }
    }
    return weight;
}

std::complex<float> Jet::get_momentum_change()
{
    std::complex<float> change;
    for(size_t i = 0; i < size(); ++i) {
        change += std::polar(pt(i)-nompt(i), phi(i));
    }
    return change;
}

void Jet::setup_jer(size_t i, GenericParticle& genJets)
{
    using namespace ROOT::Math::VectorUtil;
    float resolution = jet_resolution.evaluate({eta(i), m_pt.at(i), *rho});
    size_t g = genJetIdx.at(i);
    bool hasGenJet = genJetIdx.at(i) != -1;
    float pt_ratio = (hasGenJet) ? 1-genJets.pt(g)/m_pt.at(i) : 0.;

    for (auto var : all_vars) {
        float scale = jer_scale.evaluate({eta(i), jer_scale.name_by_var.at(var)});
        if (hasGenJet
            && deltaR(eta(i), genJets.eta(g), phi(i), genJets.phi(g)) < jet_dr/2
            && fabs(pt_ratio) < 3*resolution)
        {
            jer[var].push_back(1 + (scale-1)*pt_ratio);
        } else if (scale > 1.) {
            std::normal_distribution<> gaussian{0, resolution};
            jer[var].push_back(1 + gaussian(gen)*sqrt(pow(scale,2) - 1));
        } else {
            jer[var].push_back(1);
        }
    }
}

void Jet::setup_jes(size_t i)
{
    float delta = jes_scale.evaluate({eta(i), m_pt.at(i)});
    jes[eVar::Nominal].push_back(1);
    jes[eVar::Up].push_back(1+delta);
    jes[eVar::Down].push_back(1-delta);
}
