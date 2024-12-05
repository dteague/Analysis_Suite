#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void Jet::setup(TTreeReader& fReader, std::string data_run, std::vector<Systematic> used_jec_systs_)
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
        // Stuff for Met recalc
        t1_rawpt.setup(fReader, "CorrT1METJet_rawPt");
        t1_eta.setup(fReader, "CorrT1METJet_eta");
        t1_phi.setup(fReader, "CorrT1METJet_phi");
        t1_area.setup(fReader, "CorrT1METJet_area");
        t1_muonSub.setup(fReader, "CorrT1METJet_muonSubtrFactor");
    }
    neHEF.setup(fReader, "Jet_neHEF");
    neEmEF.setup(fReader, "Jet_neEmEF");
    muEF.setup(fReader, "Jet_muEF");
    chHEF.setup(fReader, "Jet_chHEF");
    chEmEF.setup(fReader, "Jet_chEmEF");
    nConstituents.setup(fReader, "Jet_nConstituents");
    muonSubtrFactor.setup(fReader, "Jet_muonSubtrFactor");
    area.setup(fReader, "Jet_area");

    setup_map(Level::Loose_inHEM);
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

    bjet_jec_syst[Systematic::Jet_JEC_Absolute] += yearNum.at(year_);
    bjet_jec_syst[Systematic::Jet_JEC_BBEC1] += yearNum.at(year_);
    bjet_jec_syst[Systematic::Jet_JEC_EC2] += yearNum.at(year_);
    bjet_jec_syst[Systematic::Jet_JEC_HF] += yearNum.at(year_);
    bjet_jec_syst[Systematic::Jet_JEC_RelativeSample] += yearNum.at(year_);

    used_jec_systs = used_jec_systs_;


    if (isMC) {
        // Pileup Weights
        auto jmar_set = getScaleFile("JME", "jmar");
        puid_scale = WeightHolder(jmar_set->at("PUJetID_eff"),
                                  Systematic::Jet_PUID, {"nom","up","down"});

        // BTagging Weights
        auto btag_set = getScaleFile("BTV", "btagging");
        btag_shape_scale = WeightHolder(btag_set->at("deepJet_shape"),
                                        Systematic::BJet_BTagging, {"central","up","down"});
        // Btagging Cutbased weights
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

    auto veto_set = getScaleFile("JME", "jetvetomaps");
    std::string veto_map_name = "Summer19UL"+yearNum.at(year_)+"_V1";
    veto_map = WeightHolder(veto_set->at(veto_map_name), Systematic::Nominal, {"jetvetomap", "", ""});

    m_jet_scales[Systematic::Nominal] = {{eVar::Nominal, std::vector<float>()}};
    for (auto syst : used_jec_systs) {
        m_jet_scales[syst] = {
            {eVar::Up, std::vector<float>()},
            {eVar::Down, std::vector<float>()}
        };
    }
}

void Jet::fillJet(JetOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    const auto central = m_jet_scales.at(Systematic::Nominal).at(eVar::Nominal);
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            if (btag.at(idx) > tight_bjet_cut) output.pass_btag.push_back(3);
            else if (btag.at(idx) > medium_bjet_cut) output.pass_btag.push_back(2);
            else if (btag.at(idx) > loose_bjet_cut) output.pass_btag.push_back(1);
            else output.pass_btag.push_back(0);
            output.discriminator.push_back(btag.at(idx));
            if(isMC) {
                output.flavor.push_back(hadronFlavour.at(idx));
            }
            output.pt_shift.push_back(std::vector<Float_t>());
            output.mass_shift.push_back(std::vector<Float_t>());
            for (size_t n_syst=0; n_syst < used_jec_systs.size(); ++n_syst) {
                auto syst = used_jec_systs.at(n_syst);
                output.pt_shift.back().push_back(m_pt.at(idx)*m_jet_scales[syst][eVar::Up][idx]);
                output.pt_shift.back().push_back(m_pt.at(idx)*m_jet_scales[syst][eVar::Down][idx]);
                output.mass_shift.back().push_back(m_mass.at(idx)*m_jet_scales[syst][eVar::Up][idx]);
                output.mass_shift.back().push_back(m_mass.at(idx)*m_jet_scales[syst][eVar::Down][idx]);
            }
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
            && (closeJetDr_by_index.find(i) == closeJetDr_by_index.end())
            && (pt(i) > 50 || (puId.at(i) & PU_Tight) == PU_Tight) // Issue with puid for 2016, need 111 for id now
            ) {
            if (passVeto(eta(i), phi(i))) {
                m_partList[Level::Loose]->push_back(i);
            } else if (applyHEMVeto && inHEM) {
                m_partList[Level::Loose_inHEM]->push_back(i);
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

float Jet::getMHT()
{
    auto mht = std::polar(0, 0);
    for (size_t idx = 0; idx < size(); ++idx) {
        mht -= std::polar(pt(idx), phi(idx));
    }
    return std::abs(mht);
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

float Jet::getScaleFactor()
{
    float weight = 1;
    std::string syst = systName(puid_scale);
    for (auto idx : list(Level::Loose)) {
        if (pt(idx) < 50 && genJetIdx.at(idx) != -1)
            weight *= puid_scale.evaluate({eta(idx), pt(idx), syst, "T"});
    }
    return weight;
}

float Jet::getBJetWeight(size_t idx, std::string lvl)
{
    auto scaler = (hadronFlavour.at(idx) == 0) ? btag_udsg_scale : btag_bc_scale;
    return scaler.evaluate({"central", lvl, hadronFlavour.at(idx), fabs(eta(idx)), pt(idx)});
}

float Jet::getTotalBTagWeight(std::string btag_wp)
{
    float weight = 1;
    bool charmSyst = std::find(charm_systs.begin(), charm_systs.end(), currentSyst) != charm_systs.end();
    std::string syst = "central";

    if (shape_by_syst.find(currentSyst) != shape_by_syst.end()) {
        syst = ((currentVar == eVar::Up) ? "up_" : "down_") + shape_by_syst.at(currentSyst);
    } else if (bjet_jec_syst.find(currentSyst) != bjet_jec_syst.end()) {
        syst = ((currentVar == eVar::Up) ? "up_" : "down_") + bjet_jec_syst.at(currentSyst);
    }
    for (auto i : list(Level::Loose)) {
        bool isCharm = hadronFlavour.at(i) == static_cast<Int_t>(PID::Charm);
        if (isCharm == charmSyst) {
            weight *= btag_shape_scale.evaluate({syst, hadronFlavour.at(i), fabs(eta(i)), pt(i), btag.at(i)});
        } else {
            weight *= btag_shape_scale.evaluate({"central", hadronFlavour.at(i), fabs(eta(i)), pt(i), btag.at(i)});
        }
    }
    return weight;
}

float Jet::getCutBasedBTagWeight(std::string btag_wp)
{
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
    return m_met_change[currentSyst][currentVar];
}

bool Jet::passVeto(float eta, float phi)
{
    if(phi > pi) phi = pi;
    else if (phi < -pi) phi = -pi;

    if (veto_map.evaluate({"jetvetomap", eta, phi}) < 1e-5) {
        return true;
    }  else if (!applyHEMVeto && year_ == Year::yr2018) {
        return veto_map.evaluate({"jetvetomap_hem1516", eta, phi}) > 1e-5;
    } else if (applyHEMVeto){
        inHEM = (veto_map.evaluate({"jetvetomap_hem1516", eta, phi}) > 1e-5);
    }
    return false;
}


void Jet::setupJEC(GenericParticle& genJet, size_t run)
{
    if (year_ != Year::yr2018) {
        applyHEMVeto= false;
    } else if (isMC) {
        applyHEMVeto= 1.0*std::rand()/RAND_MAX < 0.64;
    } else {
        applyHEMVeto= run >= 319077;
    }

    m_met_change[Systematic::Nominal][eVar::Nominal] = std::polar(0.f, 0.f);

    if (!isMC) {
        m_jet_scales[Systematic::Nominal][eVar::Nominal].assign(size(), 1);
    } else {
        for (auto& [syst, var_scales] : m_jet_scales ) {
            for (auto& [var, scales] : var_scales) {
                scales.resize(size());
                m_met_change[syst][var] = std::polar(0.f, 0.f);
            }
        }
    }

    for(size_t i = 0; i < size(); ++i) {
        auto jer = get_jer(m_pt.at(i), eta(i), i, genJet, false);
        float central = jer.at(0);

        m_jet_scales[Systematic::Nominal][eVar::Nominal][i] = central;
        for (auto syst: used_jec_systs) {
            if (syst == Systematic::Jet_JER) {
                m_jet_scales[syst][eVar::Up][i] = jer.at(1);
                m_jet_scales[syst][eVar::Down][i] = jer.at(2);
            } else {
                float jec_unc = get_jec_unc(m_pt.at(i), eta(i), jec_unc_vec[syst]);
                m_jet_scales[syst][eVar::Down][i] = central*(1-jec_unc);
                m_jet_scales[syst][eVar::Up][i] = central*(1+jec_unc);
            }
        }
        if ((m_pt.at(i)*(1-muonSubtrFactor.at(i)) > 15) && (neEmEF.at(i)+chEmEF.at(i) < 0.9)) {
            fix_met(rawPt(i), m_pt.at(i), eta(i), phi(i), area.at(i), muonSubtrFactor.at(i), jer);
        }
    }

    for(size_t i = 0; i < t1_rawpt.size(); ++i) {
        float jec_pt = t1_rawpt.at(i)*jec_total->evaluate({t1_area.at(i), t1_eta.at(i), t1_rawpt.at(i), *rho});
        if (jec_pt*(1-t1_muonSub.at(i)) > 15) {
            auto jer = get_jer(jec_pt, t1_eta.at(i), i, genJet, true);
            fix_met(t1_rawpt.at(i), jec_pt, t1_eta.at(i), t1_phi.at(i), t1_area.at(i), t1_muonSub.at(i), jer);

        }
    }

}

std::vector<float> Jet::get_jer(float pt, float eta, size_t idx, GenericParticle& genJets, bool isLowPt)
{

    if (!isMC) {
        return {1.0};
    }

    float resolution = jet_resolution.evaluate({eta, pt, *rho});
    float pt_ratio = 100;
    bool matched_gen = false;
    float rand = 1.;

    if (!isLowPt) {
        size_t g = genJetIdx.at(idx);
        pt_ratio = (genJetIdx.at(idx) != -1) ? 1-genJets.pt(g)/pt : 0.;
        matched_gen = (genJetIdx.at(idx) != -1)
            && (deltaR2(eta, genJets.eta(g), phi(idx), genJets.phi(g)) < 0.04) // (0.4/2)**2
            && (fabs(pt_ratio) < 3*resolution);
    } else {
        float gen_dR = 0.04;
        for (size_t i=0; i < genJets.size(); ++i) {
            float dr = deltaR2(eta, genJets.eta(i), t1_phi.at(idx), genJets.phi(i));
            if (dr < gen_dR) {
                gen_dR = dr;
                pt_ratio = 1-genJets.pt(i)/pt;
                matched_gen = (fabs(pt_ratio) < 3*resolution); // Already required dR<0.2 by def
            }
        }
    }

    if (!matched_gen) {
        std::normal_distribution<> gaussian{0, resolution};
        rand = gaussian(gen);
    }

    std::vector<float> output;
    for (auto var : all_vars) {
        float scale = jer_scale.evaluate({eta, jer_scale.name_by_var.at(var)});
        if (matched_gen) {
            output.push_back(1 + (scale-1)*pt_ratio);
        } else if (scale > 1.) {
            output.push_back(1 + rand*sqrt(pow(scale,2) - 1));
        } else {
            output.push_back(1);
        }
    }
    return output;
}

float Jet::get_jec_unc(float pt, float eta, correction::Correction::Ref& jec_unc)
{
    return jec_unc->evaluate({eta, pt});
}

void Jet::fix_met(float raw_pt, float pt, float eta, float phi, float area, float muonSub, std::vector<float> jer_)
{
    float jec_l1 = jec_l1_total->evaluate({area, eta, raw_pt, *rho});
    float l1_pt = raw_pt*(jec_l1*(1-muonSub) + muonSub);
    float met_pt = pt*(1-muonSub)+raw_pt*muonSub;

    float central = jer_.at(0);
    m_met_change[Systematic::Nominal][eVar::Nominal] += std::polar(central*met_pt-l1_pt, phi);
    for (auto syst: used_jec_systs) {
        if (syst == Systematic::Jet_JER) {
            m_met_change[syst][eVar::Up] += std::polar(jer_.at(1)*met_pt-l1_pt, phi);
            m_met_change[syst][eVar::Down] += std::polar(jer_.at(2)*met_pt-l1_pt, phi);
        } else {
            float jec_unc = get_jec_unc(met_pt, eta, jec_unc_vec[syst]);
            m_met_change[syst][eVar::Down] += std::polar((central-jec_unc)*met_pt-l1_pt, phi);
            m_met_change[syst][eVar::Up] += std::polar((central+jec_unc)*met_pt-l1_pt, phi);
        }
    }
}
