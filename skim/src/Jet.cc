#include "analysis_suite/skim/interface/Jet.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void Jet::setup(TTreeReader& fReader, std::string data_run, std::vector<Systematic> used_jec_systs_)
{
    closeLep_dr.resize(50);
    GenericParticle::setup("Jet", fReader);
    jetId.setup(fReader, "Jet_jetId");
    btag.setup(fReader, "Jet_btagDeepFlavB");
    puId.setup(fReader, "Jet_puId");
    if (isMC) {
        hadronFlavour.setup(fReader, "Jet_hadronFlavour");
        genJetIdx.setup(fReader, "Jet_genJetIdx");
    }
    neHEF.setup(fReader, "Jet_neHEF");
    neEmEF.setup(fReader, "Jet_neEmEF");
    muEF.setup(fReader, "Jet_muEF");
    chHEF.setup(fReader, "Jet_chHEF");
    chEmEF.setup(fReader, "Jet_chEmEF");
    nConstituents.setup(fReader, "Jet_nConstituents");

    setup_map(Level::Loose_inHEM);
    setup_map(Level::Loose_PU);
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

    jec.init(fReader, used_jec_systs);
}

void Jet::setup_jets(GenericParticle& genJet, size_t run)
{
    if (year_ != Year::yr2018) {
        applyHEMVeto= false;
    } else if (isMC) {
        applyHEMVeto= 1.0*std::rand()/RAND_MAX < 0.64;
    } else {
        applyHEMVeto= run >= 319077;
    }
    jec.setupJEC(*this, genJet);
    // closeLep_dr.resize(size());
    std::fill(closeLep_dr.begin(), closeLep_dr.end(), -1);

}

void Jet::setup_shift_output(TTree* tree) {
    if (o_jet_shifts.size() == 0 && used_jec_systs.size() != 0) {
        for (auto syst: used_jec_systs) {
            o_jet_shifts.push_back(new JECShiftOut());
            o_bjet_shifts.push_back(new JECShiftOut());
        }
    }
    for (size_t i = 0; i < used_jec_systs.size(); ++i) {
        std::string syst_name = get_by_val(syst_by_name, used_jec_systs.at(i));
        tree->Branch(("Jets_"+syst_name).c_str(), "JECShiftOut", &o_jet_shifts[i]);
        tree->Branch(("BJets_"+syst_name).c_str(), "JECShiftOut", &o_bjet_shifts[i]);
    }
}

void Jet::fillJet(JetOut& output, Level level, const Bitmap& event_bitmap)
{
    output.clear();
    std::vector<JECShiftOut*>& shifts = (level == Level::Bottom) ? o_bjet_shifts : o_jet_shifts;
    for (size_t n_syst=0; n_syst < used_jec_systs.size(); ++n_syst) {
        shifts[n_syst]->clear();
    }
    for (size_t idx = 0; idx < size(); ++idx) {
        bool pass = fillParticle(output, level, idx, event_bitmap);
        if (pass) {
            if (btag.at(idx) > tight_bjet_cut) output.pass_btag.push_back(3);
            else if (btag.at(idx) > medium_bjet_cut) output.pass_btag.push_back(2);
            else if (btag.at(idx) > loose_bjet_cut) output.pass_btag.push_back(1);
            else output.pass_btag.push_back(0);
            output.discriminator.push_back(btag.at(idx));
            // output.veto.push_back(passVeto(eta(idx), phi(idx)));
            if(isMC) {
                output.flavor.push_back(hadronFlavour.at(idx));
            }
            for (size_t n_syst=0; n_syst < used_jec_systs.size(); ++n_syst) {
                auto syst = used_jec_systs.at(n_syst);
                shifts[n_syst]->pt_up.push_back(pt_shift(idx, syst, eVar::Up));
                shifts[n_syst]->pt_down.push_back(pt_shift(idx, syst, eVar::Down));
                shifts[n_syst]->mass_up.push_back(mass_shift(idx, syst, eVar::Up));
                shifts[n_syst]->mass_down.push_back(mass_shift(idx, syst, eVar::Down));
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
            && closeLep_dr.at(i) < 0
            // && passVeto(eta(i), phi(i))
            ) {
            if (pt(i) < 50 && (puId.at(i) & PU_Tight) != PU_Tight) { // Issue with puid for 2016, need 111 for id now
                m_partList[Level::Loose_PU]->push_back(i);
            } else if(!passVeto(eta(i), phi(i))) {
                m_partList[Level::Loose_inHEM]->push_back(i);
            } else {
                m_partList[Level::Loose]->push_back(i);
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
    jec.set_syst(isJECSyst());
}

float Jet::getScaleFactor()
{
    float weight = 1;
    std::string syst = systName(puid_scale);
    for (auto idx : list(Level::Loose)) {
        if (pt(idx) < 50 && genJetIdx.at(idx) != -1)
            weight *= puid_scale.evaluate({eta(idx), pt(idx), syst, "T"});
    }
    for (auto idx : list(Level::Loose_PU)) {
        if (pt(idx) < 50 && genJetIdx.at(idx) != -1) {
            float eff = puid_scale.evaluate({eta(idx), pt(idx), "MCEff", "T"});
            float sf = puid_scale.evaluate({eta(idx), pt(idx), syst, "T"});
            weight *= (1-sf*eff)/(1-eff);
        }
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
        float pt_ = m_pt.at(i);
        bool isCharm = hadronFlavour.at(i) == static_cast<Int_t>(PID::Charm);
        if (isCharm == charmSyst) {
            weight *= btag_shape_scale.evaluate({syst, hadronFlavour.at(i), fabs(eta(i)), pt_, btag.at(i)});
        } else {
            weight *= btag_shape_scale.evaluate({"central", hadronFlavour.at(i), fabs(eta(i)), pt_, btag.at(i)});
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
    return jec.get_met_change();
}

bool Jet::passVeto(float eta, float phi)
{
    if(phi > pi) phi = pi;
    else if (phi < -pi) phi = -pi;
    inHEM = false;
    bool pass = (veto_map.evaluate({"jetvetomap", eta, phi}) < 1e-5);

    if (pass) {
        return true;
    }  else if (year_ == Year::yr2018) {
        inHEM = (veto_map.evaluate({"jetvetomap_hem1516", eta, phi}) > 1e-5);
        return inHEM && !applyHEMVeto;
    } else {
        return false;
    }
}
