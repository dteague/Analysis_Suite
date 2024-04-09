#include "analysis_suite/skim/interface/ThreeTop.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

enum class Channel {
    SS_Dilepton,
    SS_Multi,
    Nonprompt_Dilepton,
    Nonprompt_Multi,
    OS_MisId,
    OS_CR,
    None,
};

enum class Subchannel {
    MM,
    EM,
    EE,
    Single_E,
    Single_M,
    None,
};

void ThreeTop::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Signal_Dilepton", Channel::SS_Dilepton);
    createTree("Signal_Multi", Channel::SS_Multi);

    if (!isMC_
        || groupName_.find("ttbar") != std::string::npos
        || groupName_.find("wjet") != std::string::npos) {
        use_nonprompt = true;
        createTree("Nonprompt_Dilepton", Channel::Nonprompt_Dilepton);
        createTree("Nonprompt_Multi", Channel::Nonprompt_Multi);
    }
    if (!isMC_) { // Charge Mis-id Fake Rate
        createTree("OS_Charge_MisId", Channel::OS_MisId);
    }

    rTop.setup(fReader);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

#include "analysis_suite/skim/interface/trigger_template.h"

    LOG_FUNC << "End of Init";
}

void ThreeTop::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "ParticleOut", &o_tightMuons);
    tree->Branch("TightElectron", "ElectronOut", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("BJets", "JetOut", &o_bJets);
    tree->Branch("Tops", "TopOut", &o_top);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
    tree->Branch("NBjets_loose", &o_nb_loose);
    tree->Branch("NBjets_medium", &o_nb_medium);
    tree->Branch("NBjets_tight", &o_nb_tight);

    tree->Branch("N_loose_mu", &o_nloose_mu);
    tree->Branch("N_loose_el", &o_nloose_el);
    tree->Branch("passZVeto", &o_pass_zveto);
    tree->Branch("Zmass", &o_ZMass);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void ThreeTop::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    rTop.clear();
    LOG_FUNC << "End of clearParticles";
}

void ThreeTop::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_htb.clear();
    o_met.clear();
    o_metphi.clear();
    o_centrality.clear();
    o_nloose_mu.clear();
    o_nloose_el.clear();
    o_pass_zveto.clear();
    o_ZMass.clear();

    o_nb_loose.clear();
    o_nb_medium.clear();
    o_nb_tight.clear();

    LOG_FUNC << "End of clearOutputs";
}

void ThreeTop::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getScaleFactor();
    (*weight) *= jet.getTotalBTagWeight("M");
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    // (*weight) *= rTop.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void ThreeTop::xor_lists(Particle& part, std::vector<size_t>& outlist)
{
    outlist.clear();
    auto tight = part.begin(Level::Tight);
    for (auto fake = part.begin(Level::Fake); fake != part.end(Level::Fake); ++fake) {
        if (tight == part.end(Level::Tight) || (*tight) > (*fake)) {
            outlist.push_back(*fake);
        } else {
            tight++;
        }
    }
    part.moveLevel(Level::Fake, Level::Tight);
}

void ThreeTop::ApplyChannelScaleFactors()
{
    if (!use_nonprompt) {
        return;
    } else if ((*currentChannel_) == Channel::Nonprompt_Dilepton ||
               (*currentChannel_) == Channel::Nonprompt_Multi) {
        for (auto i: fake_mu) {
            (*weight) *= sfMaker.getNonpromptFR(muon.eta(i), muon.pt(i), PID::Muon);
        }
        for (auto i: fake_el) {
            (*weight) *= sfMaker.getNonpromptFR(elec.eta(i), elec.pt(i), PID::Electron);
        }
        (*weight) *= std::pow(-1, fake_mu.size()+fake_el.size()+1);
    } else if ((*currentChannel_) == Channel::OS_MisId) {
        float tmp_weight = 0;
        for (auto i : elec.list(Level::Tight)) {
            tmp_weight += sfMaker.getChargeMisIdFR(elec.eta(i), elec.pt(i));
        }
        (*weight) *= tmp_weight;
    }
}

void ThreeTop::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    rTop.setGoodParticles(syst);
    LOG_FUNC << "End of setOtherGoodParticles";
}

void ThreeTop::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Fake) >= 2) {
        if (muon.size(Level::Fake) * elec.size(Level::Fake) > 0) {
            subChannel_ = Subchannel::EM;
        } else if (elec.size(Level::Fake) >= 2) {
            subChannel_ = Subchannel::EE;
        }  else {
            subChannel_ = Subchannel::MM;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}

bool ThreeTop::getTriggerValue()
{
    if (subChannel_ == Subchannel::EE) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::EE)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (year_ == Year::yr2018) { // 2018 uses just EGamma not double vs single EG datasets
            return (trig_cuts.pass_cut(Subchannel::EE)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (trig_cuts.dataset_or_trig(Subchannel::EE)) {
            return trig_cuts.pass_cut(Subchannel::EE);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_E);
        }
    } else if (subChannel_ == Subchannel::MM) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::MM)
                    || trig_cuts.pass_cut(Subchannel::Single_M));
        } else if (trig_cuts.dataset_or_trig(Subchannel::MM)) {
            return trig_cuts.pass_cut(Subchannel::MM);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_M);
        }
    } else if (subChannel_ == Subchannel::EM) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::EM)
                    || trig_cuts.pass_cut(Subchannel::Single_M)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (trig_cuts.dataset_or_trig(Subchannel::EM)) {
            return trig_cuts.pass_cut(Subchannel::EM);
        } else if (trig_cuts.dataset_or_trig(Subchannel::Single_M)) {
            return trig_cuts.pass_cut(Subchannel::Single_M);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_E);
        }
    } else {
        return false;
    }
}


bool ThreeTop::isSameSign(Level level)
{
    int q_total = 0;
    for (size_t idx : muon.list(level)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(level)) {
        q_total += elec.charge(idx);
    }
    return abs(q_total) == 1 || abs(q_total) == 2;
}



bool ThreeTop::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (!baseline_cuts()) return false;

    if (signal_cuts())
        (*currentChannel_) = (nLeps(Level::Tight) == 2) ? Channel::SS_Dilepton : Channel::SS_Multi;

    if (nonprompt_cuts()) {
        xor_lists(muon, fake_mu);
        xor_lists(elec, fake_el);
        (*currentChannel_) = (nLeps(Level::Fake) == 2) ? Channel::Nonprompt_Dilepton : Channel::Nonprompt_Multi;
    }

    if (charge_misid_cuts())
        (*currentChannel_) = Channel::OS_MisId;

    if (trees.find(*currentChannel_) == trees.end()) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool ThreeTop::baseline_cuts()
{
    LOG_FUNC << "Start of baseline_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());

    passCuts &= cuts.setCut("passLeadPt", getLeadPt() > 25);
    passCuts &= cuts.setCut("passFakeLeptonNum", nLeps(Level::Fake) >= 2);

    passCuts &= cuts.setCut("passTrigger", getTriggerValue());
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 250);
    passCuts &= cuts.setCut("passLowMassCut", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));

    LOG_FUNC << "End of baseline_cuts";
    return passCuts;
}

bool ThreeTop::signal_cuts()
{
    LOG_FUNC << "Start of signal_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) == 2 || nLeps(Level::Tight) == 3);
    passCuts &= cuts.setCut("passVetoFakeLeptons", nLeps(Level::Tight) == nLeps(Level::Fake));
    passCuts &= cuts.setCut("passSameSign;", isSameSign(Level::Tight));

    // Fill Cut flow
    cuts.setCut("pass2TightLeps", nLeps(Level::Tight) == 2);
    fillCutFlow(Channel::SS_Dilepton, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("pass3TightLeps", nLeps(Level::Tight) == 3);
    fillCutFlow(Channel::SS_Multi, cuts);

    LOG_FUNC << "End of signal_cuts";
    return passCuts;
}

bool ThreeTop::nonprompt_cuts()
{
    LOG_FUNC << "Start of nonprompt_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passMoreFakeLeptons", nLeps(Level::Fake) != nLeps(Level::Tight));
    passCuts &= cuts.setCut("passSameSign;", isSameSign(Level::Fake));

    // Fill Cut flow
    cuts.setCut("pass2FakeLeps", nLeps(Level::Fake) == 2);
    fillCutFlow(Channel::Nonprompt_Dilepton, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("pass3FakeLeps", nLeps(Level::Fake) == 3);
    fillCutFlow(Channel::Nonprompt_Multi, cuts);

    LOG_FUNC << "End of nonprompt_cuts";
    return passCuts;
}

bool ThreeTop::charge_misid_cuts()
{
    LOG_FUNC << "Start of charge_misid_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passVetoFakeLeptons", nLeps(Level::Tight) == nLeps(Level::Fake));
    passCuts &= cuts.setCut("pass1Electron", elec.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passOppositeSign;", !isSameSign(Level::Tight));

    // Fill Cut flow
    fillCutFlow(Channel::OS_MisId, cuts);

    LOG_FUNC << "End of charge_misid_cuts";
    return passCuts;
}

// bool ThreeTop::opposite_sign()
// {
//     LOG_FUNC << "Start of charge_misid_cuts";
//     bool passCuts = true;
//     CutInfo cuts;

//     passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) >= 2);
//     passCuts &= cuts.setCut("passVetoFakeLeptons", nLeps(Level::Tight) == nLeps(Level::Fake));
//     passCuts &= cuts.setCut("passOppositeSign;", !isSameSign(Level::Tight));

//     // Fill Cut flow
//     fillCutFlow(Channel::OS_MisId, cuts);

//     LOG_FUNC << "End of charge_misid_cuts";
//     return passCuts;
// }

void ThreeTop::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillOutput(*o_tightMuons, Level::Tight, event_bitmap);
    elec.fillElectron(*o_tightElectrons, Level::Tight, event_bitmap);

    rTop.fillTop(*o_top, Level::Loose, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);
    jet.fillJet(*o_bJets, Level::Bottom, event_bitmap);
    float muonMass = muon.massInRange(Level::Loose);
    float elecMass = elec.massInRange(Level::Loose);
    float zmass = (muonMass < 0) ? elecMass : muonMass;

    for (size_t systNum = 0; systNum < numSystematics(); ++systNum) {
        size_t syst = syst_to_index.at(systNum);
        if (syst == 0 && systNum != 0) {
            continue;
        }
        setupSyst(systNum);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Tight, syst));

        o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
        o_nb_medium.push_back(jet.n_medium_bjet.at(syst));
        o_nb_tight.push_back(jet.n_tight_bjet.at(syst));

        o_nloose_mu.push_back(muon.size(Level::Loose)-muon.size(Level::Fake));
        o_nloose_el.push_back(elec.size(Level::Loose)-elec.size(Level::Fake));
        o_pass_zveto.push_back(zmass < 0);
        o_ZMass.push_back(zmass);
    }
    LOG_FUNC << "End of FillValues";
}

float ThreeTop::getLeadPt(size_t idx)
{
    if (idx == 0) {
        return std::max(muPt(0), elPt(0));
    } else {
        return 0;
    }

    return 0;
}

void ThreeTop::printStuff()
{
    LOG_FUNC << "Start of printStuff";
    std::cout << "Event: " << *event << std::endl;
    std::cout << "Met: " << met.pt() << std::endl;
    std::cout << "HT: " << jet.getHT(Level::Tight, 0) << std::endl;
    std::cout << "njet: " << jet.size(Level::Tight) << std::endl;
    std::cout << "nbjet: " << jet.size(Level::Bottom) << std::endl;
    std::cout << "nlep: " << muon.size(Level::Tight) << " " << elec.size(Level::Tight) << std::endl;
    std::cout << "nlep loose: " << muon.size(Level::Fake) << " " << elec.size(Level::Fake) << std::endl;
    std::cout << std::endl;
    LOG_FUNC << "End of printStuff";
}

