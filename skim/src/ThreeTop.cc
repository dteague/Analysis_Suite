#include "analysis_suite/skim/interface/ThreeTop.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

enum class Channel {
    SS_Dilepton,
    SS_Multi,
    Nonprompt_FakeRate,
    OS_MisId,
    None,
};

enum class Subchannel {
    MM,
    EM,
    ME,
    EE,
    None,
};

#define getElec(var, i) (elec.var(elec.list(Level::Tight).at(i)))
#define getMuon(var, i) (muon.var(muon.list(Level::Tight).at(i)))

void ThreeTop::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Signal_Dilepton", Channel::SS_Dilepton);
    createTree("Signal_Multi", Channel::SS_Multi);
    // Charge Mis-id Fake Rate
    if (!isMC_) {
        createTree("Nonprompt_FakeRate", Channel::Nonprompt_FakeRate);
        createTree("OS_Charge_MisId", Channel::OS_MisId);
    }

    // trigEff_leadPt.setup(fOutput, "leadPt", 4, 100, 0, 100);

    rTop.setup(fReader);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    // Dilepton triggers
    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                                      "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300"});
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300"});
    } else if(year_ == Year::yr2017) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"});
        setupTrigger(Subchannel::EM, {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"});
    } else if (year_ == Year::yr2018) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"});
        setupTrigger(Subchannel::EM, {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL"});
    }
    setupTrigger(Subchannel::None);

    LOG_FUNC << "End of Init";
}

void ThreeTop::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "ParticleOut", &o_tightMuons);
    tree->Branch("TightElectron", "ParticleOut", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("BJets", "JetOut", &o_bJets);
    tree->Branch("Tops", "TopOut", &o_top);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
    tree->Branch("N_bloose", &o_nb_loose);
    tree->Branch("N_btight", &o_nb_tight);
    tree->Branch("N_loose_muon", &o_nloose_muon);
    tree->Branch("N_loose_elec", &o_nloose_elec);
    tree->Branch("passZVeto", &o_pass_zveto);
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
    o_nb_loose.clear();
    o_nb_tight.clear();
    o_nloose_muon.clear();
    o_nloose_elec.clear();
    o_pass_zveto.clear();
    LOG_FUNC << "End of clearOutputs";
}

void ThreeTop::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    LOG_EVENT << "weight after pu scale: " << (*weight);
    (*weight) *= sfMaker.getLHESF();
    LOG_EVENT << "weight after lhe scale: " << (*weight);
    (*weight) *= sfMaker.getLHEPdf();
    LOG_EVENT << "weight after pdf scale: " << (*weight);
    (*weight) *= sfMaker.getPartonShower();
    LOG_EVENT << "weight after ps scale: " << (*weight);
    (*weight) *= jet.getScaleFactor();
    LOG_EVENT << "weight after jet scale: " << (*weight);
    (*weight) *= jet.getTotalBTagWeight();
    LOG_EVENT << "weight after bjet scale: " << (*weight);
    (*weight) *= elec.getScaleFactor();
    LOG_EVENT << "weight after elec scale: " << (*weight);
    (*weight) *= muon.getScaleFactor();
    LOG_EVENT << "weight after muon scale: " << (*weight);
    (*weight) *= rTop.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void ThreeTop::applyNonprompt(Particle& part, PID pid)
{
    auto tight = part.begin(Level::Tight);
    int parity = 1;
    for (auto fake = part.begin(Level::Fake); fake != part.end(Level::Fake); ++fake) {
        if (tight == part.end(Level::Tight) || (*tight) > (*fake)) {
            (*weight) *= parity*sfMaker.getNonpromptFR(part.eta(*fake), part.pt(*fake), pid);
            parity *= -1;
        } else {
            tight++;
        }
    }
    part.moveLevel(Level::Fake, Level::Tight);
}

void ThreeTop::ApplyDataScaleFactors()
{
    if ((*currentChannel_) == Channel::Nonprompt_FakeRate) {
        applyNonprompt(muon, PID::Muon);
        applyNonprompt(elec, PID::Electron);
    } else if ((*currentChannel_) == Channel::OS_MisId) {
        for (auto i : elec.list(Level::Tight)) {
            (*weight) *= sfMaker.getChargeMisIdFR(elec.eta(i), elec.pt(i));
        }
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

    if (nLeps(Level::Tight) == 1 && nLeps(Level::Fake) == 2) {
        if (muon.size(Level::Fake) == 2) subChannel_ = Subchannel::MM;
        else if (elec.size(Level::Fake) == 2) subChannel_ = Subchannel::EE;
        else subChannel_ = Subchannel::EM;
    } else if (nLeps(Level::Tight) == 2) {
        if (muon.size(Level::Tight) == 2) subChannel_ = Subchannel::MM;
        else if (elec.size(Level::Tight) == 2) subChannel_ = Subchannel::EE;
        else subChannel_ = Subchannel::EM;
    } else if (nLeps(Level::Tight) == 3) {
        if (muPt(0) > elPt(0) && muPt(1) > elPt(0)) subChannel_ = Subchannel::MM;
        else if (elPt(0) > muPt(0) && elPt(1) > muPt(0)) subChannel_ = Subchannel::EE;
        else subChannel_ = Subchannel::EM;
    }
    LOG_FUNC << "End of setSubChannel";
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

    if (nonprompt_cuts())
        (*currentChannel_) = Channel::Nonprompt_FakeRate;

    if (charge_misid_cuts())
        (*currentChannel_) = Channel::OS_MisId;

    if (trees.find(*currentChannel_) == trees.end()) {
        return false;
    }

    if (!isMC_)  ApplyDataScaleFactors();

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

    // if (passCuts) trigEff_leadPt.fill(getLeadPt(), trig_cuts.pass_cut(subChannel_), subChannel_, *weight);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));
    passCuts &= cuts.setCut("passLeadPt", getLeadPt() > 25);

    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passBJetNumber", jet.size(Level::Bottom) >= 1);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 250);

    LOG_FUNC << "End of baseline_cuts";
    return passCuts;
}

bool ThreeTop::signal_cuts()
{
    LOG_FUNC << "Start of signal_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) == 2 || nLeps(Level::Tight) == 3);
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

    passCuts &= cuts.setCut("pass2FakeLeptons", nLeps(Level::Fake) == 2);
    passCuts &= cuts.setCut("passSameSign;", isSameSign(Level::Fake));

    // Fill Cut flow
    cuts.setCut("pass1FakeLeps", nLeps(Level::Fake) >= 1);
    cuts.setCut("pass2FakeLeps", nLeps(Level::Fake) == 2);
    fillCutFlow(Channel::Nonprompt_FakeRate, cuts);

    LOG_FUNC << "End of nonprompt_cuts";
    return passCuts;
}

bool ThreeTop::charge_misid_cuts()
{
    LOG_FUNC << "Start of charge_misid_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("pass2Leptons", nLeps(Level::Tight) == 2);
    passCuts &= cuts.setCut("pass1Electron", elec.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passOppositeSign;", !isSameSign(Level::Tight));

    // Fill Cut flow
    fillCutFlow(Channel::OS_MisId, cuts);

    LOG_FUNC << "End of charge_misid_cuts";
    return passCuts;
}

void ThreeTop::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillOutput(*o_tightMuons, Level::Tight, event_bitmap);
    elec.fillOutput(*o_tightElectrons, Level::Tight, event_bitmap);

    rTop.fillTop(*o_top, Level::Loose, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);
    jet.fillJet(*o_bJets, Level::Bottom, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Tight, syst));
        o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
        o_nb_tight.push_back(jet.n_tight_bjet.at(syst));
        o_nloose_muon.push_back(muon.size(Level::Loose));
        o_nloose_elec.push_back(elec.size(Level::Loose));
        o_pass_zveto.push_back(muon.passZVeto() && elec.passZVeto());
    }
    LOG_FUNC << "End of FillValues";
}

float ThreeTop::getLeadPt(size_t idx)
{
    if (subChannel_ == Subchannel::None) {
        return 0.;
    } else if (idx == 0) {
        return std::max(muPt(0), elPt(0));
    } else if (idx == 1) {
        if (subChannel_ == Subchannel::MM)
            return muPt(1);
        else if (subChannel_ == Subchannel::EE)
            return elPt(1);
        else if (subChannel_ == Subchannel::EM)
            return std::min(muPt(0), elPt(0));
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
    std::cout << "lepVeto: " << muon.passZVeto() << " " << elec.passZVeto() << std::endl;
    std::cout << std::endl;
    LOG_FUNC << "End of printStuff";
}

