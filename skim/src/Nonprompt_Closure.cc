#include "analysis_suite/skim/interface/Nonprompt_Closure.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    TightTight,
    TightFake,
    FakeFake,
    None,
};

enum class Subchannel {
    MM,
    EM,
    EE,
    None,
};

void Nonprompt_Closure::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    createTree("Closure_FF", Channel::FakeFake);
    createTree("Closure_TF", Channel::TightFake);
    if (isMC_) {
        createTree("Closure_TT", Channel::TightTight);
    }

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }


    // Dilepton triggers
    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                      "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
                                      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300"});
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
                                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",});
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

void Nonprompt_Closure::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("FakeMuon", "LeptonOut_Fake", &o_fakeMuons);
    tree->Branch("TightMuon", "LeptonOut_Fake", &o_tightMuons);
    tree->Branch("FakeElectron", "LeptonOut_Fake", &o_fakeElectrons);
    tree->Branch("TightElectron", "LeptonOut_Fake", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("BJets", "JetOut", &o_bJets);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("N_bloose", &o_nb_loose);
    tree->Branch("N_btight", &o_nb_tight);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void Nonprompt_Closure::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void Nonprompt_Closure::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_htb.clear();
    o_met.clear();
    o_metphi.clear();
    o_nb_loose.clear();
    o_nb_tight.clear();
    LOG_FUNC << "End of clearOutputs";
}

void Nonprompt_Closure::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    LOG_EVENT << "weight after pu scale: " << (*weight);
    (*weight) *= sfMaker.getLHESF();
    LOG_EVENT << "weight after lhe scale: " << (*weight);
    (*weight) *= jet.getScaleFactor();
    LOG_EVENT << "weight after jet scale: " << (*weight);
    (*weight) *= elec.getScaleFactor();
    LOG_EVENT << "weight after elec scale: " << (*weight);
    (*weight) *= muon.getScaleFactor();
    LOG_EVENT << "weight after muon scale: " << (*weight);
    LOG_FUNC << "End of ApplyScaleFactors";
}

void Nonprompt_Closure::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    LOG_FUNC << "End of setOtherGoodParticles";
}

void Nonprompt_Closure::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if (nLeps(Level::Fake) == 2) {
        if (muon.size(Level::Fake) == 2) {
            subChannel_ = Subchannel::MM;
        } else if (elec.size(Level::Fake) == 2) {
            subChannel_ = Subchannel::EE;
        } else {
            subChannel_ = Subchannel::EM;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}

bool Nonprompt_Closure::isSameSign()
{
    int q_total = 0;
    for (size_t idx : muon.list(Level::Fake)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(Level::Fake)) {
        q_total += elec.charge(idx);
    }
    // if 2 leptons, SS -> +1 +1 / -1 -1 -> abs(q) == 2
    // OS cases are 0 and 3, so no overlap
    return abs(q_total) == 1 || abs(q_total) == 2;
}


bool Nonprompt_Closure::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (closure_cuts()) {
        if (nLeps(Level::Tight) == 2) (*currentChannel_) = Channel::TightTight;
        else if (nLeps(Level::Tight) == 1) (*currentChannel_) = Channel::TightFake;
        else (*currentChannel_) = Channel::FakeFake;
    }

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool Nonprompt_Closure::closure_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2FakeLep",  nLeps(Level::Fake) == 2);
    // passCuts &= cuts.setCut("passTightLep", nLeps(Level::Tight) >= 1);
    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    passCuts &= cuts.setCut("passSSLeptons", isSameSign());
    passCuts &= cuts.setCut("passZVeto", muon.passZVeto() && elec.passZVeto());
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 50);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

    // Fill Cut flow
    cuts.setCut("pass2TightLeps", nLeps(Level::Tight) == 2);
    fillCutFlow(Channel::TightTight, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("pass1TightLeps", nLeps(Level::Tight) == 1);
    fillCutFlow(Channel::TightFake, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("pass0TightLeps", nLeps(Level::Tight) == 0);
    fillCutFlow(Channel::FakeFake, cuts);

    return passCuts;
}

float Nonprompt_Closure::getLeadPt()
{
    if (subChannel_ == Subchannel::MM) {
        return muon.pt(Level::Fake, 0);
    } else if (subChannel_ == Subchannel::EE) {
        return elec.pt(Level::Fake, 0);
    } else if(subChannel_ == Subchannel::EM){
        return std::max(muon.pt(Level::Fake, 0), elec.pt(Level::Fake, 0));
    }
    return 0.;
}

void Nonprompt_Closure::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton_Iso(*o_fakeMuons, Level::FakeNotTight, event_bitmap);
    muon.fillLepton_Iso( *o_tightMuons, Level::Tight, event_bitmap);
    elec.fillLepton_Iso(*o_fakeElectrons, Level::FakeNotTight, event_bitmap);
    elec.fillLepton_Iso( *o_tightElectrons, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);
    jet.fillJet(*o_bJets, Level::Bottom, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
        o_nb_tight.push_back(jet.n_tight_bjet.at(syst));
    }
    LOG_FUNC << "End of FillValues";
}
