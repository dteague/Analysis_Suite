#include "analysis_suite/skim/interface/Lepton_MisId_Closure.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

enum class Channel {
    OS_MR,
    OS_CR,
    SS_CR,
    None,
};

enum class Subchannel {
    MM,
    EM,
    EE,
    None,
};

void Closure_MisId::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    // Charge Mis-id Fake Rate
    createTree("SS", Channel::SS_CR);
    createTree("OS_CR", Channel::OS_CR);
    createTree("OS_MR", Channel::OS_MR);
    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
        LHE_HT.setup(fReader, "LHE_HT");
        LHE_pdgId.setup(fReader, "LHEPart_pdgId");
    }

    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                                      "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                                      "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300"});
        setupTrigger(Subchannel::EE, {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300"});
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                                      "HLT_DoubleMu8_Mass8_PFHT300"});
        setupTrigger(Subchannel::EM, {"HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                                      "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
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

void Closure_MisId::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "LeptonOut", &o_tightMuons);
    tree->Branch("TightElectron", "LeptonOut", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("Nloose_Mu", &nlooseMu);
    tree->Branch("Nloose_El", &nlooseEl);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
    if (isMC_) {
        tree->Branch("LHE_HT", &o_ht_lhe);
        tree->Branch("LHE_nLeps", &o_nlhe_leps);
    }
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void Closure_MisId::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void Closure_MisId::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_ht_lhe.clear();
    o_nlhe_leps.clear();
    o_htb.clear();
    o_met.clear();
    o_metphi.clear();
    o_centrality.clear();
    LOG_FUNC << "End of clearOutputs";
}

void Closure_MisId::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    LOG_EVENT << "weight after pu scale: " << (*weight);
    (*weight) *= sfMaker.getLHESF();
    LOG_EVENT << "weight after lhe scale: " << (*weight);
    (*weight) *= sfMaker.getLHEPdf();
    LOG_EVENT << "weight after lhe pdf: " << (*weight);
    (*weight) *= jet.getScaleFactor();
    LOG_EVENT << "weight after jet scale: " << (*weight);
    (*weight) *= elec.getScaleFactor();
    LOG_EVENT << "weight after elec scale: " << (*weight);
    (*weight) *= muon.getScaleFactor();
    LOG_EVENT << "weight after muon scale: " << (*weight);

    LOG_FUNC << "End of ApplyScaleFactors";
}

bool Closure_MisId::isSameSign()
{
    int q_total = 0;
    for (size_t idx : muon.list(Level::Tight)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(Level::Tight)) {
        q_total += elec.charge(idx);
    }
    // if 2 leptons, SS -> +1 +1 / -1 -1 -> abs(q) == 2
    // if 3 leptons, SS -> +1 +1 -1 / -1 -1 +1 -> abs(q) == 1
    // OS cases are 0 and 3, so no overlap
    return abs(q_total) == 1 || abs(q_total) == 2;
}

void Closure_MisId::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Tight) == 2) {
        if (muon.size(Level::Tight) == 2) {
            subChannel_ = Subchannel::MM;
        } else if (elec.size(Level::Tight) == 2) {
            subChannel_ = Subchannel::EE;
        }  else {
            subChannel_ = Subchannel::EM;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}



bool Closure_MisId::getCutFlow()
{
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (closure_cuts()) {
        (*currentChannel_) = (isSameSign()) ? Channel::SS_CR : Channel::OS_CR;
    }
    if(measurement_cuts()) {
        (*currentChannel_) = Channel::OS_MR;
    }
    set_LHELeps();

    if (*currentChannel_ == Channel::None) {
        return false; // Pass no channels
    }

    return true;
}

bool Closure_MisId::closure_cuts() {
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2Electrons", elec.size(Level::Tight) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    float mass = get_mass();
    passCuts &= cuts.setCut("passZCut", mass > 70. && mass < 115);
    // passCuts &= cuts.setCut("passMetCut", met.pt() < 50);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) < 250);

    int charge = 0;
    if (elec.size(Level::Tight) == 2) {
        charge = elec.charge(Level::Tight, 0) * elec.charge(Level::Tight, 1);
    }

    cuts.setCut("passSSElectrons", charge > 0);
    fillCutFlow(Channel::SS_CR, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("passOSElectrons", charge < 0);
    fillCutFlow(Channel::OS_CR, cuts);

    return passCuts;
}

bool Closure_MisId::measurement_cuts() {
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2Leptons;", nLeps(Level::Tight) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    passCuts &= cuts.setCut("passZCut", get_mass() > 50);
    // passCuts &= cuts.setCut("passOppositeSign", !isSameSign());
    // passCuts &= cuts.setCut("passHasElectron", elec.size(Level::Tight) > 0);
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    // passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

    fillCutFlow(Channel::OS_MR, cuts);

    return passCuts;
}

float Closure_MisId::getLeadPt()
{
    if (subChannel_ == Subchannel::MM) {
        return muon.pt(Level::Tight, 0);
    } else if (subChannel_ == Subchannel::EE) {
        return elec.pt(Level::Tight, 0);
    } else if (subChannel_ == Subchannel::EM) {
        return std::max(muon.pt(Level::Tight, 0), elec.pt(Level::Tight, 0));
    }
    return 0.;
}

float Closure_MisId::get_mass() {
    if (subChannel_ == Subchannel::None) {
        return -1;
    } else if (subChannel_ == Subchannel::MM) {
        return (muon.p4(Level::Tight, 0) + muon.p4(Level::Tight, 1)).M();
    } else if (subChannel_ == Subchannel::EE) {
        return (elec.p4(Level::Tight, 0) + elec.p4(Level::Tight, 1)).M();
    } else {
        return (muon.p4(Level::Tight, 0) + elec.p4(Level::Tight, 0)).M();
    }
}

void Closure_MisId::set_LHELeps()
{
    nLHE_leps = 0;
    for (size_t i = 0; i < LHE_pdgId.size(); ++i) {
        if (abs(LHE_pdgId.at(i)) == 11 || abs(LHE_pdgId.at(i)) == 13 || abs(LHE_pdgId.at(i)) == 15) {
            nLHE_leps++;
        }
    }
}

void Closure_MisId::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton(*o_tightMuons, Level::Tight, event_bitmap);
    elec.fillLepton(*o_tightElectrons, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    nlooseMu = muon.size(Level::Loose);
    nlooseEl = elec.size(Level::Loose);
    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);
        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Tight, syst));
        if (isMC_) {
            o_ht_lhe.push_back(*LHE_HT);
            o_nlhe_leps.push_back(nLHE_leps);
        }
    }
    LOG_FUNC << "End of FillValues";
}

void Closure_MisId::printStuff()
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
