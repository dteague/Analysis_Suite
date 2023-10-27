#include "analysis_suite/skim/interface/Nonprompt_Closure.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    TightTight,
    TightFake,
    FakeFake,
    DY_Tight,
    DY_Fake,
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

void Nonprompt_Closure::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PF;
    BaseSelector::Init(tree);


    if (groupName_.find("ttbar") != std::string::npos ||
        groupName_.find("wjets") != std::string::npos ||
        !isMC_) {
        createTree("Closure_FF", Channel::FakeFake);
        createTree("Closure_TF", Channel::TightFake);
        if (isMC_) {
            createTree("Closure_TT", Channel::TightTight);
        }
    }

    if (groupName_.find("DY") != std::string::npos ||
        groupName_.find("zz") != std::string::npos ||
        groupName_.find("wz") != std::string::npos ||
        groupName_.find("ww") != std::string::npos ||
        !isMC_) {
        createTree("DY_Fake", Channel::DY_Fake);
        createTree("DY_Tight", Channel::DY_Tight);
    }

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);
    muon.setup_map(Level::LooseNotFake);
    elec.setup_map(Level::LooseNotFake);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }

    #include"analysis_suite/skim/interface/trigger_template.h"

    LOG_FUNC << "End of Init";
}

void Nonprompt_Closure::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    // tree->Branch("LooseMuon", "LeptonOut_Fake", &o_looseMuons);
    tree->Branch("FakeMuon", "LeptonOut_Fake", &o_fakeMuons);
    tree->Branch("TightMuon", "LeptonOut_Fake", &o_tightMuons);
    // tree->Branch("LooseElectron", "LeptonOut_Fake", &o_looseElectrons);
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
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getPileupIDWeight();
    (*weight) *= jet.getTotalBTagWeight();
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void Nonprompt_Closure::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    muon.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
    elec.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
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
    } else if (nLeps(Level::Fake) == 3) {
        if (muon.size(Level::Tight) == 2) {
            subChannel_ = Subchannel::MM;
        } else if(elec.size(Level::Tight) == 2) {
            subChannel_ = Subchannel::EE;
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

bool Nonprompt_Closure::getTriggerValue()
{
    if (subChannel_ == Subchannel::EE) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::EE)
                    || trig_cuts.pass_cut(Subchannel::Single_E));

        } else if (year_ == Year::yr2018) {
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

    if (dy_closure_cuts()) {
        if (nLeps(Level::Tight) == 3) (*currentChannel_) = Channel::DY_Tight;
        else (*currentChannel_) = Channel::DY_Fake;
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
    passCuts &= cuts.setCut("pass0LooseLep", nLeps(Level::Loose) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    passCuts &= cuts.setCut("passSSLeptons", isSameSign());
    passCuts &= cuts.setCut("passLowMassVeto", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));
    passCuts &= cuts.setCut("passZVeto", !muon.isInMassRange(Level::Loose) && !elec.isInMassRange(Level::Loose));
    // passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 50);
    // passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

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

bool Nonprompt_Closure::dy_closure_cuts()
{
    bool passCuts = true;
    CutInfo cuts;
    Lepton& tag_lep = subChannel_ == Subchannel::MM ? static_cast<Lepton&>(muon) : static_cast<Lepton&>(elec);
    Lepton& probe_lep = subChannel_ == Subchannel::MM ? static_cast<Lepton&>(elec) : static_cast<Lepton&>(muon);

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2TagLepton", tag_lep.size(Level::Tight) == 2);
    passCuts &= cuts.setCut("passProbeLepton", probe_lep.size(Level::Fake) == 1);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    int charge = 0;
    float mass = -1;
    if (subChannel_ == Subchannel::MM && tag_lep.size(Level::Tight) == 2) {
        mass = (muon.p4(Level::Tight, 0) + muon.p4(Level::Tight, 1)).M();
        charge = muon.charge(Level::Tight, 0) * muon.charge(Level::Tight, 1);
    } else if (subChannel_ == Subchannel::EE && tag_lep.size(Level::Tight) == 2) {
        mass = (elec.p4(Level::Tight, 0) + elec.p4(Level::Tight, 1)).M();
        charge = elec.charge(Level::Tight, 0) * elec.charge(Level::Tight, 1);
    }
    passCuts &= cuts.setCut("passZCut", mass > 70. && mass < 115);
    passCuts &= cuts.setCut("passOSTagLeptons", charge < 0);
    // passCuts &= cuts.setCut("passMetCut", met.pt() < 50);

    cuts.setCut("passFakeLeps", probe_lep.size(Level::Tight) == 0);
    fillCutFlow(Channel::DY_Fake, cuts);
    cuts.cuts.pop_back();
    cuts.setCut("pass1TightLeps", probe_lep.size(Level::Tight) == 1);
    fillCutFlow(Channel::DY_Tight, cuts);

    return passCuts;
}

float Nonprompt_Closure::getLeadPt()
{
    if (subChannel_ == Subchannel::MM) {
        return muon.rawpt(Level::Fake, 0);
    } else if (subChannel_ == Subchannel::EE) {
        return elec.rawpt(Level::Fake, 0);
    } else if(subChannel_ == Subchannel::EM){
        return std::max(muon.rawpt(Level::Fake, 0), elec.rawpt(Level::Fake, 0));
    }
    return 0.;
}

void Nonprompt_Closure::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    // muon.fillLepton_Iso(*o_looseMuons, jet, Level::LooseNotFake, event_bitmap);
    muon.fillLepton_Iso(*o_fakeMuons, jet, Level::FakeNotTight, event_bitmap);
    muon.fillLepton_Iso( *o_tightMuons, jet, Level::Tight, event_bitmap);
    // elec.fillLepton_Iso(*o_looseElectrons, jet, Level::LooseNotFake, event_bitmap);
    elec.fillLepton_Iso(*o_fakeElectrons, jet, Level::FakeNotTight, event_bitmap);
    elec.fillLepton_Iso( *o_tightElectrons, jet, Level::Tight, event_bitmap);
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
