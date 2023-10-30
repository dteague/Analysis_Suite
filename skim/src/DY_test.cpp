#include "analysis_suite/skim/interface/DY_test.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

enum class Channel {
    OS_FF,
    OS_TF,
    OS_TT,
    None,
};

enum class Subchannel {
    MM,
    None,
};

void DY_test::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    // Charge Mis-id Fake Rate
    createTree("OS_FF", Channel::OS_FF);
    createTree("OS_TF", Channel::OS_TF);
    createTree("OS_TT", Channel::OS_TT);
    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",});
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::MM,  Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",});
    } else if(year_ == Year::yr2017) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"});
    } else if (year_ == Year::yr2018) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",});
    }

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);

    LOG_FUNC << "End of Init";
}

void DY_test::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "LeptonOut_Fake", &o_tightMuons);
    tree->Branch("FakeMuon", "LeptonOut_Fake", &o_fakeMuons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("Nloose_Muon", &o_nlooseMu);
    tree->Branch("Nloose_Electron", &o_nlooseEl);

    tree->Branch("HT", &o_ht);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void DY_test::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void DY_test::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_nlooseEl.clear();
    o_nlooseMu.clear();
    o_met.clear();
    o_metphi.clear();
    o_centrality.clear();
    LOG_FUNC << "End of clearOutputs";
}

void DY_test::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= jet.getScaleFactor();
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void DY_test::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    LOG_FUNC << "End of setOtherGoodParticles";
}


bool DY_test::isSameSign()
{
    int q_total = 0;
    for (size_t idx : muon.list(Level::Fake)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(Level::Fake)) {
        q_total += elec.charge(idx);
    }
    // if 2 leptons, SS -> +1 +1 / -1 -1 -> abs(q) == 2
    // if 3 leptons, SS -> +1 +1 -1 / -1 -1 +1 -> abs(q) == 1
    // OS cases are 0 and 3, so no overlap
    return abs(q_total) == 1 || abs(q_total) == 2;
}

void DY_test::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Fake) == 2) {
        if (muon.size(Level::Fake) == 2) {
            subChannel_ = Subchannel::MM;
        }
    }
    LOG_FUNC << "End of setSubChannel";
}



bool DY_test::getCutFlow()
{
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (closure_cuts()) {
        if (muon.size(Level::Tight) == 0) {
            (*currentChannel_) = Channel::OS_FF;
        } else if (muon.size(Level::Tight) == 1) {
            (*currentChannel_) = Channel::OS_TF;
        } else {
            (*currentChannel_) = Channel::OS_TT;
        }

    }

    if (*currentChannel_ == Channel::None) {
        return false; // Pass no channels
    }

    return true;
}

bool DY_test::closure_cuts() {
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2Muon", muon.size(Level::Fake) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    float mass = get_mass();
    passCuts &= cuts.setCut("passZCut", mass > 70. && mass < 115);
    // passCuts &= cuts.setCut("passMetCut", met.pt() < 50);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) < 250);

    int charge = 0;
    if (muon.size(Level::Fake) == 2) {
        charge = muon.charge(Level::Fake, 0) * muon.charge(Level::Fake, 1);
    }

    cuts.setCut("passOS", charge < 0);
    fillCutFlow(Channel::OS_TT, cuts);

    return passCuts;
}

float DY_test::getLeadPt()
{
    if (subChannel_ == Subchannel::MM) {
        return muon.pt(Level::Fake, 0);
    }
    return 0.;
}

float DY_test::get_mass() {
    if (subChannel_ == Subchannel::MM) {
        return (muon.p4(Level::Fake, 0) + muon.p4(Level::Fake, 1)).M();
    } else {
        return -1;
    }

}

void DY_test::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton_Iso(*o_tightMuons, jet, Level::Tight, event_bitmap);
    muon.fillLepton_Iso(*o_fakeMuons, jet, Level::FakeNotTight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);
        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Tight, syst));
        o_nlooseMu.push_back(muon.size(Level::Loose));
        o_nlooseEl.push_back(elec.size(Level::Loose));

    }
    LOG_FUNC << "End of FillValues";
}

void DY_test::printStuff()
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
