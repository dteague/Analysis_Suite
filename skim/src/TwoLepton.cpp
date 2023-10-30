#include "analysis_suite/skim/interface/TwoLepton.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Measurement,
    None,
};

enum class Subchannel {
    MM,
    EM,
    EE,
    None,
};

void TwoLepton::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    createTree("Measurement", Channel::Measurement);

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);
    muon.setup_map(Level::LooseNotFake);
    elec.setup_map(Level::LooseNotFake);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }


    // Dilepton triggers
    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",});
        setupTrigger(Subchannel::EM, Dataset::MuonEG,
                     {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL"});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG,
                     {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",});
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::MM,  Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",});
        setupTrigger(Subchannel::EM, Dataset::MuonEG,
                     {"HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG,
                     {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",});
    } else if(year_ == Year::yr2017) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8"});
        setupTrigger(Subchannel::EM, Dataset::MuonEG,
                     {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG,
                     {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                      "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
    } else if (year_ == Year::yr2018) {
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon,
                     {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",});
        setupTrigger(Subchannel::EM, Dataset::MuonEG,
                     {"HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                      "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG,
                     {"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",});
    }

    LOG_FUNC << "End of Init";
}

void TwoLepton::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);

    tree->Branch("FakeMuon", "LeptonOut_Fake", &o_fakeMuons);
    tree->Branch("TightMuon", "LeptonOut_Fake", &o_tightMuons);

    tree->Branch("FakeElectron", "LeptonOut_Fake", &o_fakeElectrons);
    tree->Branch("TightElectron", "LeptonOut_Fake", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);

    tree->Branch("HT", &o_ht);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("N_loose_mu", &o_nloose_mu);
    tree->Branch("N_loose_el", &o_nloose_el);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void TwoLepton::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void TwoLepton::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_met.clear();
    o_metphi.clear();
    o_nloose_mu.clear();
    o_nloose_el.clear();
    LOG_FUNC << "End of clearOutputs";
}

void TwoLepton::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getLHESF();
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPrefire();
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPartonShower();
    // std::cout << (*weight) << " ";
    (*weight) *= jet.getPileupIDWeight();
    // std::cout << (*weight) << " ";
    (*weight) *= elec.getScaleFactor();
    // std::cout << (*weight) << " ";
    (*weight) *= muon.getScaleFactor();
    // std::cout << (*weight) << " " << "end" << std::endl;
    LOG_FUNC << "End of ApplyScaleFactors";
}

void TwoLepton::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    muon.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
    elec.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
    LOG_FUNC << "End of setOtherGoodParticles";
}

void TwoLepton::setSubChannel()
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

bool TwoLepton::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (measurement_cuts()) {
        (*currentChannel_) = Channel::Measurement;
    }

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool TwoLepton::measurement_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2FakeLep",  nLeps(Level::Fake) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    return passCuts;
}

float TwoLepton::getLeadPt()
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

void TwoLepton::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton_Iso(*o_fakeMuons, jet, Level::FakeNotTight, event_bitmap);
    muon.fillLepton_Iso( *o_tightMuons, jet, Level::Tight, event_bitmap);
    elec.fillLepton_Iso(*o_fakeElectrons, jet, Level::FakeNotTight, event_bitmap);
    elec.fillLepton_Iso( *o_tightElectrons, jet, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_nloose_mu.push_back(muon.size(Level::LooseNotFake));
        o_nloose_el.push_back(elec.size(Level::LooseNotFake));
    }
    LOG_FUNC << "End of FillValues";
}
