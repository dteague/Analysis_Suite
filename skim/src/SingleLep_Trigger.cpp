#include "analysis_suite/skim/interface/SingleLep_Trigger.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Measurement,
    None,
};

enum class Subchannel {
    E, M,
    None,
};

void SingleLep_Trigger::Init(TTree* tree)
{
    met_type = MET_Type::PUPPI;
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Measurement", Channel::Measurement);

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }
    hlt_lo_mu.setup(fReader, "HLT_Mu8_TrkIsoVVL");
    hlt_hi_mu.setup(fReader, "HLT_Mu17_TrkIsoVVL");
    hlt_lo_el.setup(fReader, "HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30");
    hlt_hi_el.setup(fReader, "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30");


    // Single Lepton Triggers
    setupTrigger(Subchannel::M, Dataset::DoubleMuon,
                 {"HLT_Mu8_TrkIsoVVL",
                  "HLT_Mu17_TrkIsoVVL",
                  // "HLT_Mu8",
                  // "HLT_Mu17"
                 });
    setupTrigger(Subchannel::E, Dataset::DoubleEG,
                 {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", // Was 12, changed to 8
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", // was 23, changed to 17
                  // "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
                  // "HLT_Ele17_CaloIdM_TrackIdM_PFJet30"
                 });

    // // Single Lepton Triggers
    // setupTrigger(Subchannel::Met, Dataset::Met,
    //              {"HLT_PFMETNoMu120_PFMHTNoMu120_IDTight",});

    LOG_FUNC << "End of Init";
}

void SingleLep_Trigger::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("FakeMuon", "LeptonOut_Fake", &o_fakeMuons);
    tree->Branch("TightMuon", "LeptonOut_Fake", &o_tightMuons);
    tree->Branch("FakeElectron", "ElectronOut_Fake", &o_fakeElectrons);
    tree->Branch("TightElectron", "ElectronOut_Fake", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);

    tree->Branch("HT", &o_ht);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("hlt_low_mu", &o_hlt_lo_mu);
    tree->Branch("hlt_high_mu", &o_hlt_hi_mu);
    tree->Branch("hlt_low_el", &o_hlt_lo_el);
    tree->Branch("hlt_high_el", &o_hlt_hi_el);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void SingleLep_Trigger::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void SingleLep_Trigger::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_met.clear();
    o_metphi.clear();
    LOG_FUNC << "End of clearOutputs";
}

void SingleLep_Trigger::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= jet.getScaleFactor();
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

void SingleLep_Trigger::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    LOG_FUNC << "End of setOtherGoodParticles";
}

void SingleLep_Trigger::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if (nLeps(Level::Fake) == 1) {
        subChannel_ = muon.size(Level::Fake) == 1 ? Subchannel::M : Subchannel::E;
    }
    LOG_FUNC << "End of setSubChannel";
}


bool SingleLep_Trigger::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;

    if (measurement_cuts()) (*currentChannel_) = Channel::Measurement;

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool SingleLep_Trigger::measurement_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass1FakeLep", nLeps(Level::Fake) == 1);
    passCuts &= cuts.setCut("pass0LooseLep", nLeps(Level::Loose) == 1);
    // passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    // auto met_vec = std::polar(met.pt(), met.phi());
    // float ht = jet.getHT(Level::Tight);
    // for (auto i : muon.list(Level::Fake)) {
    //     met_vec += std::polar(muon.rawpt(i), muon.phi(i));
    //     ht -= muon.rawpt(i);
    // }
    // passCuts &= cuts.setCut("passMetCut", std::abs(met_vec) > 110);
    // passCuts &= cuts.setCut("passHTCut", ht > 110);

    fillCutFlow(Channel::Measurement, cuts);

    return passCuts;
}

void SingleLep_Trigger::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton_Iso(*o_fakeMuons, jet, Level::FakeNotTight, event_bitmap);
    muon.fillLepton_Iso( *o_tightMuons, jet, Level::Tight, event_bitmap);
    elec.fillElectron_Iso(*o_fakeElectrons, jet, Level::FakeNotTight, event_bitmap);
    elec.fillElectron_Iso( *o_tightElectrons, jet, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    o_hlt_lo_mu = *trig_cuts.trigs[Subchannel::M].at(0);
    o_hlt_hi_mu = *trig_cuts.trigs[Subchannel::M].at(1);
    o_hlt_lo_el = *trig_cuts.trigs[Subchannel::E].at(0);
    o_hlt_hi_el = *trig_cuts.trigs[Subchannel::E].at(1);
    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
    }
    LOG_FUNC << "End of FillValues";
}
