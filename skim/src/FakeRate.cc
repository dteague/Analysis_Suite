#include "analysis_suite/skim/interface/FakeRate.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Measurement,
    SideBand,
    None,
};

enum class Subchannel {
    M,
    E,
    None,
};

void FakeRate::Init(TTree* tree)
{
    met_type = MET_Type::PUPPI;
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Measurement", Channel::Measurement);
    createTree("SideBand", Channel::SideBand);

    muon.setup_map(Level::FakeNotTight);
    elec.setup_map(Level::FakeNotTight);
    muon.setup_map(Level::LooseNotFake);
    elec.setup_map(Level::LooseNotFake);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }


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


    LOG_FUNC << "End of Init";
}

void FakeRate::SetupOutTreeBranches(TTree* tree)
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

    tree->Branch("HT", &o_ht);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("bjet_scale", &bjet_scales);
    tree->Branch("hlt_loPt_prescale", &hlt_loPt_prescale);
    tree->Branch("hlt_hiPt_prescale", &hlt_hiPt_prescale);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void FakeRate::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void FakeRate::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_met.clear();
    o_metphi.clear();
    bjet_scales.clear();
    LOG_FUNC << "End of clearOutputs";
}

void FakeRate::ApplyScaleFactors()
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

void FakeRate::ApplyDataSpecifics()
{
    if ((*currentChannel_) == Channel::Measurement) {
        size_t trigIdx = 0;
        if (subChannel_ == Subchannel::E && elec.rawpt(Level::Fake, 0) > 25) trigIdx = 1;
        else if (subChannel_ == Subchannel::M && muon.rawpt(Level::Fake, 0) > 20) trigIdx = 1;

        (*weight) *= sfMaker.getPrescale(*run, *lumiblock, trig_cuts.trig_name(subChannel_, trigIdx));
    }
}

void FakeRate::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    muon.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    elec.xorLevel(Level::Fake, Level::Tight, Level::FakeNotTight);
    muon.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
    elec.xorLevel(Level::Loose, Level::Fake, Level::LooseNotFake);
    LOG_FUNC << "End of setOtherGoodParticles";
}

void FakeRate::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if (nLeps(Level::Fake) == 1) {
        subChannel_ = muon.size(Level::Fake) == 1 ? Subchannel::M : Subchannel::E;
    }
    LOG_FUNC << "End of setSubChannel";
}

bool FakeRate::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (measurement_cuts()) (*currentChannel_) = Channel::Measurement;

    if (sideband_cuts()) (*currentChannel_) = Channel::SideBand;

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    if (!isMC_) {
        ApplyDataSpecifics();
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool FakeRate::single_lep_cuts(CutInfo& cuts)
{
    bool passCuts = true;
    bool haveOneFake = nLeps(Level::Fake) == 1;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass1FakeLep", haveOneFake);
    passCuts &= cuts.setCut("pass0LooseLep", nLeps(Level::Loose) == 1);
    // TriggerCuts

    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));
    passCuts &= cuts.setCut("pass1Jet", jet.size(Level::Tight) == 1);

    bool hasFarJet = false;
    if (haveOneFake && jet.size(Level::Tight) == 1) {
        if (subChannel_ == Subchannel::M) {
            hasFarJet = deltaR(muon.eta(Level::Fake, 0), jet.eta(Level::Tight, 0),
                               muon.phi(Level::Fake, 0), jet.phi(Level::Tight, 0)) > 1;
        } else {
            hasFarJet = deltaR(elec.eta(Level::Fake, 0), jet.eta(Level::Tight, 0),
                               elec.phi(Level::Fake, 0), jet.phi(Level::Tight, 0)) > 1;
        }
    }
    passCuts &= cuts.setCut("passHasFarJet", hasFarJet);

    return passCuts;
}

bool FakeRate::measurement_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= single_lep_cuts(cuts);
    // Split the trigger cuts
    bool trig_match = false;
    if (subChannel_ != Subchannel::None) {
        size_t trig_idx = 0;
        if (subChannel_ == Subchannel::E && elec.rawpt(Level::Fake, 0) > 25) trig_idx = 1;
        else if (subChannel_ == Subchannel::M && muon.rawpt(Level::Fake, 0) > 20) trig_idx = 1;
        trig_match = trig_cuts.pass_cut(subChannel_, trig_idx);
    }
    passCuts &= cuts.setCut("passSplitTrigger", trig_match);
    passCuts &= cuts.setCut("passMetCut", met.pt() < 30);

    // Fill Cut flow
    fillCutFlow(Channel::Measurement, cuts);

    return passCuts;
}

bool FakeRate::sideband_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= single_lep_cuts(cuts);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 30);
    if (groupName_.find("qcd") == std::string::npos) {
        passCuts &= cuts.setCut("passTightLep", nLeps(Level::Tight) == 1);
    }

    // Fill Cut flow
    fillCutFlow(Channel::SideBand, cuts);

    return passCuts;
}

void FakeRate::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton_Iso(*o_fakeMuons, jet, Level::FakeNotTight, event_bitmap);
    muon.fillLepton_Iso( *o_tightMuons, jet, Level::Tight, event_bitmap);
    elec.fillLepton_Iso(*o_fakeElectrons, jet, Level::FakeNotTight, event_bitmap);
    elec.fillLepton_Iso( *o_tightElectrons, jet, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    if (isMC_) {
        bjet_scales.push_back(jet.getBJetWeight(jet.idx(Level::Tight, 0), "L"));
        bjet_scales.push_back(jet.getBJetWeight(jet.idx(Level::Tight, 0), "M"));
        bjet_scales.push_back(jet.getBJetWeight(jet.idx(Level::Tight, 0), "T"));
        hlt_loPt_prescale = static_cast<int>(trig_cuts.pass_cut(subChannel_, 0));
        hlt_hiPt_prescale = static_cast<int>(trig_cuts.pass_cut(subChannel_, 1));
    } else {
        bjet_scales.push_back(1.);
        bjet_scales.push_back(1.);
        bjet_scales.push_back(1.);
        hlt_loPt_prescale = sfMaker.getPrescale(*run, *lumiblock, trig_cuts.trig_name(subChannel_, 0)) * static_cast<int>(trig_cuts.pass_cut(subChannel_, 0));
        hlt_hiPt_prescale = sfMaker.getPrescale(*run, *lumiblock, trig_cuts.trig_name(subChannel_, 1)) * static_cast<int>(trig_cuts.pass_cut(subChannel_, 1));
    }


    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);

        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());

    }
    LOG_FUNC << "End of FillValues";
}
