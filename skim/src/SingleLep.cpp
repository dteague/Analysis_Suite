#include "analysis_suite/skim/interface/SingleLep.h"

#include "analysis_suite/skim/interface/logging.h"

namespace Channel {
    enum  {
        Measurement,
        None,
    };
}

void SingleLep::Init(TTree* tree)
{
    met_type = MET_Type::PUPPI;
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Measurement", Channel::Measurement);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }

    // Single Lepton Triggers
    setupTrigger(Subchannel::Single_M, Dataset::DoubleMuon,
                 {"HLT_Mu8_TrkIsoVVL",
                  "HLT_Mu17_TrkIsoVVL",
                 });
    setupTrigger(Subchannel::Single_E, Dataset::DoubleEG,
                 {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", // Was 12, changed to 8
                  "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", // was 23, changed to 17
                 });

    muon.mvaCut = 2;
    elec.mvaCut = 2;
    LOG_FUNC << "End of Init";
}

void SingleLep::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("FakeMuon", "LeptonOut_small", &o_fakeMuons);
    tree->Branch("FakeElectron", "LeptonOut_small", &o_fakeElectrons);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void SingleLep::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    BaseSelector::clearParticles();
    LOG_FUNC << "End of clearParticles";
}

void SingleLep::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_met.clear();
    o_metphi.clear();
    LOG_FUNC << "End of clearOutputs";
}

void SingleLep::ApplyScaleFactors()
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

void SingleLep::setOtherGoodParticles(size_t syst)
{
    LOG_FUNC << "Start of setOtherGoodParticles";
    LOG_FUNC << "End of setOtherGoodParticles";
}

void SingleLep::setSubChannel()
{
    LOG_FUNC << "Start of setSubChannel";
    subChannel_ = Subchannel::None;

    if (nLeps(Level::Fake) == 1) {
        subChannel_ = muon.size(Level::Fake) == 1 ? Subchannel::Single_M : Subchannel::Single_E;
    }
    LOG_FUNC << "End of setSubChannel";
}


bool SingleLep::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (measurement_cuts()) (*currentChannel_) = Channel::Measurement;

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool SingleLep::measurement_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass1FakeLep", nLeps(Level::Fake) == 1);
    passCuts &= cuts.setCut("pass0LooseLep", nLeps(Level::Loose) == 1);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    fillCutFlow(Channel::Measurement, cuts);

    return passCuts;
}

void SingleLep::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillMuon_small(*o_fakeMuons, Level::Fake, event_bitmap);
    elec.fillElectron_small(*o_fakeElectrons, Level::Fake, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);
    }
    LOG_FUNC << "End of FillValues";
}
