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

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    } else {
        sfMaker.setup_prescale();
    }


    // Single Lepton Triggers
    setupTrigger(Subchannel::M, {"HLT_Mu8_TrkIsoVVL",
                                 "HLT_Mu17_TrkIsoVVL",
                                 // "HLT_Mu8",
                                 // "HLT_Mu17"
        });
    setupTrigger(Subchannel::E, {"HLT_Ele8_CaloIdL_TrackIdL_IsoVL_PFJet30", // Was 12, changed to 8
                                 "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", // was 23, changed to 17
                                 // "HLT_Ele8_CaloIdM_TrackIdM_PFJet30",
                                 // "HLT_Ele17_CaloIdM_TrackIdM_PFJet30"
        });

    setupTrigger(Subchannel::None);

    LOG_FUNC << "End of Init";
}

void FakeRate::SetupOutTreeBranches(TTree* tree)
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
    o_htb.clear();
    o_met.clear();
    o_metphi.clear();
    o_nb_loose.clear();
    o_nb_tight.clear();
    LOG_FUNC << "End of clearOutputs";
}

void FakeRate::ApplyScaleFactors()
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

void FakeRate::ApplyDataSpecifics()
{
    if ((*currentChannel_) == Channel::Measurement || (*currentChannel_) == Channel::SideBand) {
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
    set_leadlep();

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
    // TriggerCuts
    bool trig_match = false;
    if (haveOneFake) {
        size_t trig_idx = 0;
        if (subChannel_ == Subchannel::E && elec.rawpt(Level::Fake, 0) > 25) trig_idx = 1;
        else if (subChannel_ == Subchannel::M && muon.rawpt(Level::Fake, 0) > 20) trig_idx = 1;
        trig_match = trig_cuts.pass_cut(subChannel_, trig_idx); // || trig_cuts.pass_cut(subChannel_, trig_idx+2);
    }
    passCuts &= cuts.setCut("passLeadPtCut", trig_match);
    passCuts &= cuts.setCut("passTrigger", trig_cuts.pass_cut(subChannel_));

    bool hasFarJet = false;
    if (haveOneFake) {
        for (auto jidx: jet.list(Level::Tight)) {
            if (deltaR(lead_lep.Eta(), jet.eta(jidx), lead_lep.Phi(), jet.phi(jidx)) > 1) {
                hasFarJet = true;
                break;
            }
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
    passCuts &= cuts.setCut("passMetCut", met.pt() < 20);
    passCuts &= cuts.setCut("passMtCut", met.mt(lead_lep.Pt(), lead_lep.Phi()) < 20);
    // Fill Cut flow
    fillCutFlow(Channel::Measurement, cuts);

    return passCuts;
}

bool FakeRate::sideband_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= single_lep_cuts(cuts);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 20);
    passCuts &= cuts.setCut("passTightLep", nLeps(Level::Fake) == 1);
    passCuts &= cuts.setCut("passTightLep", nLeps(Level::Tight) == 1);
    passCuts &= cuts.setCut("passLeadLepPt", lead_lep.Pt() > 20);

    // Fill Cut flow
    fillCutFlow(Channel::SideBand, cuts);

    return passCuts;
}

void FakeRate::set_leadlep()
{
    if (nLeps(Level::Fake) != 1) {
        lead_lep *= 0;
        return;
    } else if (muon.size(Level::Fake) == 1) {
        lead_lep = muon.p4(Level::Fake, 0);
    } else {
        lead_lep = elec.p4(Level::Fake, 0);
    }
}

void FakeRate::FillValues(const Bitmap& event_bitmap)
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
