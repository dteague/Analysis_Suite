#include "analysis_suite/skim/interface/FakeRate.h"

#include "analysis_suite/skim/interface/logging.h"

namespace Channel {
    enum {
        Measurement,
        SideBand,
        None,
    };
}

void FakeRate::Init(TTree* tree)
{
    met_type = MET_Type::PF;
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

    auto egamma_dataset = (year_ == Year::yr2017) ? Dataset::Single_E : Dataset::DoubleEG;
    highpt_e_dataset = (year_ == Year::yr2018) ? Dataset::DoubleEG : Dataset::Single_E;

    // Single Lepton Triggers
    std::vector<std::string> mu_list = {"HLT_Mu8_TrkIsoVVL",
                                        "HLT_Mu17_TrkIsoVVL",};
    std::vector<std::string> el_list = {"HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30",
                                        "HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30",};

    setupTrigger(Subchannel::Single_M, Dataset::DoubleMuon, mu_list);
    setupTrigger(Subchannel::Single_E, egamma_dataset, el_list);

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

    tree->Branch("HT", &o_ht);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    // tree->Branch("bjet_scale", &bjet_scales);
    if (isMC_) {
        tree->Branch("wgt_nobtag", &o_wgt_nobtag);
        tree->Branch("bwgt_cb_l", &o_bwgt_loose);
        tree->Branch("bwgt_cb_m", &o_bwgt_medium);
        tree->Branch("bwgt_cb_t", &o_bwgt_tight);
    }
    tree->Branch("N_bloose", &o_nb_loose);
    tree->Branch("N_bmedium", &o_nb_medium);
    tree->Branch("N_btight", &o_nb_tight);
    LOG_FUNC << "End of SetupOutTreeBranches";
}


void FakeRate::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_met.clear();
    o_metphi.clear();
    o_nb_loose.clear();
    o_nb_medium.clear();
    o_nb_tight.clear();
    o_bwgt_loose.clear();
    o_bwgt_medium.clear();
    o_bwgt_tight.clear();
    o_wgt_nobtag.clear();
    // bjet_scales.clear();
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
    (*weight) *= jet.getTotalBTagWeight("M");
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

size_t FakeRate::get_trigger()
{
    if (subChannel_ == Subchannel::Single_E) {
        float pt = elec.rawpt(Level::Fake, 0);
        if (pt > 40) return 1;
        else if (pt > 25) return 1;
        else return 0;
    } else {
        float pt = muon.rawpt(Level::Fake, 0);
        if (pt > 30) return 1;
        else if (pt > 20) return 1;
        else return 0;
    }
}

bool FakeRate::apply_trigger() {
    if (subChannel_ == Subchannel::None) {
        return false;
    }
    size_t trig_idx = get_trigger();
    if (subChannel_ == Subchannel::Single_M && trig_idx == 2) {
        return ((Dataset::Single_M == dataset_ ||
                 dataset_ == Dataset::None)
                && *trig_cuts.trigs[subChannel_].at(trig_idx));
    } else if (subChannel_ == Subchannel::Single_E && trig_idx == 2) {
        return ((highpt_e_dataset == dataset_ ||
                 dataset_ == Dataset::None)
                && *trig_cuts.trigs[subChannel_].at(trig_idx));
    } else {
        return trig_cuts.pass_cut(subChannel_, trig_idx);
    }
}

void FakeRate::ApplyDataSpecifics()
{
    if ((*currentChannel_) != Channel::None) {
        size_t trigIdx = get_trigger();
        if (trigIdx > 1)
            return;
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
        subChannel_ = muon.size(Level::Fake) == 1 ? Subchannel::Single_M : Subchannel::Single_E;
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
    passCuts &= cuts.setCut("passTrigger", apply_trigger());
    passCuts &= cuts.setCut("pass1Jet", jet.size(Level::Loose) >= 1);

    bool hasFarJet = false;
    if (haveOneFake && subChannel_ == Subchannel::Single_E) {
        float e_eta = elec.eta(Level::Fake, 0);
        float e_phi = elec.phi(Level::Fake, 0);
        for (auto i : jet.list(Level::Loose)) {
            hasFarJet |= deltaR2(e_eta, jet.eta(i), e_phi, jet.phi(i)) > 1;
        }
    } else {
        hasFarJet = true;
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
    muon.fillLepton_Iso(*o_tightMuons, jet, Level::Tight, event_bitmap);
    elec.fillLepton_Iso(*o_fakeElectrons, jet, Level::FakeNotTight, event_bitmap);
    elec.fillLepton_Iso(*o_tightElectrons, jet, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Loose, event_bitmap);

    for (size_t systNum = 0; systNum < numSystematics(); ++systNum) {
        size_t syst = syst_to_index.at(systNum);
        if (syst == 0 && systNum != 0) {
            continue;
        }
        setupSyst(systNum);

        o_ht.push_back(jet.getHT(Level::Loose, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        if (isMC_) {
            o_wgt_nobtag.push_back(o_weight[systNum]/jet.getTotalBTagWeight("M"));
            o_bwgt_loose.push_back(jet.getCutBasedBTagWeight("L"));
            o_bwgt_medium.push_back(jet.getCutBasedBTagWeight("M"));
            o_bwgt_tight.push_back(jet.getCutBasedBTagWeight("T"));
        }
        o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
        o_nb_medium.push_back(jet.n_medium_bjet.at(syst));
        o_nb_tight.push_back(jet.n_tight_bjet.at(syst));
    }
    LOG_FUNC << "End of FillValues";
}
