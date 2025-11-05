#include "analysis_suite/skim/interface/ThreeTop.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

// Uncomment and recompile if want to run preselection
// #define PRESELECTION

namespace Channel {
    enum {
        SS_Dilepton,
        SS_Multi,
        Nonprompt_Dilepton,
        Nonprompt_Multi,
        OS_MisId,
        None,
    };
}

void ThreeTop::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    DileptonBase::Init(tree);
    bool nonprompt_mc = (groupName_.find("ttbar") != std::string::npos
                         || groupName_.find("wjet") != std::string::npos
                         || groupName_.find("DY") != std::string::npos);

#ifdef PRESELECTION
    std::cout << "RUNNING PRESELECTION" << std::endl;
#endif
    createTree("Signal_Dilepton", Channel::SS_Dilepton);
    createTree("Signal_Multi", Channel::SS_Multi);

    if (!isMC_ || nonprompt_mc) {
        use_nonprompt = true;
        createTree("Nonprompt_Dilepton", Channel::Nonprompt_Dilepton);
        createTree("Nonprompt_Multi", Channel::Nonprompt_Multi);
    }
    if (!isMC_) { // Charge Mis-id Fake Rate
        createTree("OS_Charge_MisId", Channel::OS_MisId);
    }

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    LOG_FUNC << "End of Init";
}

void ThreeTop::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);

    tree->Branch("TightMuon", "LeptonOut", &o_tightMuons);
    tree->Branch("TightElectron", "LeptonOut", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("BJets", "JetOut", &o_bJets);

    if (isMC_) {
        tree->Branch("wgt_nobtag", &o_wgt_nobtag);
        tree->Branch("bwgt_cb_l", &o_bwgt_loose);
        tree->Branch("bwgt_cb_m", &o_bwgt_medium);
        tree->Branch("bwgt_cb_t", &o_bwgt_tight);
        jet.setup_shift_output(tree);
    }
    tree->Branch("hasVetoJet", &o_hasVetoJet);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
    tree->Branch("NBjets_loose", &o_nb_loose);
    tree->Branch("NBjets_medium", &o_nb_medium);
    tree->Branch("NBjets_tight", &o_nb_tight);

    tree->Branch("N_loose_mu", &o_nloose_mu);
    tree->Branch("N_loose_el", &o_nloose_el);
    tree->Branch("passZVeto", &o_pass_zveto);
    tree->Branch("Zmass", &o_ZMass);
    LOG_FUNC << "End of SetupOutTreeBranches";
}


void ThreeTop::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_ht.clear();
    o_htb.clear();
    o_met.clear();
    o_metphi.clear();
    o_centrality.clear();
    o_nloose_mu.clear();
    o_nloose_el.clear();
    o_pass_zveto.clear();
    o_ZMass.clear();

    o_nb_loose.clear();
    o_nb_medium.clear();
    o_nb_tight.clear();
    o_bwgt_loose.clear();
    o_bwgt_medium.clear();
    o_bwgt_tight.clear();
    o_wgt_nobtag.clear();

    LOG_FUNC << "End of clearOutputs";
}

void ThreeTop::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getScaleFactor();
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    (*weight) *= jet.getTotalBTagWeight("M");
    LOG_FUNC << "End of ApplyScaleFactors";
}

void ThreeTop::ApplyChannelScaleFactors()
{
    if (!use_nonprompt) {
        return;
    } else if ((*currentChannel_) == Channel::Nonprompt_Dilepton ||
               (*currentChannel_) == Channel::Nonprompt_Multi) {
        // Code for nonprompt scale factor
        // sign is to account for sign flip for every additional fake lepton in the event
        int sign = -1;
        auto muon_tight = muon.list(Level::Tight);
        for (auto fake: muon.list(Level::Fake)) {
            if (std::find(muon_tight.begin(), muon_tight.end(), fake) == muon_tight.end()) {
                (*weight) *= sfMaker.getNonpromptFR(muon.eta(fake), muon.pt(fake), PID::Muon);
                sign *= -1;
            }
        }
        auto elec_tight = elec.list(Level::Tight);
        for (auto fake: elec.list(Level::Fake)) {
            if (std::find(elec_tight.begin(), elec_tight.end(), fake) == elec_tight.end()) {
                (*weight) *= sfMaker.getNonpromptFR(elec.eta(fake), elec.pt(fake), PID::Electron);
                sign *= -1;
            }
        }
        (*weight) *= sign;
    } else if ((*currentChannel_) == Channel::OS_MisId) {
        // Code for applying charge-misid scale factor. Scales are added instead of multiplied
        float tmp_weight = 0;
        for (auto i : elec.list(Level::Tight)) {
            tmp_weight += sfMaker.getChargeMisIdFR(elec.eta(i), elec.pt(i));
        }
        (*weight) *= tmp_weight;
    }
}


bool ThreeTop::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (!baseline_cuts()) {
        return false;
    }

    if (signal_cuts()) {
        (*currentChannel_) = (nLeps(Level::Tight) == 2) ? Channel::SS_Dilepton : Channel::SS_Multi;
        return true;
    } else if(isMC_ && !use_nonprompt) {
        return false;
    } else if (nonprompt_cuts()) {
        (*currentChannel_) = (nLeps(Level::Fake) == 2) ? Channel::Nonprompt_Dilepton : Channel::Nonprompt_Multi;
        return true;
    } else if (!isMC_ && charge_misid_cuts()) {
        (*currentChannel_) = Channel::OS_MisId;
        return true;
    }

    return false;
}

bool ThreeTop::passPreselection()
{
    return (metfilters.pass()
            && trig_cuts.pass_any_trig()
            && (muon.size() + elec.size() >= 2)
            && (jet.size() >= 1)
            );
}

bool ThreeTop::baseline_cuts()
{
    LOG_FUNC << "Start of baseline_cuts";
    bool passCuts = true;

    passCuts &= getLeadPt() > 25;
    passCuts &= nLeps(Level::Fake) >= 2;
    passCuts &= getTriggerValue();
    passCuts &= jet.size(Level::Loose) >= 2;
    passCuts &= !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.);

    #ifdef PRESELECTION
    // Preselection
    passCuts &= met.pt() > 25;
    passCuts &= jet.getHT(Level::Loose) > 150;
    # else
    // Signal Region
    passCuts &= jet.size(Level::Bottom) >= 1;
    passCuts &= met.pt() > 50;
    passCuts &= jet.getHT(Level::Loose) > 250;
    #endif

    LOG_FUNC << "End of baseline_cuts";
    return passCuts;
}

bool ThreeTop::signal_cuts()
{
    LOG_FUNC << "Start of signal_cuts";
    bool passCuts = true;
    passCuts &= nLeps(Level::Tight) >= 2;
    passCuts &= nLeps(Level::Tight) == nLeps(Level::Fake);
    passCuts &= isSameSign(Level::Tight);

    LOG_FUNC << "End of signal_cuts";
    return passCuts;
}

bool ThreeTop::nonprompt_cuts()
{
    LOG_FUNC << "Start of nonprompt_cuts";
    bool passCuts = true;
    passCuts &= nLeps(Level::Fake) >= 2;
    passCuts &= nLeps(Level::Fake) != nLeps(Level::Tight);
    passCuts &= isSameSign(Level::Fake);

    LOG_FUNC << "End of nonprompt_cuts";
    return passCuts;
}

bool ThreeTop::charge_misid_cuts()
{
    LOG_FUNC << "Start of charge_misid_cuts";
    bool passCuts = true;
    passCuts &= nLeps(Level::Tight) == 2;
    passCuts &= nLeps(Level::Tight) == nLeps(Level::Fake);
    passCuts &= elec.size(Level::Tight) >= 1;
    passCuts &= !isSameSign(Level::Tight);
    passCuts &= !elec.isInMassRange(Level::Loose);

    LOG_FUNC << "End of charge_misid_cuts";
    return passCuts;
}



void ThreeTop::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton(*o_tightMuons, Level::Fake, event_bitmap);
    elec.fillLepton(*o_tightElectrons, Level::Fake, event_bitmap);

    o_hasVetoJet = jet.size(Level::Loose_inHEM) > 0;

    jet.fillJet(*o_jets, Level::Loose, event_bitmap);
    jet.fillJet(*o_bJets, Level::Bottom, event_bitmap);
    float muonMass = muon.massInRange(Level::Loose);
    float elecMass = elec.massInRange(Level::Loose);
    float zmass = (muonMass < 0) ? elecMass : muonMass;

    for (size_t systNum = 0; systNum < numSystematics(); ++systNum) {
        size_t syst = syst_to_index.at(systNum);
        if (syst == 0 && systNum != 0) {
            continue;
        }
        setupSyst(systNum);

        o_ht.push_back(jet.getHT(Level::Loose, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Loose, syst));

        o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
        o_nb_medium.push_back(jet.n_medium_bjet.at(syst));
        o_nb_tight.push_back(jet.n_tight_bjet.at(syst));

        if (isMC_) {
            o_wgt_nobtag.push_back(o_weight[systNum]/jet.getTotalBTagWeight("M"));
            o_bwgt_loose.push_back(jet.getCutBasedBTagWeight("L"));
            o_bwgt_medium.push_back(jet.getCutBasedBTagWeight("M"));
            o_bwgt_tight.push_back(jet.getCutBasedBTagWeight("T"));
        }

        o_nloose_mu.push_back(muon.size(Level::Loose)-muon.size(Level::Fake));
        o_nloose_el.push_back(elec.size(Level::Loose)-elec.size(Level::Fake));
        o_pass_zveto.push_back(zmass < 0);
        o_ZMass.push_back(zmass);
    }
    LOG_FUNC << "End of FillValues";
}

void ThreeTop::printStuff()
{
    LOG_FUNC << "Start of printStuff";
    std::cout << "Event: " << *event << std::endl;
    std::cout << "Met: " << met.pt() << std::endl;
    std::cout << "HT: " << jet.getHT(Level::Loose, 0) << std::endl;
    std::cout << "njet: " << jet.size(Level::Loose) << std::endl;
    std::cout << "nbjet: " << jet.size(Level::Bottom) << std::endl;
    std::cout << "nlep: " << muon.size(Level::Tight) << " " << elec.size(Level::Tight) << std::endl;
    std::cout << "nlep loose: " << muon.size(Level::Loose) << " " << elec.size(Level::Loose) << std::endl;
    std::cout << std::endl;
    LOG_FUNC << "End of printStuff";
}
