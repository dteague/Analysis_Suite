#include "analysis_suite/skim/interface/Lepton_MisId_Closure.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

namespace Channel {
    enum {
        OS_MR,
        OS_CR,
        SS_CR,
        None,
    };
}

void Closure_MisId::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    DileptonBase::Init(tree);

    // Charge Mis-id Fake Rate
    createTree("SS", Channel::SS_CR);
    createTree("OS_CR", Channel::OS_CR);
    createTree("OS_MR", Channel::OS_MR);
    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    LOG_FUNC << "End of Init";
}

void Closure_MisId::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("TightMuon", "LeptonOut", &o_tightMuons);
    tree->Branch("TightElectron", "LeptonOut", &o_tightElectrons);
    tree->Branch("Jets", "JetOut", &o_jets);

    tree->Branch("HT", &o_ht);
    tree->Branch("HT_b", &o_htb);
    tree->Branch("Met", &o_met);
    tree->Branch("Met_phi", &o_metphi);
    tree->Branch("Centrality", &o_centrality);
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
    (*weight) *= sfMaker.getLHESF();
    (*weight) *= sfMaker.getLHEPdf();
    (*weight) *= sfMaker.getPrefire();
    (*weight) *= sfMaker.getPartonShower();
    (*weight) *= sfMaker.getTriggerSF(elec, muon);
    (*weight) *= jet.getScaleFactor();
    (*weight) *= jet.getTotalBTagWeight("M");
    (*weight) *= elec.getScaleFactor();
    (*weight) *= muon.getScaleFactor();
    LOG_FUNC << "End of ApplyScaleFactors";
}

bool Closure_MisId::getCutFlow()
{
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (closure_cuts()) {
        (*currentChannel_) = (isSameSign(Level::Tight)) ? Channel::SS_CR : Channel::OS_CR;
    }
    if(measurement_cuts()) {
        (*currentChannel_) = Channel::OS_MR;
    }

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
    passCuts &= cuts.setCut("pass2LooseLep", nLeps(Level::Loose) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", getTriggerValue());

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
    passCuts &= cuts.setCut("pass2LooseLep", nLeps(Level::Loose) == 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", getTriggerValue());

    passCuts &= cuts.setCut("passZCut", get_mass() > 50);
    // passCuts &= cuts.setCut("passOppositeSign", !isSameSign());
    // passCuts &= cuts.setCut("passHasElectron", elec.size(Level::Tight) > 0);
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    // passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

    fillCutFlow(Channel::OS_MR, cuts);

    return passCuts;
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

void Closure_MisId::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    muon.fillLepton(*o_tightMuons, Level::Tight, event_bitmap);
    elec.fillLepton(*o_tightElectrons, Level::Tight, event_bitmap);
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);

    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        setupSyst(syst);
        o_ht.push_back(jet.getHT(Level::Tight, syst));
        o_htb.push_back(jet.getHT(Level::Bottom, syst));
        o_met.push_back(met.pt());
        o_metphi.push_back(met.phi());
        o_centrality.push_back(jet.getCentrality(Level::Tight, syst));
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
    std::cout << std::endl;
    LOG_FUNC << "End of printStuff";
}
