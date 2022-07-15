#include "analysis_suite/Analyzer/interface/ZZAnalysis.h"

#include"analysis_suite/Analyzer/interface/logging.h"
#include"analysis_suite/Analyzer/interface/CommonFuncs.h"

enum class Channel {
    Signal,
    None,
};

enum class Subchannel {
    EEEE,
    EEMM,
    MMMM,
    None,
};

void ZZAnalysis::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Signal", Channel::Signal);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    // setupTrigger(Subchannel::MM, {"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
    //                               "HLT_DoubleMu8_Mass8_PFHT300"});
    setupTrigger(Subchannel::None);

    LOG_FUNC << "End of Init";
}

void ZZAnalysis::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    // tree->Branch("LooseMuon", "ParticleOut", &o_looseMuons);
    // tree->Branch("LooseElectron", "ParticleOut", &o_looseElectrons);
    // tree->Branch("TightMuon", "LeptonOut", &o_tightMuons);
    // tree->Branch("TightElectron", "LeptonOut", &o_tightElectrons);
    // tree->Branch("TightLeptons", "ParticleOut", &o_tightLeptons);
    // tree->Branch("Jets", "JetOut", &o_jets);
    // tree->Branch("BJets", "JetOut", &o_bJets);
    // tree->Branch("ResolvedTops", "TopOut", &o_resolvedTop);

    // tree->Branch("HT", &o_ht);
    // tree->Branch("HT_b", &o_htb);
    // tree->Branch("Met", &o_met);
    // tree->Branch("Met_phi", &o_metphi);
    // tree->Branch("Centrality", &o_centrality);
    // tree->Branch("N_bloose", &o_nb_loose);
    // tree->Branch("N_btight", &o_nb_tight);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void ZZAnalysis::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    // o_ht.clear();
    // o_htb.clear();
    // o_met.clear();
    // o_metphi.clear();
    // o_centrality.clear();
    // o_nb_loose.clear();
    // o_nb_tight.clear();
    LOG_FUNC << "End of clearOutputs";
}

void ZZAnalysis::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    // LOG_EVENT << "weight: " << (*weight);
    // (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    // LOG_EVENT << "weight after pu scale: " << (*weight);
    // (*weight) *= sfMaker.getLHESF();
    // LOG_EVENT << "weight after lhe scale: " << (*weight);
    // (*weight) *= jet.getScaleFactor();
    // LOG_EVENT << "weight after jet scale: " << (*weight);
    // (*weight) *= elec.getScaleFactor();
    // LOG_EVENT << "weight after elec scale: " << (*weight);
    // (*weight) *= muon.getScaleFactor();
    // LOG_EVENT << "weight after muon scale: " << (*weight);
    // (*weight) *= rTop.getScaleFactor(rGen);
    LOG_FUNC << "End of ApplyScaleFactors";
}

bool ZZAnalysis::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (signal_cuts()) {
        (*currentChannel_) = Channel::Signal;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool ZZAnalysis::signal_cuts()
{
    LOG_FUNC << "Start of signal_cuts";
    bool passCuts = true;
    CutInfo cuts;

    // passCuts = baseline_cuts(cuts);
    // passCuts &= cuts.setCut("passZVeto", muon.passZVeto() && elec.passZVeto());
    // passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) == 2 || nLeps(Level::Tight) == 3);
    // passCuts &= cuts.setCut("passSameSign;", isSameSign(Level::Tight));

    // // Fill Cut flow
    // cuts.setCut("pass2TightLeps", nLeps(Level::Tight) == 2);
    // fillCutFlow(Channel::SS_Dilepton, cuts);
    // cuts.cuts.pop_back();
    // cuts.setCut("pass3TightLeps", nLeps(Level::Tight) == 3);
    // fillCutFlow(Channel::SS_Multi, cuts);

    LOG_FUNC << "End of signal_cuts";
    return passCuts;
}

void ZZAnalysis::FillValues(const std::vector<bool>& passVec)
{
    LOG_FUNC << "Start of FillValues";
    BaseSelector::FillValues(passVec);
    size_t pass_bitmap = 0;
    for (size_t i = 0; i < passVec.size(); ++i) {
        pass_bitmap += passVec.at(i) << i;
    }

    // fillParticle(muon, Level::Loose, *o_looseMuons, pass_bitmap);
    // fillParticle(elec, Level::Loose, *o_looseElectrons, pass_bitmap);
    // fillLepton(muon, Level::Tight, *o_tightMuons, pass_bitmap);
    // fillLepton(elec, Level::Tight, *o_tightElectrons, pass_bitmap);
    // fillJet(jet, Level::Tight, *o_jets, pass_bitmap);
    // fillJet(jet, Level::Bottom, *o_bJets, pass_bitmap);
    // fillParticle(rTop, Level::Loose, *o_resolvedTop, pass_bitmap);
    // fillAllLeptons(muon, elec, *o_tightLeptons, pass_bitmap);

    // for (size_t syst = 0; syst < numSystematics(); ++syst) {
    //     o_ht.push_back(jet.getHT(Level::Tight, syst));
    //     o_htb.push_back(jet.getHT(Level::Bottom, syst));
    //     o_met.push_back(met.pt());
    //     o_metphi.push_back(met.phi());
    //     o_centrality.push_back(jet.getCentrality(Level::Tight, syst));
    //     o_nb_loose.push_back(jet.n_loose_bjet.at(syst));
    //     o_nb_tight.push_back(jet.n_tight_bjet.at(syst));
    // }
    LOG_FUNC << "End of FillValues";
}

void ZZAnalysis::printStuff()
{
    LOG_FUNC << "Start of printStuff";
    std::cout << "Event: " << *event << std::endl;
    std::cout << "Met: " << met.pt() << std::endl;
    std::cout << "HT: " << jet.getHT(Level::Tight, 0) << std::endl;
    std::cout << "njet: " << jet.size(Level::Tight) << std::endl;
    std::cout << "nbjet: " << jet.size(Level::Bottom) << std::endl;
    std::cout << "nlep: " << muon.size(Level::Tight) << " " << elec.size(Level::Tight) << std::endl;
    std::cout << "nlep loose: " << muon.size(Level::Fake) << " " << elec.size(Level::Fake) << std::endl;
    std::cout << "lepVeto: " << muon.passZVeto() << " " << elec.passZVeto() << std::endl;
    std::cout << std::endl;
    LOG_FUNC << "End of printStuff";
}

