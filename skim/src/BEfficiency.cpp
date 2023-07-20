#include "analysis_suite/skim/interface/BEfficiency.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Signal,
    None,
};

enum class Subchannel {
    All_Dilep_HT,
    All_Dilep,
    All_Multi,
    None,
};

void BEfficiency::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    createTree("Signal", Channel::Signal);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    if (year_ == Year::yr2016pre) {
        setupTrigger(Subchannel::All_Dilep_HT, Dataset::None, {
                "HLT_DoubleMu8_Mass8_PFHT300",
                "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
                "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
            });
        setupTrigger(Subchannel::All_Dilep, Dataset::None, {
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            });
        setupTrigger(Subchannel::All_Multi, Dataset::None, {
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
                "HLT_TripleMu_12_10_5",
            });
    } else if (year_ == Year::yr2016post) {
        setupTrigger(Subchannel::All_Dilep_HT, Dataset::None, {
                "HLT_DoubleMu8_Mass8_PFHT300",
                "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
                "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
            });
        setupTrigger(Subchannel::All_Dilep, Dataset::None, {
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                "HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DoubleMu8_Mass8_PFHT300",
                "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT300",
                "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT300",
            });
        setupTrigger(Subchannel::All_Multi, Dataset::None, {
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL",
                "HLT_TripleMu_12_10_5"
            });
    } else if (year_ == Year::yr2017) {
        setupTrigger(Subchannel::All_Dilep_HT, Dataset::None, {
                "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
                "HLT_DoubleMu8_Mass8_PFHT350",
                "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ",
            });
        setupTrigger(Subchannel::All_Dilep, Dataset::None, {
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
            });
        setupTrigger(Subchannel::All_Multi, Dataset::None, {
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "HLT_TripleMu_12_10_5",
                "HLT_TripleMu_10_5_5_DZ",
            });
    } else if (year_ == Year::yr2018) {
        setupTrigger(Subchannel::All_Dilep_HT, Dataset::None, {
                "HLT_DoubleEle8_CaloIdM_TrackIdM_Mass8_PFHT350",
                "HLT_Mu8_Ele8_CaloIdM_TrackIdM_Mass8_PFHT350_DZ",
                "HLT_DoubleMu4_Mass3p8_DZ_PFHT350",
            });
        setupTrigger(Subchannel::All_Dilep, Dataset::None, {
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
            });
        setupTrigger(Subchannel::All_Multi, Dataset::None, {
                "HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL",
                "HLT_Mu8_DiEle12_CaloIdL_TrackIdL",
                "HLT_DiMu9_Ele9_CaloIdL_TrackIdL_DZ",
                "HLT_TripleMu_12_10_5",
                "HLT_TripleMu_10_5_5_DZ",
            });
    }


    LOG_FUNC << "End of Init";
}

void BEfficiency::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("BJets", "BEffOut", &o_beff);
    tree->Branch("HLT_dilepton_HT", &o_hlt_dilepton_ht);
    tree->Branch("HLT_dilepton", &o_hlt_dilepton);
    tree->Branch("HLT_trilepton", &o_hlt_trilepton);
    tree->Branch("passZVeto_loose", &o_pass_zveto_loose);
    tree->Branch("passZVeto_fake", &o_pass_zveto_fake);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void BEfficiency::clearOutputs()
{
    o_hlt_dilepton_ht.clear();
    o_hlt_dilepton.clear();
    o_hlt_trilepton.clear();
    o_pass_zveto_loose.clear();
    o_pass_zveto_fake.clear();
}

void BEfficiency::ApplyScaleFactors()
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

bool BEfficiency::isSameSign()
{
    int q_total = 0;
    for (size_t idx : muon.list(Level::Tight)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(Level::Tight)) {
        q_total += elec.charge(idx);
    }
    // if 2 leptons, SS -> +1 +1 / -1 -1 -> abs(q) == 2
    // OS cases are 0 and 3, so no overlap
    return abs(q_total) == 1 || abs(q_total) == 2;
}


bool BEfficiency::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";

    if (signal_cuts()) {
        (*currentChannel_) = Channel::Signal;
        return true;
    } else {
        (*currentChannel_) = Channel::None;
        return false;
    }
}


bool BEfficiency::signal_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("passLeadPt", std::max(muPt(0), elPt(0)) > 25);
    passCuts &= cuts.setCut("pass2FakeLep",  nLeps(Level::Tight) >= 2);

    passCuts &= cuts.setCut("passSSLeptons", isSameSign());
    passCuts &= cuts.setCut("passLowMassVeto", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

    return passCuts;
}

void BEfficiency::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    jet.fillJetEff(*o_beff, Level::Loose, event_bitmap);
    for (size_t syst = 0; syst < numSystematics(); ++syst) {
        o_hlt_dilepton_ht.push_back(trig_cuts.pass_cut_any(Subchannel::All_Dilep_HT));
        o_hlt_dilepton.push_back(trig_cuts.pass_cut_any(Subchannel::All_Dilep));
        o_hlt_trilepton.push_back(trig_cuts.pass_cut_any(Subchannel::All_Multi));
        o_pass_zveto_loose.push_back(!muon.isInMassRange(Level::Loose) && !elec.isInMassRange(Level::Loose));
        o_pass_zveto_fake.push_back(!muon.isInMassRange(Level::Fake) && !elec.isInMassRange(Level::Fake));
    }
    LOG_FUNC << "End of FillValues";
}
