#include "analysis_suite/skim/interface/BEfficiency.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Signal,
    None,
};


void BEfficiency::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    BaseSelector::Init(tree);

    createTree("Signal", Channel::Signal);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    LOG_FUNC << "End of Init";
}

void BEfficiency::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("BJets", "BEffOut", &o_beff);
    LOG_FUNC << "End of SetupOutTreeBranches";
}

void BEfficiency::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";

    LOG_EVENT << "weight: " << (*weight);
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    LOG_EVENT << "weight after pu scale: " << (*weight);
    (*weight) *= sfMaker.getLHESF();
    LOG_EVENT << "weight after lhe scale: " << (*weight);
    (*weight) *= sfMaker.getLHEPdf();
    LOG_EVENT << "weight after lhe pdf: " << (*weight);
    (*weight) *= jet.getScaleFactor();
    LOG_EVENT << "weight after jet scale: " << (*weight);
    (*weight) *= elec.getScaleFactor();
    LOG_EVENT << "weight after elec scale: " << (*weight);
    (*weight) *= muon.getScaleFactor();
    LOG_EVENT << "weight after muon scale: " << (*weight);

    LOG_FUNC << "End of ApplyScaleFactors";
}

bool BEfficiency::isSameSign()
{
    int q_total = 0;
    for (size_t idx : muon.list(Level::Fake)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(Level::Fake)) {
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
    passCuts &= cuts.setCut("pass2FakeLep",  nLeps(Level::Fake) == 2);

    passCuts &= cuts.setCut("passSSLeptons", isSameSign());
    passCuts &= cuts.setCut("passZVeto", muon.passZVeto() && elec.passZVeto());
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 2);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 100);

    return passCuts;
}

void BEfficiency::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    jet.fillJetEff(*o_beff, Level::Loose, event_bitmap);
    LOG_FUNC << "End of FillValues";
}
