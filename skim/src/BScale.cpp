#include "analysis_suite/skim/interface/BScale.h"

#include"analysis_suite/skim/interface/logging.h"
#include"analysis_suite/skim/interface/CommonFuncs.h"

namespace Channel {
    enum {
        Preselection,
        None,
    };
}

void BScale::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    DileptonBase::Init(tree);

    createTree("Preselection", Channel::Preselection);
    Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");

    LOG_FUNC << "End of Init";
}

void BScale::SetupOutTreeBranches(TTree* tree)
{
    LOG_FUNC << "Start of SetupOutTreeBranches";
    BaseSelector::SetupOutTreeBranches(tree);
    tree->Branch("Jets", "JetOut", &o_jets);
    tree->Branch("wgt_nobtag", &o_wgt_nobtag);

    LOG_FUNC << "End of SetupOutTreeBranches";
}

void BScale::clearOutputs()
{
    LOG_FUNC << "Start of clearOutputs";
    o_wgt_nobtag.clear();

    LOG_FUNC << "End of clearOutputs";
}

void BScale::ApplyScaleFactors()
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

bool BScale::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (!baseline_cuts()) return false;

    if (signal_cuts())
        (*currentChannel_) = Channel::Preselection;

    if (trees.find(*currentChannel_) == trees.end()) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool BScale::baseline_cuts()
{
    LOG_FUNC << "Start of baseline_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());

    passCuts &= cuts.setCut("passLeadPt", getLeadPt() > 25);
    passCuts &= cuts.setCut("passFakeLeptonNum", nLeps(Level::Fake) >= 2);

    passCuts &= cuts.setCut("passTrigger", getTriggerValue());
    passCuts &= cuts.setCut("passJetNumber", jet.size(Level::Tight) >= 1);
    passCuts &= cuts.setCut("passMetCut", met.pt() > 25);
    passCuts &= cuts.setCut("passHTCut", jet.getHT(Level::Tight) > 250);
    passCuts &= cuts.setCut("passLowMassCut", !muon.isInMassRange(Level::Loose, 0., 12.) && !elec.isInMassRange(Level::Loose, 0., 12.));

    LOG_FUNC << "End of baseline_cuts";
    return passCuts;
}

bool BScale::signal_cuts()
{
    LOG_FUNC << "Start of signal_cuts";
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("pass2Or3Leptons", nLeps(Level::Tight) == 2 || nLeps(Level::Tight) == 3);
    passCuts &= cuts.setCut("passVetoFakeLeptons", nLeps(Level::Tight) == nLeps(Level::Fake));
    passCuts &= cuts.setCut("passSameSign;", isSameSign(Level::Tight));

    // Fill Cut flow
    fillCutFlow(Channel::Preselection, cuts);

    LOG_FUNC << "End of signal_cuts";
    return passCuts;
}

void BScale::FillValues(const Bitmap& event_bitmap)
{
    LOG_FUNC << "Start of FillValues";
    jet.fillJet(*o_jets, Level::Tight, event_bitmap);
    for (size_t systNum = 0; systNum < numSystematics(); ++systNum) {
        size_t syst = syst_to_index.at(systNum);
        if (syst == 0 && systNum != 0) {
            continue;
        }
        setupSyst(systNum);
        o_wgt_nobtag.push_back(o_weight[systNum]/jet.getTotalBTagWeight("M"));
    }
    LOG_FUNC << "End of FillValues";
}
