#include "analysis_suite/skim/interface/TwoLepton.h"

#include "analysis_suite/skim/interface/logging.h"

enum class Channel {
    Measurement,
    None,
};

void TwoLepton::Init(TTree* tree)
{
    LOG_FUNC << "Start of Init";
    met_type = MET_Type::PUPPI;
    BaseSelector::Init(tree);

    createTree("Measurement", Channel::Measurement);

    if (isMC_) {
        Pileup_nTrueInt.setup(fReader, "Pileup_nTrueInt");
    }

    setup_trigger();

    LOG_FUNC << "End of Init";
}

void TwoLepton::ApplyScaleFactors()
{
    LOG_FUNC << "Start of ApplyScaleFactors";
    LOG_EVENT << "weight: " << (*weight);
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPileupSF(*Pileup_nTrueInt);
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getLHESF();
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPrefire();
    // std::cout << (*weight) << " ";
    (*weight) *= sfMaker.getPartonShower();
    // std::cout << (*weight) << " ";
    (*weight) *= jet.getScaleFactor();
    // std::cout << (*weight) << " ";
    (*weight) *= elec.getScaleFactor();
    // std::cout << (*weight) << " ";
    (*weight) *= muon.getScaleFactor();
    // std::cout << (*weight) << " " << "end" << std::endl;
    LOG_FUNC << "End of ApplyScaleFactors";
}

bool TwoLepton::getCutFlow()
{
    LOG_FUNC << "Start of passSelection";
    (*currentChannel_) = Channel::None;
    setSubChannel();

    if (measurement_cuts()) {
        (*currentChannel_) = Channel::Measurement;
    }

    if (*currentChannel_ == Channel::None) {
        return false;
    }

    LOG_FUNC << "End of passSelection";
    return true;
}

bool TwoLepton::measurement_cuts()
{
    bool passCuts = true;
    CutInfo cuts;

    passCuts &= cuts.setCut("passPreselection", true);
    passCuts &= cuts.setCut("passMETFilter", metfilters.pass());
    passCuts &= cuts.setCut("pass2FakeLep",  nLeps(Level::Fake) >= 2);

    // Trigger Cuts
    passCuts &= cuts.setCut("passLeadPtCut", getLeadPt() > 25);
    passCuts &= cuts.setCut("passTrigger", getTriggerValue());
    return passCuts;
}

