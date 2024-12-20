#include "analysis_suite/skim/interface/DileptonBase.h"

void DileptonBase::Init(TTree* tree)
{
    BaseSelector::Init(tree);
    std::string scaleDir = getenv("CMSSW_BASE");
    scaleDir += "/src/analysis_suite/data";
    // std::ifstream trigger_file(scaleDir + "/golden_json/golden_json_" + yearMap.at(year_).substr(0,4) + ".json");
    // golden_json_file >> dilep_trigs;
    setup_trigger();

}

void DileptonBase::setup_trigger()
{
    if (year_ == Year::yr2016pre || year_ == Year::yr2016post) {
        setupTrigger(Subchannel::Single_M, Dataset::Single_M,
                     {"HLT_IsoMu24", "HLT_IsoTkMu24"});
        setupTrigger(Subchannel::Single_E, Dataset::Single_E,
                     {"HLT_Ele27_WPTight_Gsf",});

        if (data_run != "H") {
            setupTrigger(Subchannel::EM, Dataset::MuonEG, {
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL"});
            setupTrigger(Subchannel::MM, Dataset::DoubleMuon, {
                    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL",
                    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"});
        } else {
            setupTrigger(Subchannel::EM, Dataset::MuonEG, {
                    "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                    "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"});
            setupTrigger(Subchannel::MM, Dataset::DoubleMuon, {
                    "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ",
                    "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"});
        }

        setupTrigger(Subchannel::EE, Dataset::DoubleEG, {
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DoubleEle33_CaloIdL_MW",
                "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL",
                "HLT_Ele27_WPTight_Gsf",
            });
    } else if (year_ == Year::yr2017) {
        setupTrigger(Subchannel::Single_M, Dataset::Single_M,
                     {"HLT_IsoMu27",});
        setupTrigger(Subchannel::Single_E, Dataset::Single_E,
                     {"HLT_Ele35_WPTight_Gsf",});

        std::vector<std::string> mm_list;
        if (isMC_ || data_run == "B") {
            mm_list.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
        }
        if (isMC_ || data_run != "B") {
            mm_list.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8");
            mm_list.push_back("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
        }
        setupTrigger(Subchannel::MM, Dataset::DoubleMuon, mm_list);
        setupTrigger(Subchannel::EM, Dataset::MuonEG, {
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG, {
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL",
                "HLT_DoubleEle33_CaloIdL_MW",});
    } else if (year_ == Year::yr2018) {
        setupTrigger(Subchannel::Single_M, Dataset::Single_M, {
                "HLT_IsoMu24"});
        setupTrigger(Subchannel::Single_E, Dataset::DoubleEG, {
                "HLT_Ele32_WPTight_Gsf",});

        setupTrigger(Subchannel::MM, Dataset::DoubleMuon, {
                "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8",});
        setupTrigger(Subchannel::EM, Dataset::MuonEG, {
                "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",});
        setupTrigger(Subchannel::EE, Dataset::DoubleEG, {
                "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ",
                "HLT_DoubleEle25_CaloIdL_MW",
                "HLT_Ele32_WPTight_Gsf",});

    }
    setupTrigger(Subchannel::None, Dataset::None, {});
}
bool DileptonBase::getTriggerValue()
{
    if (subChannel_ == Subchannel::EE) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::EE)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (year_ == Year::yr2018) { // 2018 uses just EGamma not double vs single EG datasets
            return (trig_cuts.pass_cut(Subchannel::EE)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (trig_cuts.dataset_or_trig(Subchannel::EE)) {
            return trig_cuts.pass_cut(Subchannel::EE);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_E);
        }
    } else if (subChannel_ == Subchannel::MM) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::MM)
                    || trig_cuts.pass_cut(Subchannel::Single_M));
        } else if (trig_cuts.dataset_or_trig(Subchannel::MM)) {
            return trig_cuts.pass_cut(Subchannel::MM);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_M);
        }
    } else if (subChannel_ == Subchannel::EM) {
        if (isMC_) {
            return (trig_cuts.pass_cut(Subchannel::EM)
                    || trig_cuts.pass_cut(Subchannel::Single_M)
                    || trig_cuts.pass_cut(Subchannel::Single_E));
        } else if (trig_cuts.dataset_or_trig(Subchannel::EM)) {
            return trig_cuts.pass_cut(Subchannel::EM);
        } else if (trig_cuts.dataset_or_trig(Subchannel::Single_M)) {
            return trig_cuts.pass_cut(Subchannel::Single_M);
        } else {
            return trig_cuts.pass_cut(Subchannel::Single_E);
        }
    } else {
        return false;
    }
}
float DileptonBase::getLeadPt()
{
    if (subChannel_ == Subchannel::MM) {
        return muon.rawpt(Level::Fake, 0);
    } else if (subChannel_ == Subchannel::EE) {
        return elec.rawpt(Level::Fake, 0);
    } else if(subChannel_ == Subchannel::EM){
        return std::max(muon.rawpt(Level::Fake, 0), elec.rawpt(Level::Fake, 0));
    }
    return 0.;
}


void DileptonBase::setSubChannel()
{
    subChannel_ = Subchannel::None;

    if(nLeps(Level::Fake) >= 2) {
        if (muon.size(Level::Fake) * elec.size(Level::Fake) > 0) {
            subChannel_ = Subchannel::EM;
        } else if (muon.size(Level::Fake) == 0) {
            subChannel_ = Subchannel::EE;
        }  else {
            subChannel_ = Subchannel::MM;
        }
    }
}

bool DileptonBase::isSameSign(Level level)
{
    int q_total = 0;
    for (size_t idx : muon.list(level)) {
        q_total += muon.charge(idx);
    }
    for (size_t idx : elec.list(level)) {
        q_total += elec.charge(idx);
    }
    return abs(q_total) == 1 || abs(q_total) == 2;
}
