#ifndef TRIGGER_TEMPLATE_H_
#define TRIGGER_TEMPLATE_H_


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
            "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW",
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
            "HLT_DoubleEle25_CaloIdL_MW",});

 }
setupTrigger(Subchannel::None, Dataset::None, {});

#endif // TRIGGER_TEMPLATE_H_
