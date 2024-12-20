#include "analysis_suite/skim/interface/BaseSelector.h"

#include"analysis_suite/skim/interface/logging.h"

#include "analysis_suite/skim/interface/Systematic.h"
#include "analysis_suite/skim/interface/ScaleFactors.h"
#include "analysis_suite/skim/interface/CommonFuncs.h"

void BaseSelector::SetupOutTreeBranches(TTree* tree)
{
    tree->Branch("weight", &o_weight);
    tree->Branch("PassEvent", &o_pass_event);
    run.setupBranch(tree);
    lumiblock.setupBranch(tree);
    event.setupBranch(tree);
 
}

void BaseSelector::Init(TTree* tree)
{
    // Set verbosity
    for (auto item : *fInput) {
        std::string itemName = item->GetName();
        if (itemName == "Verbosity") {
            loguru::g_stderr_verbosity = std::stoi(item->GetTitle());
        }
    }

    LOG_FUNC << "Start of Init";
    loguru::g_preamble_thread = false;
    loguru::g_preamble_time = false;
    if (!tree)
        return;

    // Read in python inputs and set them in the Selector
    TList* rootSystList = new TList();
    rootSystList->SetName("Systematics");
    rootSystList->Add(new TNamed("Nominal", "Nominal"));

    TList* syst_index_list = new TList();
    syst_index_list->SetName("Syst_Index");
    syst_index_list->Add(new TNamed("0", "0"));

    std::vector<Systematic> used_jec_systs;

    LOG_POST <<  "Start Reading python inputs";
    for (auto item : *fInput) {
        std::string itemName = item->GetName();
        if (itemName == "MetaData") {
            fOutput->Add(item);
            for (auto data : *static_cast<TList*>(item)) {
                std::string dataName = data->GetName();
                if (dataName == "Year") {
                    std::string year = data->GetTitle();
                    year_ = get_by_val(yearMap, year);
                } else if (dataName == "isData") {
                    isMC_ = static_cast<std::string>(data->GetTitle()) == "False";
                } else if (dataName == "Group") {
                    groupName_ = data->GetTitle();
                } else if (dataName == "Dataset") {
                    dataset_ = dataset_name_to_enum.at(data->GetTitle());
                } else if (dataName == "DataRun") {
                    data_run = data->GetTitle();
                }
            }
        } else if (itemName == "Systematics") {
            for (auto systNamed : *static_cast<TList*>(item)) {
                std::string systName = systNamed->GetName();
                bool data_syst = std::find(data_systs.begin(), data_systs.end(), systName) != data_systs.end();
                if (isMC_ != !data_syst)
                    continue;
                // Add systematic to list used by selector as well as TList for writing out
                auto cur_syst = syst_by_name.at(systName);
                if (std::find(jec_systs.begin(), jec_systs.end(), cur_syst) != jec_systs.end()) {
                    used_jec_systs.push_back(cur_syst);
                }
                systematics_.push_back(cur_syst);
                size_t syst_start = 2*systematics_.size() - 3;
                if (std::find(systs_that_change.begin(), systs_that_change.end(), systematics_.back()) != systs_that_change.end()) {
                    novel_systs.insert(novel_systs.end(), {syst_to_index.size(), syst_to_index.size()+1});
                    size_t start = novel_systs.size()-2;
                    syst_to_index.insert(syst_to_index.end(), {start, start+1});
                    syst_index_list->Add(new TNamed(std::to_string(syst_start), std::to_string(start)));
                    syst_index_list->Add(new TNamed(std::to_string(syst_start+1), std::to_string(start+1)));
                } else {
                    syst_index_list->Add(new TNamed(std::to_string(syst_start), "0"));
                    syst_index_list->Add(new TNamed(std::to_string(syst_start+1), "0"));
                    syst_to_index.insert(syst_to_index.end(), {0, 0});
                }
                rootSystList->Add(new TNamed((systName + "_up").c_str(), systName.c_str()));
                rootSystList->Add(new TNamed((systName + "_down").c_str(), systName.c_str()));
            }
        } else if (itemName == "NEvents") {
            stop_point = std::stoi(item->GetTitle());
        }
    }

    for (auto syst : systematics_) {
        std::vector<eVar> vars = (syst != Systematic::Nominal) ? syst_vars : nominal_var;
        for (auto var : vars) {
            syst_var_pair.push_back(std::make_pair(syst, var));
        }
    }

    outdir = outfile->mkdir(groupName_.c_str());
    fOutput->Add(rootSystList);
    fOutput->Add(syst_index_list);
    LOG_POST << "Finished setting python inputs";
    setupSystematicInfo();

    fReader.SetTree(tree);
    if (isMC_) {
        genWeight.setup(fReader, "genWeight");
    }

    // Setup basic variables used and needed by the base selector
    run.setup(fReader, "run");
    lumiblock.setup(fReader, "luminosityBlock");
    event.setup(fReader, "event");
    PV_npvs.setup(fReader, "PV_npvs");

    metfilters.setup(fReader, {"Flag_goodVertices", "Flag_globalSuperTightHalo2016Filter", "Flag_HBHENoiseFilter", "Flag_HBHENoiseIsoFilter",
                               "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter","Flag_ecalBadCalibFilter",
                               "Flag_eeBadScFilter", "Flag_BadPFMuonDzFilter",
        });

    // Set output vector. Fixed to size of systematics (one per variation)
    o_weight.resize(numSystematics());
    o_channels.resize(novel_systs.size());
    o_pass_event.resize(novel_systs.size());

    // Setup particle branches
    sfMaker.init(fReader);
    // met.setup(fReader, used_jec_systs, met_type);
    met.setup(fReader, met_type);
    muon.setup(fReader);
    elec.setup(fReader);
    jet.setup(fReader, data_run, used_jec_systs);
    if (isMC_) {
        rGen.setup(fReader);
        rGenJet.setup(fReader);
    }
    trig_cuts.set_dataset(dataset_);

    if (loguru::g_stderr_verbosity > 0) {
        bar.set_total(tree->GetEntries());
        bar.print_header();
    }
    LOG_FUNC << "End of Init";
}

Bool_t BaseSelector::Process(Long64_t entry)
{
    if (stop_point > 0 && entry > stop_point) return false;
    // if (loguru::g_stderr_verbosity > 8 && entry > 2) return false;
    // else if (loguru::g_stderr_verbosity > 2 && entry > 10000) return false;

    LOG_FUNC << "Start of Process";
    if (loguru::g_stderr_verbosity > 0) {
        bar++;
        bar.print_bar();
    }
    fReader.SetLocalEntry(entry);
    clearParticles();

    // Remove non-golden lumi stuff
    if (!isMC_ && !sfMaker.inGoldenLumi(*run, *lumiblock)) {
        return false;
    } else if (!passPreselection()) {
        return false;
    }

    std::vector<bool> systPassSelection;
    bool passAny = false;
    jet.setup_jets(rGenJet, *run);
    for (auto it = syst_var_pair.begin(); it != syst_var_pair.end(); ++it) {
        Systematic syst = it->first;
        eVar var = it->second;
        size_t systNum = std::distance(syst_var_pair.begin(), it);
        LOG_EVENT_IF(syst != Systematic::Nominal) << "Systematic is " << get_by_val(syst_by_name, syst)
                                                  << "| Variation is: " << varName_by_var.at(var);
        SetupEvent(systNum);
        // Add result of setting channel to vector used for writing out results
        if (syst_to_index.at(systNum) > 0 || systNum == 0) {
            systPassSelection.push_back(getCutFlow());
            passAny |= systPassSelection.back();
        }
        if (systPassSelection.back()) {
            ApplyChannelScaleFactors();
        }
    }

    if (systPassSelection.size() > 0 && systPassSelection[0]) {
        bar.pass(); // Currenly, only events in first channel are put in the bar
    }
    // Don't write out if nothing passed
    if (!passAny) return false;

    // Fill up variables then fill up the tree
    setupSyst(0);
    for (auto& [chan, tree] : trees) {
        bool passedChannel = false;
        Bitmap event_bitmap;
        for (size_t syst = 0; syst < novel_systs.size(); ++syst) {
            o_pass_event[syst] = systPassSelection[syst] && chan == o_channels[syst];
            event_bitmap[syst] = o_pass_event[syst];
            passedChannel |= o_pass_event[syst];
        }

        if (passedChannel) {
            clearOutputs();

            run.fill();
            lumiblock.fill();
            event.fill();
            FillValues(event_bitmap);
            tree.tree->Fill();
        }
    }
    LOG_FUNC << "End of Process";
    return kTRUE;
}

void BaseSelector::SlaveTerminate()
{
    if (loguru::g_stderr_verbosity > 0) bar.print_trailer();
    LOG_POST << "End job";
}

void BaseSelector::setupSyst(size_t systNum)
{
    size_t novelNum = syst_to_index.at(systNum);
    SystematicWeights::currentSyst = syst_var_pair.at(systNum).first;
    SystematicWeights::currentVar = syst_var_pair.at(systNum).second;

    muon.setSyst(novelNum);
    elec.setSyst(novelNum);
    jet.setSyst(novelNum);
    met.setSyst();
}

void BaseSelector::SetupEvent(size_t systNum)
{
    LOG_FUNC << "Start of SetupEvent";
    size_t novelNum = syst_to_index.at(systNum);
    setupSyst(systNum);
    weight = &o_weight[systNum];
    currentChannel_ = &o_channels[novelNum];

    (*weight) = 1.0;

    if (novelNum > 0 || systNum == 0) {
        met.setupMet(jet, *run, *PV_npvs);
        // Setup good particle lists
        muon.setGoodParticles(novelNum, jet, rGen);
        elec.setGoodParticles(novelNum, jet, rGen);
        jet.setGoodParticles(novelNum);
        // Class dependent setting up
        setOtherGoodParticles(novelNum);
    }

    (*weight) *= muon.getRocCorrection(rGen);
    if (isMC_) {
        (*weight) *= *genWeight;
        ApplyScaleFactors();
    }

    LOG_FUNC << "End of SetupEvent";
}

void BaseSelector::clearParticles()
{
    LOG_FUNC << "Start of clearParticles";
    muon.clear();
    elec.clear();
    jet.clear();
    rGen.clear();
    rGenJet.clear();

    std::fill(o_weight.begin(), o_weight.end(), 1.);
    LOG_FUNC << "End of clearParticles";
}

void BaseSelector::setupSystematicInfo()
{
    LOG_FUNC << "Start of setupSystematicInfo";
    std::string scaleDir = getenv("CMSSW_BASE");
    scaleDir += "/src/analysis_suite/data";

    SystematicWeights::nSyst = novel_systs.size();
    SystematicWeights::year_ = year_;
    SystematicWeights::scaleDir_ = scaleDir;
    SystematicWeights::isMC = isMC_;

    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_Absolute] = "jesAbsolute_20"+yearNum.at(year_);
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_Absolute_corr] = "jesAbsolute";
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_BBEC1] = "jesBBEC1_20"+yearNum.at(year_);
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_BBEC1_corr] = "jesBBEC1";
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_EC2] = "jesEC2_20"+yearNum.at(year_);
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_EC2_corr] = "jesEC2";
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_HF] = "jesHF_20"+yearNum.at(year_);
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_HF_corr] = "jesHF";
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_RelativeBal] = "jesRelativeBal";
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_RelativeSample] = "jesRelativeSample_20"+yearNum.at(year_);
    SystematicWeights::jecSyst_to_string[Systematic::Jet_JEC_FlavorQCD] = "jesFlavorQCD";
    LOG_FUNC << "End of setupSystematicInfo";
}
