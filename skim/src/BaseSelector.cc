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
    LOG_FUNC << "Start of Init";
    loguru::g_preamble_thread = false;
    loguru::g_preamble_time = false;
    if (!tree)
        return;

    // Read in python inputs and set them in the Selector
    TList* rootSystList = new TList();
    rootSystList->SetName("Systematics");
    rootSystList->Add(new TNamed("Nominal", "Nominal"));
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
                }
            }
        } else if (itemName == "Systematics") {
            for (auto systNamed : *static_cast<TList*>(item)) {
                std::string systName = systNamed->GetName();
                bool data_syst = std::find(data_systs.begin(), data_systs.end(), systName) != data_systs.end();
                if (isMC_ != !data_syst)
                    continue;
                // Add systematic to list used by selector as well as TList for writing out
                systematics_.push_back(syst_by_name.at(systName));
                rootSystList->Add(new TNamed((systName + "_up").c_str(), systName.c_str()));
                rootSystList->Add(new TNamed((systName + "_down").c_str(), systName.c_str()));
            }
        } else if (itemName == "Verbosity") {
            loguru::g_stderr_verbosity = std::stoi(item->GetTitle());
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
                               "Flag_EcalDeadCellTriggerPrimitiveFilter", "Flag_BadPFMuonFilter","Flag_ecalBadCalibFilter"});

    // Set output vector. Fixed to size of systematics (one per variation)
    o_weight.resize(numSystematics());
    o_channels.resize(numSystematics());
    o_pass_event.resize(numSystematics());

    // Setup particle branches
    sfMaker.init(fReader);
    met.setup(fReader, met_type);
    muon.setup(fReader);
    elec.setup(fReader);
    jet.setup(fReader);
    if (isMC_) {
        rGen.setup(fReader);
        rGenJet.setup(fReader);
    }

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
    }

    std::vector<bool> systPassSelection;
    bool passAny = false;
    jet.setupJEC(rGenJet);
    for (auto it = syst_var_pair.begin(); it != syst_var_pair.end(); ++it) {
        Systematic syst = it->first;
        eVar var = it->second;
        LOG_EVENT_IF(syst != Systematic::Nominal) << "Systematic is " << get_by_val(syst_by_name, syst)
                                                  << "| Variation is: " << varName_by_var.at(var);
        SetupEvent(std::distance(syst_var_pair.begin(), it));
        // Add result of setting channel to vector used for writing out results
        systPassSelection.push_back(getCutFlow());
        passAny |= systPassSelection.back();
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
        for (size_t syst = 0; syst < numSystematics(); ++syst) {
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
    SystematicWeights::currentSyst = syst_var_pair.at(systNum).first;
    SystematicWeights::currentVar = syst_var_pair.at(systNum).second;

    muon.setSyst(systNum);
    elec.setSyst(systNum);
    jet.setSyst(systNum);
    met.setSyst();
}

void BaseSelector::SetupEvent(size_t systNum)
{
    LOG_FUNC << "Start of SetupEvent";
    setupSyst(systNum);
    weight = &o_weight[systNum];
    currentChannel_ = &o_channels[systNum];

    (*weight) = isMC_ ? *genWeight : 1.0;

    met.setupMet(jet, *run, *PV_npvs);
    // Setup good particle lists
    muon.setGoodParticles(systNum, jet, rGen);
    elec.setGoodParticles(systNum, jet, rGen);
    jet.setGoodParticles(systNum);
    // Class dependent setting up
    setOtherGoodParticles(systNum);

    if (isMC_) {
        ApplyScaleFactors();
    }
    // (*weight) *= muon.getRocCorrection(rGen, isMC_);

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

    SystematicWeights::nSyst = numSystematics();
    SystematicWeights::year_ = year_;
    SystematicWeights::scaleDir_ = scaleDir;
    SystematicWeights::isMC = isMC_;
    LOG_FUNC << "End of setupSystematicInfo";
}
