#include "analysis_suite/skim/interface/ScaleFactors.h"
#include "analysis_suite/skim/interface/CommonEnums.h"

#include <sstream>
#include <fstream>

void ScaleFactors::init(TTreeReader& fReader)
{
    LOG_FUNC << "Start init";
    if (isMC) {
        // Generator Weights
        LHEScaleWeight.setup(fReader, "LHEScaleWeight", false);
        LHEPdfWeight.setup(fReader, "LHEPdfWeight", false);
        PSWeight.setup(fReader, "PSWeight", false);

        PrefireWeight.setup(fReader, "L1PreFiringWeight_Nom");
        PrefireWeight_up.setup(fReader, "L1PreFiringWeight_Up");
        PrefireWeight_down.setup(fReader, "L1PreFiringWeight_Dn");
        // Pileup Weights
        auto corr_set = getScaleFile("LUM", "puWeights");
        pu_scale = WeightHolder(corr_set->at("Collisions" + yearNum.at(year_)+ "_UltraLegacy_goldenJSON"),
                                Systematic::Pileup, {"nominal", "up", "down"});
        // Trigger Scale Factor
        // auto trig_set = getScaleFile("USER", "trigger_sf");
        // ee_scale = WeightHolder(trig_set->at("EE_TriggerSF"), Systematic::TriggerSF,
        //                         {"nom", "up", "down"});
        // em_scale = WeightHolder(trig_set->at("EM_TriggerSF"), Systematic::TriggerSF,
        //                         {"nom", "up", "down"});
        // mm_scale = WeightHolder(trig_set->at("MM_TriggerSF"), Systematic::TriggerSF,
        //                         {"nom", "up", "down"});
        auto trig_set = getScaleFile("USER", "my_trigger_sf");
        ee_scale = WeightHolder(trig_set->at("EE"), Systematic::TriggerSF,
                                {"nom", "up", "down"});
        em_scale = WeightHolder(trig_set->at("EM"), Systematic::TriggerSF,
                                {"nom", "up", "down"});
        mm_scale = WeightHolder(trig_set->at("MM"), Systematic::TriggerSF,
                                {"nom", "up", "down"});
    } else if (!isMC) {
        // Golden Json
        std::ifstream golden_json_file(scaleDir_ + "/golden_json/golden_json_" + yearMap.at(year_).substr(0,4) + ".json");
        golden_json_file >> golden_json;
    }

    // Fake Rates
    try {
        auto corr_set = getScaleFile("USER", "fake_rates");
        charge_misId = WeightHolder(corr_set->at("Charge_MisId"), Systematic::ChargeMisId_stat,
                                    {"nom", "up", "down"});
        nonprompt_muon = WeightHolder(corr_set->at("Nonprompt_muon"), Systematic::Nonprompt_Mu_stat,
                                      {"nom", "up", "down"});
        nonprompt_elec = WeightHolder(corr_set->at("Nonprompt_electron"), Systematic::Nonprompt_El_stat,
                                      {"nom", "up", "down"});
    } catch (...) {
        LOG_WARN << "Fake Rates not found for this year. May not be necessary, will continue";
    }

    LOG_FUNC << "End init";
}

void ScaleFactors::setup_prescale()
{
    std::ifstream prescale_json_file(scaleDir_ + "/prescale/prescales_" + yearMap.at(year_).substr(0,4) + ".json");
    prescale_json_file >> prescale_json;
    for (size_t i = 0; i < prescale_json["trigs"].size(); ++i) {
        trigger_idx[prescale_json["trigs"][i]] = i;
    }
    for (auto& [run, info]: prescale_json.items()) {
        if (run == "trigs") continue;
        std::istringstream iss(run);
        size_t run_num;
        iss >> run_num;
        auto lumis = info["lumi"].get<std::vector<size_t>>();
        std::vector<std::vector<size_t>> ps_info;
        for (auto& val: lumis) {
            ps_info.push_back(info[std::to_string(val)].get<std::vector<size_t>>());
        }
        prescale_info[run_num] = {lumis, ps_info};
    }
}

size_t ScaleFactors::getPrescale(size_t run, size_t lumi, std::string trig)
{
    auto run_info = prescale_info[run];
    size_t lumi_idx = run_info.distance(lumi) - 1;

    return run_info.prescales[lumi_idx][trigger_idx[trig]];
}

float ScaleFactors::getPrefire()
{
    if (currentSyst == Systematic::Prefire) {
        return (currentVar == eVar::Up) ? *PrefireWeight_up : *PrefireWeight_down;
    } else {
        return *PrefireWeight;
    }
}

float ScaleFactors::getPileupSF(int nPU)
{
    return pu_scale.evaluate({static_cast<float>(nPU), systName(pu_scale)});
}

float ScaleFactors::getTriggerSF(Particle& elec, Particle& muon)
{
    if (elec.size(Level::Fake) + muon.size(Level::Fake) < 2) return 1.;

    if (muon.size(Level::Fake) == 0)  {
        float pt1 = (elec.pt(Level::Fake, 0) < 20) ? 20 : elec.pt(Level::Fake, 0);
        float pt2 = (elec.pt(Level::Fake, 1) < 20) ? 20 : elec.pt(Level::Fake, 1);
        return ee_scale.evaluate({pt1, pt2});
    } else if (elec.size(Level::Fake) == 0) {
        float pt1 = (muon.pt(Level::Fake, 0) < 20) ? 20 : muon.pt(Level::Fake, 0);
        float pt2 = (muon.pt(Level::Fake, 1) < 20) ? 20 : muon.pt(Level::Fake, 1);
        return mm_scale.evaluate({pt1, pt2});
    } else {
        float pte = (elec.pt(Level::Fake, 0) < 20) ? 20 : elec.pt(Level::Fake, 0);
        float ptm = (muon.pt(Level::Fake, 0) < 20) ? 20 : muon.pt(Level::Fake, 0);
        return em_scale.evaluate({pte, ptm});
    }
}

float ScaleFactors::getLHESF()
{
    // [0] is muR=0.5 muF=0.5 ; [1] is muR=0.5 muF=1.0 ; [2] is muR=0.5 muF=2.0 ;
    // [3] is muR=0.1 muF=0.5 ; [4] is muR=1.0 muF=1.0 ; [5] is muR=1.0 muF=2.0 ;
    // [6] is muR=2.0 muF=0.5 ; [7] is muR=2.0 muF=1.0 ; [8] is muR=2.0 muF=2.0 ;

    if (LHEScaleWeight.size() != 9) {
        return 1.;
    }
    if (isMC) {
        int varIdx = 4;
        if (currentSyst == Systematic::LHE_muF) {
            varIdx = (currentVar == eVar::Up) ? 5 : 3;
        } else if (currentSyst == Systematic::LHE_muR) {
            varIdx = (currentVar == eVar::Up) ? 7 : 1;
        }
        return LHEScaleWeight.at(varIdx);
    }
    return 1.;
}

float ScaleFactors::getPartonShower()
{
    // [0] is ISR=2 FSR=1; [1] is ISR=1 FSR=2
    // [2] is ISR=0.5 FSR=1; [3] is ISR=1 FSR=0.5;
    if (!isMC || PSWeight.size() != 4) return 1.;

    if (currentSyst == Systematic::PS_ISR) {
        return (currentVar == eVar::Up) ? PSWeight.at(0) : PSWeight.at(2);
    } else if (currentSyst == Systematic::PS_FSR) {
        return (currentVar == eVar::Up) ? PSWeight.at(1) : PSWeight.at(3);
    }
    return 1.;
}

float ScaleFactors::getLHEPdf()
{
    // [0] central value
    // [1-100] PDF replica variations
    // [101] alphaZ down ; [102] alphaZ up
    // http://nnpdf.mi.infn.it/wp-content/uploads/2019/03/NNPDFfits_Positivity_PhysicalCrossSections_v2.pdf
    if (currentSyst == Systematic::PDF_unc) {
        if (LHEPdfWeight.size() < 101) {
            return 1.;
        }
        std::sort(LHEPdfWeight.begin()+1, LHEPdfWeight.begin()+100);
        float err = (LHEPdfWeight.at(84) - LHEPdfWeight.at(17))/2;
        return LHEPdfWeight.at(0) + ((currentVar == eVar::Up) ? err : -err);
    } else if (currentSyst == Systematic::PDF_alphaZ) {
        if (LHEPdfWeight.size() != 103) {
            return 1.;
        } else  {
            return (currentVar == eVar::Up) ? LHEPdfWeight.at(102) : LHEPdfWeight.at(101);
        }
    } else {
        if (LHEPdfWeight.size() > 1) {
            return LHEPdfWeight.at(0);
        } else {
            return 1.;
        }
    }
}


bool ScaleFactors::inGoldenLumi(UInt_t run, UInt_t lumi)
{
    if (golden_json.contains(std::to_string(run))) {
        for (auto lumi_pair : golden_json[std::to_string(run)]) {
            if (lumi < lumi_pair[0]) {
                return false;
            } else if (lumi <= lumi_pair[1]) {
                return true;
            }
        }
    }
    return false;
}

float ScaleFactors::getChargeMisIdFR(float eta, float pt)
{
    std::string syst = systName(charge_misId);
    if (pt < 20) pt = 20;
    return charge_misId.evaluate({syst, fabs(eta), pt});
}

float ScaleFactors::getNonpromptFR(float eta, float pt, PID pid)
{
    if (pt < 20) pt = 20;
    if (pid == PID::Muon) {
        std::string syst = systName(nonprompt_muon);
        return nonprompt_muon.evaluate({syst, fabs(eta), pt});
    } else if (pid == PID::Electron) {
        std::string syst = systName(nonprompt_elec);
        return nonprompt_elec.evaluate({syst, fabs(eta), pt});
    }
    return 1.0;
}
