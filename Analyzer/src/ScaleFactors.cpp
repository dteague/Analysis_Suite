#include "analysis_suite/Analyzer/interface/ScaleFactors.h"
#include "analysis_suite/Analyzer/interface/CommonEnums.h"

#include <sstream>


void ScaleFactors::init(bool isMC_, TTreeReader& fReader)
{
    isMC = isMC_;
    if (isMC) {
        // Generator Weights
        LHEScaleWeight.setup(fReader, "scaleWeights");
        // LHEPdfWeight.setup(fReader, "pdfWeights");
        // PSWeight.setup(fReader, "PSWeight");

        // Pileup Weights
        auto corr_set = getScaleFile("LUM", "puWeights");
        pu_scale = WeightHolder(corr_set->at("Collisions" + yearNum.at(year_)+ "_UltraLegacy_goldenJSON"),
                                Systematic::Pileup, {"nominal", "up", "down"});
    } else if (!isMC) {
        // Golden Json
        std::ifstream golden_json_file(scaleDir_ + "/golden_json/golden_json_" + yearMap.at(year_).substr(0,4) + ".json");
        golden_json_file >> golden_json;
    }
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

float ScaleFactors::getPileupSF(int nPU)
{
    return pu_scale.evaluate({static_cast<float>(nPU), systName(pu_scale)});
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
    // if (!isMC || PSWeight.size() != 4) return 1.;

    // if (currentSyst == Systematic::PS_ISR) {
    //     return (currentVar == eVar::Up) ? PSWeight.at(0) : PSWeight.at(2);
    // } else if (currentSyst == Systematic::PS_FSR) {
    //     return (currentVar == eVar::Up) ? PSWeight.at(1) : PSWeight.at(3);
    // }
    return 1.;
}

float ScaleFactors::getLHEPdf()
{
    // [0] central value
    // [1-100] PDF replica variations
    // [101] alphaZ down ; [102] alphaZ up
    // http://nnpdf.mi.infn.it/wp-content/uploads/2019/03/NNPDFfits_Positivity_PhysicalCrossSections_v2.pdf
    // if (currentSyst == Systematic::PDF_unc) {
    //     if (LHEPdfWeight.size() < 101) {
    //         return 1.;
    //     }
    //     std::sort(LHEPdfWeight.begin()+1, LHEPdfWeight.begin()+100);
    //     float err = (LHEPdfWeight.at(17) - LHEPdfWeight.at(84))/2;
    //     return LHEPdfWeight.at(0) + ((currentVar == eVar::Up) ? err : -err);
    // } else if (currentSyst == Systematic::PDF_alphaZ) {
    //     if (LHEPdfWeight.size() != 103) {
    //         return 1.;
    //     } else  {
    //         return (currentVar == eVar::Up) ? LHEPdfWeight.at(102) : LHEPdfWeight.at(101);
    //     }
    // } else {
    //     if (LHEPdfWeight.size() > 1) {
    //         return LHEPdfWeight.at(0);
    //     } else {
    //         return 1.;
    //     }
    // }
    return 1.;
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
