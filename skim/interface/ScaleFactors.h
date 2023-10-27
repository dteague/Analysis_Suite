#ifndef __SCALEFACTORS_H_
#define __SCALEFACTORS_H_

#include <TTreeReader.h>

#include <nlohmann/json.hpp>
#include <algorithm>

#include "analysis_suite/skim/interface/Systematic.h"
#include "analysis_suite/skim/interface/Variable.h"
#include "analysis_suite/skim/interface/Particle.h"

class ScaleFactors : public SystematicWeights {
public:
    ScaleFactors() {}

    void init(TTreeReader& fReader);

    void setup_prescale();

    float getPileupSF(int nPU);

    float getLHESF();
    float getLHEPdf();
    float getPartonShower();
    float getPrefire();
    float getTriggerSF(Particle& elec, Particle& muon);

    size_t getPrescale(size_t run, size_t lumi, std::string trig);

    bool inGoldenLumi(UInt_t run, UInt_t lumi);

    float getChargeMisIdFR(float eta, float pt);
    float getNonpromptFR(float eta, float pt, PID pid);

private:

    TRArray<Float_t> LHEScaleWeight;
    TRArray<Float_t> LHEPdfWeight;
    TRArray<Float_t> PSWeight;
    TRVariable<Float_t> PrefireWeight, PrefireWeight_up, PrefireWeight_down;

    nlohmann::json golden_json, prescale_json;

    std::unordered_map<std::string, size_t> trigger_idx;

    struct Prescale_Info {
        std::vector<size_t> lumis;
        std::vector<std::vector<size_t>> prescales;

        size_t distance(size_t lumi) {
            return std::distance(lumis.begin(), std::upper_bound(lumis.begin(), lumis.end(), lumi));
        }
    };

    WeightHolder pu_scale;
    WeightHolder ee_scale, em_scale, mm_scale;
    WeightHolder charge_misId, nonprompt_elec, nonprompt_muon;

    std::unordered_map<size_t, Prescale_Info> prescale_info;
};

#endif // __SCALEFACTORS_H_
