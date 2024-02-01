#ifndef __COMMONENUMS_H_
#define __COMMONENUMS_H_

#include <unordered_map>
#include <vector>
#include <string>

const size_t MAX_PARTICLES = 65535;

enum class Year {
    yr2016pre,
    yr2016post,
    yr2017,
    yr2018,
    yrDefault,
};

static const std::unordered_map<Year, std::string> yearMap = {
    { Year::yr2016pre, "2016preVFP" },
    { Year::yr2016post, "2016postVFP" },
    { Year::yr2017, "2017" },
    { Year::yr2018, "2018" }
};

static const std::unordered_map<Year, std::string> yearNum = {
    { Year::yr2016pre, "16" },
    { Year::yr2016post, "16" },
    { Year::yr2017, "17" },
    { Year::yr2018, "18" }
};

enum class PID {
    Muon = 13,
    Electron = 11,
    Top = 6,
    Bottom = 5,
    Charm = 4,
    Jet = 0
};

enum class Level {
    Loose_NoPUID,
    Loose,
    Fake,
    Tight,
    FakeNotTight,
    LooseNotFake,
    Top,
    Bottom,
    Jet,
};

enum class Systematic {
    Nominal,
    LHE_muF,
    LHE_muR,
    PS_ISR,
    PS_FSR,
    PDF_alphaZ,
    PDF_unc,

    Prefire,
    TriggerSF,

    BJet_BTagging,
    BJet_Eff,
    BJet_Shape_hf,
    BJet_Shape_hfstats1,
    BJet_Shape_hfstats2,
    BJet_Shape_lf,
    BJet_Shape_lfstats1,
    BJet_Shape_lfstats2,
    BJet_Shape_cferr1,
    BJet_Shape_cferr2,

    Muon_Scale,
    Electron_Scale,
    Muon_tthMVA,
    Electron_tthMVA,
    Top_SF,
    Pileup,
    Jet_JER,
    Jet_JES,
    Jet_PUID,

    Jet_JEC_Absolute,
    Jet_JEC_Absolute_corr,
    Jet_JEC_BBEC1,
    Jet_JEC_BBEC1_corr,
    Jet_JEC_EC2,
    Jet_JEC_EC2_corr,
    Jet_JEC_HF,
    Jet_JEC_HF_corr,
    Jet_JEC_RelativeBal,
    Jet_JEC_RelativeSample,
    Jet_JEC_FlavorQCD,

    Electron_EScale,
    Electron_ESigma,

    ChargeMisId_stat,
    Nonprompt_Mu_stat,
    Nonprompt_El_stat,
};

enum class eVar {
    Nominal,
    Up,
    Down,
};

static const std::unordered_map<std::string, Systematic> syst_by_name = {
    { "LHE_muF", Systematic::LHE_muF },
    { "LHE_muR", Systematic::LHE_muR },
    { "PDF_unc", Systematic::PDF_unc },
    { "PDF_alphaZ", Systematic::PDF_alphaZ},
    { "PS_ISR", Systematic::PS_ISR },
    { "PS_FSR", Systematic::PS_FSR },
    { "Prefire", Systematic::Prefire },
    { "TriggerSF" , Systematic::TriggerSF },

    { "BJet_BTagging", Systematic::BJet_BTagging },
    { "BJet_Eff", Systematic::BJet_Eff },
    { "BJet_Shape_lf", Systematic::BJet_Shape_lf },
    { "BJet_Shape_lfstats1", Systematic::BJet_Shape_lfstats1 },
    { "BJet_Shape_lfstats2", Systematic::BJet_Shape_lfstats2 },
    { "BJet_Shape_hf", Systematic::BJet_Shape_hf },
    { "BJet_Shape_hfstats1", Systematic::BJet_Shape_hfstats1 },
    { "BJet_Shape_hfstats2", Systematic::BJet_Shape_hfstats2 },
    { "BJet_Shape_cferr1", Systematic::BJet_Shape_cferr1 },
    { "BJet_Shape_cferr2", Systematic::BJet_Shape_cferr2 },

    { "Muon_Scale", Systematic::Muon_Scale },
    { "Electron_Scale", Systematic::Electron_Scale },
    { "Muon_tthMVA", Systematic::Muon_tthMVA},
    { "Electron_tthMVA", Systematic::Electron_tthMVA},
    { "Top_SF", Systematic::Top_SF },
    { "Pileup", Systematic::Pileup },
    { "Jet_PUID", Systematic::Jet_PUID },

    { "Jet_JER", Systematic::Jet_JER },
    { "Jet_JES", Systematic::Jet_JES },
    { "Jet_JEC_Absolute", Systematic::Jet_JEC_Absolute },
    { "Jet_JEC_Absolute_corr", Systematic::Jet_JEC_Absolute_corr },
    { "Jet_JEC_BBEC1", Systematic::Jet_JEC_BBEC1},
    { "Jet_JEC_BBEC1_corr", Systematic::Jet_JEC_BBEC1_corr },
    { "Jet_JEC_EC2", Systematic::Jet_JEC_EC2},
    { "Jet_JEC_EC2_corr", Systematic::Jet_JEC_EC2_corr },
    { "Jet_JEC_HF", Systematic::Jet_JEC_HF},
    { "Jet_JEC_HF_corr", Systematic::Jet_JEC_HF_corr },
    { "Jet_JEC_AbsoluteBal", Systematic::Jet_JEC_RelativeBal},
    { "Jet_JEC_AbsoluteSample", Systematic::Jet_JEC_RelativeSample},
    { "Jet_JEC_FlavorQCD", Systematic::Jet_JEC_FlavorQCD},

    { "ChargeMisId_stat", Systematic::ChargeMisId_stat },
    { "Nonprompt_Mu_stat", Systematic::Nonprompt_Mu_stat },
    { "Nonprompt_El_stat", Systematic::Nonprompt_El_stat },
};

static const std::vector<Systematic> systs_that_change = {
    Systematic::Jet_JER,
    Systematic::Jet_JES,
    Systematic::Jet_JEC_Absolute,
    Systematic::Jet_JEC_Absolute_corr,
    Systematic::Jet_JEC_BBEC1,
    Systematic::Jet_JEC_BBEC1_corr,
    Systematic::Jet_JEC_EC2,
    Systematic::Jet_JEC_EC2_corr,
    Systematic::Jet_JEC_HF,
    Systematic::Jet_JEC_HF_corr,
    Systematic::Jet_JEC_RelativeBal,
    Systematic::Jet_JEC_RelativeSample,
    Systematic::Jet_JEC_FlavorQCD,
};

const std::vector<Systematic> jec_systs = {
    Systematic::Jet_JER,
    Systematic::Jet_JES,
    Systematic::Jet_JEC_Absolute,
    Systematic::Jet_JEC_Absolute_corr,
    Systematic::Jet_JEC_BBEC1,
    Systematic::Jet_JEC_BBEC1_corr,
    Systematic::Jet_JEC_EC2,
    Systematic::Jet_JEC_EC2_corr,
    Systematic::Jet_JEC_HF,
    Systematic::Jet_JEC_HF_corr,
    Systematic::Jet_JEC_RelativeBal,
    Systematic::Jet_JEC_RelativeSample,
    Systematic::Jet_JEC_FlavorQCD,
};

const std::vector<std::string> data_systs = {
    "ChargeMisId_stat",
    "Nonprompt_Mu_stat",
    "Nonprompt_El_stat",
};

static const std::vector<eVar> syst_vars = { eVar::Up, eVar::Down };
static const std::vector<eVar> nominal_var = { eVar::Nominal };
static const std::vector<eVar> all_vars = { eVar::Nominal, eVar::Up, eVar::Down };

static const std::unordered_map<eVar, std::string> varName_by_var = {
    { eVar::Nominal, "central" },
    { eVar::Up, "up" },
    { eVar::Down, "down" },
};

enum class Dataset {
    DoubleMuon,
    MuonEG,
    DoubleEG,
    Single_E,
    Single_M,
    MET,
    MHT,
    None,
};

static const std::unordered_map<std::string, Dataset> dataset_name_to_enum {
    {"DoubleMuon", Dataset::DoubleMuon},
    {"MuonEG", Dataset::MuonEG},
    {"DoubleEG", Dataset::DoubleEG},
    {"SingleElectron", Dataset::Single_E},
    {"SingleMuon", Dataset::Single_M},
    {"MET", Dataset::MET},
    {"MHT", Dataset::MHT},
    {"None", Dataset::None},
};

static const float ZMASS = 91.188;
static const float ZWINDOW = 15;
static const float LOW_ENERGY_CUT = 12;

#endif // __COMMONENUMS_H_
