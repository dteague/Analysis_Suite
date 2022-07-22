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
    Loose,
    Tight,
    Iso,
    FakeNotTight,
    Top,
    Bottom,
    Jet,
    Electron,
    Muon,
};

enum class Systematic {
    Nominal,
    LHE_muF,
    LHE_muR,
    PS_ISR,
    PS_FSR,
    PDF_alphaZ,
    PDF_unc,

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
    Top_SF,
    Pileup,
    Jet_JER,
    Jet_JES,
    Jet_PUID,

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
    { "Top_SF", Systematic::Top_SF },
    { "Pileup", Systematic::Pileup },
    { "Jet_JER", Systematic::Jet_JER },
    { "Jet_JES", Systematic::Jet_JES },
    { "Jet_PUID", Systematic::Jet_PUID },
};

const std::vector<Systematic> jec_systs = {
    Systematic::Jet_JER,
    Systematic::Jet_JES,
};

static const std::vector<eVar> syst_vars = { eVar::Up, eVar::Down };
static const std::vector<eVar> nominal_var = { eVar::Nominal };

static const std::unordered_map<eVar, std::string> varName_by_var = {
    { eVar::Nominal, "central" },
    { eVar::Up, "up" },
    { eVar::Down, "down" },
};

#endif // __COMMONENUMS_H_
