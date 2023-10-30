#!/usr/bin/env python3
from analysis_suite.combine.systematics import Systematic
from analysis_suite.commons.constants import lumi
import numpy as np

ratio = {
    "2016": round(lumi['2016']/lumi['all'], 3),
    "2017": round(lumi['2017']/lumi['all'], 3),
    "2018": round(lumi['2018']/lumi['all'], 3),
}

mc_samples = ['ttt', 'xg', 'ttw', 'tth', 'ttz', 'ttXY', 'rare', '4top', 'tttj', 'tttw']

systematics = [
    Systematic("LUMI_RUN2", "lnN", True).add(1+0.006, groups=mc_samples, year="2016pre") \
                                        .add(1+0.006, groups=mc_samples, year="2016post") \
                                        .add(1+0.009, groups=mc_samples, year="2017")    \
                                        .add(1+0.020, groups=mc_samples, year="2018"),
    Systematic("LUMI_17_18", "lnN", True).add(1+0.006, groups=mc_samples, year="2017") \
                                         .add(1+0.002, groups=mc_samples, year="2018"),
    Systematic("LUMI_", "lnN").add(1+0.010, groups=mc_samples, year="2016pre") \
                             .add(1+0.010, groups=mc_samples, year="2016post") \
                             .add(1+0.020, groups=mc_samples, year="2017")    \
                             .add(1+0.015, groups=mc_samples, year="2018"),

    #### Theory systs
    Systematic("LHE_muF", "shape", True).add(1, groups=mc_samples) \
                                        .dname("MUF"),
    Systematic("LHE_muR", "shape", True).add(1, groups=mc_samples) \
                                        .dname("MUR"),
    Systematic("PDF_unc", "shape", True).add(1, groups=mc_samples) \
                                        .dname("PDF"),
    Systematic("PDF_alphaZ", "shape", True).add(1, groups=mc_samples)
                                           .dname("ALPHAS"),
    Systematic("PS_ISR", "shape", True).add(1, groups=mc_samples)
                                       .dname("ISR"),
    Systematic("PS_FSR", "shape", True).add(1, groups=mc_samples)
                                       .dname("FSR"),

    ### BTagging scales
    Systematic("BJet_Shape_lf", "shape", True).add(1, groups=mc_samples) \
                                        .dname("LF"),
    Systematic("BJet_Shape_lfstats1", "shape").add(1, groups=mc_samples) \
                                                    .dname("LFSTATS1"),
    Systematic("BJet_Shape_lfstats2", "shape").add(1, groups=mc_samples) \
                                                    .dname("LFSTATS2"),
    Systematic("BJet_Shape_hf", "shape", True).add(1, groups=mc_samples) \
                                        .dname("HF"),
    Systematic("BJet_Shape_hfstats1", "shape").add(1, groups=mc_samples) \
                                                    .dname("HFSTATS1"),
    Systematic("BJet_Shape_hfstats2", "shape").add(1, groups=mc_samples) \
                                                    .dname("HFSTATS2"),
    Systematic("BJet_Shape_cferr1", "shape", True).add(1, groups=mc_samples) \
                                            .dname("CFERR1"),
    Systematic("BJet_Shape_cferr2", "shape", True).add(1, groups=mc_samples)
                                            .dname("CFERR2"),

    #### Lepton/Jet ID stuff
    Systematic("Muon_Scale", "shape").add(1, groups=mc_samples) \
                                     .dname("ID_MU_"),
    Systematic("Electron_Scale", "shape").add(1, groups=mc_samples) \
                                         .dname("ID_EL_"),
    Systematic("Muon_tthMVA", "shape").add(1, groups=mc_samples),
    Systematic("Electron_tthMVA", "shape").add(1, groups=mc_samples),
    Systematic("Jet_PUID", "shape", True).add(1, groups=mc_samples) \
                                         .dname("PILEUPJETID"),

    #### JEC stuff
    Systematic("Jet_JER", "shape").add(1, groups=mc_samples) \
                                  .dname("JER"),
    Systematic("Jet_JEC_Absolute", "shape").add(1, groups=mc_samples) \
                                           .dname("JECABSOLUTE"),
    Systematic("Jet_JEC_Absolute_corr", "shape", True).add(1, groups=mc_samples) \
                                                      .dname("JECABSOLUTE"),
    Systematic("Jet_JEC_BBEC1", "shape").add(1, groups=mc_samples) \
                                        .dname("JECBBEC1"),
    Systematic("Jet_JEC_BBEC1_corr", "shape", True).add(1, groups=mc_samples) \
                                                   .dname("JECBBEC1"),
    Systematic("Jet_JEC_EC2", "shape").add(1, groups=mc_samples) \
                                      .dname("JECEC2"),
    Systematic("Jet_JEC_EC2_corr", "shape", True).add(1, groups=mc_samples) \
                                                 .dname("JECEC2"),
    Systematic("Jet_JEC_HF", "shape").add(1, groups=mc_samples) \
                                     .dname("JECHF"),
    Systematic("Jet_JEC_HF_corr", "shape", True).add(1, groups=mc_samples) \
                                                .dname("JECHF"),
    Systematic("Jet_JEC_AbsoluteBal", "shape", True).add(1, groups=mc_samples) \
                                                    .dname("JECRELATIVEBAL"),
    Systematic("Jet_JEC_AbsoluteSample", "shape").add(1, groups=mc_samples) \
                                                 .dname("JECRELATIVESAMPLE"),
    Systematic("Jet_JEC_FlavorQCD", "shape", True).add(1, groups=mc_samples) \
                                                  .dname("JECFLAVORQCD"),

    ##### Other
    Systematic("Prefire", "shape").add(1, groups=mc_samples) \
                                  .dname("PREFIRE"),
    Systematic("TriggerSF", "shape").add(1, groups=mc_samples),
    Systematic("Pileup", "shape", True).add(1, groups=mc_samples) \
                                       .dname("PILEUP"),

    ## Fake Rate stuff
    Systematic("ChargeMisId_stat", "shape").add(1, 'charge_flip'),
    Systematic("Nonprompt_Mu_stat", "shape").add(1, 'nonprompt'),
    Systematic("Nonprompt_El_stat", "shape").add(1, 'nonprompt'),
    Systematic("Nonprompt_closure", 'lnN').add(1.3, groups="nonprompt"),
    Systematic("ChargeMisId_closure", 'lnN').add(1.2, groups="charge_flip"),

    # Normalization stuff
    # Systematic("CMS_norm_tttt", "lnN").add(1.2, groups="tttt"),
    # Systematic("CMS_norm_ttw", "lnN").add(1.5, groups="ttw"),
    # Systematic("CMS_norm_ttz", "lnN").add(1.5, groups="ttz"),
    # Systematic("CMS_norm_tth", "lnN").add(1.5, groups="tth"),
    # Systematic("CMS_norm_xg", "lnN").add(1.5, groups="xg"),
    # Systematic("CMS_norm_rare", "lnN").add(1.5, groups="rare"),

    ##### Not used any more
    # # Systematic("Top_SF", "shape").add(1),
    # # Systematic("Jet_JES", "shape").add(1, groups=mc_samples),
    # # Systematic("BJet_BTagging", "shape").add(1, groups=mc_samples)
    # # Systematic("BJet_Eff", "shape").add(1, groups=mc_samples),
]

dummy = [Systematic("dummy", "lnN").add(1.0001, groups="rare"),]

def get_shape_systs(year):
    systs = [(syst.name, syst.get_name(year)) for syst in systematics if syst.syst_type == "shape"]
    systs = np.vstack((np.char.add(systs, "_up"), np.char.add(systs, "_down")))
    return np.concatenate((systs, [["Nominal", "Nominal"]]))
