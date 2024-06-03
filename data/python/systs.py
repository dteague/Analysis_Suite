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
mc_samples_norp = ['ttt', 'xg', 'ttw', 'tth', 'ttXY', 'rare', '4top', 'tttj', 'tttw']
mc_samples_nosig = ['xg', 'ttw', 'tth', 'ttz', 'ttXY', 'rare', '4top']
qcd = ['ttw', 'tth', 'ttz', '4top', 'ttXY', 'xg']
ewk = ['tttj', 'tttw', 'rare']
Systematic.default_groups = mc_samples


systematics = [
    Systematic("LUMI_RUN2", "lnN", True).add(1.006, groups=mc_samples_norp, year="2016pre") \
                                        .add(1.006, groups=mc_samples_norp, year="2016post") \
                                        .add(1.009, groups=mc_samples_norp, year="2017")    \
                                        .add(1.020, groups=mc_samples_norp, year="2018"),
    Systematic("LUMI_17_18", "lnN", True).add(1.006, groups=mc_samples_norp, year="2017") \
                                         .add(1.002, groups=mc_samples_norp, year="2018"),
    Systematic("LUMI_", "lnN").add(1.010, groups=mc_samples_norp, year="2016pre") \
                              .add(1.010, groups=mc_samples_norp, year="2016post") \
                              .add(1.020, groups=mc_samples_norp, year="2017")    \
                              .add(1.015, groups=mc_samples_norp, year="2018"),

    #### Theory systs
    Systematic("PDF_unc", corr=True).add(groups=mc_samples_nosig).dname("PDF"),
    Systematic("PDF_alphaZ", corr=True).add(groups=mc_samples_nosig).dname("ALPHAS"),
    Systematic("PS_FSR", corr=True).add().dname("FSR"),

    Systematic("LHE_muF", corr=True).add(groups=qcd).dname("MUFQCD"),
    Systematic("LHE_muF", corr=True).add(groups=ewk).dname("MUFEWK"),
    Systematic("LHE_muR", corr=True).add(groups=qcd).dname("MURQCD"),
    Systematic("LHE_muR", corr=True).add(groups=ewk).dname("MUREWK"),
    Systematic("PS_ISR", corr=True).add(groups=qcd).dname("ISRQCD"),
    Systematic("PS_ISR", corr=True).add(groups=ewk).dname("ISREWK"),

    ### BTagging scales
    Systematic("BJet_Shape_lf", corr=True).add().dname("LF"),
    Systematic("BJet_Shape_lfstats1").add().dname("LFSTATS1"),
    Systematic("BJet_Shape_lfstats2").add().dname("LFSTATS2"),
    Systematic("BJet_Shape_hf", corr=True).add().dname("HF"),
    Systematic("BJet_Shape_hfstats1").add().dname("HFSTATS1"),
    Systematic("BJet_Shape_hfstats2").add().dname("HFSTATS2"),
    Systematic("BJet_Shape_cferr1", corr=True).add().dname("CFERR1"),
    Systematic("BJet_Shape_cferr2", corr=True).add().dname("CFERR2"),

    #### Lepton/Jet ID stuff
    Systematic("Muon_Scale", corr=True).add().dname("ID_SS_MU"),
    Systematic("Electron_Scale", corr=True).add().dname("ID_SS_EL"),
    Systematic("Lepton_tthMVA_SS", "lnN", corr=True).add(1.04, chan='Dilepton') # 1.02**2 \
                                                    .add(1.061, chan=['Multi', 'ttzCR', 'ttttCR']), # 1.02**3
    Systematic("Trigger_SS", 'lnN').add(1.05),
    Systematic("Jet_PUID", corr=True).add().dname("PILEUPJETID"),

    #### JEC stuff
    Systematic("Jet_JER").add().dname("JER"),
    Systematic("Jet_JEC_Absolute").add().dname("JECABSOLUTE"),
    Systematic("Jet_JEC_Absolute_corr", corr=True).add().dname("JECABSOLUTE"),
    Systematic("Jet_JEC_BBEC1").add().dname("JECBBEC1"),
    Systematic("Jet_JEC_BBEC1_corr", corr=True).add().dname("JECBBEC1"),
    Systematic("Jet_JEC_EC2").add().dname("JECEC2"),
    Systematic("Jet_JEC_EC2_corr", corr=True).add().dname("JECEC2"),
    Systematic("Jet_JEC_HF").add().dname("JECHF"),
    Systematic("Jet_JEC_HF_corr", corr=True).add().dname("JECHF"),
    Systematic("Jet_JEC_AbsoluteBal", corr=True).add().dname("JECRELATIVEBAL"),
    Systematic("Jet_JEC_AbsoluteSample").add().dname("JECRELATIVESAMPLE"),
    Systematic("Jet_JEC_FlavorQCD", corr=True).add().dname("JECFLAVORQCD"),

    ##### Other
    Systematic("Prefire").add().dname("PREFIRE"),
    Systematic("Pileup", corr=True).add().dname("PILEUP"),

    ## Fake Rate stuff
    Systematic("ChargeMisId_Stat").add(groups='charge_flip'),
    Systematic("ChargeMisId_Closure", 'lnN').add(1.2, groups="charge_flip"),

    Systematic("Nonprompt_Mu_Stat").add(groups='nonprompt'),
    Systematic("Nonprompt_El_Stat").add(groups='nonprompt'),
    Systematic("Nonprompt_Closure", 'lnN', corr=True).add(1.3, groups="nonprompt"),

    # Normalization stuff
    Systematic("XSEC_TTTT", "lnN", corr=True).add('0.86/1.08', groups="4top"),
    # Systematic("XSEC_TTTW", "lnN", corr=True).add('0.85/1.16', groups="tttw"),
    # Systematic("XSEC_TTTJ", "lnN", corr=True).add('0.9/1.12', groups="tttj"),
    Systematic("XSEC_TTW", "lnN", corr=True).add(1.15, groups="ttw"),
    # Systematic("CMS_norm_ttz", "lnN").add(1.20, groups="ttz"),
    Systematic("XSEC_TTH", "lnN", corr=True).add(1.20, groups="tth"),
    Systematic("XSEC_XG", "lnN", corr=True).add(1.5, groups="xg"),
    Systematic("XSEC_RARE", "lnN", corr=True).add(1.5, groups="rare"),

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

