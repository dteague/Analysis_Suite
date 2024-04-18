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
qcd = ['tth', 'ttz', '4top', 'ttXY']
ewk = ['ttw', 'xg', 'rare', 'tttw', 'rare']
Systematic.default_groups = mc_samples


systematics = [
    Systematic("LUMI_RUN2", "lnN", True).add(1.006, groups=mc_samples, year="2016pre") \
                                        .add(1.006, groups=mc_samples, year="2016post") \
                                        .add(1.009, groups=mc_samples, year="2017")    \
                                        .add(1.020, groups=mc_samples, year="2018"),
    Systematic("LUMI_17_18", "lnN", True).add(1.006, groups=mc_samples, year="2017") \
                                         .add(1.002, groups=mc_samples, year="2018"),
    Systematic("LUMI_", "lnN").add(1.010, groups=mc_samples, year="2016pre") \
                              .add(1.010, groups=mc_samples, year="2016post") \
                              .add(1.020, groups=mc_samples, year="2017")    \
                              .add(1.015, groups=mc_samples, year="2018"),

    #### Theory systs
    Systematic("PDF_unc", corr=True).add().dname("PDF"),
    Systematic("PDF_alphaZ", corr=True).add().dname("ALPHAS"),
    Systematic("PS_FSR", corr=True).add().dname("FSR"),

    Systematic("LHE_muF", corr=True).add(groups=qcd).dname("MUF_QCD"),
    Systematic("LHE_muF", corr=True).add(groups=ewk).dname("MUF_EWK"),
    Systematic("LHE_muR", corr=True).add(groups=qcd).dname("MUR_QCD"),
    Systematic("LHE_muR", corr=True).add(groups=ewk).dname("MUR_EWK"),
    Systematic("PS_ISR", corr=True).add(groups=qcd).dname("ISR_QCD"),
    Systematic("PS_ISR", corr=True).add(groups=ewk).dname("ISR_EWK"),

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
    Systematic("Muon_Scale").add().dname("ID_MU_"),
    Systematic("Electron_Scale").add().dname("ID_EL_"),
    Systematic("Lepton_tthMVA", "lnN").add(1.04, chan='Dilepton') # 1.02**2 \
                                      .add(1.061, chan=['Multi', 'ttzCR']), # 1.02**3
    Systematic("TriggerSF", 'lnN').add(1.03),
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
    Systematic("ChargeMisId_stat").add(groups='charge_flip'),
    Systematic("ChargeMisId_closure", 'lnN').add(1.2, groups="charge_flip"),

    Systematic("Nonprompt_Mu_stat").add(groups='nonprompt'),
    Systematic("Nonprompt_El_stat").add(groups='nonprompt'),
    Systematic("Nonprompt_closure", 'lnN').add(1.3, groups="nonprompt"),

    # Normalization stuff
    Systematic("CMS_norm_tttt", "lnN").add(1.2, groups="tttt"),
    Systematic("CMS_norm_ttw", "lnN").add(1.5, groups="ttw"),
    Systematic("CMS_norm_ttz", "lnN").add(1.5, groups="ttz"),
    Systematic("CMS_norm_tth", "lnN").add(1.5, groups="tth"),
    Systematic("CMS_norm_xg", "lnN").add(1.5, groups="xg"),
    Systematic("CMS_norm_rare", "lnN").add(1.5, groups="rare"),

    ##### Not used any more
    # # Systematic("Top_SF", "shape").add(1),
    # # Systematic("Jet_JES", "shape").add(1, groups=mc_samples),
    # # Systematic("BJet_BTagging", "shape").add(1, groups=mc_samples)
    # # Systematic("BJet_Eff", "shape").add(1, groups=mc_samples),
]

dummy = [Systematic("dummy", "lnN").add(1.0001, groups="rare"),]

def get_syst_name(systname, year):
    if systname == 'Nominal':
        return 'Nominal'
    for syst in systematics:
        if syst.name in systname:
            updown = systname[systname.rfind('_'):]
            return syst.get_name(year)+updown
    print('ERROR')
    return None

def get_shape_systs(year):
    systs = [(syst.name, syst.get_name(year)) for syst in systematics if syst.syst_type == "shape"]
    systs = np.vstack((np.char.add(systs, "_up"), np.char.add(systs, "_down")))
    return np.concatenate((systs, [["Nominal", "Nominal"]]))

