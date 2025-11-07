#!/usr/bin/env python3

info = {
    # TTT
    "tttj" : {
        "cross_section" : 0.0007,
        "Source of cross section" : "AN2018-062",
        "DAS": 'TTTJ_Tune.*',
    },
    "tttw" : {
        "cross_section" : 0.0013,
        "Source of cross section" : "AN2018-062",
        "DAS": 'TTTW_Tune.*',
    },
    "tttj_m" : {
        "cross_section" : 0.00049,
        "Source of cross section" : "https://github.com/gdurieux/triple-top-nlo",
        "DAS": 'TTTJminus-DR1_Tune.*',
    },
    "tttj_p" : {
        "cross_section" : 0.00021,
        "Source of cross section" : "https://github.com/gdurieux/triple-top-nlo",
        "DAS": 'TTTJplus-DR1_Tune.*',
    },
    "tttw_m" : {
        "cross_section" : 0.00066,
        "Source of cross section" : "https://github.com/gdurieux/triple-top-nlo",
        "DAS": 'TTTWminus-DR1_Tune.*',
    },
    "tttw_p" : {
        "cross_section" : 0.00066,
        "Source of cross section" : "https://github.com/gdurieux/triple-top-nlo",
        "DAS": 'TTTWplus-DR1_Tune.*',
    },

    # TTTT
    "tttt" : {
        "cross_section" : 13.37e-3,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1711.02116",
        "DAS": 'TTTT_Tune.*',
    },

    # TTXY
    "ttzh" : {
        "cross_section" : 0.001535,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'TTZH_Tune.*',
    },
    "tthh" : {
        "cross_section" : 0.000757,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'TTHH_Tune.*',
        "ISR": "QCD"
    },
    "ttwh" : {
        "cross_section" : 0.001582,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'TTWH_Tune.*',
    },
    "ttwz" : {
        "cross_section" : 0.002453,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'TTWZ_Tune.*',
    },
    "ttzz" : {
        "cross_section" : 0.001982,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'TTZZ_Tune.*',
    },
    "ttww" : {
        "cross_section" : 0.01150,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1610.0792",
        "DAS": 'TTWW_Tune.*',
    },

    # TTX
    "ttw" : {
        "cross_section" : 0.2316,
        "Level": "NLO+NNLL",
        "Source of cross section" : "arXiv:1812.08622",
        "DAS": 'TTWJetsToLNu_Tune.*'
    },
    "ttw_ewk" : {
        "cross_section" : 0.01127,
        "Source of cross section" : "arXiv:1812.08622",
        "DAS": 'ttWJetsToLNu_5f_EWK_Tune.*',
    },
    "ttz": {
        "cross_section" : 0.281,
        "Level": "NLO+NNLL",
        "Source of cross section" : "arXiv:1812.08622",
        "DAS": 'TTZToLLNuNu_M-10_Tune.*',
    },
    "ttz_m1-10": {
        "cross_section" : 0.0532,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1812.08622",
        "DAS": "TTZToLL_M-1to10_Tune.*",
    },
    "tth" : {
        "cross_section" : 0.209,
        "Level": "NLO+NNLL",
        "Source of cross section" : "LHCHXSWG",
        "DAS" : "ttHToNonbb_M125_Tune.*",
    },
    "tth_bb" : {
        "cross_section" : 0.291,
        "Level": "NLO+NNLL",
        "Source of cross section" : "LHCHXSWG",
        "DAS" : "ttHTobb_M125_Tune.*",
    },


    # VVV
    "www": {
        "cross_section" : 0.2109,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1405.0301",
        "DAS": 'WWW_4F_Tune.*',
    },
    "wwz" : {
        "cross_section" : 0.1679,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1405.0301",
        "DAS": 'WWZ_4F_Tune.*',
    },
    "wzz": {
        "cross_section" : 0.05550,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1405.0301",
        "DAS": 'WZZ_Tune.*'
    },
    "zzz" : {
        "cross_section" : 0.01417,
        "Level": "NLO",
        "Source of cross section" : "arXiv:1405.0301",
        "DAS": 'ZZZ_Tune.*',
    },
    "wwg" : {
        "cross_section" : 0.3369,
        "Level": "NLO",
        "Source of cross section" : "XSDB",
        "DAS": 'WWG_Tune.*',
    },
    "wzg" : {
        "cross_section" : 0.07876,
        "Level": "NLO",
        "Source of cross section" : "XSDB",
        "DAS": 'WZG_Tune.*',
    },

    # VV
    "wzTo3lnu" : {
        "cross_section" : 4.9173,
        "Level": "NNLO",
        "Source of cross section" : "https://indico.cern.ch/event/963958/contributions/4058869/attachments/2119607/3566972/2020_10_09_WZ_cross_section.pdf",
        "DAS": 'WZTo3LNu_Tune.*',
    },
    "zz4l": {
        "cross_section" : 1.256,
        "Level": "NNLO",
        "Source of cross section" : "arxiv:1405.2219",
        "DAS": 'ZZTo4L_Tune.*',
    },
    "ssww" : {
        "cross_section" : 0.02794,
        "Level": "NLO",
        "Source of cross section" : "arxiv:2005.01173",
        "DAS": 'SSWW_Tune.*',
    },
    "vh2nonbb" : {
        "cross_section" : 0.9561,
        "Level": "NLO",
        "Source of cross section" : "AN2023-082",
        "DAS": 'VHToNonbb_M125.*',
    },

    "ww_doubleScatter" : {
        "cross_section" : 0.2094,
        "Source of cross section" : "AN2018-062",
        "DAS": 'WWTo2L2Nu_.*DoubleScattering_.*',
    },
    "wpwpjj_ewk" : {
        "cross_section" : 0.03711,
        "Source of cross section" : "AN2018-062",
        "DAS": 'WpWpJJ_EWK-QCD_Tune.*',
    },

    ## VV Inclusive
    "zz" : {
        "cross_section": 16.19,
        "Level": "NNLO",
        "Source of cross section": "arXiv:1405.2219",
        "DAS": "ZZ_TuneC.*",
    },
    "wz" : {
        "cross_section": 51.11,
        "Level": "NNLO",
        "Source of cross section": "arxiv:1604.08576",
        "DAS": "WZ_TuneC.*",
    },
    "ww" : {
        "cross_section": 118.7,
        "Level": "NNLO",
        "Source of cross section": "arxiv:1408.5243",
        "DAS": "WW_Tune.*",
    },

    # X+g
    "wg" : {
        "cross_section" : 412.7,
        "Source of cross section" : "AN2018-062",
        "DAS": 'WGToLNuG_Tune.*',
    },
    "zg" : {
        "cross_section" : 51.53,
        'Level': 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'ZGToLLG_01J_5f_Tune.*',
    },
    "tg" : {
        "cross_section" : 2.997,
        'Level': 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'TGJets_Tune.*',
    },

    "ttg_dilep" : {
        "cross_section" : 2.22,
        "Level": 'NLO',
        'kfactor': 1.4852,
        "Source of cross section" : "AN2019_227",
        "DAS": 'TTGamma_Dilept_Tune.*',
    },
    "ttg_hadronic" : {
        "cross_section" : 4.178,
        "Level": 'NLO',
        'kfactor': 1.4852,
        "Source of cross section" : "AN2019_227",
        "DAS": 'TTGamma_Hadronic_Tune.*',
    },
    "ttg_singleLept" : {
        "cross_section" : 5.095,
        "Level": 'NLO',
        'kfactor': 1.4852,
        "Source of cross section" : "AN2019_227",
        "DAS": 'TTGamma_SingleLept_Tune.*',
    },

    # other
    "tzq" : {
        "cross_section" : 0.07561,
        "Level": 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'tZq_ll_4f_ckm_NLO_Tune.*',
    },
    "twz_wl" : {
        "cross_section" : 0.003004 ,
        "Level": 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'TWZToLL_thad_Wlept_5f_DR_Tune.*',
    },
    "twz_tl" : {
        "cross_section" : 0.003004 ,
        "Level": 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'TWZToLL_tlept_Whad_5f_DR_Tune.*',
    },
    "twz_twl" : {
        "cross_section" : 0.0015,
        "Level": 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'TWZToLL_tlept_Wlept_5f_DR_Tune.*',
    },
    "thq" : {
        "cross_section" : 0.0743,
        "Level": "QCD NLO",
        "Source of cross section" : "arXiv:1610.07922",
        "DAS": 'THQ_ctcvcp_4f_Hincl_Tune.*',
    },


    "wjets" : {
        # "cross_section" : 61334.9,
        "cross_section" : 61526.7,
        "Source of cross section" : "FEWZ",
        "DAS": 'WJetsToLNu_Tune.*',
    },
    "ggh2zz" : {
        "cross_section" : 0.01212,
        "Level": "NLO",
        "Source of cross section" : "XSDB",
        "DAS": 'GluGluHToZZTo4L_M125_Tune.*',
    },
    "st_twll" : {
        "cross_section" : 0.01123,
        "Level": 'NLO',
        "Source of cross section" : "XSDB",
        "DAS": 'ST_tW_Dilept_5f.*',
    },

    "data" : {
        "cross_section" : 1,
        "DAS": "dummy",
    },

    # TTBar samples
    "ttbar" : {
        "cross_section" : 833.9,
        "Level": "NNLO+NLL",
        "Source of cross section": "XSDB",
        "DAS": 'TT_Tune.*',
    },
    "ttbar_2l2n": {
        "cross_section" : 87.56,
        "Level": "NNLO+NLL",
        "Source of cross section": "XSDB",
        "DAS" : "TTTo2L2Nu_Tune.*",
    },
    "ttbar_semilep": {
        "cross_section" : 362.99,
        "Level": "NNLO+NLL",
        "Source of cross section": "XSDB",
        "DAS" : "TTToSemiLeptonic_Tune.*",
    },
    "ttbar_hadronic": {
        "cross_section" : 381.09,
        "Level": "NNLO+NLL",
        "Source of cross section": "XSDB",
        "DAS" : "TTToHadronic_Tune.*",
    },

    # Drell-Yan
    'DYm10-50': {
        "cross_section" : 15810.0,
        "Source of cross section" : "1G25ns",
        "Level": "NNLO",
        "DAS": 'DYJetsToLL_M-10to50_Tune.*',
    },
    'DYm50': {
        "cross_section" : 6077.22,
        "Source of cross section" : "FEWZ",
        "Level": "NNLO",
        "DAS": 'DYJetsToLL_M-50_Tune.*madgraphMLM-pythia8',
    },
    'DYm50_amc': {
        "cross_section" : 6077.22,
        "Source of cross section" : "FEWZ",
        "Level": "NNLO",
        "DAS": 'DYJetsToLL_M-50_Tune.*amcatnloFXFX-pythia8',
    },
    "DYm50_ht40-70" : {
        "cross_section": 311.4,
        "DAS" : "DYJetsToLL_M-50_HT-40to70_Tune.*",
    },
    "DYm50_ht70-100" : {
        "cross_section": 169.9, #145.5,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-70to100_Tune.*",
    },
    "DYm50_ht100-200" : {
        "cross_section": 147.4, #160.7,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-100to200_Tune.*",
    },
    "DYm50_ht200-400" : {
        "cross_section": 40.99, #48.63,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-200to400_Tune.*",
    },
    "DYm50_ht400-600" : {
        "cross_section": 5.678, #6.993,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-400to600_Tune.*",
    },
    "DYm50_ht600-800" : {
        "cross_section": 1.367,  #1.761,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-600to800_Tune.*",
    },
    "DYm50_ht800-1200" : {
        "cross_section": 0.6304, #0.8021,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-800to1200_Tune.*",
    },
    "DYm50_ht1200-2500" : {
        "cross_section": 0.1514, #0.1937,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-1200to2500_Tune.*",
    },
    "DYm50_ht2500-Inf" : {
        "cross_section": 0.003565, #0.003514,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_M-50_HT-2500toInf_Tune.*",
    },
    "DY_0J" : {
        "cross_section": 5090.0,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_0J_Tune.*",
    },
    "DY_1J" : {
        "cross_section": 983.5,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_1J_Tune.*",
    },
    "DY_2J" : {
        "cross_section": 353.6,
        "kfactor": 1.23,
        "DAS" : "DYJetsToLL_2J_Tune.*",
    },
    "DY_tautau" : {
        "cross_section": 353.6,
        "DAS" : "DYJetsToTauTau_.*",
    },


    # W + Jets
    "wjets_ht70-100" : {
        "cross_section": 1563,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-70To100_TuneC*",
    },
    "wjets_ht100-200": {
        "cross_section": 1627,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-100To200_TuneC*",
    },
    "wjets_ht200-400": {
        "cross_section": 435.237,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-200To400_TuneC*",
    },
    "wjets_ht400-600": {
        "cross_section": 59.181,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-400To600_TuneC*",
    },
    "wjets_ht600-800": {
        "cross_section": 14.581,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-600To800_TuneC*",
    },
    "wjets_ht800-1200": {
        "cross_section": 6.656,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-800To1200_TuneC*",
    },
    "wjets_ht1200-2500": {
        "cross_section": 1.681,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-1200To2500_TuneC*",
    },
    "wjets_ht2500-Inf": {
        "cross_section": 0.039,
        "Level": 'NLO',
        "Source of cross section": "1G25ns",
        "DAS" : "WJetsToLNu_HT-2500ToInf_TuneC*",
    },

    # QCD Mu enriched
    "qcd_mu15_pt20-Inf" : {
        "cross_section" : 239000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-20_MuEnrichedPt15_Tune.*',
    },
    "qcd_mu_pt15-20" : {
        "cross_section" : 2797000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-15To20_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt20-30" : {
        "cross_section" : 2518000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-20To30_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt30-50" : {
        "cross_section" : 1361000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-30To50_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt50-80" : {
        "cross_section" : 377800,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-50To80_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt80-120" : {
        "cross_section" : 88620,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-80To120_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt120-170" : {
        "cross_section" : 21070,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-120To170_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt170-300" : {
        "cross_section" : 7019,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-170To300_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt300-470" : {
        "cross_section" : 622.4,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-300To470_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt470-600" : {
        "cross_section" : 58.86,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-470To600_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt600-800" : {
        "cross_section" : 18.22,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-600To800_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt800-1000" : {
        "cross_section" : 3.25,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-800To1000_MuEnrichedPt5_Tune.*',
    },
    "qcd_mu_pt1000-Inf" : {
        "cross_section" : 1.078,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS": 'QCD_Pt-1000_MuEnrichedPt5_Tune.*',
    },

    # electron enriched
    "qcd_em_pt15-20" : {
        "cross_section" : 1324000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-15to20_EMEnriched_Tune.*",
    },
    "qcd_em_pt20-30" : {
        "cross_section" : 4896000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-20to30_EMEnriched_Tune.*",
    },
    "qcd_em_pt30-50" : {
        "cross_section" : 6447000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-30to50_EMEnriched_Tune.*",
    },
    "qcd_em_pt50-80" : {
        "cross_section" : 1988000,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-50to80_EMEnriched_Tune.*",
    },
    "qcd_em_pt80-120" : {
        "cross_section" : 367500,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-80to120_EMEnriched_Tune.*",
    },
    "qcd_em_pt120-170" : {
        "cross_section" : 66590,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-120to170_EMEnriched_Tune.*",
    },
    "qcd_em_pt170-300" : {
        "cross_section" : 16620,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-170to300_EMEnriched_Tune.*",
    },
    "qcd_em_pt300-Inf" : {
        "cross_section" : 1104,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt-300toInf_EMEnriched_Tune.*",
    },

    # QCD bc
    "qcd_bcToE_pt15-20" : {
        "cross_section" : 186500,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_15to20_bcToE_Tune.*",
    },
    "qcd_bcToE_pt20-30" : {
        "cross_section" : 308900,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_20to30_bcToE_Tune.*",
    },
    "qcd_bcToE_pt30-80" : {
        "cross_section" : 361800,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_30to80_bcToE_Tune.*",
    },
    "qcd_bcToE_pt80-170" : {
        "cross_section" : 34180.0,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_80to170_bcToE_Tune.*",
    },
    "qcd_bcToE_pt170-250" : {
        "cross_section" : 2109.0,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_170to250_bcToE_Tune.*",
    },
    "qcd_bcToE_pt250-Inf" : {
        "cross_section" : 568.2,
        "Level": "LO",
        "Source of cross section" : "XSDB",
        "DAS" : "QCD_Pt_250toInf_bcToE_Tune.*",
    },

    # Daniel stuff
    "ABCDNN" : {"cross_section" : 1,},
    "TOP" : {"cross_section" : 1,},
    "TTH" : {"cross_section" : 1,},
    "TTTT" : {"cross_section" : 1,},
}
