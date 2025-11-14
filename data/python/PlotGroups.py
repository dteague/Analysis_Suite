# -*- coding: utf-8 -*-
info ={
    "xg": {
        "Name": r"X+\gamma",
        "Members": [
            "ttg_hadronic",
	        "ttg_singleLept",
	        "ttg_dilep",
	        "zg",
            # "wg",
            "tg"
        ],
    },

    "tttt": {
        "Name": r"t\bar{t}t\bar{t}",
        "Members": ["tttt",]
    },
    "4top": {
        "Name": r"t\bar{t}t\bar{t}",
        "Members": ["tttt",],
        "Combine": 'TTTT',
    },
    "tttj_nlo": {
        "Name": r"t\bar{t}tj",
        "Members": ['tttj_m', 'tttj_p'],
        "Combine": 'TTTJ',
    },
    "tttj_p": {
        "Name": r"t\bar{t}tj^{+}",
        "Members": ['tttj_p'],
        "Combine": 'TTTJp',
    },
    "tttj_m": {
        "Name": r"t\bar{t}tj^{-}",
        "Members": ['tttj_m'],
        "Combine": 'TTTJm',
    },
    "tttj_lo": {
        "Name": r"t\bar{t}tj",
        "Members": ['tttj'],
        "Combine": 'TTTJ',
    },
    "tttw_nlo": {
        "Name": r"t\bar{t}tW",
        "Members": ['tttw_m', 'tttw_p'],
        "Combine": 'TTTW',
    },
    "tttw_p": {
        "Name": r"t\bar{t}tW^{+}",
        "Members": ['tttw_p'],
        "Combine": 'TTTWp',
    },
    "tttw_m": {
        "Name": r"t\bar{t}tW^{-}",
        "Members": ['tttw_m'],
        "Combine": 'TTTWm',
    },
    "tttw_lo": {
        "Name": r"t\bar{t}tW",
        "Members": ['tttw'],
        "Combine": 'TTTW',
    },
    "ttt_lo": {
        "Composite": True,
        "Name": "ttt",
        "Members": ["tttj_lo", "tttw_lo"],
        "Combine": 'TTT',
    },
    "ttt_nlo": {
        "Composite": True,
        "Name": "ttt",
        "Members": ["tttj_nlo", "tttw_nlo"],
        "Combine": 'SIG',
    },

    "ttw": {
        "Name": r"t\bar{t}W", 
        "Members": ["ttw", 'ttw_ewk']
    },
    "ttz": {
        "Name": r"t\bar{t}Z/\gamma*", 
        "Members": ["ttz", "ttz_m1-10"]
    },
    "tth": {
        "Name": r"t\bar{t}H", 
        "Members": ["tth", 'tth_bb']
    },

    "ttXY": {
        "Name": r"t\bar{t}VV",
        "Members": [
            "ttww",
	        "ttwz",
	        "ttzz",
	        "tthh",
            "ttzh",
	        "ttwh"
        ]
    },
    "vv_inc": {
        "Name": "VV",
        "Members": [
            "zz",
            "wz",
            "ww"
        ],
    },
    "vv": {
        "Name": "VV",
        "Members": [
            "vh2nonbb",
            "zz4l",
            "wpwpjj_ewk",
            "ww_doubleScatter",
            "wzTo3lnu",
            "ssww",
        ],
    },
    "vv_nowz": {
        "Name": "vv",
        "Members": [
            "ssww",
            "vh2nonbb",
            "zz4l",
            "wpwpjj_ewk",
            "ww_doubleScatter",
        ],
    },
    "vvv": {
        "Name": "vv",
        "Members": [
            "wwz",
            "wzz",
	        "www",
            "zzz",
            "wzg",
	        "wwg",
        ],
    },
    "vv+vvv": {
          "Composite": True,
        "Name": r"VV(V)",
        "Members": [
            "vvv",
            "vv_nowz",
        ]
    },
    "ttX": {
        "Composite": True,
        "Name": r"t\bar{t}X",
        "Members": [
            "ttw",
            "ttz",
            "tth"
        ]
    },
    'other': {
        "Composite": True,
        "Name": r'Other',
        "Members": ['xg', 'ttXY', 'rare'],
    },
    'twz': {
        "Members": [
            'twz_wl',
            'twz_tl',
            'twz_twl',
        ]
    },
    "rare": {
        "Composite": True,
        "Name": r"Rare",
        "Members": [
            'wz',
            'rare_nowz'
        ]
    },
    "rare_nowz": {
        "Composite": True,
        "Name": r"Rare",
        "Members": [
            "vvv",
            "vv_nowz",
            "st_twll",
	        "ggh2zz",
            'tzq',
            'twz',
            'thq',
            'ttXY',
        ]
    },
    "misc": {
        "Composite": True,
        "Name": r"Rare",
        "Members": [
            "st_twll",
	        "ggh2zz",
            'tzq',
            'twz',
            'thq',
        ]
    },
    'wz': {
        "Name": r"WZ",
        "Members": [
            "wzTo3lnu",
        ]
    },

    "qcd" : {
        "Composite": True,
        "Name": r"QCD",
        "Members": ["qcd_mu", "qcd_em"],
    },

    "qcd_mu" : {
        "Name": r"QCD",
        "Members": [
            'qcd_mu15_pt20-Inf',
            "qcd_mu_pt15-20",
            "qcd_mu_pt20-30",
            "qcd_mu_pt30-50",
            "qcd_mu_pt50-80",
            "qcd_mu_pt80-120",
            "qcd_mu_pt120-170",
            "qcd_mu_pt170-300",
            "qcd_mu_pt300-470",
            "qcd_mu_pt470-600",
            "qcd_mu_pt600-800",
            "qcd_mu_pt800-1000",
            "qcd_mu_pt1000-Inf",
        ]
    },
    "qcd_em" : {
        "Name": r"QCD",
        "Members": [
            "qcd_em_pt15-20",
            "qcd_em_pt20-30",
            "qcd_em_pt30-50",
            "qcd_em_pt50-80",
            "qcd_em_pt80-120",
            "qcd_em_pt120-170",
            "qcd_em_pt170-300",
            "qcd_em_pt300-Inf",
            "qcd_bcToE_pt15-20",
            "qcd_bcToE_pt20-30",
            "qcd_bcToE_pt30-80",
            "qcd_bcToE_pt80-170",
            "qcd_bcToE_pt170-250",
            "qcd_bcToE_pt250-Inf",
        ]
    },
    "ttjets_lep": {
        "Name": r"t\bar{t}",
        "Members": [
            "ttjets_dilep",
            "ttjets_single_t",
            "ttjets_single_tbar",
        ],
    },
    "ttbar_lep": {
        "Name": r"t\bar{t}",
        "Members" : [
            "ttbar_2l2n",
            "ttbar_semilep",
            "ttbar_hadronic",
        ],
    },
    "ttbar": {
        "Name": r"t\bar{t}",
        "Members" : [
            "ttbar",
        ],
    },
    "wjet_ht": {
        "Name": r"W+jets",
        "Members": [
            "wjets_ht70-100",
            "wjets_ht100-200",
            "wjets_ht200-400",
            "wjets_ht400-600",
            "wjets_ht600-800",
            "wjets_ht800-1200",
            "wjets_ht1200-2500",
            "wjets_ht2500-Inf"
        ],
    },
    "wjets" : {
        "Name": r"W+jets",
        "Members": [
            "wjets",
        ],
    },
    "ewk" : {
        "Name": r"EWK",
        "Members": [
            # "ttjet",
            'ttbar',
            "DYm50",
            "DYm50_amc",
            "DYm10-50",
            "wjets",
        ]
    },
    "ewk_lep" : {
        "Name": r"EWK",
        "Composite": True,
        "Members": [
            'ttbar_lep',
            "DYm50",
            "DYm50_amc",
            "DYm10-50",
            "wjets",
        ]
    },
    "DY_ht": {
        "Name": r"DY",
        "Members" : [
            # "DYm50_amc",
            "DYm50_ht40-70",
            "DYm50_ht70-100",
            "DYm50_ht100-200",
            "DYm50_ht200-400",
            "DYm50_ht400-600",
            "DYm50_ht600-800",
            "DYm50_ht800-1200",
            "DYm50_ht1200-2500",
            "DYm50_ht2500-Inf",
            "DYm10-50"
        ]
    },
    "DY": {
        "Name": r"DY",
        "Members" : [
            "DYm50_amc",
            "DYm50",
            "DYm10-50"
        ]
    },
    "DY_J": {
        "Name": r"DY",
        "Members" : [
            "DY_0J",
            "DY_1J",
            "DY_2J",
        ]
    },

    "data" : {
        "Name": r"Data",
        "Members": [
            "data"
        ],
    },
    'charge_flip': {
        "Name": r"Charge Misid",
        "Members": ["data"],
        "DataDriven": True,
    },
    'nonprompt': {
        "Name": r"Nonprompt",
        "Members": ["data"],
        "DataDriven": True,
    },
    'nonprompt_mc': {
        "Composite": True,
        "Name": r'Nonprompt MC',
        "Members": ['ttbar_lep', "wjet_ht", 'DY', 'DY_J', 'DY_ht'],
    },
    'nonprompt_nofr': {
        "Composite": True,
        "Name": r'Nonprompt MC raw',
        "Members": ['ttbar_lep', "wjet_ht", 'DY', 'DY_J', 'DY_ht'],
    },
    'data_driven': {
        "Composite": True,
        "Name": r'Data-Driven',
        "Members": ['nonprompt', "charge_flip",
                    'ttbar_lep', "wjet_ht", 'DY', 'DY_J', "DY_ht"], # mc nonprompt
    },

    'large_Xsec': {
        "Composite": True,
        "Name": r'Nonprompt',
        "Members": ['ttbar_lep', "wjet_ht", 'DY_ht', 'DY', 'DY_J'],
    },
###### Daniel stuff
    'abcdnn': {
         "Name": r'ABCDnn',
         'Combine': "ABCDNN",
         'Members': ['ABCDNN'],
    },
    'top': {
         "Name": r'ttW+ttZ+tt+VV',
         'Combine': "TOP",
         'Members': ['TOP'],
    },

    "all": {
        "Composite": True,
        "Name": "All",
        "Members": [
            "ttt",
            "tttt",
            "ttXY",
            "ttX",
            "xg",
            "vv",
            "vvv",
        ],
    }
}

def expand_composite(group):
    output = []
    if group not in info:
        output.append(group)
    elif 'DataDriven' in info[group]:
        output.append(group)
    else:
        is_composite = 'Composite' in info[group]
        for subgroup in info[group]['Members']:
            if is_composite and subgroup in info:
                output += expand_composite(subgroup)
            else:
                output.append(subgroup)
    return output

for group, subinfo in info.items():
    if 'Composite' in subinfo:
        subinfo['Members'] = expand_composite(group)

