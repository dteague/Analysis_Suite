#!/usr/bin/env python3

# First Training
params_first = {
    'colsample_bytree': 0.839967031847491, 'eta': 0.3571422057598092,
    'eval_metric': 'logloss', 'gamma': 3.442922836043457,
    'max_depth': 2, 'min_child_weight': 9.0,
    'n_estimators': 400, 'subsample': 0.563142636246572

    # 'colsample_bytree': 0.9171238285708372,
    # 'eta': 0.0246275336378838,
    # 'eval_metric': 'auc',
    # 'gamma': 5.9083936869037315,
    # 'max_depth': 4,
    # 'min_child_weight': 6,
    # 'n_estimators': 600,
    # 'subsample': 0.8307763305715123
}
# params_first = {
#     'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
#     'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
#     'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
# }
signal_first = "4top_sig"
groups_first = {
    "Signal": ['4top'], 'Background': ['tttj_nlo', 'tttw_nlo', 'ttt_nlo']
}


# Second Training
# cut = 0.865
cut = {
    '2018': 0.88,
    '2017': 0.86,
    '2016post': 0.80,
    '2016pre': 0.80,
}
params_second = {
    # 'colsample_bytree': 0.8004529026415877, 'eta': 0.03874070011624177,
    'colsample_bytree': 0.8004529026415877, 'eta': 0.01,
    'eval_metric': 'auc', 'gamma': 3.8585624071563256,
    'max_depth': 3, 'min_child_weight': 2.0,
    'n_estimators': 800, 'subsample': 0.5058560953790723
    # 'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
    # 'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
    # 'min_child_weight': 1.0, 'n_estimators': 800, 'subsample': 0.9
}
signal_second = "3top_sig"
groups_second = {
    "Signal": ['tttj_nlo', 'tttw_nlo', 'ttt_nlo'], "NotTrained": ['4top'],
}
