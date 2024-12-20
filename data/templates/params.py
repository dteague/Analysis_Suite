#!/usr/bin/env python3

# First Training
params_first = {
    'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
    'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
    'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
}
signal_first = "4top_sig"
groups_first = {
    "Signal": ['4top'], 'Background': ['tttj_nlo', 'tttw_nlo']
}


# Second Training
cut = 0.93
params_second = {
    'colsample_bytree': 0.5283131479037887, 'eta': 0.037489501191992,
    'eval_metric': 'logloss', 'gamma': 5.439613578244588, 'max_depth': 3,
    'min_child_weight': 1.0, 'n_estimators': 500, 'subsample': 0.9
}
signal_second = "3top_sig"
groups_second = {
    "Signal": ['tttj_nlo', 'tttw_nlo'], "NotTrained": ['4top'],
}
