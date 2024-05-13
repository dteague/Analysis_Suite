#!/usr/bin/env python3
from skopt import BayesSearchCV
import pandas as pd
import numpy as np
from copy import copy
import uproot
from matplotlib import colors as clr
import boost_histogram as bh
from skopt.space import Real, Categorical, Integer
import xgboost as xgb
from sklearn import metrics

from analysis_suite.commons.configs import get_list_systs, get_inputs, get_ntuple_info
import analysis_suite.commons.user as user

def get_set(typ, indir, name='Signal', year='2018'):
    extra = "_signal" if typ == "test" else ""

    final_df = pd.DataFrame()
    filename = f"{typ}_Nominal.root"
    file_dir = indir/'train_files'
    with uproot.open(file_dir/filename) as f:
        for key, cls in f.classnames().items():
            if "/" in key or cls != "TTree":
                continue
            df = f[key].arrays(library='pd')
            classID = np.unique(df['classID'])[0]
            # print(key, np.sum(df.split_weight), len(df), np.unique(df.split_weight), np.sum(df.scale_factor), np.sum(df.split_weight)/np.sum(df.scale_factor))
            final_df = pd.concat([final_df, df])
    return final_df

output_base = user.analysis_area/"corr_test"
indir = output_base / "4top_first_Dec"
all_data = get_set("train", indir)
validate = get_set("validation", indir)

usevars = get_inputs(user.workspace_area/"October_run").usevars

X_train = all_data[usevars]
y_train = all_data.classID
X_test = validate[usevars]
y_test = validate.classID

from hyperopt import fmin, tpe, hp, STATUS_OK

# Define the hyperparameter space
space = {
    'max_depth': hp.quniform('max_depth', 2, 8, 1),
    'learning_rate': hp.loguniform('learning_rate', -5, -2),
    'subsample': hp.uniform('subsample', 0.5, 1)
}

# Define the objective function to minimize
def objective(params):
    default_params = {
        "eta": 0.09,
        'gamma': 0,
        'min_child_weight': 1e-6,
        'n_estimators': 500,
        'subsample': 0.7,
        'colsample_bytree': 0.75,
        'max_depth': 5,
        'objective': 'binary:logistic',
        'eval_metric': 'logloss',
        'use_label_encoder': False,
    }
    params.update(default_params)

    xgb_model = xgb.XGBClassifier(**params)
    xgb_model.fit(X_train, y_train)
    y_pred = xgb_model.predict(X_test)
    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred)
    auc = metrics.auc(fpr, tpr)
    print(auc)
    return {'loss': -auc, 'status': STATUS_OK}

# Perform the optimization
best_params = fmin(objective, space, algo=tpe.suggest, max_evals=100)
print("Best set of hyperparameters: ", best_params)
