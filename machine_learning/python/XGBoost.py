"""
.. module:: XGBoostMaker
   :synopsis: Takes in ROOT file to run a BDT training over it using XGBoost
.. moduleauthor:: Dylan Teague
"""
import numpy as np
import xgboost as xgb
import warnings
from dataclasses import dataclass, InitVar, asdict
from analysis_suite.commons.plot_utils import plot, cms_label

from .dataholder import MLHolder

formatter = {'extra_format': 'pdf',}

def fom_metric(y_pred, dtrain):
    # prev = fom_metric.prev
    nbins = 50
    y_true = dtrain.get_label()
    if len(y_pred.shape) > 1:
        y_pred = 1/(1+np.exp(-y_pred.T[1]))
    else:
        y_pred = 1/(1+np.exp(-y_pred))

    bins = np.linspace(np.min(y_pred), np.max(y_pred), nbins+1)
    weight = dtrain.get_weight()

    b = np.histogram(y_pred[y_true==0], bins, weights=weight[y_true==0])[0]
    s = np.histogram(y_pred[y_true==1], bins, weights=weight[y_true==1])[0]
    with np.errstate(invalid='ignore'):
        fom = -np.sqrt(2*np.sum((s+b)*np.log(np.where(b > 1e-5, 1+s/b, 1))-s))

    # diff = abs((fom-prev)/fom)
    # fom_metric.prev = fom
    return 'fom', fom
# fom_metric.prev = 0

@dataclass
class Params:
    eta: float = 0.09
    gamma: float = 0
    # reg_alpha: float = 0.0
    min_child_weight: float = 1e-6
    n_estimators: float = 500
    # reg_lambda: float = 0.05
    subsample: float = 1
    # base_score: float = 0.5
    # colsample_bylevel: int = 1
    colsample_bytree: int = 0.75
    # learning_rate: float = 1
    max_depth: int = 10
    # max_delta_step: int = 0
    objective: str = 'binary:logistic'
    eval_metric: str = "logloss"
    # use_label_encoder: bool = False
    # num_class: int = 2

    params: InitVar = None

    def __post_init__(self, params):
        if params is not None:
            for key, val in params.items():
                self.__setattr__(key, val)

    def __getitem__(self, args):
        if isinstance(args, str):
            return {args: self.__getattribute__(args)}
        else:
            return {attr: self.__getattribute__(attr) for attr in args}

class XGBoostMaker(MLHolder):
    def __init__(self, *args, **kwargs):
        """Constructor method
        """
        super().__init__(*args, **kwargs)
        self.param = Params(params=kwargs.get("params"))

    def update_params(self, params):
        self.param = Params(params=params)

    def csv(self):
        train_matrix = xgb.DMatrix(data=x_train,label=y_train, weight=w_train)
        params = asdict(self.param)
        n_rounds = params["n_estimators"]
        del params["n_estimators"]
        cv_train = xgb.cv(dtrain=train_matrix, params=params, nfold=5,
                          num_boost_round=n_rounds, early_stopping_rounds=20, stratified=True,
                          as_pandas=True, seed=123, verbose_eval=25,
                          feval=fom_metric)
        print(cv_train.to_string())

    def get_sets(self, workset):
        for sample, value in self.sample_map.items():
            if sample in self.nonTrained:
                workset = workset[workset.sampleName != value]
        x = workset[self.use_vars]
        y = workset.classID
        split_wgt = workset.split_weight.copy()
        scale = abs(workset.scale_factor.to_numpy())
        return x, y, split_wgt, scale

    def train(self, verbose=20):
        """**Train for multiclass BDT**

        Does final weighting of the data (normalize all groups total
        weight to the same value), train the BDT, and fill the
        predictions ie the BDT values for each group.

        Returns:
          xgboost.XGBClassifer: XGBoost model that was just trained

        """
        x_train, y_train, split_train, scale_train = self.get_sets(self.train_set)
        x_test, y_test, split_test, scale_test = self.get_sets(self.validation_set)

        sig_total = np.sum(split_train[y_train == 1])
        bkg_total = np.sum(split_train[y_train != 1])
        if sig_total > bkg_total:
            split_train[y_train == 1] *= bkg_total/sig_total
        else:
            split_train[y_train != 1] *= sig_total/bkg_total

        if len(np.unique(y_test)) > 2:
            self.param.objective = 'multi:softprob'
            self.param.num_class = len(np.unique(y_test))
            self.param.eval_metric = 'mlogloss'

        fit_model = xgb.XGBClassifier(**asdict(self.param))
        fit_model.fit(x_train, y_train, sample_weight=split_train,
                      eval_set=[(x_train, y_train), (x_test, y_test)],
                      sample_weight_eval_set=[scale_train, scale_test],
                      early_stopping_rounds=1500, verbose=verbose)
        self.results = fit_model.evals_result()
        self.best_iter = fit_model.get_booster().best_iteration
        fit_model.save_model(f'{self.outdir}/model.bin')

    def predict(self, use_set, directory):
        fit_model = xgb.XGBClassifier()  # init model
        fit_model.load_model(str(directory / "model.bin"))  # load data
        prediction = fit_model.predict_proba(use_set[self.use_vars])
        self.minmax[0] = np.min([np.min(prediction), self.minmax[0]])
        self.minmax[1] = np.max([np.max(prediction), self.minmax[1]])
        return prediction

    def get_importance(self, typ='total_gain'):
        x_train = self.train_set[self.use_vars]
        fit_model = xgb.XGBClassifier()  # init model
        fit_model.load_model(str(self.outdir / "model.bin"))  # load data
        impor = fit_model.get_booster().get_score(importance_type=typ)
        sorted_import = {x_train.columns[int(k[1:])]: v for k, v in sorted(impor.items(), key=lambda item: item[1]) }

        with plot(f"{self.outdir}/importance_{typ}.png") as ax:
            ax.barh(range(len(sorted_import)), list(sorted_import.values()),
                    align='center',
                    height=0.5,)
            ax.set_yticks(range(len(sorted_import)))
            ax.set_yticklabels(sorted_import.keys(), fontsize=12)
            ax.set_xscale("log")
            ax.set_xlabel(typ)
            ax.set_title("Variable Importance")

    def plot_training_progress(self):
        for typ in self.results['validation_0'].keys():
            with plot(self.outdir/f"{typ}_training.png", **formatter) as ax:
                ax.plot(self.results['validation_0'][typ], label='Training Set')
                ax.plot(self.results['validation_1'][typ], label='Validation Set')
                ax.set_xlabel("Training Epoch")
                ax.set_ylabel(typ)
                ax.legend()
                cms_label(ax, year='all')
