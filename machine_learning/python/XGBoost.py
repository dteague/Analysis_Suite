"""
.. module:: XGBoostMaker
   :synopsis: Takes in ROOT file to run a BDT training over it using XGBoost
.. moduleauthor:: Dylan Teague
"""
import numpy as np
import xgboost as xgb
from dataclasses import dataclass, InitVar, asdict

from .dataholder import MLHolder


def fom_metric(y_pred, dtrain):
    prev = fom_metric.prev
    nbins = 20
    bins = np.linspace(0, 1, nbins+1)

    y_true = dtrain.get_label()
    y_pred = 1/(1+np.exp(-y_pred))
    weight = dtrain.get_weight()
    b = np.histogram(y_pred[y_true==0], bins, weights=weight[y_true==0])[0]
    s = np.histogram(y_pred[y_true==1], bins, weights=weight[y_true==1])[0]
    fom = -np.sqrt(2*np.sum((s+b)*np.log(1+s/(b+1e-5))-s))
    # print(min(y_pred), max(y_pred))

    diff = abs((fom-prev)/fom)
    # print(fom, diff)
    # if diff > 0.35 and prev < -0.09:
    #     fom = prev
    fom_metric.prev = fom
    return 'fom', fom
fom_metric.prev = 0

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

    def train(self, outdir):
        """**Train for multiclass BDT**

        Does final weighting of the data (normalize all groups total
        weight to the same value), train the BDT, and fill the
        predictions ie the BDT values for each group.

        Returns:
          xgboost.XGBClassifer: XGBoost model that was just trained

        """

        x_train = self.train_set.drop(self._drop_vars, axis=1)
        # w_train = self.train_set.train_weight.copy()
        w_train = self.train_set.split_weight.copy()
        w_train2 = abs(self.train_set.scale_factor.to_numpy())
        y_train = self.train_set.classID

        # train_matrix = xgb.DMatrix(data=x_train,label=y_train, weight=w_train)
        # params = asdict(self.param)
        # n_rounds = params["n_estimators"]
        # del params["n_estimators"]
        # cv_train = xgb.cv(dtrain=train_matrix, params=params, nfold=5,
        #                   num_boost_round=n_rounds, early_stopping_rounds=20, stratified=True,
        #                   as_pandas=True, seed=123, verbose_eval=25,
        #                   feval=fom_metric)
        # print(cv_train.to_string())
        # exit()

        x_test = self.validation_set.drop(self._drop_vars, axis=1)
        y_test = self.validation_set.classID
        w_test = abs(self.validation_set.scale_factor.to_numpy())
        # w_test = self.validation_set.scale_factor

        _, group_tot = np.unique(y_train, return_counts=True)
        # w_train[y_train == 0] /= np.sum(w_train[y_train == 0])
        # w_test = self.validation_set.scale_factor.copy()
        # w_train[y_train == 1] /= np.sum(w_train[y_train == 1])

        # w_test[y_test == 0] /= np.sum(w_test[y_test == 0])
        # w_test[y_test == 1] /= np.sum(w_test[y_test == 1])

        # w_train[self.train_set["classID"] == 0] *= max(group_tot)/group_tot[0]
        # w_train[self.train_set["classID"] == 1] *= max(group_tot)/group_tot[1]

        fit_model = xgb.XGBClassifier(**asdict(self.param))
        fit_model.fit(x_train, y_train, sample_weight=w_train,
                      eval_metric=fom_metric,
                      eval_set=[(x_train, y_train), (x_test, y_test)],
                      sample_weight_eval_set=[w_train2, w_test],
                      early_stopping_rounds=1500, verbose=20)
        self.results = fit_model.evals_result()

        self.best_iter = fit_model.get_booster().best_iteration
        fit_model.save_model(f'{outdir}/model.bin')

    def predict(self, use_set, directory):
        fit_model = xgb.XGBClassifier({'nthread': 4})  # init model
        fit_model.load_model(str(directory / "model.bin"))  # load data
        return fit_model.predict_proba(use_set.drop(self._drop_vars, axis=1))

    def get_importance(self, directory):
        x_train = self.train_set.drop(self._drop_vars, axis=1)
        fit_model = xgb.XGBClassifier({'nthread': 4})  # init model
        fit_model.load_model(str(directory / "model.bin"))  # load data
        impor = fit_model.get_booster().get_score(importance_type= "total_gain")
        sorted_import = {x_train.columns[int(k[1:])]: v for k, v in sorted(impor.items(), key=lambda item: item[1]) }

        from analysis_suite.commons.plot_utils import plot, color_options
        with plot("{}/importance.png".format(directory)) as ax:
            ax.barh(range(len(sorted_import)), list(sorted_import.values()),
                    align='center',
                    height=0.5,)
            ax.set_yticks(range(len(sorted_import)))
            ax.set_yticklabels(sorted_import.keys(), fontsize=12)
            ax.set_xscale("log")
            ax.set_xlabel("Total Gain")
            ax.set_title("Variable Importance")

    def plot_training_progress(self, outdir):
        for typ in self.results['validation_0'].keys():
            from analysis_suite.commons.plot_utils import plot, color_options
            with plot(f"{outdir}/{typ}_training.png") as ax:
                ax.plot(self.results['validation_0'][typ], label='Training Set')
                ax.plot(self.results['validation_1'][typ], label='Validation Set')
                ax.legend()

    def approx_likelihood(self, var, bins, year, comb_bkg=True):
        use_set = self.test_sets[year]

        sig_mask = use_set.classID.astype(int) == 1
        sig_wgt = use_set.scale_factor[sig_mask]
        sig_var = use_set[var][sig_mask]
        sig_pred = self.pred_test[year]["Signal"][sig_mask]

        bkg_mask = use_set.classID.astype(int) != 1
        bkg_wgt = use_set.scale_factor[bkg_mask]
        bkg_var = use_set[var][bkg_mask]
        bkg_pred = self.pred_test[year]["Signal"][bkg_mask]

        max_fom = [0, -1]
        for i in np.linspace(0, 1, 21):
            sig = np.histogram(sig_var[sig_pred > i], bins=bins, weights=sig_wgt[sig_pred > i])[0]
            bkg = np.histogram(bkg_var[bkg_pred > i], bins=bins, weights=bkg_wgt[bkg_pred > i])[0]
            value = np.sqrt(np.sum(2*np.nan_to_num((sig+bkg)*np.log(1+sig/bkg) - sig)))
            if max_fom[0] < value:
                max_fom = [value, i]
        return max_fom
