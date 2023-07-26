#!/usr/bin/env python3

import os
import logging
import numpy as np
from dataclasses import dataclass, InitVar
from sklearn.metrics import roc_auc_score, roc_curve, auc
from sklearn.model_selection import ShuffleSplit, cross_val_score

# os.environ["KERAS_BACKEND"] = "tensorflow"
# import keras
import tensorflow.keras as keras

from .dataholder import MLHolder

@dataclass
class Params:
    shuffle: bool = True
    initial_nodes: int = 60
    hidden_layers: int = 2
    node_pattern: str = "static"
    learning_rate: float = 0.002
    regulator: str = "dropout"
    activation: str = "elu"
    monitor: str = "loss"
    patience: int = 5
    mode: str = "auto"
    period: int = 1
    save_best_only: bool = True
    save_weights_only: bool = False
    epochs: int = 200
    batch_power: int = 14
    batch_size: int = 2**14
    validation_split: float = 0.25
    verbose: bool = False

    params: InitVar = None

    def __post_init__(self, params):
        if params is not None:
            for key, val in params.items():
                self.__setattr__(key, val)
            if "batch_power" in params:
                self.__setattr__("batch_size", params["batch_power"])

    def __getitem__(self, args):
        if isinstance(args, str):
            return {args: self.__getattribute__(args)}
        else:
            return {attr: self.__getattribute__(attr) for attr in args}


class KerasMaker(MLHolder):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.update_params(kwargs.get("params"))

    def update_params(self, params):
        self.build = Params(params=params)
        self.params = self.build["epochs", "batch_size", "shuffle", "validation_split"]
        self.early_stop = self.build["monitor", "patience"]
        self.checkpoint = self.build["save_best_only", "save_weights_only", "mode", "period"]


    def build_model(self, input_size="auto"):
        input_dim = len(self.all_vars) - len(self._drop_vars)
        node_lengths = self.build.initial_nodes * np.ones(
            self.build.hidden_layers, dtype=int
        )
        if self.build.node_pattern == "dynamic":
            partition = self.build.initial_nodes // self.build.hidden_layers
            node_lengths -= partition * np.arange(self.build.hidden_layers)

        model = keras.models.Sequential()
        model.add(keras.Input(input_dim))

        model.add(
            keras.layers.Dense(
                self.build.initial_nodes,
                kernel_initializer="glorot_normal",
                input_dim=input_dim,
                activation=self.build.activation,
            )
        )

        for nodes in node_lengths:
            model.add(keras.layers.BatchNormalization())
            if self.build.regulator == "dropout":
                model.add(keras.layers.Dropout(0.5))
            model.add(
                keras.layers.Dense(
                    nodes,
                    kernel_initializer="glorot_normal",
                    activation=self.build.activation,
                )
            )

        # Final classification node
        model.add(keras.layers.Dense(1, activation="sigmoid"))
        model.compile(
            optimizer=keras.optimizers.Adam(lr=self.build.learning_rate),
            loss="binary_crossentropy",
            metrics=["accuracy", "AUC"],
        )
        if self.build.verbose:
            model.summary()
        return model

    def train(self, outdir):
        shuf_train = self.train_set.sample(frac=1, random_state=0).reset_index(drop=True)

        x_train = shuf_train.drop(self._drop_vars, axis=1)
        w_train = shuf_train["train_weight"].to_numpy()
        y_train = shuf_train.classID

        x_test = self.validation_set.drop(self._drop_vars, axis=1)
        y_test = self.validation_set.classID

        _, group_tot = np.unique(y_train, return_counts=True)
        w_train[shuf_train["classID"] == 0] *= max(group_tot) / group_tot[0]
        w_train[shuf_train["classID"] == 1] *= max(group_tot) / group_tot[1]

        classW = {i: max(group_tot)/tot for i, tot in enumerate(group_tot)}
        callback = [
            keras.callbacks.EarlyStopping(**self.early_stop),
            # keras.callbacks.ModelCheckpoint("test.h5", verbose=0, **self.checkpoint),
        ]
        # Train
        logging.info("\n>> Training.")
        fit_model = self.build_model()
        history = fit_model.fit(
            x_train, y_train.astype(bool),
            # sample_weight=w_train,
            class_weight= classW,
            callbacks=callback,
            verbose=self.build.verbose,
            **self.params
        )

        # Test
        logging.info("\n>> Testing.")
        loss, accuracy, auc = fit_model.evaluate(x_test, y_test, verbose=1)
        logging.debug(f'loss: {loss}')
        logging.debug(f'accuracy: {accuracy}')
        logging.debug(f'auc: {auc}')

        fit_model.save(outdir / 'model.h5')

    def predict(self, use_set, directory):
        fit_model = keras.models.load_model(directory / 'model.h5')
        pred = fit_model.predict(use_set.drop(self._drop_vars, axis=1))
        return np.column_stack((1-pred, pred))
