#!/usr/bin/env python3
import awkward as ak
import numpy as np
from copy import copy

from .basegetter import BaseGetter


class FlatGetter(BaseGetter):
    """ """

    def __init__(self, upfile, member):
        super().__init__()
        if member not in upfile or "TTree" not in repr(upfile[member]):
            return
        self.arr = upfile[member].arrays()
        self._base_mask = np.ones(len(self.arr), dtype=bool)
        self._mask = copy(self._base_mask)
        self._scale = ak.to_numpy(self.arr["scale_factor"])
        self.branches = self.arr.fields
        self.syst_name = "Nominal"
        self.correct_syst = True
        if "weights" in upfile:
            self.syst_weights = upfile[f'weights/{member}'].arrays()
            systs = self.list_systs()
            if len(systs) == 1:
                self.syst_name = systs[0]

    def __getitem__(self, key):
        if key not in self.branches:
            raise AttributeError(f"{key} not found")
        return np.nan_to_num(self.arr[key][self.mask], nan=-10000)

    def includes_syst(self, syst):
        return syst in self.syst_weights.fields

    def list_systs(self, input_systs=None):
        vals = self.syst_weights.fields
        vals.remove("index")
        if input_systs is None:
            return vals
        else:
            return np.array(input_systs)[np.in1d(input_systs, vals)]

    def set_systematic(self, syst):
        if syst in self.syst_weights.fields:
            self._scale = ak.to_numpy(self.syst_weights[syst])
            self.syst_name = syst
            self.correct_syst = True
        else:
            self.correct_syst = False
            self.syst_name = None

    @BaseGetter.mask.setter
    def mask(self, mask):
        """ """
        if isinstance(mask, str):
            mask = self.get_cut(mask)
        super(FlatGetter, type(self)).mask.fset(self, mask)

    def cut(self, cut):
        """Cut the dataframe base on some input cut

        Parameter
        ---------
        cut : string or array
            Cut string or a mask used for masking the whole dataframe
        """
        if cut is None:
            return
        elif isinstance(cut, str):
            super().cut(self.get_cut(cut))
        else:
            super().cut(cut)

    def get_cut(self, cut):
        """Get mask associated with a cut string

        Parameters
        ----------
        cut : string
            Cut string used to apply cuts to the dataframe (of form `var >/</== cutval`)

        Returns
        -------
        array
            Returns mask corresponding to the cut in the given cut string
        """
        if cut is None:
            return self.mask
        if len(cutter := cut.split("<")) > 1:
            return self[cutter[0]] < float(cutter[1])
        elif len(cutter := cut.split(">")) > 1:
            return self[cutter[0]] > float(cutter[1])
        elif len(cutter := cut.split("==")) > 1:
            return self[cutter[0]] == float(cutter[1])
        else:
            raise Exception(f"{cut} is not formatted correctly")
