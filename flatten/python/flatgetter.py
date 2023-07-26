#!/usr/bin/env python3
import awkward as ak
import numpy as np
from copy import copy

from .basegetter import BaseGetter


class FlatGetter(BaseGetter):
    """ """

    def __init__(self, upfile, member):
        super().__init__()
        if member not in upfile:
            return
        self.arr = upfile[member].arrays()
        self.syst_weights = upfile[f'weights/{member}'].arrays()
        self._base_mask = np.ones(len(self.arr), dtype=bool)
        self._mask = copy(self._base_mask)
        self._scale = ak.to_numpy(self.arr["scale_factor"])
        self.branches = self.arr.fields

    def __getitem__(self, key):
        if key not in self.branches:
            raise AttributeError(f"{key} not found")
        return np.nan_to_num(self.arr[key][self.mask], nan=-10000)

    def list_systs(self):
        vals = self.syst_weights.fields
        vals.remove("index")
        return vals

    def set_syst(self, syst):
        if syst in self.syst_weights.fields:
            self._scale = self.syst_weights[syst]

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
