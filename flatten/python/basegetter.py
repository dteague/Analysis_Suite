#!/usr/bin/env python3
import numpy as np
from copy import copy
import awkward as ak


class BaseGetter:
    """ """

    def __init__(self):
        self._mask = None
        self._scale = None

    def __bool__(self):
        return self._mask is not None

    def __len__(self):
        return np.count_nonzero(self.mask)

    def __getattr__(self, key):
        return self.__getitem__(key)

    @property
    def scale(self):
        """ """
        return self._scale[self.mask]

    @scale.setter
    def scale(self, scale):
        """Way of changing the event weights multiplicatively

        Parameters
        ----------
        scale : int/array or tuple of int/array and array
            int/array that will be multiplied to the event scale or tuple with int/array to multiply to scale
            as well as the mask used to apply scale to. Mask must have as many True values as weights in the
            scaling array
        """
        if isinstance(scale, tuple):
            scale, mask = scale
            self._scale[self.get_submask(mask)] = (
                scale * self._scale[self.get_submask(mask)]
            )
        else:
            self._scale[self.mask] = scale * self._scale[self.mask]

    def get_submask(self, mask):
        """Helper function for getting a submask.
        Needed since the normal mask will have False values that reduce the overall size
        of the output variables meaning masks derived from variables will be smaller that
        the 'self.mask' length

        Parameters
        ----------
        mask : array
            Mask used to further mask the normal `self.mask` output

        Returns
        -------
        array
            Mask of length of `self.mask` but with correct masking needed
        """
        submask = copy(self.mask)
        submask[submask] = mask
        return submask

    @property
    def mask(self):
        """ """
        return self._mask

    @mask.setter
    def mask(self, mask):
        """ """
        if callable(mask):
            mask = ak.to_numpy(mask(self))
        self._mask[self._mask] *= mask

    def clear_mask(self):
        """ """
        self._mask = copy(self._base_mask)

    def cut(self, mask):
        """ """
        if callable(mask):
            mask = ak.to_numpy(mask(self))
        self._base_mask[self._base_mask] = mask
        self.reset()

    def get_graph(self, graph, *args):
        """
        Parameters
        ----------
        graph : GraphInfo object
            GraphInfo object that graphs the variable information
        """
        if graph.cuts is not None:
            self.mask = graph.cuts
        if callable(graph.func):
            return graph.func(self, *args)
        else:
            return self.get_hist(graph.func)

    def reset(self):
        self.clear_mask()

    def get_hist(self, var):
        """

        Parameters
        ----------
        var : string
            Name of variable that is to be grabbed

        Returns
        -------
        tuple of (array, array)
        """
        if not hasattr(self, var):
            raise AttributeError
        return self[var], self.scale
