#!/usr/bin/env python3
import awkward as ak
import uproot as uproot
import numpy as np
from copy import copy

from .basegetter import BaseGetter


class NtupleGetter(BaseGetter):
    """ """

    single_branch = ["run", "event", "lumiBlock"]

    def __init__(self, filename, treename, group, xsec, syst=0, systName=""):
        super().__init__()
        self.syst = syst
        self.syst_bit = 2**self.syst
        self.systName = ""
        self.part_name = []
        self.arr = dict()
        self.parts = dict()
        self.branches = []

        f = uproot.open(filename)
        if group not in f or treename not in f[group]:
            return
        self.tree = f[group][treename]
        self.branches = [
            key for key, array in self.tree.items() if len(array.keys()) == 0
        ]
        self.part_name = np.unique(
            [br.split("/")[0] for br in self.branches if "/" in br]
        )
        self._base_mask = ak.to_numpy(self._get_var("PassEvent"))
        self._mask = copy(self._base_mask)
        self._scale = ak.to_numpy(self._get_var("weight"))
        self._setSyst(systName)

        if group == "data":
            _, unique_idx = np.unique(ak.to_numpy(self["event"]), return_index=True)
            self.mask = np.in1d(np.arange(len(self)), unique_idx)
            self._base_mask = copy(self._mask)
        elif "sumweight" in f[group]:
            sumweight, _ = f[group]["sumweight"].to_numpy()
            self._scale = xsec / sumweight[0] * self._scale
        else:
            print(f"Problem with group {group}")
            self._mask = None

    def _get_var(self, name):
        return self.tree[name].array()[:, self.syst]

    def _get_var_nosyst(self, name):
        return self.tree[name].array()

    def __getitem__(self, key):
        if key in self.part_name:
            if key not in self.parts:
                self.parts[key] = Particle(key, self)
            return self.parts[key]
        elif not self.exists(key):
            raise AttributeError(f"{key} not found")
        elif key not in self.arr:
            if "/" in key or key in self.single_branch:
                self.arr[key] = self._get_var_nosyst(key)
            else:
                self.arr[key] = self._get_var(key)
        return self.arr[key][self.mask]

    def _setSyst(self, systName):
        if systName in ["JES", "JER"]:
            jec = systName.lower().split("_")
            updown = {"down": "first", "up": "second"}
            self.jec_var = f"jec/{jec[1]}.{updown[jec[2]]}"
        else:
            self.jec_var = ""

    def reset(self):
        self.clear_mask()
        for part in self.parts.values():
            part.reset()

    def mergeParticles(self, merge, *parts):
        self.part_name = np.append(self.part_name, merge)
        self.parts[merge] = MergeParticle([self[part] for part in parts])

    def exists(self, branch):
        return branch in self.branches

    def close(self):
        self.tree.close()

    def dr(self, part1, idx1, part2, idx2):
        """Calculates the DeltaR between two particles in the event

        Parameters
        ----------
        part1 : string
            Name of the first particle
        idx1 : int
            Index of the first particle
        part2 : string
            Name of the second particle
        idx2 : int
            Index of the second particle

        Returns
        -------
        array
            Array of DeltaR between part1 and idx1 and part2 at idx2
        """
        deta = self[part1]["eta", idx1] - self[part2]["eta", idx2]
        dphi = self.dphi(part1, idx1, part2, idx2)
        return np.sqrt(deta**2 + dphi**2)

    def dphi(self, part1, idx1, part2, idx2):
        """Calculates the DeltaPhi (normalized) between two particles in the event

        Parameters
        ----------
        part1 : string
            Name of the first particle
        idx1 : int
            Index of the first particle
        part2 : string
            Name of the second particle
        idx2 : int
            Index of the second particle

        Returns
        -------
        array
            Array of DeltaPhi between part1 and idx1 and part2 at idx2
        """
        dphi = ak.to_numpy(self[part1]["phi", idx1] - self[part2]["phi", idx2])
        dphi[dphi > np.pi] = dphi[dphi > np.pi] - np.pi
        dphi[dphi < -np.pi] = dphi[dphi < -np.pi] + np.pi
        return dphi

    def cosDtheta(self, part1, idx1, part2, idx2):
        """Calculates the cos(DeltaTheta) between two particles in the event

        Parameters
        ----------
        part1 : string
            Name of the first particle
        idx1 : int
            Index of the first particle
        part2 : string
            Name of the second particle
        idx2 : int
            Index of the second particle

        Returns
        -------
        array
            Array of cos(DeltaTheta) between part1 and idx1 and part2 at idx2
        """
        cosh_eta = np.cosh(self[part1]["eta", idx1]) * np.cosh(self[part2]["eta", idx2])
        sinh_eta = np.sinh(self[part1]["eta", idx1]) * np.sinh(self[part2]["eta", idx2])
        cos_dphi = np.cos(self[part1]["phi", idx1] - self[part2]["phi", idx2])
        return (cos_dphi + sinh_eta) / cosh_eta

    def mass(self, part1, idx1, part2, idx2):
        """Calculates the mass between two particles in the event

        Parameters
        ----------
        part1 : string
            Name of the first particle
        idx1 : int
            Index of the first particle
        part2 : string
            Name of the second particle
        idx2 : int
            Index of the second particle

        Returns
        -------
        array
            Array of mass between part1 and idx1 and part2 at idx2
        """
        energy = self[part1]["energy", idx1] + self[part2]["energy", idx2]
        px = self[part1]["px", idx1] + self[part2]["px", idx2]
        py = self[part1]["py", idx1] + self[part2]["py", idx2]
        pz = self[part1]["pz", idx1] + self[part2]["pz", idx2]
        return np.sqrt(energy**2 - px**2 - py**2 - pz**2)

    # def top_mass(self, part1, part2):
    #     energy = self.combine(self[part1].energy(), self[part2].energy()) + self["Met"]
    #     px = self.combine(self[part1].px(), self[part2].px()) + self["Met"]*np.cos(self["Met_phi"])
    #     py = self.combine(self[part1].py(), self[part2].py()) + self["Met"]*np.sin(self["Met_phi"])
    #     pz = self.combine(self[part1].pz(), self[part2].pz())
    #     top = np.sqrt(energy**2 - px**2 - py**2 - pz**2)
    #     return top[ak.argmin(np.abs(top - 172.76), axis=-1, keepdims=True)][:, 0]


class ParticleBase:
    def __init__(self, vg):
        self.vg = vg

    def reset(self):
        pass

    def num(self):
        return ak.num(self.shape(), axis=-1)

    def scale(self, idx):
        """Get the event weight for the collection

        Parameters
        ----------
        idx : int
            Index of particle needed for scale. Masks the scale based on number of particles needed

        Returns
        -------
        array
            Scale masked to the number of particles needed by the index
        """
        if idx == -1:
            return ak.broadcast_arrays(self.vg.scale, self.shape())[0]
        else:
            return self.vg.scale[self.num() > idx]


    def get_hist(self, var, idx):
        """Get the values and scales to make a histogram for a given variable

        Parameters
        ----------
        var : string
            Variable or function name used to get a variable
        idx : int or (int, bool)
            Index of variable per event, e.g. 0 for 0th particle information. Can include bool for padding

        Returns
        -------
        tuple of (array, array)
            Gives values and scales packaged for use in histogram creation
        """
        if self.__len__() == 0:
            return ak.Array([]), ak.Array([])
        if idx == -1:
            return ak.flatten(self[var, idx]), ak.flatten(self.scale(idx))
        else:
            return self[var, idx], self.scale(idx)

    def get_hist2d(self, var1, var2, idx):
        """Get the values and scales to make a 2d histogram for the given variables

        Parameters
        ----------
        var1 : string
            Variable or function name used to get a variable
        var2 : string
            Variable or function name used to get a variable
        idx : int or (int, bool)
            Index of variable per event, e.g. 0 for 0th particle information. Can include bool for padding

        Returns
        -------
        tuple of form ((array, array), array)
            Gives values (tuple of the arrays of each variable) and scales packaged for use in histogram creation
        """
        if idx == -1:
            return (ak.flatten(self[var1, idx]), ak.flatten(self[var2, idx])), ak.flatten(self.scale(idx))
        else:
            return (self[var1, idx], self[var2, idx]), self.scale(idx)


class Particle(ParticleBase):
    """ """

    def __init__(self, name, vg):
        super().__init__(vg)
        self.name = name
        self._mask = (
            np.bitwise_and(vg._get_var_nosyst(f"{self.name}/syst_bitMap"), vg.syst_bit)
            != 0
        )

    def __getattr__(self, var):
        return self.vg[f"{self.name}/{var}"][self.mask]

    def __getitem__(self, idx):
        pad = False
        if len(idx) == 2:
            var, idx = idx
        else:
            var, idx, pad = idx

        if callable(getattr(self, var)):
            return getattr(self, var)(idx)
        else:
            return self._get_val(var, idx, pad)

    def __len__(self):
        return len(self.mask)

    def _get_val(self, var, idx, pad=False):
        if pad:
            vals = ak.pad_none(self.vg[f"{self.name}/{var}"][self.mask], idx + 1)
        elif idx == -1:
            return self.vg[f"{self.name}/{var}"][self.mask]
        else:
            vals = self.vg[f"{self.name}/{var}"][self.mask][self.num() > idx]
        return ak.to_numpy(vals[:, idx])

    def shape(self):
        return self.pt()

    # Special function for allowing jec in pt
    def pt(self, idx=-1):
        if self.vg.jec_var and "Jet" in self.name:
            return self._get_val("pt", idx) * self[f"{self.vg.jec_var}"]
        else:
            return self._get_val("pt", idx)

    @property
    def mask(self):
        return self._mask[self.vg.mask]

    # Functions for a particle

    def abseta(self, idx=-1):
        return np.abs(self["eta", idx])

    def num(self):
        return ak.to_numpy(ak.count(self.mask, axis=1))

    def px(self, idx=-1):
        return self["pt", idx] * np.cos(self["phi", idx])

    def py(self, idx=-1):
        return self["pt", idx] * np.sin(self["phi", idx])

    def pz(self, idx=-1):
        return self["pt", idx] * np.sinh(self["eta", idx])

    def energy(self, idx=-1):
        return np.sqrt(
            self["mass", idx] ** 2 + (self["pt", idx] * np.cosh(self["eta", idx])) ** 2
        )

    def mt(self, idx=-1):
        mask = self.num() > idx
        angle_part = 1 - np.cos(self["phi", idx] - self.vg["Met_phi"][mask])
        return np.sqrt(2 * self["pt", idx] * self.vg["Met"][mask] * angle_part)

    def pt_ratio(self, idx=-1):
        return 1 / (1 + self["ptRatio", idx])


class MergeParticle(ParticleBase):
    """ """

    def __init__(self, parts):
        super().__init__(parts[0].vg)
        self.parts = parts
        self.reset()

    def reset(self):
        self._idx_sort = ak.argsort(self._get_combined_item("pt"), ascending=False)

    def __getattr__(self, var):
        return self._get_combined_item(var)[self._sort]

    def __getitem__(self, idx):
        var, *idx = idx
        vals = self._get_combined_item(var)[self._sort]
        if len(idx) == 1:
            if idx[0] == -1:
                return vals
            else:
                return vals[ak.num(vals) > idx[0], idx[0]]
        else:
            idx, pad = idx
            return ak.to_numpy(ak.pad_none(vals, idx + 1)[:, idx])

    def __len__(self):
        return len(self._sort)

    def _get_combined_item(self, var):
        if callable(getattr(self.parts[0], var)):
            return ak.concatenate(
                (getattr(part, var)() for part in self.parts), axis=-1
            )
        else:
            return ak.concatenate(
                (part.__getattr__(var) for part in self.parts), axis=-1
            )
    def shape(self):
        return self._sort

    @property
    def _sort(self):
        return self._idx_sort[self.mask]

    @property
    def mask(self):
        """ """
        return self.vg.mask[self.vg._base_mask]
