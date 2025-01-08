#!/usr/bin/env python3
import awkward as ak
import uproot as uproot
import numpy as np
from copy import copy
import multiprocessing as mp

from .basegetter import BaseGetter

class NtupleGetter(BaseGetter):
    """ """

    single_branch = ["run", "event", "luminosityBlock", "bjet_scale"]

    def __init__(self, root_file, treename, group, xsec, systName="", **kwargs):
        super().__init__()
        self.part_name = []
        self.arr = dict()
        self.parts = dict()
        self.branches = []
        self.xsec = xsec
        self.isData = (group == "data")

        if group not in root_file or treename not in root_file[group]:
            return
        self.all_systs = []
        for name in root_file[group]["Systematics"]:
            mem = name.member("fName")
            if mem in self.all_systs:
                break
            else:
                self.all_systs.append(mem)

        self.syst_indices = {}
        if "Syst_Index" in root_file[group]:
            for name in root_file[group]["Syst_Index"]:
                i, j = int(name.member("fName")), int(name.member("fTitle"))
                if i in self.syst_indices:
                    break
                else:
                    self.syst_indices[i] = j

        self.sumw_hist, _ = root_file[group]["sumweight"].to_numpy()
        self.tree = root_file[group][treename]
        self.branches = [key for key, arr in self.tree.items() if len(arr.keys()) == 0]
        self.part_name = np.unique(
            [br.split("/")[0] for br in self.branches if "/" in br]
        )
        self.set_systematic(systName)

    def _get_var(self, name):
        return self.tree[name].array()[:, self.syst]

    def _get_var_nosyst(self, name):
        return self.tree[name].array()

    def _get_arr(self, name):
        if name not in self.arr:
            self[name]
        return self.arr[name]

    def __getitem__(self, key):
        if key in self.part_name:
            if key not in self.parts:
                self.parts[key] = Particle(key, self)
            return self.parts[key]
        elif not self.exists(key):
            raise AttributeError(f"{key} not found")
        elif key not in self.arr:
            if "/" in key or 'vector' not in self.tree[key].typename:
                self.arr[key] = self._get_var_nosyst(key)
            else:
                self.arr[key] = self._get_var(key)
        return self.arr[key][self.mask]

    def get_sf(self, systName):
        if self.xsec == 1:
            return 1.
        sumweight_change = {"LHE_muR_down": 2, "LHE_muF_down": 4, "LHE_muF_up": 6, "LHE_muR_up": 8,
                            "PDF_unc_down": 10, "PDF_unc_up": 11, "PDF_alphaZ_down": 12, "PDF_alphaZ_up": 13}
        sw_idx = sumweight_change.get(systName, 0)
        sumw = self.sumw_hist[sw_idx] if self.sumw_hist[sw_idx]/self.sumw_hist[0] > 0.1 else self.sumw_hist[0]
        return self.xsec/sumw

    def set_systematic(self, systname):
        if systname not in self.all_systs:
            print("HERE, setting nominal")
            systname = "Nominal"

        self.syst_unique = self.all_systs.index(systname)
        if self.syst_indices:
            self.syst = self.syst_indices[self.syst_unique]
            self.systNames = [(s, self.all_systs[s]) for s, n in self.syst_indices.items() if n == self.syst]
        else:
            self.systNames = []
            self.syst = self.syst_unique
        self.syst_bit = 2**self.syst

        # Redo scales and masks
        self._base_mask = ak.to_numpy(self._get_var("PassEvent"))
        self._mask = copy(self._base_mask)
        self._scale = ak.to_numpy(self.tree['weight'].array()[:, self.syst_unique])
        self.is_jec_unc = "JER" in systname or "JEC" in systname

        if not self.isData:
            self._scale = self.get_sf(systname) * self._scale

        for part in self.parts:
            self[part].reset_mask()

        for key in list(self.arr.keys()):
            if "/" not in key and 'vector' in self.tree[key].typename:
                del self.arr[key]

    def get_all_weights(self):
        all_weights = {}
        for i, systName in self.systNames:
            base_wgt = ak.to_numpy(self.tree['weight'].array()[:, i])
            scale = self.get_sf(systName)*base_wgt
            all_weights[systName] = scale[self.mask]
        return all_weights

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
        dphi[dphi > np.pi] = dphi[dphi > np.pi] - 2*np.pi
        dphi[dphi < -np.pi] = dphi[dphi < -np.pi] + 2*np.pi
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

    def dimass(self, part1, idx1, part2, idx2):
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

    def dipart_mt(self, part1, idx1, part2, idx2):
        ang_part = 1 - np.cos(self[part1]['phi', idx1] - self[part2]['phi', idx2])
        return np.sqrt(2*self[part1]['pt', idx1]*self[part2]['pt', idx2]*ang_part)

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

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d

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


    def get_hist(self, var, idx=-1):
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
        if var == "num":
            return self.num(), self.vg.scale
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
        self.reset_mask()

    def __getattr__(self, var):
        return self.vg[f"{self.name}/{var}"][self.mask]

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d

    def __call__(self, *args):
        return self[args]

    def __getitem__(self, args):
        pad = False
        idx = -1
        if len(args) == 1:
            var = args[0]
        elif len(args) == 2:
            var, idx = args
        else:
            var, idx, pad = args
        if callable(getattr(self, var)):
            return getattr(self, var)(idx, pad)
        else:
            return self._get_val(var, idx, pad)

    def __len__(self):
        return len(self.mask)

    def _get_val(self, var, idx=-1, pad=False):
        if pad:
            # vals = ak.where(self.num()>idx, self.vg[f"{self.name}/{var}"][self.mask, idx:idx+1]  pad)
            vals = ak.fill_none(ak.firsts(self.vg[f"{self.name}/{var}"][self.mask][:,idx:idx+1]), pad)
            return ak.to_numpy(vals)
            # vals = ak.fill_none(ak.pad_none(self.vg[f"{self.name}/{var}"][self.mask], idx + 1), pad)
        elif idx == -1:
            return self.vg[f"{self.name}/{var}"][self.mask]
        else:
            vals = self.vg[f"{self.name}/{var}"][self.mask][self.num() > idx]
        return ak.to_numpy(vals[:, idx])

    def reset_mask(self):
        self._base_mask = (
            np.bitwise_and(self.vg._get_var_nosyst(f"{self.name}/syst_bitMap"), self.vg.syst_bit)
            != 0
        )
        self.clear_mask()

    def clear_mask(self):
        self._mask = copy(self._base_mask)

    def shape(self):
        return self.pt()

    # Special function for allowing jec in pt
    def pt(self, *args):
        if self.vg.is_jec_unc and "Jet" in self.name:
            if f'{self.name}/pt_fix' not in  self.vg.arr:
                self.vg.arr[f'{self.name}/pt_fix'] = self.vg.tree[f"{self.name}/pt_shift"].array()[:, :, self.vg.syst-1]
                self.vg.branches.append(f'{self.name}/pt_fix')
            return self._get_val('pt_fix', *args)
        else:
            return self._get_val("pt", *args)

    # Special function for allowing jec in pt
    def mass(self, *args):
        if self.vg.is_jec_unc and "Jet" in self.name:
            if f'{self.name}/mass_fix' not in  self.vg.arr:
                self.vg.arr[f'{self.name}/mass_fix'] = self.vg.tree[f"{self.name}/mass_shift"].array()[:, :, self.vg.syst-1]
                self.vg.branches.append(f'{self.name}/mass_fix')
            return self._get_val('mass_fix', *args)
        else:
            return self._get_val("mass", *args)

    @property
    def mask(self):
        return self._mask[self.vg.mask]

    def mask_part(self, var, func):
        self._mask = func(self.vg._get_arr(f'{self.name}/{var}')) * self._mask

    # Functions for a particle

    def abseta(self, *args):
        return np.abs(self("eta", *args))

    def num(self):
        return ak.to_numpy(ak.count_nonzero(self.mask, axis=1))

    def px(self, *args):
        return self("pt", *args) * np.cos(self("phi", *args))

    def py(self, *args):
        return self("pt", *args) * np.sin(self("phi", *args))

    def pz(self, *args):
        return self('pt', *args) * np.sinh(self("eta", *args))

    def energy(self, *args):
        return np.sqrt(
            self("mass", *args) ** 2 + (self("pt", *args) * np.cosh(self("eta", *args))) ** 2
        )

    def mt(self, idx=-1, *args):
        mask = self.num() > idx
        angle_part = 1 - np.cos(self("phi", idx, *args) - self.vg["Met_phi"][mask])
        return np.sqrt(2 * self("pt", idx, *args) * self.vg["Met"][mask] * angle_part)

    def mt_fix(self, idx=-1, *args):
        mask = self.num() > idx
        angle_part = 1 - np.cos(self("phi", idx, *args) - self.vg["Met_phi"][mask])
        return np.sqrt(2 * 35 * self.vg["Met"][mask] * angle_part)

    def mt_lin(self, idx=-1, *args):
        mask = self.num() > idx
        met_lin = self.vg["Met"][mask]*np.cos(self("phi", idx, *args) - self.vg["Met_phi"][mask])
        return np.sqrt(2 * met_lin *(self("energy", idx, *args)-self("pt", idx, *args)))


class MergeParticle(ParticleBase):
    """ """

    def __init__(self, parts):
        super().__init__(parts[0].vg)
        self.parts = parts
        self.reset()

    def reset_mask(self):
       self.reset()

    def reset(self):
        self._idx_sort = ak.argsort(self._get_combined_item("pt"), ascending=False)

    def __getstate__(self):
        return self.__dict__

    def __setstate__(self, d):
        self.__dict__ = d

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
            return ak.to_numpy(ak.fill_none(ak.pad_none(vals, idx + 1)[:, idx], pad))

    def __len__(self):
        return len(self._sort)

    def _get_combined_item(self, var):
        if callable(getattr(self.parts[0], var)):
            return ak.concatenate(
                [getattr(part, var)() for part in self.parts], axis=-1
            )
        else:
            return ak.concatenate(
                [part.__getattr__(var) for part in self.parts], axis=-1
            )

    def shape(self):
        return self._sort

    def mask_part(self, var, func):
        for part in self.parts:
            part.mask_part(var, func)
        self.reset()

    @property
    def _sort(self):
        return self._idx_sort[self.mask]

    @property
    def mask(self):
        """ """
        return self.vg.mask[self.vg._base_mask]
