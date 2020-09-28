#########################################################################
#    Copyright 2020 Anton Douginets
#    This file is part of solensim.
#
#    solensim is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    solensim is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with solensim.  If not, see <https://www.gnu.org/licenses/>.
#########################################################################

import scipy.constants as const
import numpy as np

from solensim.units import *

class TrackModule():
    def __init__(self, astra_interface):
        self.astra = astra_interface  # use Astra Frontend, same as the one provided in main script
        self.calc_phi_v = np.vectorize(self.calc_phi)
        self.calc_dphi_v = np.vectorize(self.calc_dphi)


    def calc_phi(self, cos, y):
        """
         Get true angle in [0, 2pi] from cosine and y-coordinate, in uniits of pi (!)
        """
        if y >=0: return np.arccos(cos)/np.pi
        else: return np.arcsin(cos)/np.pi + 3/2

    def calc_dphi(self, phi, phi0):
        """
        Get phase shift, accounting for cyclicity of phi; assume phis in units of pi
        """
        if phi0 > phi:
            return phi - phi0 + 2
        else:
            return phi - phi0

    def process_states(self, s):
        """
        take astra output and add polar coordinates and respective impulses;
        returns new data, reference particle and screen positions (indices)
        (drops q, t, type and flags - disabled)
        """
        refs = s.query("particle==0")
        s.loc[:, "r"] = np.sqrt(s["x"].values**2 + s["y"].values**2)
        todrop = s.query("r==0").index
        s = s.drop(todrop)
        s.loc[:, "cosphi"] = s["x"].values/s["r"].values
        s.loc[:, "sinphi"] = s["y"].values/s["r"].values
        s.loc[:, "phi"] = self.calc_phi_v(s["cosphi"].values, s["y"].values)
        s.loc[0, "dphi"] = 0
        zpos = s.index.levels[0]
        for z in zpos[1:]:
            s.loc[z, "dphi"] = self.calc_dphi_v(s.loc[z, "phi"].values, s.loc[0, "phi"].values)
        s.loc[:, "pr"] = s["cosphi"].values*s["px"].values + s["sinphi"].values*s["py"].values
        s.loc[:, "pphi"] = - s["sinphi"].values*s["px"].values + s["cosphi"].values*s["py"].values
#        refs = refs.drop(columns = ["q", "t", "type", "flag"])
#        s = s.drop(columns = ["sinphi", "cosphi", "q", "t", "type", "flag"])
        s = s.drop(columns = ["sinphi", "cosphi"])
        zpos = s.index.levels[0]
        return zpos, s, refs

    def monochrome_sweep(self, Emin, Emax, n=10):
        Es = np.linspace(Emin, Emax, num=n)
        pass
