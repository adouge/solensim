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


    def screens_trafo(self):
        """
        Read astra output and add polar coordinates and respective impulses;
        returns new data, reference particle and screen positions (indices)
        drops q, t, type and flags
        """
        s = self.astra.read_screens()
        refs = s.query("particle==0")
        s.loc[:, "r"] = abs(s["x"].values**2 + s["y"].values**2)**0.5
        todrop = s.query("r==0").index
        s = s.drop(todrop)
        def get_phi(cos, y):
            if y >=0: return np.arccos(cos)
            else: return np.arcsin(cos) + np.pi*3/2
        get_phi_v = np.vectorize(get_phi)
        s.loc[:, "cosphi"] = s["x"].values/s["r"].values
        s.loc[:, "sinphi"] = s["y"].values/s["r"].values
        s.loc[:, "phi"] = get_phi_v(s["cosphi"].values, s["y"].values)
        s.loc[:, "pr"] = s["cosphi"].values*s["px"].values + s["sinphi"].values*s["py"].values
        s.loc[:, "pphi"] = - s["sinphi"].values*s["px"].values + s["cosphi"].values*s["py"].values
        refs = refs.drop(columns = ["q", "t", "type", "flag"])
        s = s.drop(columns = ["sinphi", "cosphi", "q", "t", "type", "flag"])
        zpos = s.index.levels[0]
        return zpos, s, refs

    def monochrome_sweep(self, Emin, Emax, n=10):
        Es = np.linspace(Emin, Emax, num=n)
        pass
