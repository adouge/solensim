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

# debug
#import pysnooper

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

    def calc_dphi(self, phi2, phi1):
        """
        Get phase shift, accounting for cyclicity of phi; assume phis in units of pi
        """
        if np.round(phi1) > np.round(phi2):
            return phi2 - phi1 + 2
        else:
            return phi2 - phi1

    #@pysnooper.snoop()
    def process_states(self, s):
        """
        take astra output and add polar coordinates (r, phi) and respective impulses (pr, pphi);
        calculates phase change relative to start (turn) and previous state (dphi)
        updates particle z, pz, t with ref particle values
            --> transform into system static relative to field/"lab"; eliminates need for reference particle
        returns z positions,  new data
        """
        refs = s.query("particle==0")
        s = s.query("particle>0").copy()
        s.loc[:, "r"] = np.sqrt(s["x"].values**2 + s["y"].values**2)
# Here, assume only ref particles are particles with r=0
#        todrop = s.query("r==0").index
#        s = s.drop(todrop)
        s.loc[:, "cosphi"] = s["x"].values/s["r"].values
        s.loc[:, "sinphi"] = s["y"].values/s["r"].values
        s.loc[:, "phi"] = self.calc_phi_v(s["cosphi"].values, s["y"].values)
        zpos = s.index.levels[0]
        s.loc[zpos[0], "turn"] = 0
        s.loc[zpos[0], "dphi"] = 0
        i = 1
        while i < len(zpos):
            z = zpos[i]
            s.loc[z, "z"] = s.loc[z, "z"].values + refs.loc[z].get("z").values
            s.loc[z, "pz"] = s.loc[z, "pz"].values + refs.loc[z].get("pz").values
            s.loc[z, "t"] = s.loc[z, "t"].values + refs.loc[z].get("t").values
            s.loc[z, "turn"] = self.calc_dphi_v(s.loc[z, "phi"].values, s.loc[0, "phi"].values)
            s.loc[z, "dphi"] = self.calc_dphi_v(s.loc[z, "phi"].values, s.loc[zpos[i-1], "phi"].values)
            i += 1
        s.loc[:, "pr"] = s["cosphi"].values*s["px"].values + s["sinphi"].values*s["py"].values
        s.loc[:, "pphi"] = - s["sinphi"].values*s["px"].values + s["cosphi"].values*s["py"].values
#        refs = refs.drop(columns = ["q", "t", "type", "flag"])
#        s = s.drop(columns = ["sinphi", "cosphi", "q", "t", "type", "flag"])
        s = s.drop(columns = ["sinphi", "cosphi"])
        return zpos, s

    def monochrome_sweep(self, Emin, Emax, n=10):
        Es = np.linspace(Emin, Emax, num=n)
        pass
