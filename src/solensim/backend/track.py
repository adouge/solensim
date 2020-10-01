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
import pandas as pd

from solensim.units import *

# debug
#import pysnooper

class TrackModule():
    def __init__(self, astra_interface):
        self.astra = astra_interface  # use Astra Frontend, same as the one provided in main script
        self.calc_phi_v = np.vectorize(self.calc_phi)
        self.calc_dphi_v = np.vectorize(self.calc_dphi)
        self.round_phase_v = np.vectorize(self.round_phase)

    def calc_phi(self, cos, y):
        """
         Get true angle in [0, 2pi] from cosine and y-coordinate, in uniits of pi (!)
        """

        if y >=0: return np.arccos(cos)/np.pi
        else: return np.arcsin(cos)/np.pi + 3/2

    def round_phase(self, phase):
        rounded = np.round(phase)
        if rounded == 2 or rounded == 0: return 0
        else: return 1

    def calc_dphi(self, phi2, phi1):
        """
        Get phase shift, accounting for cyclicity of phi; assume phis in units of pi
        """
        if self.round_phase_v(phi2) < self.round_phase_v(phi1):
            return phi2 - phi1 + 2
        else:
            return phi2 - phi1
        #return phi2 - phi1


    def process_states(self, s):
        """
        Process ASTRA beam state output
            s, zpos, parts, pref = process_states(s)
        returns:
            s - new dataframe (state-sorted)
            zpos - state z-coordinate list,
            parts - particle index list,
            pref - reference particle states
        Does:
            - add polar coordinates (r, phi) and respective impulses (pr, pphi);
            - calculates phase change relative to start (turn) and previous state (dphi)
            -updates particle z, pz, t with ref particle values
                --> transform into system static relative to field/"lab"; eliminates need for reference particle
            - drops charge, flag, type columns, but leaves them in pref for redundancy
        """

        pref = s.query("particle==0")
        N = len(s.index.levels[1])
        s = s.query("particle>0").copy()

        pref_broadcast = pd.concat([pref]*(N-1)).sort_index() # replicate each ref row N-1 times, since the ref particle is out
        s.loc[:, "z"] = s.loc[:, "z"].values + pref_broadcast.get("z").values
        s.loc[:, "pz"] = s.loc[:, "pz"].values + pref_broadcast.get("pz").values
        s.loc[:, "t"] = s.loc[:, "t"].values + pref_broadcast.get("t").values

        s.loc[:, "r"] = np.sqrt(s["x"].values**2 + s["y"].values**2)
        s.loc[:, "onaxis"] = s["r"].values==0
        modr = s.loc[:, "r"].values
        modr[s.loc[:, "onaxis"].values] = 1  # avoid 0-division
        cosphi = s.loc[:, "x"].values/modr
        sinphi = s.loc[:, "y"].values/modr

        s.loc[:, "phi"] = self.calc_phi_v(cosphi, s.loc[:, "y"].values)
        s.loc[:, "pr"] = cosphi*s.loc[:, "px"].values + sinphi*s.loc[:, "py"].values
        s.loc[:, "pphi"] = - sinphi*s["px"].values + cosphi*s["py"].values

        zpos = s.index.levels[0]
        z0 = zpos[0]
        s.loc[z0, "turn"] = 0
        s.loc[z0, "dphi"] = 0
        beam_start = s.loc[z0].copy()
        beam_start_broadcast = pd.concat([beam_start]*(len(zpos)-1))
        zs = zpos[1:]
        z0s = zpos[:-1]
        s.loc[zs, "turn"] = self.calc_dphi_v(s.loc[zs, "phi"].values, beam_start_broadcast.get("phi").values)
        s.loc[zs, "dphi"] = self.calc_dphi_v(s.loc[zs, "phi"].values, s.loc[z0s, "phi"].values)

        parts = np.arange(1, N, 1)
        todrop = ["type", "flag", "q"]
        s = s.drop(columns=todrop)


        return s, zpos, parts, pref

    def make_heads(self, s, zpos):
        """
            makes some "statistics" on each state
        """
        headers = {}
        for z in zpos:
            r_avg = s.loc[z, "r"].mean()
            r_std = s.loc[z, "r"].std()
            pr_avg = s.loc[z, "pr"].mean()
            pr_std = s.loc[z, "pr"].std()
            pphi_avg = s.loc[z, "pphi"].mean()
            pphi_std = s.loc[z, "pphi"].std()
            dphi_avg = s.loc[z, "dphi"].mean()
            dphi_std = s.loc[z, "dphi"].std()
            turn_avg = s.loc[z, "turn"].mean()
            turn_std = s.loc[z, "turn"].std()
            header = {
                "r_avg" : r_avg,
                "r_std" : r_std,
                "pr_avg" : pr_avg,
                "pr_std" : pr_std,
                "pphi_avg" : pphi_avg,
                "pphi_std" : pphi_std,
                "dphi_avg" : dphi_avg,
                "dphi_std" : dphi_std,
                "turn_avg" : turn_avg,
                "turn_std" : turn_std
            }
        headers[z] = header
        heads = pd.DataFrame.from_dict(headers).transpose()
        return heads

    def monochrome_sweep(self, Emin, Emax, n=10):
        Es = np.linspace(Emin, Emax, num=n)
        pass
