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
import scipy.signal as signal
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

# State processing
    def calc_phi(self, x, y):
        """
         Get true angle in [0, 2pi] from cosine and y-coordinate, in uniits of pi (!)
        """
        phi = np.arctan2(y, x)/np.pi
        return 2 + phi if phi < 0 else phi

    def calc_dphi(self, phi2, phi1):
        """
        Get phase shift, accounting for cyclicity of phi; assume phis in units of pi, and only rotate in one direction
        """
        if np.sign(phi2-1) < np.sign(phi1-1):
            return phi2 - phi1 + 2
        else:
            return phi2 - phi1

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

        phi = self.calc_phi_v(s.get("x").values, s.get("y").values)
        s.loc[:, "phi"] = phi
        s.loc[:, "pr"] = np.cos(phi*np.pi)*s.loc[:, "px"].values + np.sin(phi*np.pi)*s.loc[:, "py"].values
        s.loc[:, "pphi"] = - np.sin(phi*np.pi)*s["px"].values + np.cos(phi*np.pi)*s["py"].values

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
            header = {
                "z" : z,
                "r_avg" : s.loc[z, "r"].mean(),
                "r_std" : s.loc[z, "r"].std(),
                "r_min" : s.loc[z, "r"].min(),
                "r_max" : s.loc[z, "r"].max(),
                "pr_avg" : s.loc[z, "pr"].mean(),
                "pr_std" : s.loc[z, "pr"].std(),
                "pr_min" : s.loc[z, "pr"].min(),
                "pr_max" : s.loc[z, "pr"].max(),
                "pphi_avg" : s.loc[z, "pphi"].mean(),
                "pphi_std" : s.loc[z, "pphi"].std(),
                "pphi_min" : s.loc[z, "pphi"].min(),
                "pphi_max" : s.loc[z, "pphi"].max(),
                "dphi_avg" : s.loc[z, "dphi"].mean(),
                "dphi_std" : s.loc[z, "dphi"].std(),
                "dphi_min" : s.loc[z, "dphi"].min(),
                "dphi_max" : s.loc[z, "dphi"].max(),
                "turn_avg" : s.loc[z, "turn"].mean(),
                "turn_std" : s.loc[z, "turn"].std(),
                "turn_min" : s.loc[z, "turn"].min(),
                "turn_max" : s.loc[z, "turn"].max()
            }
            headers[z] = header
        heads = pd.DataFrame.from_dict(headers).transpose()
        return heads

### Run variants & corresponding subroutines
# Overview run:
    def overview_run(self, E, label=None, beam_type="line"):
        """
            Do an overview run at a particular energy
        """
        self.msg("Doing overview ASTRA run for %.2f MeV electrons."%E)
        kinE_MeV = E - const.m_e/const.e*const.c**2/MeV
        self.astra.clean()
        self.astra.verbose = False
        self.astra.gen_preset = beam_type
        self.astra.genfile["input"]["ref_ekin"] = kinE_MeV
        self.astra.update_genfile()
        self.astra.generate()
        self.astra.track_preset = "overview"
        z_solenoid = self.astra.runfile["solenoid"]["s_pos"]
        self.astra.verbose = False
        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s = self.astra.read_states()
        s, zpos, parts, pref = self.process_states(s)
        heads = self.make_heads(s, zpos)
        z_neck, dz_neck, r_neck, sigma_neck = self.get_focal_region(heads)
        self.msg(("Beam neck approx. at z: (%.3f +- %.3f)m (+- focal region)"%(z_neck, dz_neck)))
        self.msg("Mean radius @ neck: %.3f mm, sigma %.3f mm"%(r_neck/mm, sigma_neck/mm))

        prefocal = s.query("z>@z_solenoid and z<(@z_neck - @dz_neck)")
        larmor = prefocal.get("turn").max()
        self.msg("Maximum Larmor angle: %.3f rad/pi"%larmor)

        result = {"E":E, "beam_preset":self.astra.beam_preset, "z_solenoid":z_solenoid, "z_neck":z_neck, "dz_neck":dz_neck, "max_larmor":larmor}
        data = {"s":s, "heads":heads, "zpos":zpos, "parts":parts, "pref":pref}

        if label==None:
            lbl = self._run_ticker
            self._run_ticker += 1
        else:
            lbl = label
            self._run_ticker += 1
        self.msg("Saving results under \"%s\"."%lbl)
        result = pd.DataFrame(result, index=[lbl])
        self.results = self.results.append(result)
        self.data[lbl] = data


    def get_focal_region(self, heads, method="pphi"):
        """
            todo
        """
        if method == "pr":  # Broken?
            z_solenoid = self.astra.runfile["solenoid"]["s_pos"]
            left = heads.query("z>@z_solenoid and pr_max == 0").get("z").max()
            right =  heads.query("z>@z_solenoid and pr_min == 0").get("z").min()

        elif method == "pphi":
            pphi = heads["pphi_avg"].values
            width_info = signal.peak_widths(np.abs(pphi), np.where(pphi==pphi.min())[0])
            lind = np.floor(width_info[2])
            rind = np.floor(width_info[3]) + 1
            left = heads["z"].values[int(lind)]
            right = heads["z"].values[int(rind)]

        z_neck = heads.get("r_avg").idxmin()
        r_neck = heads.loc[z_neck, "r_avg"]
        sigma_neck = heads.loc[z_neck, "r_std"]

        if not ((left <= z_neck) and (z_neck <= right)): raise(ValueError("Beam neck outside approximated focal region."))

        dz_neck = np.abs(np.array((left, right))-z_neck).max()

        return z_neck, dz_neck, r_neck, sigma_neck


    def focal_run(self):
        pass
