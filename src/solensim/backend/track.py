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
    def generate_beam(self, E, beam_type):
        kinE_MeV = E - const.m_e/const.e*const.c**2/MeV
        self.astra.gen_preset = beam_type
        self.astra.genfile["input"]["ref_ekin"] = kinE_MeV
        self.astra.update_genfile()
        self.astra.generate()

    def overview_run(self, E, label=None, beam_type="line", use_heads=False):
        """
            Do an overview run at a particular energy
        """
        self.msg("Doing overview ASTRA run for %.2f MeV electrons."%E)
        self.astra.clean()
        self.astra.verbose = False
        self.generate_beam(E, beam_type)
        self.astra.track_preset = "overview"
        z_solenoid = self.astra.runfile["solenoid"]["s_pos"]
        z, Bz = self.astra.read_field()
        field_width = z[-1]-z[0]
        self.astra.verbose = False
        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s = self.astra.read_states()
        s, zpos, parts, pref = self.process_states(s)

        if use_heads:
            self.msg("Calculating state headers...")
            heads = self.make_heads(s, zpos)
            left, right = self.get_focal_region(heads, from_heads=True)
        else:
            heads = None
            left, right = self.get_focal_region(s, from_heads=False)

        if left==right:
            step = (self.astra.runfile["output"]["zstop"]-self.astra.runfile["output"]["zstart"])/self.astra.runfile["output"]["zphase"]
            self.msg("WARN: left/right margins overlap, expanding by stepsize = %.3f m"%step)
            left-=step
            right+=step

        self.msg("Focal region at z: [%.3f; %.3f] m "%(left, right))

        if use_heads:
            prefocal = heads.query("z>@z_solenoid and  z < @left")
            postfocal = heads.query("z>@z_solenoid and  z > @right")
            larmor = prefocal.get("turn_avg").max()
            larmor2 = postfocal.get("turn_avg").max() - 1
        else:
            prefocal = s.query("z>@z_solenoid and  z < @left")
            postfocal = s.query("z>@z_solenoid and  z > @right")
            larmor = prefocal.get("turn").max()
            larmor2 = postfocal.get("turn").max() - 1

        self.msg("Maximum Larmor angle (pre-focus): %.3f rad/pi"%larmor)
        self.msg("Maximum Larmor angle (post-focus, minus focal flip): %.3f rad/pi"%larmor2)

        result = {
            "E":E, "beam_type":beam_type, "z_solenoid":z_solenoid, "field_width":field_width,
            "use_heads": use_heads,
            "f_left":left, "f_right":right,
            "larmor_left":larmor, "larmor_right":larmor2
            }
        data = {"s":s, "heads":heads, "zpos":zpos, "parts":parts, "pref":pref}

        if type(label)==type(None):
            self._run_ticker += 1
            lbl = self._run_ticker
        else:
            lbl = label
            self._run_ticker += 1
        self.msg("Saving results under \"%s\"."%lbl)
        result = pd.DataFrame(result, index=[lbl])
        self.results = self.results.append(result)
        self.data[lbl] = data
        self.astra.verbose=True


    def get_focal_region(self, data, from_heads=False):
        """
            todo
        """
        z_solenoid = self.astra.runfile["solenoid"]["s_pos"]
        step = (self.astra.runfile["output"]["zstop"]-self.astra.runfile["output"]["zstart"])/self.astra.runfile["output"]["zphase"]
        decimals = np.int(np.round(-np.log10(step)))

        if not from_heads:
            # ACTHUNG works only with nice beams; averaging would help to avoid spontaneous "oops I crossed the axis" moments,
            # but takes more time.
            paraxial_prefocal = data.query("z>@z_solenoid and pr < 0").get("z").max()
            offaxis_postfocal =  data.query("z>@z_solenoid and pr > 0").get("z").min()
            choose = np.round(np.array([paraxial_prefocal, offaxis_postfocal]), decimals=decimals)
            return choose.min(), choose.max()#, (paraxial_prefocal, offaxis_postfocal)
        else:
            return self.get_focal_region_from_heads(data)

    def get_focal_region_from_heads(self, heads, method="pr"):
        """
            todo
        """
        #if method == "pr":
        z_solenoid = self.astra.runfile["solenoid"]["s_pos"]
        left = heads.query("z>@z_solenoid and pr_avg < 0").get("z").max()
        right =  heads.query("z>@z_solenoid and pr_avg > 0").get("z").min()

        #elif method == "pphi":
        #    pphi = heads["pphi_avg"].values
        #    width_info = signal.peak_widths(np.abs(pphi), np.where(pphi==pphi.min())[0])
        #    lind = np.floor(width_info[2])
        #    rind = np.floor(width_info[3]) + 1
        #    left = heads["z"].values[int(lind)]
        #    right = heads["z"].values[int(rind)]

        #z_neck = heads["r_avg"].idxmax()
        #if not left <= z_neck <= right:
        #    if method != "pr":
        #        self.msg("WARN: Approximate beam neck not in determined focal region, falling back on pr method!")
        #        return self.get_focal_region_from_heads(heads, method="pr")
        #    else:
        #        self.msg("WARN: still not ok, returning potentially wrong results, be advised.")

        return left, right

#Focal run:
    def ray_model(self, z, f, drdz):  # move to core.Model?
        return np.abs((z-f)*drdz)

    def focal_run(self, label=None):
        self.astra.verbose=False
        self.astra.clean()
        if label == None:
            lbl = self._run_ticker
        else:
            lbl = label
        self.msg("Performing focal sweep at label: %s"%lbl)

        left = self.results.loc[lbl, "f_left"]
        right = self.results.loc[lbl, "f_right"]
        self.msg("Sweep region: [%.2f, %.2f]"%(left, right))

        E = self.results.loc[lbl, "E"]
        z_solenoid = self.results.loc[lbl, "z_solenoid"]
        beam_type = self.results.loc[lbl, "beam_type"]
        self.astra.clean()
        self.generate_beam(E, beam_type)

        self.astra.track_preset="overview"
        self.astra.runfile["solenoid"]["s_pos"] = z_solenoid
        self.astra.runfile["output"]["zstart"] = left
        self.astra.runfile["output"]["zstop"] = right
        zphase = np.round((right-left)/0.0001)
        self.astra.runfile["output"]["zphase"] = zphase
        self.astra.runfile["output"]["zemit"] = zphase
        self.astra.update_runfile()
        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s_f = self.astra.read_states()
        s_f, zpos_f, parts, pref = self.process_states(s_f)

        self.data[lbl]["s_focal"] = s_f
        self.data[lbl]["zpos_focal"] = zpos_f

        self.msg("Fitting particle trajectories (linear approximation!)...")

        foci = pd.DataFrame()
        p_f = s_f.swaplevel()
        f_guess = np.mean((left, right))
        drdz_guess = s_f.loc[zpos_f[-1], "z"].mean()/(zpos_f[-1]-f_guess)
        for part in parts:
            p, dp = self.linked_core.fit_to_model(self.ray_model,
                x=p_f.query("z>@z_solenoid").loc[part, "z"].values,
                y=p_f.query("z>@z_solenoid").loc[part, "r"].values,
                p0=(f_guess, drdz_guess))
            foci.loc[part, "r"] =  p_f.query("zpos==0").loc[part,"r"].values
            foci.loc[part, "f"] = p[0]
            foci.loc[part, "df"] = dp[0]
            foci.loc[part, "drdz"] = p[1]
            foci.loc[part, "ddrdz"] = dp[1]
        f_prelim = foci.get("f").max()-z_solenoid
        self.results.loc[lbl, "f_prelim"] = f_prelim
        self.data[lbl]["foci"] = foci
        self.msg("Preliminary focal length value: %.2f cm"%(f_prelim*100))
        self.astra.verbose=True
