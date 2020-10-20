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
import scipy.optimize as opt
import numpy as np
import pandas as pd

from solensim.aux import *

# debug
#import pysnooper

class TrackModule():
    def __init__(self, astra_interface):
        self.astra = astra_interface  # use Astra Frontend, same as the one provided in main script

        self.calc_phi_v = np.vectorize(self.calc_phi)
        self.calc_dphi_v = np.vectorize(self.calc_dphi)

        self.E = 3.5  # default enrgy, placeholder
        self.N = 100
        self.sig_r = 10

# Astra setup

    def setup_tracking(self, bounds=[0,3], step=10):
        """
            Configure ASTRA to track within bounds=[left, right] with step = step [mm].
            Defaults: 0, 3m; step 1 cm
            Intended to follow track.use_field() and precede track.generate_beam(), followed by track.astra.run()
        """

        self.msg("Setting up tracking...")
        self.astra.clean()
        unit = step*mm
        self.z_solenoid = ceil(self.field_width/2, unit)
        self.astra.runfile["solenoid"]["s_pos"] = self.z_solenoid
        left = floor(bounds[0], unit)
        right = ceil(bounds[1], unit)
        self.astra.runfile["output"]["zstart"] = left
        self.astra.runfile["output"]["zstop"] = right
        zphase = np.ceil((right-left)/unit)
        self.astra.runfile["output"]["zphase"] = zphase
        self.astra.runfile["output"]["zemit"] = zphase
        self.msg("Tracking in: [%.2f, %.2f] m, %.2f mm step."%(left, right, step))
        self.astra.update_runfile()

    def generate_beam(self, E, N=100, sig_r=10):
        self.msg("Generating %.2f MeV, N=%d particle beam..."%(E, N))
        kinE_MeV = E - const.m_e/const.e*const.c**2/MeV
        self.astra.gen_preset = "line"
        self.astra.genfile["input"]["ref_ekin"] = kinE_MeV
        self.astra.genfile["input"]["sig_x"] = sig_r
        self.astra.genfile["input"]["ipart"] = N
        self.astra.genfile["input"]["q_total"] = N*const.e*10**(-9)  # nC charge units
        self.astra.update_genfile()
        self.astra.generate()

# Astra output processing
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
#                "r_min" : s.loc[z, "r"].min(),
#                "r_max" : s.loc[z, "r"].max(),
                "pr_avg" : s.loc[z, "pr"].mean(),
                "pr_std" : s.loc[z, "pr"].std(),
#                "pr_min" : s.loc[z, "pr"].min(),
#                "pr_max" : s.loc[z, "pr"].max(),
                "pphi_avg" : s.loc[z, "pphi"].mean(),
                "pphi_std" : s.loc[z, "pphi"].std(),
#                "pphi_min" : s.loc[z, "pphi"].min(),
#                "pphi_max" : s.loc[z, "pphi"].max(),
                "dphi_avg" : s.loc[z, "dphi"].mean(),
                "dphi_std" : s.loc[z, "dphi"].std(),
#                "dphi_min" : s.loc[z, "dphi"].min(),
#                "dphi_max" : s.loc[z, "dphi"].max(),
                "turn_avg" : s.loc[z, "turn"].mean(),
                "turn_std" : s.loc[z, "turn"].std(),
#                "turn_min" : s.loc[z, "turn"].min(),
#                "turn_max" : s.loc[z, "turn"].max()
            }
            headers[z] = header
        heads = pd.DataFrame.from_dict(headers).transpose()
        return heads

### the "Algorithm" components
## Init:
    def init_run(self, label, newrun):
        if newrun:
            message = "Init'ed new run at label: %s"
        else:
            message = "Loaded existing run at label: %s"
        self.run_label = label
        self.msg(message%label)
        self.linked_core.sample_field(self.field_z, self.field_Bz)
        self.linked_core.FM = "interpol"
        if newrun:
            F2 = self.linked_core.fint("lolkek", 2)
            f_predict = self.linked_core.get_f("lolkek", self.E)
            l = self.linked_core.get_fwhm("lolkek")
            self.runs.loc[label, "E"] = self.E
            self.runs.loc[label, "F2"] = F2
            self.runs.loc[label, "FWHM"] = l
            self.runs.loc[label, "f_predict"] = f_predict
            if f_predict < self.field_width/2:
                self.msg("[WARN] Predicted focus inside field definition boundary!")
                self.runs.loc[label, "warn_focus_in_field"] = True
            else:
                self.runs.loc[label, "warn_focus_in_field"] = False
        else:
            F2 = self.runs.loc[label, "F2"]
            f_predict = self.runs.loc[label, "f_predict"]
            l = self.runs.loc[label, "FWHM"]

        self.msg("F2: %.3f mT^2"%(F2/mm**2))
        self.msg("FWHM: %.2f cm"%(l/cm))
        self.msg("Predicted focal length for %.2f MeV: %.2f cm"%(self.E, f_predict/cm))

## Initial data collection:
    def overview_run(self, step=10):
        """
            TBA
        """
        label = self.run_label
        self.msg("%s : Performing overview tracking of %d x %.2f MeV electrons."%(label, self.N, self.E))
        self.setup_tracking(step=step)
        self.runs.loc[label, "z_solenoid"] = self.z_solenoid
        self.runs.loc[label, "field_width"] = self.field_width
        self.runs.loc[label, "step_overview"] = step
        self.generate_beam(E=self.E, N=self.N)
        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s = self.astra.read_states()
        s, zpos, parts, pref = self.process_states(s)

        self.msg("Saving results...")
        data = {"s":s, "zpos":zpos, "parts":parts, "pref":pref}
        self.data[label] = data

## Focal region estimation
    def get_focal_region(self, average=False):
        """
            todo
        """
        label = self.run_label
        s = self.data[label]["s"]
        z_solenoid = self.runs.loc[label, "z_solenoid"]
        self.msg("%s: estimating focal region. Use averaging : %s"%(label, average))
        if not average:
            # ACTHUNG works only with nice beams; averaging would help to avoid spontaneous "oops I crossed the axis" moments,
            # but takes more time.
            paraxial_prefocal = s.query("z>@z_solenoid and pr < 0").get("z").max()
            offaxis_postfocal =  s.query("z>@z_solenoid and pr > 0").get("z").min()
            unit = self.runs.loc[label, "step_overview"]*mm
            left = floor(offaxis_postfocal, unit)
            right = ceil(paraxial_prefocal, unit)
            if left > right:
                self.msg("[WARN] Paraxial trajectories seem to focus before off-axis ones. Reversing bound order.")
                self.runs.loc[label, "z_focal_left"] = right
                self.runs.loc[label, "z_focal_right"] = left
                self.runs.loc[label, "warn_bad_focal_bound_order"] = True
                self.msg("Focal region at ~ z: [%.2f, %.2f] cm; z_solenoid %.2f m"%(right/cm, left/cm, self.z_solenoid))
            else:
                self.runs.loc[label, "z_focal_left"] = left
                self.runs.loc[label, "z_focal_right"] = right
                self.runs.loc[label, "warn_bad_focal_bound_order"] = False
                self.msg("Focal region at ~ z: [%.2f, %.2f] cm; z_solenoid %.2f m"%(left/cm, right/cm, self.z_solenoid))
        else:
            wip()  # TODO, if needed

    # deprecated method
    def get_focal_region_from_heads(self, heads, method="pr"):
        """
            useless???
        """
        if method == "pr":
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
    def focal_run(self, step=0.1):
        """
            TBA
        """
        label = self.run_label
        self.msg("%s: Performing focal sweep..."%label)
        left = self.runs.loc[label, "z_focal_left"]
        right = self.runs.loc[label, "z_focal_right"]
        self.setup_tracking(bounds=[left, right], step=step)
        self.runs.loc[label, "step_focal"] = step


        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s_f = self.astra.read_states()
        s_f, zpos_f, parts, pref = self.process_states(s_f)

        self.msg("Saving results...")
        self.data[label]["s_f"] = s_f
        self.data[label]["zpos_f"] = zpos_f

    def ray_model(self, z, f, r_min, dxdz, dydz, cos0):
        return np.sqrt(r_min**2 + 2*r_min*(z-f)*(cos0*dxdz + np.sqrt(1-cos0**2)*dydz) + (z-f)**2*(dxdz**2 + dydz**2))

#Fitting trajectories @ focus:
    def focal_fitting(self):
        label = self.run_label

        self.msg("%s: Analyzing focal region data..."%label)
        s_f = self.data[label]["s_f"]
        parts = self.data[label]["parts"]
        z_solenoid = self.runs.loc[label, "z_solenoid"]

        fits = pd.DataFrame()
        p_f = s_f.swaplevel()

        for part in parts:
            r_min_idx = p_f.loc[part, "r"].idxmin()
            f_guess = p_f.loc[(part, r_min_idx), "z"]
            r_min_guess = p_f.loc[(part, r_min_idx), "r"]
            dxdz_guess = np.mean(p_f.loc[part, "x"].diff()/p_f.loc[part, "z"].diff())
            dydz_guess = np.mean(p_f.loc[part, "y"].diff()/p_f.loc[part, "z"].diff())
            cos0_guess = np.cos(p_f.loc[(part, r_min_idx), "phi"]*np.pi)

            p, pcov = opt.curve_fit(self.ray_model,
                xdata=p_f.query("z>@z_solenoid").loc[part, "z"].values,
                ydata=p_f.query("z>@z_solenoid").loc[part, "r"].values,
                p0=[f_guess, r_min_guess, dxdz_guess, dydz_guess, cos0_guess],
                bounds=([z_solenoid, 0, -np.inf, -np.inf, -1], [np.inf, self.sig_r*mm, np.inf, np.inf, 1]))
            dp = np.sqrt(np.diag(pcov))
            fits.loc[part, "r0"] =  p_f.query("zpos==0").loc[part,"r"].values
            columns = ["z_f", "r_min", "dxdz", "dydz", "cos0", "dz_f", "dr_min", "ddxdz", "ddydz", "dcos0"]
            fits.loc[part, columns] = [*p, *dp]
        f_prelim = fits.get("z_f").max()-self.runs.loc[label, "z_solenoid"]
        self.runs.loc[label, "f_prelim"] = f_prelim
        self.data[label]["fits"] = fits
        self.msg("Preliminary focal length value: %.2f cm"%(f_prelim*100))
