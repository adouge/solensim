"""solensim track module."""
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
import scipy.optimize as opt
import numpy as np
import pandas as pd

from solensim.aux import ceil, mm, cm, MeV, floor

# debug
# import pysnooper


class Model():
    """Model subblock."""

    def __init__(self):
        """Init model subblock."""
        pass

    def off_axis_trajectory(self, z, f, r_min, dxdz, dydz, cos0):
        """Model off-axis particle."""
        return np.sqrt(r_min**2 + 2*r_min*(z-f)*(cos0*dxdz + np.sqrt(1-cos0**2)*dydz) + (z-f)**2*(dxdz**2 + dydz**2))

    def axial_trajectory(self, z, f, drdz):
        """Model on-axis particle."""
        return np.abs(drdz*(z-f))

    def f_expansion(self, r, f, *args):
        """Model polynomial f expansion."""
        f_real = f
        order = len(args)
        for i in np.arange(order):
            f_real -= r**(2*(order-i))*args[-(i+1)]
        return f_real


class TrackModule():
    """Backend track block."""

    def __init__(self, astra_interface):
        """Init track block."""
        self.astra = astra_interface  # use Astra Frontend, same as the one provided in main script
        self.model = Model()

        self.calc_phi_v = np.vectorize(self.calc_phi)
        self.calc_dphi_v = np.vectorize(self.calc_dphi)

        self.E = 3.5  # default energy, placeholder
        self.N = 250
        self.sig_r = 10
        self.baseline_f = 1.25

# Astra setup

    def setup_tracking(self, bounds=[0, 3], step=10):
        """
        Configure tracking.

        Track within bounds=[left, right] with step = step [mm].
        Defaults: 0, 3m; step 1 cm
        Intended to follow track.use_field() and precede track.generate_beam(), followed by track.astra.run()
        """
        self.msg("Setting up tracking...")
        self.astra.clean()
        unit = step*mm
        self.z_solenoid = ceil(self.field_width/2, unit) + 0.25  # INFO: add a safety margin of 25 cm on the left
        self.astra.runfile["solenoid"]["s_pos"] = self.z_solenoid
        left = floor(bounds[0], unit)
        right = ceil(bounds[1], unit)
        self.astra.runfile["output"]["zstart"] = left
        self.astra.runfile["output"]["zstop"] = right
        zphase = np.ceil((right-left)/unit)
        self.astra.runfile["output"]["zphase"] = zphase
        self.astra.runfile["output"]["zemit"] = zphase
        self.msg("Tracking in: [%.2f, %.2f] m, %.2f mm step." % (left, right, step))
        self.astra.update_runfile()

    def generate_beam(self, E, N, sig_r, distribution="gauss", twodim=True):
        """Generate beam."""
        self.msg("Generating %.2f MeV, N=%d %s beam; sigma: %.1f mm, 2D - %s" % (E, N, distribution, sig_r, twodim))
        kinE_MeV = E - const.m_e/const.e*const.c**2/MeV
        self.astra.gen_preset = distribution
        self.astra.genfile["input"]["ref_ekin"] = kinE_MeV
        self.astra.genfile["input"]["sig_x"] = sig_r
        if twodim:
            self.astra.genfile["input"]["sig_y"] = sig_r
        self.astra.genfile["input"]["ipart"] = N
        self.astra.genfile["input"]["q_total"] = N*const.e*10**(-9)  # nC charge units
        self.astra.update_genfile()
        self.astra.generate()

# Astra output processing
    def calc_phi(self, x, y):
        """Get true angle in [0, 2pi] from cosine and y-coordinate, in uniits of pi (!)."""
        phi = np.arctan2(y, x)/np.pi
        return 2 + phi if phi < 0 else phi

    def calc_dphi(self, phi2, phi1):
        """Get phase shift, accounting for cyclicity of phi; assume phis in units of pi, and only rotate in one direction."""
        if np.sign(phi2-1) < np.sign(phi1-1):
            return phi2 - phi1 + 2
        else:
            return phi2 - phi1

    def process_states(self, s):
        """
        Process ASTRA beam state output.

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

        pref_broadcast = pd.concat([pref] * (N - 1)).sort_index()  # replicate each ref row N-1 times, since the ref particle is out
        s.loc[:, ["z_rel", "pz_rel"]] = s.loc[:, ["z", "pz"]].values
        s.loc[:, "z"] = s.loc[:, "z"].values + pref_broadcast.get("z").values
        s.loc[:, "pz"] = s.loc[:, "pz"].values + pref_broadcast.get("pz").values
        s.loc[:, "t"] = s.loc[:, "t"].values + pref_broadcast.get("t").values

        s.loc[:, "r"] = np.sqrt(s["x"].values**2 + s["y"].values**2)
        s.loc[:, "onaxis"] = s["r"].values == 0

        phi = self.calc_phi_v(s.get("x").values, s.get("y").values)
        s.loc[:, "phi"] = phi
        s.loc[:, "pr"] = np.cos(phi * np.pi) * s.loc[:, "px"].values + np.sin(
            phi * np.pi) * s.loc[:, "py"].values
        s.loc[:, "pphi"] = - np.sin(phi * np.pi) * s["px"].values + np.cos(
            phi * np.pi) * s["py"].values

        zpos = s.index.levels[0]
        z0 = zpos[0]
        s.loc[z0, "turn"] = 0
        s.loc[z0, "dphi"] = 0
        beam_start = s.loc[z0].copy()
        beam_start_broadcast = pd.concat([beam_start] * (len(zpos) - 1))
        zs = zpos[1:]
        z0s = zpos[:-1]
        s.loc[zs, "turn"] = self.calc_dphi_v(s.loc[zs, "phi"].values, beam_start_broadcast.get("phi").values)
        s.loc[zs, "dphi"] = self.calc_dphi_v(s.loc[zs, "phi"].values, s.loc[z0s, "phi"].values)

        parts = np.arange(1, N, 1)
        todrop = ["type", "flag", "q"]
        s = s.drop(columns=todrop)

        return s, zpos, parts, pref

    def make_heads(self, s, zpos):
        """make some "statistics" on each state."""
        headers = {}
        for z in zpos:
            header = {
                "z": z,
                "r": s.loc[z, "r"].mean(),
                # "r_std" : s.loc[z, "r"].std(),
                "pr": s.loc[z, "pr"].mean(),
                # "pr_std" : s.loc[z, "pr"].std(),
                "pphi": s.loc[z, "pphi"].mean(),
                # "pphi_std" : s.loc[z, "pphi"].std(),
                "dphi": s.loc[z, "dphi"].mean(),
                # "dphi_std" : s.loc[z, "dphi"].std(),
                "turn": s.loc[z, "turn"].mean(),
                # "turn_std" : s.loc[z, "turn"].std(),
                # "r_min" : s.loc[z, "r"].min(),
                # "r_max" : s.loc[z, "r"].max(),
                # "pr_min" : s.loc[z, "pr"].min(),
                # "pr_max" : s.loc[z, "pr"].max(),
                # "pphi_min" : s.loc[z, "pphi"].min(),
                # "pphi_max" : s.loc[z, "pphi"].max(),
                # "dphi_min" : s.loc[z, "dphi"].min(),
                # "dphi_max" : s.loc[z, "dphi"].max(),
                # "turn_min" : s.loc[z, "turn"].min(),
                # "turn_max" : s.loc[z, "turn"].max()
            }
            headers[z] = header
        heads = pd.DataFrame.from_dict(headers).transpose()
        return heads

# the "Algorithm" components
# Init:
    def init_run(self, label, rel_decrement=0):
        """Initialize new run."""
        message = ">>> Initiating new run at label: %s"
        self.run_label = label
        self.msg(message % label)
        self.linked_core.sample_field(self.field_z, self.field_Bz)
        F2 = self.linked_core.fint(2)
        f_predict = self.linked_core.get_f(self.E)
        fwhm = self.linked_core.get_fwhm()
        self.runs.loc[label, "E"] = self.E
        self.runs.loc[label, "N"] = self.N
        self.runs.loc[label, "sigma_r"] = self.sig_r
        self.runs.loc[label, "rel_decrement"] = rel_decrement
        self.runs.loc[label, "F2"] = F2
        self.runs.loc[label, "FWHM"] = fwhm
        self.runs.loc[label, "f_predict"] = f_predict
        if f_predict < self.field_width/2:
            self.msg("[WARN] Predicted focus inside field definition boundary!")
            self.runs.loc[label, "warn_focus_in_field"] = True
        self.msg("F2: %.3f mT^2"%(F2/mm**2))
        self.msg("FWHM: %.2f cm" % (fwhm/cm))
        self.msg("Predicted focal length for %.2f MeV: %.2f cm"%(self.E, f_predict/cm))


## Initial data collection:
    def overview_run(self, step=10, beam="gauss", beam_2d=True):
        """
            TBA
        """
        label = self.run_label
        E = self.runs.loc[label, "E"]
        N = self.runs.loc[label, "N"]
        sig_r = self.runs.loc[label, "sigma_r"]
        self.msg(">>> at %s: Performing overview tracking of %d electrons, E=%.2f"%(label, N, E))
        self.runs.loc[label, "field_width"] = self.field_width
        right_margin = np.max([self.runs.loc[label, "f_predict"], self.field_width/2+0.25])
        bounds = [0, self.field_width/2+0.25 + right_margin+0.5]
        self.setup_tracking(bounds, step=step)
        self.runs.loc[label, "z_solenoid"] = self.z_solenoid
        self.runs.loc[label, "step_overview"] = step
        self.generate_beam(E=E, N=N, sig_r=sig_r, distribution=beam, twodim=beam_2d)
        self.runs.loc[label, "beam_distr"] = beam
        self.runs.loc[label, "beam_dim"] = 2 if beam_2d else 1
        self.msg("Running ASTRA...")
        self.astra.run()

        self.msg("Processing ASTRA output...")
        s = self.astra.read_states()
        s, zpos, parts, pref = self.process_states(s)

        self.msg("Saving results...")
        data = {"s":s, "zpos":zpos, "parts":parts, "pref":pref}
        for key in data.keys():
            self.data[label][key] = data[key]

    def add_heads(self):
        label = self.run_label
        s = self.data[label]["s"]
        zpos = self.data[label]["zpos"]
        heads = self.make_heads(s, zpos)
        self.data[label]["heads"] = heads

## Focal region estimation
    def get_focal_region(self):
        # using edge case method.
        label = self.run_label
        self.msg(">>> at %s: estimating focal region."%label)
        s = self.data[label]["s"]
        z_solenoid = self.runs.loc[label, "z_solenoid"]
        paraxial = s.loc[0, "r"].idxmin()
        offaxis = s.loc[0, "r"].idxmax()
        p = s.swaplevel()
        step = self.runs.loc[label, "step_overview"]  # add step, b/c the indices represet beam states before/after the respective focal events
        right = p.loc[paraxial].query("z>@z_solenoid and pr < 0").get("z").idxmax() + step*mm
        left =  p.loc[offaxis].query("z>@z_solenoid and pr > 0").get("z").idxmin() - step*mm
        self.runs.loc[label, "z_focal_left"] = left
        self.runs.loc[label, "z_focal_right"] = right
        self.msg("Focal region at ~ z: [%.2f, %.2f] cm; z_solenoid %.2f m"%(left/cm, right/cm, self.z_solenoid))

#Focal run:
    def focal_run(self, step=0.1):
        """
            TBA
        """
        label = self.run_label
        self.msg(">>> at %s: Performing focal sweep..."%label)
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

#Fitting trajectories @ focus:
    def fit_focal_traj(self, model="offset"):  #TODO
        label = self.run_label

        self.msg(">>> at %s: Analyzing focal region data..."%label)
        self.msg("Trajectory model: %s"%model)
        self.runs.loc[label, "ray_model"] = model
        s_f = self.data[label]["s_f"]
        parts = self.data[label]["parts"]
        z_solenoid = self.runs.loc[label, "z_solenoid"]

        fits = pd.DataFrame()
        p_f = s_f.swaplevel()

        if model=="offset":
            for part in parts:
                r_min_idx = p_f.loc[part, "r"].idxmin()
                f_guess = p_f.loc[(part, r_min_idx), "z"]
                r_min_guess = p_f.loc[(part, r_min_idx), "r"]
                dxdz_guess = np.mean(p_f.loc[part, "x"].diff()/p_f.loc[part, "z"].diff())
                dydz_guess = np.mean(p_f.loc[part, "y"].diff()/p_f.loc[part, "z"].diff())
                cos0_guess = np.cos(p_f.loc[(part, r_min_idx), "phi"]*np.pi)

                p, pcov = opt.curve_fit(self.model.off_axis_trajectory,
                    xdata=p_f.query("z>@z_solenoid").loc[part, "z"].values,
                    ydata=p_f.query("z>@z_solenoid").loc[part, "r"].values,
                    p0=[f_guess, r_min_guess, dxdz_guess, dydz_guess, cos0_guess],
                    bounds=([z_solenoid, 0, -np.inf, -np.inf, -1], [np.inf, self.sig_r*mm, np.inf, np.inf, 1]))
                dp = np.sqrt(np.diag(pcov))
                fits.loc[part, "r0"] =  p_f.query("zpos==0").loc[part,"r"].values
                columns = ["z_f", "r_min", "dxdz", "dydz", "cos0", "dz_f", "dr_min", "ddxdz", "ddydz", "dcos0"]
                fits.loc[part, columns] = [*p, *dp]

        if model=="axial":
            for part in parts:
                r_min_idx = p_f.loc[part, "r"].idxmin()
                f_guess = p_f.loc[(part, r_min_idx), "z"]
                drdz_guess = np.mean(p_f.loc[part, "r"].diff()/p_f.loc[part, "z"].diff())

                p, pcov = opt.curve_fit(self.model.axial_trajectory,
                    xdata=p_f.query("z>@z_solenoid").loc[part, "z"].values,
                    ydata=p_f.query("z>@z_solenoid").loc[part, "r"].values,
                    p0=[f_guess, drdz_guess])
                    #bounds=([z_solenoid, 0, -np.inf, -np.inf, -1], [np.inf, self.sig_r*mm, np.inf, np.inf, 1]))
                dp = np.sqrt(np.diag(pcov))
                fits.loc[part, "r0"] =  p_f.query("zpos==0").loc[part,"r"].values
                columns = ["z_f", "drdz", "dz_f", "ddrdz"]
                fits.loc[part, columns] = [*p, *dp]

        f_max_observed = fits.get("z_f").max()-self.runs.loc[label, "z_solenoid"]
        self.runs.loc[label, "f_max_observed"] = f_max_observed
        f_min_observed = fits.get("z_f").min()-self.runs.loc[label, "z_solenoid"]
        self.runs.loc[label, "f_min_observed"] = f_min_observed
        self.data[label]["fits"] = fits
        self.msg("Minimum focal length value observed: %.2f cm"%(f_min_observed*100))
        self.msg("Maximum focal length value observed: %.2f cm"%(f_max_observed*100))

#Assessing aberration:
    def fit_cs_expansion(self, order=1, sigma="radius2", sigma_abs=False):
        label = self.run_label
        self.msg(">>> at %s: Fitting foci to cs expansion model."%label)
        fits = self.data[label]["fits"].sort_values("r0")
        z_solenoid = self.runs.loc[label, "z_solenoid"]
        f = fits.get("z_f").values - z_solenoid
        r = fits.get("r0").values
        if not sigma:
            dr = np.ones(len(r))
        elif sigma=="offset":
            dr = fits.get("r_min").values
            self.msg("Weighing against axis offset.")
        elif sigma=="radius":
            dr = fits.get("r0").values
            self.msg("Weighing against initial radial position.")
        elif sigma=="radius2":
            dr = fits.get("r0").values**2
            self.msg("Weighing against squared initial radial position.")
        else: raise ValueError("Incorrect weighing option: %s"%sigma)
        self.runs.loc[label, "f_expansion_sigma"] = sigma
        f_guess = self.runs.loc[label, "f_max_observed"]
        c2_guess = 1
        #self.msg("f guess: %.2f cm"%(f_guess/cm))
        #self.msg("c1 guess: %.2e m"%c1_guess)
        self.msg("Expansion order: %d"%order)
        p0_1 = np.array((f_guess, c2_guess))
        p0 = np.concatenate((p0_1, np.zeros(order-1)))
        bounds_lower = np.concatenate(((f_guess*0.9, 0), np.zeros(order-1)))
        bounds_upper = np.concatenate(((f_guess*1.1, np.inf), np.ones(order-1)*np.inf))
        c, ccov = opt.curve_fit(self.model.f_expansion,
        xdata=r,
        ydata=f,
        p0=p0,
        sigma=dr,
        absolute_sigma=sigma_abs,
        bounds=(bounds_lower, bounds_upper))

        dc = np.sqrt(np.diag(ccov))
        self.msg("Resulting focus: (%.3f +/- %.3f) cm"%(c[0]/cm, dc[0]/cm))
        cs = c[1]*c[0]**2
        quality = "Effective" if order==1 else "Minimal"
        self.msg("%s cs: %.3e m"%(quality, cs))
        self.runs.loc[label, "cs"] = cs
        self.runs.loc[label, "f"] = c[0]
        self.runs.loc[label, "df"] = dc[0]
        self.runs.loc[label, "f_expansion_order"] = order
        self.data[label]["exp_coeff"] = c[1:]
        self.data[label]["d_exp_coeff"] = dc[1:]
