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

# # # # # parameters:
# geometry = [r, a, b] # mean readius, radial width, axial width, Windungsdichte [mm]
# E = energy [MeV]
# R = beam radius

from pycode.methods import *
import scipy.constants as const
import scipy.optimize as opt
import numpy as np

mm = 10**(-3)

class Core():
    def update_settings(self):
        self.P = impuls(self.E)
        self.R = self.R_mm*mm

    def __init__(self, E, R):
        self.E = E
        self.R_mm = R
        self.update_settings()

    # descriptive method:
    def calc(self, scaling, geometry):
        geomp = parse_geometry(geometry)

        f1,df1 = F1(scaling, geomp)
        f2,df2 = F2(scaling, geomp)
        f3,df3 = F3(scaling, geomp)
        f4,df4 = F4(scaling, geomp)

        f = focal(f2, self.P)
        cs = aberr(f3, f4, self.P, self.R)
        l = l_eff(scaling, geomp)
        B0 = peak_B(scaling, geomp)

        result = (B0, l, f, cs)
        return result

    # constraint evaluation:

    def get_Bpeak(self, params):
        scaling, r, a, b = params
        geomp = parse_geometry((r,a,b))
        return peak_B(scaling, geomp)

    def get_f(self, params):
        scaling, r, a, b = params
        geomp = parse_geometry((r,a,b))
        f2,df2 = F2(scaling, geomp)
        return focal(f2, self.P)

    def get_l(self, params):
        scaling, r, a, b = params
        geomp = parse_geometry((r,a,b))
        return l_eff(scaling, geomp)

    def get_cs(self, params):
        scaling, r, a, b = params
        geomp = parse_geometry((r,a,b))
        f3,df3 = F3(scaling, geomp)
        f4,df4 = F4(scaling, geomp)
        return aberr(f3, f4, self.P, self.R)

    # optimization methods

    def opt_cs(self, params):  # function for use in minimization
        NI, r, a, b = params
        geomp = parse_geometry((r,a,b))
        f3,df3 = F3(NI, geomp)
        f4,df4 = F4(NI, geomp)
        return aberr(f3, f4, self.P, self.R)

    #####
        # constrained trust region algorithm:

    def define_ctr_constraints(self, margin=5, t_Bpeak="None", t_l="None", t_f="None", t_p="None"):
        """
        Define constraints. Defaults to unconstrained.
        B, l: [lower, upper] or target (margin of X% (def. 5%) assumed)
        f: target - lower fixed, upper unbounded; or [lower, upper]
        Geometry: defaults to unconstrained;
            supply (lb, ub) as [Rmean, a, b] in m;
            or target list +- margin
        """

        t_margin = margin/100
        # target constraints:
        constraints = []
        # peak B:
        if t_Bpeak != "None":
            if type(t_Bpeak) in [np.float64, float]:
                con_Bpeak = opt.NonlinearConstraint(self.get_Bpeak, t_Bpeak*(1-t_margin), t_Bpeak*(1+t_margin))
            elif type(t_Bpeak) in [np.array, list]:
                con_Bpeak = opt.NonlinearConstraint(self.get_Bpeak, t_Bpeak[0], t_Bpeak[1])
            else: raise ValueError("Incorrect maxB constraint provided.")
            constraints.append(con_Bpeak)
        # FWHM:
        if t_l != "None":
            if type(t_l) in [np.float64, float]:
                con_l = opt.NonlinearConstraint(self.get_l, t_l*(1-t_margin), t_l*(1+t_margin))
            elif type(t_l) in [np.array, list]:
                con_l = opt.NonlinearConstraint(self.get_l, t_l[0], t_l[1])
            else: raise ValueError("Incorrect FWHM constraint provided.")
            constraints.append(con_l)
        # focal length:
        if t_f != "None":
            if type(t_f) in [np.float64, float]:
                con_f = opt.NonlinearConstraint(self.get_f, t_f, np.inf)
            elif type(t_f) in [np.array, list]:
                con_f = opt.NonlinearConstraint(self.get_f, t_f[0], t_f[1])
            else: raise ValueError("Incorrect f constraint provided.")
            constraints.append(con_f)

        # parameter bounds:
        A = np.array([[1,0,0,0],[0,1,-1/2,0],[0,0,1,0],[0,0,0,1]])  # lin abb to verify p
        lower_bound = np.array([0,self.R*5,0,0])
        if t_p != "None":
            if len(t_p)==4:  # anticipating a target parameter setting
                if type(t_p) != np.array: t_p = np.array(t_p)
                if not (A.dot(t_p)>lower_bound).all(): raise ValueError("Negative a,b, or inner radius below 5x beam radius provided")
                con_p = opt.LinearConstraint(np.identity(4),t_p*(1-t_margin),t_p*(1+t_margin))
            elif len(t_p)==2:
                if not (A.dot(t_p[0])>lower_bound).all(): raise ValueError("Negative a,b, or inner radius below 5x beam radius provided")
                if not (A.dot(t_p[1])>lower_bound).all(): raise ValueError("Negative a,b, or inner radius below 5x beam radius provided")
                con_p = opt.LinearConstraint(np.identity(4),t_p[0],t_p[1])
            else: raise ValueError("Improper geometry bounds provided.")
            constraints.append(con_p)
        else:
            con_p = opt.LinearConstraint(A, lower_bound, np.inf)
            constraints.appennd(con_p)

        return constraints

    def ctr_minimize(self, start_p, constraints, max_iter=1000, ptol=6, verbose=2):
        opt_out = opt.minimize(self.opt_cs, start_p,
            constraints=constraints,
            options={"maxiter":max_iter, "verbose":verbose, "xtol":np.power(10.,-ptol)},
            method="trust-constr")
        return opt_out
