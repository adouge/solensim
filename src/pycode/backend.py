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

from pycode.methods import *
import scipy.constants as const
import scipy.optimize as opt
import numpy as np

mm = 10**(-3)
cm = 10**(-2)

class Core():
    def __init__(self, E, R):
        self.E = E
        self.R_mm = R

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

    def get_B(self, s, g, grain=3):
        geomp = parse_geometry(g)
        z = np.linspace(-1,1,num=2*10**grain+1)
        return get_Bz(z, s, *geomp)

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

    def get_l(self, params, grain=3):
        scaling, r, a, b = params
        geomp = parse_geometry((r,a,b))
        return l_eff(scaling, geomp, decimal_places=grain)

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

    def define_ctr_constraints(self, margin=5):
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
        if str(self.target_Bpeak) != "None":
            t_Bpeak = np.array(self.target_Bpeak)*mm
            if type(t_Bpeak) in [np.float64, float]:
                con_Bpeak = opt.NonlinearConstraint(self.get_Bpeak, t_Bpeak*(1-t_margin), t_Bpeak*(1+t_margin))
            elif type(t_Bpeak) in [np.ndarray, list]:
                con_Bpeak = opt.NonlinearConstraint(self.get_Bpeak, t_Bpeak[0], t_Bpeak[1])
            else: raise ValueError("Incorrect maxB constraint provided.")
            constraints.append(con_Bpeak)
        # FWHM:
        if str(self.target_l) != "None":
            t_l = np.array(self.target_l)*mm
            if type(t_l) in [np.float64, float]:
                con_l = opt.NonlinearConstraint(self.get_l, t_l*(1-t_margin), t_l*(1+t_margin))
            elif type(t_l) in [np.ndarray, list]:
                con_l = opt.NonlinearConstraint(self.get_l, t_l[0], t_l[1])
            else: raise ValueError("Incorrect FWHM constraint provided.")
            constraints.append(con_l)
        # focal length:
        if str(self.target_f) != "None":
            t_f = np.array(self.target_f)*cm
            if type(t_f) in [np.float64, float]:
                con_f = opt.NonlinearConstraint(self.get_f, t_f, np.inf)
            elif type(t_f) in [np.ndarray, list]:
                con_f = opt.NonlinearConstraint(self.get_f, t_f[0], t_f[1])
            else: raise ValueError("Incorrect f constraint provided.")
            constraints.append(con_f)

        # geometry, scaling bounds:
        A = np.array([[1,0,0,0],[0,1,-1/2,0],[0,0,1,0],[0,0,0,1]])  # lin abb to verify p, general case
        Ag = np.array([[0,0,0,0],[0,1,-1/2,0],[0,0,1,0],[0,0,0,1]])  # verifying geometry
        As = np.array([[1,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]])  # making sure scaling is positive
        lower_bound = np.array([0,self.R*5,0,0])

        if self.target_s != "None":
            t_s = np.array(self.target_s)
            if (t_s.shape == ()):
                lb = np.array(((1-t_margin)*t_s,0,0,0))
                ub = np.array(((1+t_margin)*t_s,0,0,0))
                con_s = opt.LinearConstraint(As, lb, ub)
            else:
                con_s = opt.LinearConstraint(As, np.array((t_s[0],0,0,0)), np.array((t_s[1],0,0,0)))
            constraints.append(con_s)

        if str(self.target_g) != "None":
            t_g = np.array(self.target_g)
            if len(t_g) == 3:
                con_g = opt.LinearConstraint(Ag, np.array((0,*t_g*(1-t_margin))),np.array((0,*t_g*(1+t_margin))))
            elif len(t_g) == 2:
                con_g = opt.LinearConstraint(Ag, np.array((0,*t_g[0])),np.array((0,*t_g[1])))
            else: raise ValueError("Improper geometry bounds provided.")
            constraints.append(con_g)

        if (str(self.target_s) == "None") and (str(self.target_g) == "None"):
            con_validity = opt.LinearConstraint(A, lower_bound, np.inf)
            constraints.append(con_validity)

        return constraints

    def ctr_minimize(self, constraints, max_iter=1000, ptol=6, gtol=6, verbose=2, penalty=0):
        opt_out = opt.minimize(self.opt_cs, (self.s, *self.g),
            constraints=constraints,
            options={"maxiter":max_iter,
                "verbose":verbose,
                "xtol": np.power(10.,-ptol),
                "gtol": np.power(10.,-gtol),
                "initial_constr_penalty": np.power(10.,penalty)},
            method="trust-constr")
        return opt_out
