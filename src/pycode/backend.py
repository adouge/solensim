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
        self.target_l = 0.05
        self.target_f = 0.5
        self.target_Bpeak = 0.1
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

    def verify_geometry(self, params):
        scaling, r, a, b = params
        inner_radius = r-a/2
        return r, a, b, inner_radius

    # optimization methods

    def opt_cs(self, params):  # function for use in minimization
        NI, r, a, b = params
        geomp = parse_geometry((r,a,b))
        f3,df3 = F3(NI, geomp)
        f4,df4 = F4(NI, geomp)
        return aberr(f3, f4, self.P, self.R)

    #####
        # constrained trust region algorithm:

    def define_ctr_constraints(self, target_Bpeak, target_l, target_f, target_margin=0.0499, bounded_f=False, geom_lb=[0,0,0], geom_ub=[1000, 1000,1000]):
        # target constraints:
        l_con = opt.NonlinearConstraint(self.get_l, target_l*(1-target_margin), target_l*(1+target_margin))
        Bpeak_con = opt.NonlinearConstraint(self.get_Bpeak, target_Bpeak*(1-target_margin), target_Bpeak*(1+target_margin))
        if bounded_f: f_con = opt.NonlinearConstraint(self.get_f, target_f, target_f*(1+target_margin))
        else: f_con = opt.NonlinearConstraint(self.get_f, target_f, np.inf)
        # geometric constraints:
        ir = lambda r,a : r-a/2
        glb = geom_lb.copy()
        gub = geom_ub.copy()
        glb.append(ir(*glb[:2]))
        gub.append(ir(*gub[:2]))
        glb[-1] += self.R*5  # inserting 5 sigmas of beam radius
        geometry_con = opt.NonlinearConstraint(self.verify_geometry, glb, gub)
        constraints = [Bpeak_con, f_con, l_con, geometry_con]
        return constraints

    def ctr_minimize(self, start_scaling, start_geometry, constraints, max_iter=1000, ptol=6, verbose=2):
        start_params = (start_scaling, *start_geometry)
        opt_out = opt.minimize(self.opt_cs, start_params,
            constraints=constraints,
            options={"maxiter":max_iter, "verbose":verbose, "xtol":np.power(10.,-ptol)},
            method="trust-constr")
        return opt_out
