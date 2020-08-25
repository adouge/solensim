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
from scipy.integrate import quad as integral
from scipy.misc import derivative
import numpy as np
import plugins.model as model


mm = 10**(-3)
cm = 10**(-2)

class Core():
    """
    TODO
    """
    def __init__(self):
        pass

    field = {
    "twoloop" : model.twoloop
    }

# E, P relationship:
    def get_E(self):
        return self._E
    def set_E(self, E):
        self._E = E
        if str(E) != "None": self.P = model.impuls(E)
        else: self.P = 0
    E = property(get_E, set_E)

# Field description:
    def FN(self, s, g, n):
        if n == 3:
            integrand = lambda z: -1/2*self.field[self.M](z, s, g)*derivative(self.field[self.M], z, n=2, args=(s,g))
        else:
            integrand = lambda z: self.field[self.M](z, s, g)**n
        I, dI = integral(integrand, -np.inf, np.inf)
        return I

    def get_B(self, s, g, grain=3):
        z = np.linspace(-1,1,num=2*10**grain+1)
        return self.field[self.M](z, s, g)

    def get_l(self, s, g):
        """
        Assuming a symmetrical field with max at 0
        """
        tol = 4  # (0.1 milimeter precision)
        Bhalb = self.field[self.M](0,s,g)/2
        f = lambda x: self.field[self.M](x,s,g) - Bhalb
        return opt.root_scalar(f, bracket=[0,1], xtol=10**(-tol)).root*2

    def get_f(self, s, g):
        f2 = self.FN(s, g, 2)
        return 1/((const.e/2/self.P)**2*f2)

    def get_Bpeak(self, s, g):
        return self.field[self.M](0,s,g)

# Aberrations and the like:
    def get_cs(self, s, g):  # current opt function
        f3 = self.FN(s, g, 3)
        f4 = self.FN(s, g, 4)
        rad = self.R*mm
        return const.e**2*rad**4/4/self.P**2*f3 + const.e**4*rad**4/12/self.P**4*f4

### OPT section:

    def opt(self, p):
        return self.get_cs(p[0], p[1:])

    def char(self, p):  # characteristic vector for constraints, order Bmax FWHM Focal
        Bpeak = self.get_Bpeak(p[0], p[1:])
        l = self.get_l(p[0], p[1:])
        f = self.get_f(p[0], p[1:])
        return (Bpeak, l, f)


    def define_ctr_constraints(self):
        """
        Define constraints. Defaults to unconstrained.
        B, l: [lower, upper] or target (margin of X% (def. 5%) assumed)
        f: target - lower fixed, upper unbounded; or [lower, upper]
        Geometry: defaults to unconstrained;
            supply (lb, ub) as [Rmean, a, b] in m;
            or target list +- margin
        """

        t_margin = self.margin/100
        constraints = []

        # characteristics: Bpeak, fwhm, f
        lower_bound = [0,0,0]
        upper_bound = [np.inf,np.inf,np.inf]

        constrained_B = str(self.target_Bpeak) != "None"
        constrained_FWHM = str(self.target_l) != "None"
        constrained_f = str(self.target_f) != "None"

        if constrained_B or constrained_FWHM or constrained_f:
            # Bpeak
            if constrained_B:
                t_Bpeak = np.array(self.target_Bpeak)*mm
                if type(t_Bpeak) in [np.float64, float]:
                    lower_bound[0] = t_Bpeak*(1-t_margin)
                    upper_bound[0] = t_Bpeak*(1+t_margin)
                elif type(t_Bpeak) in [np.ndarray, list]:
                    lower_bound[0] = t_Bpeak[0]
                    upper_bound[0] = t_Bpeak[1]
                else: raise ValueError("Incorrect Bpeak constraint provided.")

            if constrained_FWHM:
                t_l = np.array(self.target_l)*mm
                if type(t_l) in [np.float64, float]:
                    lower_bound[1] = t_l*(1-t_margin)
                    upper_bound[1] = t_l*(1+t_margin)
                elif type(t_l) in [np.ndarray, list]:
                    lower_bound[1] = t_l[0]
                    upper_bound[1] = t_l[1]
                else: raise ValueError("Incorrect FWHM constraint provided.")

            if constrained_f:
                t_f = np.array(self.target_f)*cm
                if type(t_f) in [np.float64, float]:
                    lower_bound[2] = t_f*(1-t_margin)
                    upper_bound[2] = t_f*(1+t_margin)
                elif type(t_f) in [np.ndarray, list]:
                    lower_bound[2] = t_f[0]
                    upper_bound[2] = t_f[1]
                else: raise ValueError("Incorrect f constraint provided.")

            con_char = opt.NonlinearConstraint(self.char, lower_bound, upper_bound)
            constraints.append(con_char)

        # parameter bounds:
        lower_bound = [0,self.minRin,0,0]
        upper_bound = [np.inf, np.inf, np.inf, np.inf]

        if str(self.target_s) != "None":
            t_s = np.array(self.target_s)
            if (t_s.shape == ()):
                lower_bound[0] = (1-t_margin)*t_s
                upper_bound[0] = (1+t_margin)*t_s
            else:
                lower_bound[0] = t_s[0]
                upper_bound[0] = t_s[1]

        if str(self.target_g) != "None":
            t_g = np.array(self.target_g)
            if len(t_g) == 3:
                lower_bound = [lower_bound[0], *t_g*(1-t_margin)]
                upper_bound = [upper_bound[0], *t_g*(1+t_margin)]
            elif len(t_g) == 2:
                lower_bound = [lower_bound[0], *t_g[0]]
                upper_bound = [upper_bound[0], *t_g[1]]
            else: raise ValueError("Error handling geometry bounds.")

        con_p = opt.LinearConstraint(np.identity(4), lower_bound, upper_bound)
        constraints.append(con_p)

        return constraints

    def ctr_minimize(self, constraints, max_iter=100, ptol=9, gtol=9, verbose=2, penalty=0):
        opt_out = opt.minimize(self.opt, (self.s, *self.g),
            constraints=constraints,
            options={"maxiter":max_iter,
                "verbose":verbose,
                "xtol": np.power(10.,-ptol),
                "gtol": np.power(10.,-gtol),
                "initial_constr_penalty": np.power(10.,penalty)},
            method="trust-constr")
        return opt_out
