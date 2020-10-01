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
from numpy.lib.scimath import sqrt as csqrt
from numpy.lib.scimath import power as cpow
import scipy.integrate as integrate
from scipy.misc import derivative
import scipy.optimize as opt
import scipy.interpolate as interpolate
import numpy as np

from solensim.units import *
import solensim.backend.track as track

class Model():
    """
    TODO - fiedl calculation & beam model subblock
    """
    def __init__(self, linked_core):
        self.field = {
            "twoloop" : self.twoloop,
            "interpol" : self.data_interpol
        }
        self.linked_core = linked_core

#### Native field calculation
    def twoloop(self, z, p):
        s = p[0]
        g = p[1:]
        r_in, a, b = np.array(g, dtype=np.complex128)*mm
        r = r_in + a/2
        Rsq = r + a**2/24/r
        c = csqrt((b**2 - a**2)/12)
        pterm = (Rsq + c)**2/(z**2+(Rsq+c)**2)**(1.5)
        mterm = (Rsq - c)**2/(z**2+(Rsq-c)**2)**(1.5)
        B = 1/4*const.mu_0*s*(pterm + mterm)
        return np.real(B)  # complex part is 0 anyways

    def data_interpol(self, z, p):
        return self.linked_core.interpol_field(z)

#### Beam properties for characterization
    def impuls(self, E):  # [MeV] relativistic impulse
        return 1/const.c*np.sqrt((E*const.e*MeV)**2 - (const.m_e*const.c**2)**2)



class Core():
    """
    TODO - main class
    """
    def __init__(self):
        self.Model = Model(self)
        self.FM = "twoloop"  # default field model
        self.zmax = 1
        self.grain = 4  # 0.1 mm precision
        # Beam:
        self.E = "None"
        self.R = 1  # 1 mm beam "radius"
        self.sample_field(np.zeros(10), np.zeros(10))

# E, P relationship:
    def get_E(self):
        return self._E
    def set_E(self, E):
        self._E = E
        if str(E) != "None": self.P = self.Model.impuls(E)
        else: self.P = 0
    E = property(get_E, set_E)

    def sample_field(self, z, Bz):
        """
        enter z, Bz to create an interpolator.
        Use core.FM = "interpol" to calculate field based on samples
        (the p parameters are then irrelevant)
        """
        self.interpol_field = interpolate.interp1d(z, Bz, fill_value="extrapolate")

    def fint(self, p, n):
        """
        Compute nth field integral
        """
        if n == 3:
            integrand = lambda z: -1/2*self.Model.field[self.FM](z, p)*derivative(self.Model.field[self.FM], z, n=2, args=[p])
        else:
            integrand = lambda z: self.Model.field[self.FM](z, p)**n
        I, dI = integrate.quad(integrand, -np.inf, np.inf)
        return I

    def get_Bz(self, p):
        z = np.linspace(-self.zmax, self.zmax, num=2*10**self.grain+1)
        return self.Model.field[self.FM](z, p)

    def get_fwhm(self, p):
        """
        Get z(FWHM), assuming a field symmetrical around 0, small enough to contain FWHM within 1 meter
        """
        Bhalb = self.get_Bmax(p)/2
        f = lambda z: self.Model.field[self.FM](z,p) - Bhalb
        return opt.root_scalar(f, bracket=[0,self.zmax], xtol=10**(-self.grain)).root*2

    def get_f(self, p):
        f2 = self.fint(p, 2)
        return 1/((const.e/2/self.P)**2*f2)

    def get_Bmax(self, p):
        z = np.linspace(-self.zmax, self.zmax, num=2*10**self.grain+1)
        return np.max(self.Model.field[self.FM](z,p))

# Aberrations and the like:
# Spherical aberrations:
    def get_cs(self, p):
        f3 = self.fint(p, 3)
        f4 = self.fint(p, 4)
        rad = self.R*mm
        return const.e**2*rad**4/4/self.P**2*f3 + const.e**4*rad**4/12/self.P**4*f4
