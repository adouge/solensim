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

from solensim.aux import *
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
    def impuls_SI(self, E):  # SI relativistic impulse from MeV total energy
        return 1/const.c*np.sqrt((E*const.e*MeV)**2 - (const.m_e*const.c**2)**2)

class Core():
    """
    TODO - main class
    """
    def __init__(self):
        self.model = Model(self)
        self.FM = "twoloop"  # default field model
        self.bcalc_zmax = 1
        self.bcalc_zgrain = 4  # 0.1 mm precision

## Field model switching:
    def get_FM(self):
        return self._FM
    def set_FM(self, FM):
        if FM not in self.model.field.keys():
            raise(ValueError("Unknown field model: %s"%FM))
        else:
            self._FM = FM
            self.msg("Setting field model to \"%s\"."%FM)
    FM = property(get_FM, set_FM)

## Field data handling:
    def sample_field(self, z, Bz, extrapolate=False):
        """
        enter z, Bz to create an interpolator.
        Use core.FM = "interpol" to calculate field based on samples
        (the p parameters are then irrelevant)
        """
        if extrapolate: fill = "extrapolate"
        else: fill = 0
        self.interpol_field = interpolate.interp1d(z, Bz, fill_value=fill, bounds_error=False)
        self.msg("Sampled field for interpolation.")
        self.FM = "interpol"

    def get_scale_factor(self, f_N, E_N):  # [m], [MeV]
        """
            Assuming the to-be-scaled field is already sampled
        """
        p_N = self.model.impuls_SI(E_N)
        F2 = self.fint(2)
        k = 2*p_N/const.e/np.sqrt(f_N*F2)
        self.msg("Factor to scale sampled field to %.2f m @ %.2f MeV: k=%.6f"%(f_N, E_N, k))
        return k

## Field characterization:
    def fint(self, n, p="numeric"):
        """
        Compute nth field integral
        """
        if n == 3:
            integrand = lambda z: -1/2*self.model.field[self.FM](z, p)*derivative(self.model.field[self.FM], z, n=2, args=[p])
        else:
            integrand = lambda z: self.model.field[self.FM](z, p)**n
        I, dI = integrate.quad(integrand, -np.inf, np.inf)
        return I

    def get_z(self):
        return np.linspace(-self.bcalc_zmax, self.bcalc_zmax, num=2*10**self.bcalc_zgrain+1)

    def get_Bz(self, p="numeric"):
        z = self.get_z()
        return self.model.field[self.FM](z, p)

    def get_fwhm(self, p="numeric"):
        """
        Get z(FWHM), assuming a field symmetrical around 0, small enough to contain FWHM within 1 meter
        """
        Bhalb = self.get_Bmax(p)/2
        f = lambda z: self.model.field[self.FM](z,p) - Bhalb
        return opt.root_scalar(f, bracket=[0,self.bcalc_zmax], xtol=10**(-self.bcalc_zgrain)).root*2

    def get_f(self, E, p="numeric"):
        f2 = self.fint(2, p)
        P = self.model.impuls_SI(E)
        return 1/((const.e/2/P)**2*f2)

    def get_Bmax(self, p):
        z = np.linspace(-self.bcalc_zmax, self.bcalc_zmax, num=2*10**self.bcalc_zgrain+1)
        return np.max(self.model.field[self.FM](z,p))
