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

# OLD FILE; TO BE REPLACED

import scipy.constants as const
import numpy as np
from numpy.lib.scimath import sqrt as csqrt
from numpy.lib.scimath import power as cpow
from scipy.integrate import quad as integral
from scipy.misc import derivative

mm = 10**(-3)
MeV = 10**6

def parse_geometry(geometry):
    r_in, a, b = np.array(geometry, dtype=np.complex128)*mm
    r = r_in + a/2
    Rsq = r + a**2/24/r
    c = csqrt((b**2 - a**2)/12)
    return [Rsq, c]

def impuls(E):  # [MeV] relativistic impulse
    return 1/const.c*np.sqrt((E*const.e*MeV)**2 - (const.m_e*const.c**2)**2)

### Models
def get_Bz_2loop(z, scaling, Rsq, c):
    pterm = (Rsq + c)**2/(z**2+(Rsq+c)**2)**(1.5)
    mterm = (Rsq - c)**2/(z**2+(Rsq-c)**2)**(1.5)
    B = 1/4*const.mu_0*scaling*(pterm + mterm)
    return np.real(B)  # complex part is 0 anyways

### Old Field integrals:
def F1(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz_2loop(z, scaling, Rsq, c)**1
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I, dI

def F2(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz_2loop(z, scaling, Rsq, c)**2
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I, dI

def F3(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        Bz = get_Bz_2loop(z, scaling, Rsq, c)
        ddBz = derivative(get_Bz_2loop, z, n=2, args=(scaling, Rsq, c))
        return -ddBz*Bz/2
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I, dI

def F4(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz_2loop(z, scaling, Rsq, c)**4
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I, dI

### Old Resulting values:
def focal(f2, p):
    inverted = (const.e/2/p)**2*f2
    return 1/inverted

def aberr_s(f3, f4, p, R):
    rad = R*mm
    return const.e**2*rad**4/4/p**2*f3 + const.e**4*rad**4/12/p**4*f4

def l_eff(scaling, geomp, decimal_places=3):
    z = np.linspace(-1,1,2*10**decimal_places+1)
    Bz = get_Bz_2loop(z, scaling, *geomp)
    ober = z[Bz>=np.max(Bz)/2]
    fwhm = (2*abs(ober[0]))
    return fwhm  # l_eff at fwhm +- 1/10^decimal_places

def peak_B(scaling, geomp):
    return get_Bz_2loop(0, scaling, *geomp)
