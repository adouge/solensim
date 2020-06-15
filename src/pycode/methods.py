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
# geometry = [r, a, b] # mean readius, radial width, axial width, Windungsdichte
# E = energy
# R = beam radius

import scipy.constants as const
import numpy as np
from scipy.integrate import quad as integral
from scipy.misc import derivative


def parse_geometry(geometry):
    [r, a, b] = geometry
    Rsq = r*(1 + a**2/(24*r**2))
    c = np.sqrt((a**2 - b**2)/12)
    return [Rsq, c]

def impuls(E):
    return np.sqrt(2*const.m_e*(E*const.e - const.m_e*const.c**2))

# # # # #

# axial field:
def get_Bz(z, NI, Rsq, c):
    pterm = (Rsq + c)**2/(z**2+(Rsq+c)**2)**(3/2)
    mterm = (Rsq - c)**2/(z**2+(Rsq-c)**2)**(3/2)
    B = const.mu_0*NI*(pterm + mterm)
    return B

def unscaled_Bz(z, Rsq, c):
    pterm = (Rsq + c)**2/(z**2+(Rsq+c)**2)**(3/2)
    mterm = (Rsq - c)**2/(z**2+(Rsq-c)**2)**(3/2)
    B = const.mu_0*(pterm + mterm)
    return B

# Field integrals:
def F1(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz(z, scaling, Rsq, c)**1
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I

def F2(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz(z, scaling, Rsq, c)**2
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I

def F3(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        Bz = get_Bz(z, scaling, Rsq, c)
        ddBz = derivative(get_Bz, z, n=2, args=(scaling, Rsq, c))
        return -ddBz*Bz/2
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I

def F4(scaling, geomp):
    def integrand(z, scaling, Rsq, c):
        return get_Bz(z, scaling, Rsq, c)**4
    I, dI = integral(integrand, -np.inf, np.inf, args=(scaling, *geomp))
    return I


# Resulting values:
def focal(f2, p):
    inverted = (const.e/2/p)**2*f2
    return 1/inverted

def aberr(f3, f4, p, rad):
    return const.e**2*rad**4/4/p**2*f3 + const.e**4*rad**4/12/p**4*f4

def l_eff(scaling, geomp, decimal_places=3):
    z = np.linspace(-1,1,2*10**decimal_places+1)
    Bz = get_Bz(z, scaling, *geomp)
    ober = z[Bz>=np.max(Bz)/2]
    fwhm = np.round((2*abs(ober[0]) + 1/10**decimal_places), decimal_places)
    return fwhm  # l_eff at fwhm +- 1/10^decimal_places

def peak_B(scaling, geomp):
    return np.max(get_Bz(0, scaling, *geomp))
