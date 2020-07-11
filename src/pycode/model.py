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
import numpy as np
from numpy.lib.scimath import sqrt as csqrt
from numpy.lib.scimath import power as cpow

mm = 10**(-3)
MeV = 10**6

### Field models
def twoloop(z, s, g):
    r_in, a, b = np.array(g, dtype=np.complex128)*mm
    r = r_in + a/2
    Rsq = r + a**2/24/r
    c = csqrt((b**2 - a**2)/12)
    pterm = (Rsq + c)**2/(z**2+(Rsq+c)**2)**(1.5)
    mterm = (Rsq - c)**2/(z**2+(Rsq-c)**2)**(1.5)
    B = 1/4*const.mu_0*s*(pterm + mterm)
    return np.real(B)  # complex part is 0 anyways


### Beam properties
def impuls(E):  # [MeV] relativistic impulse
    return 1/const.c*np.sqrt((E*const.e*MeV)**2 - (const.m_e*const.c**2)**2)
