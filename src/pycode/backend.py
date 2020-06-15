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

from pycode.methods import *
import scipy.constants as const
import numpy as np

class Core():
    def __init__(self):
        pass

    def scalc(self):
        geometry = self.p["g"]
        scaling = self.p["s"]
        E = self.p["E"]
        R = self.p["R"]

        geomp = parse_geometry(geometry)
        p = impuls(E)

        f1 = F1(scaling, geomp)
        f2 = F2(scaling, geomp)
        f3 = F3(scaling, geomp)
        f4 = F4(scaling, geomp)

        f = focal(f2, p)
        cs = aberr(f3, f4, p, R)
        l = l_eff(scaling, geomp)
        B0 = peak_B(scaling, geomp)

        returns = (B0, l, f, cs)
        self.result = returns
        return returns

    def calc(self, geometry, scaling, E, R):
        geomp = parse_geometry(geometry)
        p = impuls(E)

        f1 = F1(scaling, geomp)
        f2 = F2(scaling, geomp)
        f3 = F3(scaling, geomp)
        f4 = F4(scaling, geomp)

        f = focal(f2, p)
        cs = aberr(f3, f4, p, R)
        l = l_eff(scaling, geomp)
        B0 = peak_B(scaling, geomp)

        print("Peak axial field:", B0*1000, "mT")
        print("Effective field length:", l*1000,"mm")
        print("Focal distance for given E:", f*100,"cm")
        print("Spherical aberration for given E:", cs)


        return (B0, l, f, cs)
