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

import pycode.wrapper as wrap
import numpy as np
from pycode.methods import impuls

mm = 10**(-3)
fm = 10**(-15)
cm = 10**(-2)

class API_iPython(wrap.PWrapper):

    def process_E_R(self):
        self.P = impuls(self.E)
        self.R = self.R_mm*mm

    def __init__(self, E, R):
        wrap.PWrapper.__init__(self, E, R)
        self.g = "None"
        self.s = "None"
        self.target_Bpeak = "None"
        self.target_l = "None"
        self.target_f = "None"
        self.target_g = "None"
        self.target_s = "None"
        self.process_E_R()

    def help(self):
        """
        Print some helpful info !WIP!
        """
        print("To view settings: handle.settings()")
        print("To view set targets: handle.targets()")
        print("To set targets, initial parameter X:")
        print("handle.target_x = new value / interval")
        print("handle.x = new value")
        print("Parameters: g (Rmean, a, b) [mm], s [Ampere-Turns]")
        print("Targets: Bpeak [mT], l [mm], f [cm]; \n g, s as above.")
        print("To disable constraints in a parameter, set to \"None\".")

    def settings(self):
        print("E: Electron energy [MeV]:", self.E)
        print("R: RMS Beam radius [mm]:", self.R_mm)
        print("g: Geomtery [mm]:", self.g)
        print("s: Ampere-turns [A*N]:", self.s)

    def targets(self):
        print("Target peak B [mT]:", self.target_Bpeak)
        print("Target FWHM [mm]:", self.target_l)
        print("Target f [cm]:", self.target_f)
        print("Target g [mm]:", self.target_g)
        print("Target s [N*A]:", self.target_s)

    def get_spot(self, f, cs):
        rspot = cs*(self.R/(f-self.R**2*cs/f**2))**3
        return rspot

    def describe(self, s, g):
        (B0, l, f, cs) = self.calc(s,g)
        spotsize = self.get_spot(f, cs)
        print("Peak axial field: %.3f mT"%(B0/mm))
        print("Effective field length: %.3f mm"%(l/mm))
        print("Focal distance for given E: %.3f cm"%(f/cm))
        print("Spherical aberration for given E: %.3E m"%cs)
        print("Focal spot radius: %.3E fm"%(spotsize/fm))

    def run_ctr(self, margin=5, maxiter=1000, ptol=8, gtol=8, verbose=2, penalty=0):
        wrap.PWrapper.run_ctr(self, margin=margin, verbose=verbose,
        ptol=ptol, gtol=gtol, penalty=penalty, maxiter=maxiter)

        print("Arrived at parameters:\n - s: %.3f"%self.s_opt)
        print(" - R_mean: %.3f mm, a: %.3f mm, b: %.3f mm"%(self.g_opt[0], self.g_opt[1], self.g_opt[2]))
        self.describe(s,g)
