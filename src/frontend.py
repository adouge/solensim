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
        self.results = []

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
        print("Parameters:\n - s: %.3f"%s)
        print(" - R_mean: %.3f mm, a: %.3f mm, b: %.3f mm"%(g[0], g[1], g[2]))
        print("Resulting characteristics:")
        print(" - Peak axial field: %.3f mT"%(B0/mm))
        print(" - Effective field length: %.3f mm"%(l/mm))
        print(" - Focal distance for given E: %.3f cm"%(f/cm))
        print(" - Spherical aberration for given E: %.3E m"%cs)
        print(" - Focal spot radius: %.3E fm"%(spotsize/fm))

    def illustrate(self, s, g):
        B = e.get_B(s, g)*1000
        z = np.linspace(-1,1,len(B))
        plt.figure()
        plt.plot(z, B)
        plt.axis([-0.2,0.2, 0, max(B)])
        plt.xlabel("Position on axis [m]")
        plt.ylabel("B_z [mT]")
        plt.show()

    def result(self, n):
        result = self.results[n]
        print("Settings:")
        print(" - g [mm]:", result["g"])
        print(" - s [A*N]:", result["s"])
        print("Targets:")
        print(" - peak B [mT]:", result["t_B"])
        print(" - FWHM [mm]:", result["t_l"])
        print(" - f [cm]:", result["t_f"])
        print(" - g [mm]:", result["t_g"])
        print(" - s [N*A]:", result["t_s"])
        print("Result:")
        self.describe(result["sopt"],result["gopt"])

    def append_result(self):
        result = {
        "g" : self.g,
        "s" : self.s,
        "t_g" : self.target_g,
        "t_s" : self.target_s,
        "t_B" : self.target_Bpeak,
        "t_l" : self.target_l,
        "t_f" : self.target_f,
        "gopt": self.g_opt,
        "sopt": self.s_opt
        }
        self.results.append(result)

    def run_ctr(self, margin=5, maxiter=1000, ptol=8, gtol=8, verbose=2, penalty=0):
        wrap.PWrapper.run_ctr(self, margin=margin, verbose=verbose,
        ptol=ptol, gtol=gtol, penalty=penalty, maxiter=maxiter)
        self.append_result()
        self.result(-1)
        self.illustrate(self.s_opt, self.g_opt)
