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
import matplotlib.pyplot as plt

mm = 10**(-3)
fm = 10**(-15)
cm = 10**(-2)

class API_iPython(wrap.PWrapper):

    def __init__(self, E, R):
        wrap.PWrapper.__init__(self, E, R)
        self.g = "None"
        self.s = "None"
        self.target_Bpeak = "None"
        self.target_l = "None"
        self.target_f = "None"
        self.target_g = "None"
        self.target_s = "None"
        self.results = []

    def help(self):
        """
        Display placeholder help text.
        """
        helptext = """
            To view settings:
                handle.settings()
            To view set targets:
                handle.targets()
            To set targets, initial parameter x:
                handle.target_x = new value / interval
                handle.x = new value
            Permanent settings:
                E [MeV] - electron energy (currently only monochrome beams handled)
                R_mm [mm] - "beam radius"
              call handle.process_E_R() to apply changes - currently a non-developed feature;
              simpler would be to create a new handle, using "handle_name = front.API_iPython([E in MeV], [R im mm]).
            Starting optimization settings:
                g (Rmean, a, b) [mm]
                s [Ampere-Turns]
            Target values - single values or tuples (lower bound, upper bound):
                Bpeak [mT] - peak field on axis
                l [mm] - FWHM
                f [cm] - focal length for given E
            To disable constraints in a parameter, set to \"None\".

            Main optimization routine:
            handle.run_ctr(margin=5, maxiter=1000, ptol=8, gtol=8, verbose=2)
                margin: maximum tolerable percent deviation from target values, for non-interval settings
                maxiter: maximum iteration number
                ptol, gtol: convergence tolerance (10 to negative power of), not really tested yet
        """
        print(helptext)

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

    def describe(self, s, g):
        """
        Calculate characteristic values from given parameters s, g
        """
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
        """
        Draw a Bz(x) plot from given s, g
        """
        B = self.get_B(s, g)*1000
        z = np.linspace(-1,1,len(B))
        plt.figure()
        plt.plot(z, B)
        plt.axis([-0.2,0.2, 0, max(B)])
        plt.xlabel("Position on axis [m]")
        plt.ylabel("B_z [mT]")
        plt.show()

    def result(self, n):
        """
        Show n-th (0, 1, .... -1) result
        """
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
        """
        CTR Optimization routine wrapper.
        Set parameters as attributes of the handle;
        Keyword arguments:
            maxiter - maximum number of algorithm iterations (default 1000)
            margin - percent margin for non-interval target values (default 5%)
        """
        wrap.PWrapper.run_ctr(self, margin=margin, verbose=verbose,
        ptol=ptol, gtol=gtol, penalty=penalty, maxiter=maxiter)
        self.append_result()
        self.result(-1)
        self.illustrate(self.s_opt, self.g_opt)
