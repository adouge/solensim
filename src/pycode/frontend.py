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

    def __init__(self):
        wrap.PWrapper.__init__(self)
        self.M = "twoloop"
        self.E = "None"
        self.R = "None"
        self.minRin = "None"
        self.g = "None"
        self.s = "None"
        self.target_Bpeak = "None"
        self.target_l = "None"
        self.target_f = "None"
        self.target_g = "None"
        self.target_s = "None"
        self.margin = 10
        self.results = []

### Help text
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
            Beam settings:
                E [MeV] - electron energy (currently only monochrome beams handled)
                R [mm] - "beam radius"
                minRin [mm] - Lower bound on solenoid inner radius
            Starting optimization settings:
                g (Inner radius, a, b) [mm]
                s [Ampere-Turns]
            Target values - single values or tuples (lower bound, upper bound):
                Bpeak [mT] - peak field on axis
                l [mm] - FWHM
                f [cm] - focal length for given E

            handle.margin [%] - target margin for single-value settings
            To disable constraints in a parameter, set to \"None\".

            Main optimization routine:
            handle.run_ctr(margin=5, maxiter=1000, ptol=9, gtol=9, verbose=2)
                margin: maximum tolerable percent deviation from target values, for non-interval settings
                maxiter: maximum iteration number
                ptol, gtol: convergence tolerance (10 to the negative power of),
                    should be set according to desired result precision and scale of the end solution (WIP)
        """
        print(helptext)

### Current settings/targets readout
    def settings(self):
        print("E: Electron energy [MeV]:", self.E)
        print("R: Beam radius [mm]:", self.R)
        print("minRin: Minimal inner radius [mm]:", self.minRin)
        print("g: Starting geomtery [mm]:", self.g)
        print("s: Starting scaling factor [A*N]:", self.s)

    def targets(self):
        print("Target peak B [mT]:", self.target_Bpeak)
        print("Target FWHM [mm]:", self.target_l)
        print("Target f [cm]:", self.target_f)
        print("Target g [mm]:", self.target_g)
        print("Target s [N*A]:", self.target_s)
        print("Target margin for non-interval bounds [%]:", self.margin)

### descriptive methods
    def calc(self, scaling, geometry):
        f2 = self.FN(scaling, geometry, 2)
        f3 = self.F3(scaling, geometry)
        f4 = self.FN(scaling, geometry, 4)

        f = self.get_f(scaling, geometry)
        cs = self.get_cs(scaling, geometry)
        l = self.get_l(scaling,geometry)
        B0 = self.get_Bpeak(scaling, geometry)

        result = (B0, l, f, cs)
        return result

    def get_spot(self, f, cs):
        R = self.R*mm
        return cs*(R/(f-R**2*cs/f**2))**3

    def describe(self, s, g):
        """
        Calculate characteristic values from given parameters s, g
        """
        (B0, l, f, cs) = self.calc(s,g)
        spotsize = self.get_spot(f,cs)
        print("Parameters:\n - s: %.3f"%s)
        print(" - R: %.3f mm, a: %.3f mm, b: %.3f mm"%(g[0], g[1], g[2]))
        print("Resulting characteristics:")
        print(" - Peak axial field: %.3f mT"%(B0/mm))
        print(" - Effective field length: %.3f mm"%(l/mm))
        print(" - Focal distance for given E: %.3f cm"%(f/cm))
        print(" - Spherical aberration for given E, R: %.3E m"%cs)
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



### main routine
    def run_ctr(self, maxiter=100, ptol=9, gtol=9, verbose=2, penalty=0):
        """
        CTR Optimization routine wrapper.
        Set parameters as attributes of the handle;
        Keyword arguments:
            maxiter - maximum number of algorithm iterations (default 1000)
            margin - percent margin for non-interval target values (default 5%)
        """
        wrap.PWrapper.run_ctr(self, verbose=verbose,
        ptol=ptol, gtol=gtol, penalty=penalty, maxiter=maxiter)
        self.append_result()
        self.result(-1)
