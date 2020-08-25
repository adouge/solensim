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

import wrapper
import numpy as np
import matplotlib.pyplot as plt

mm = 10**(-3)
fm = 10**(-15)
cm = 10**(-2)

class API_iPython(wrapper.Wrapper):

    def __init__(self):
        wrapper.Wrapper.__init__(self)
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

            handle.describe(s, g):
                characterize a given parameter set
            handle.illustrate(s,g):
                draw a plot

            handle.result(n):
                show result #n; follows standard python indexing:
                    - n -1: previous result, -2 the one before etc.
                    - n 0: first result, 1 the second etc.
                A WIP feature, will be improved.

            handle.get_B(s,g,model):
                return axial field from -1 to 1 m, using specified field model (default - "twoloop")
            handle.FN(s, g, n):
                get the nth field integral of field with parameters s, g

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
        f3 = self.FN(scaling, geometry, 3)
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
        print(" - Margin [%]:", result["t_margin"])
        print("Result:")
        self.describe(result["sopt"],result["gopt"])
        self.illustrate(result["sopt"],result["gopt"])

### main routine
    def run_ctr(self, maxiter=100, ptol=9, gtol=9, verbose=2, penalty=0):
        """
        CTR Optimization routine wrapper.
        Set parameters as attributes of the handle;
        Keyword arguments:
            maxiter - maximum number of algorithm iterations (default 1000)
            margin - percent margin for non-interval target values (default 5%)
        """
        wrapper.Wrapper.run_ctr(self, verbose=verbose,
        ptol=ptol, gtol=gtol, penalty=penalty, maxiter=maxiter)
        self.append_result()
        self.result(-1)
