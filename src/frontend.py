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

# frontend code segment

import wrapper as wrap
import numpy as np
from pycode.methods import impuls

mm = 10**(-3)

def load_conf():
    """
    Rudimentary configuration, possibly replace with regexp
    """
    work_dir = wrap.workdir()
    config = open(wrap.workdir()+"/solensim.cfg", "r")
    config.readline()
    vstring = config.readline().split("=")[1].strip()
    config.close()
    return [vstring]

def test_matlab():
    print("Testing Matlab interface:")
    M = wrap.mWrapper()
    print("\nRunning Kamp's script:")
    M.run("doNewCoil.m", nargout=0)
    print("done. \n Testing mwraptest.m:")
    A = M.magic(5)
    B = M.magic(5)
    C = M.mwraptest(A,B)
    print(A,B)
    print("product:")
    print(C)
    input("Press Enter to continue...")
    wrap.stop(M)

def test_octave():
    O = wrap.OWrapper()
    print("Testing Oct2Py Wrapper:\nRunning Kamps's script:")
    O.run("doNewCoil")
    input("Press Enter to continue...")
    O.close()
    print("done. \n Testing wraptest.m:")
    A = np.random.rand(5,5)
    B = np.random.rand(5,5)
    C = O.wraptest(A,B)
    print(A,B)
    print("product:")
    print(C)
    wrap.stop(O)

class API_iPython(wrap.PWrapper):

    def process_E_R(self):
        self.P = impuls(self.E)
        self.R = self.R_mm*mm

    def __init__(self, E, R):
        wrap.PWrapper.__init__(self, E, R)
        self.g = "None"
        self.s = "None"
        self.target_Bpeak = 0.1
        self.target_l = 0.05
        self.target_f = 0.5
        self.target_g = "None"
        self.target_s = "None"
        self.process_E_R()

    def help(self):
        """
        Print some helpful info !WIP!
        """
        print("To view settings: handle.settings()")
        print("To view set targets: handle.targets()")
        print("To set targets, initial parameter X:"")
        print("handle.target_x = new value / interval")
        print("handle.x = new value")
        print("Parameters: g (Rmean, a, b) [mm], s [Ampere-Turns]")
        print("Targets: Bpeak [mT], l [mm], f [mm], g, s as above.")
        print("To disable constraints in a parameter, set to \"None\".")

    def settings(self):
        print("E: Electron energy [MeV]:", self.E)
        print("R: RMS Beam radius [mm]:", self.R_mm)
        print("g: Geomtery [mm]:", self.g)
        print("s: Ampere-turns [A*N]:", self.s)

    def targets(self):
        print("Target peak B [mT]:",self.target_Bpeak)
        print("Target FWHM [mm]:", self.target_l)
        print("Target f [mm]:",self.target_f)

    def describe(self, result):
        (B0, l, f, cs) = result
        print("Peak axial field:", B0/mm, "mT")
        print("Effective field length:", l/mm,"mm")
        print("Focal distance for given E:", f/mm,"mm")
        print("Spherical aberration for given E:", cs)

    def run_ctr(self, margin=5, maxiter=1000, ptol=6, verbose=2):
        wrap.PWrapper.run_ctr(self, margin=margin, maxiter=maxiter, ptol=ptol, verbose=verbose)
        result = self.calc(self.s_opt, self.g_opt)
        self.describe(result)
