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
    def __init__(self, E, R):
        wrap.PWrapper.__init__(self, E, R)
        self.g = "Not set",
        self.s = "Not set"
        self.sg = self.g
        self.ss = self.s
        self.target_Bpeak = 100*mm
        self.target_l = 50*mm
        self.target_f = 500*mm

    def help(self):
        """
        Print some helpful info !WIP!
        """
        pass

    def show_settings(self):
        print("g: Geomtery [mm]:", self.g)
        print("s: Ampere-turns [A*N]:", self.s)
        print("E: Electron energy [MeV]:", self.E)
        print("R: RMS Beam radius [mm]:", self.R_mm)

    def show_targets(self):
        print("Target peak B [mT]:",self.target_Bpeak*1000)
        print("Target FWHM [mm]:", self.target_l*1000)
        print("Target f [cm]:",self.target_f*100)

    def describe(self, result):
        (B0, l, f, cs) = result
        print("Peak axial field:", B0*1000, "mT")
        print("Effective field length:", l*1000,"mm")
        print("Focal distance for given E:", f*100,"cm")
        print("Spherical aberration for given E:", cs)

    def run_ctr(self, margin=5, maxiter=1000, ptol=6, verbose=2,
        target_Bpeak = "None",
        target_l = "None",
        target_f = "None",
        target_p = "None"
        ):
        if type(target_p) != np.array: target_p = np.array(target_p)
        wrap.PWrapper.run_ctr(self,
            margin=margin, maxiter=maxiter, ptol=ptol, verbose=verbose,
            # units conversion: mT, mm, cm, mm
            target_Bpeak = target_Bpeak*mm,
            target_l = target_l*mm,
            target_f = target_f*mm*10,
            target_p = target_p*mm
            )
        result = self.calc(self.p_opt[0], self.p_opt[1:])
        self.describe(result)
