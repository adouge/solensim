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

# Wrapper code segment

import oct2py
import os
import matlab.engine
import pycode.backend
import numpy as np

def workdir():
    work_dir = os.path.dirname(os.path.realpath(__file__))
    return work_dir

class OWrapper(oct2py.Oct2Py):
    """
    wrapper class for Oct2Py
    Octave instances start with mcode in PATH
    """
    work_dir = workdir()
    mcode_path = work_dir  + "/mcode"
    def __init__(self):
        oct2py.Oct2Py.__init__(self)
        self.addpath(self.mcode_path)

    def restart(self):
        oct2py.Oct2Py.restart(self)
        self.addpath(self.mcode_path)

def mWrapper():
    """
    Returns a MATLAB engine instance, with mcode in PATH
    """
    engine = matlab.engine.start_matlab()
    work_dir = workdir()
    mcode_path = work_dir  + "/mcode"
    engine.addpath(mcode_path)
    return engine

def stop(Wrapper):
    """
    Stop the engine, delete handle
    """
    Wrapper.exit()
    del(Wrapper)

class PWrapper(pycode.backend.Core):
    """
    User-facing methods of the python backend
    """
    def __init__(self, E, R):
        pycode.backend.Core.__init__(self, E, R)
        self.g = "Not set",
        self.s = "Not set"
        self.sg = self.g
        self.ss = self.s

    def show_settings(self):
        print("g: Geomtery [mm]:", self.g)
        print("s: Ampere-turns [A*N]:", self.s)
        print("E: Electron energy [MeV]:", self.E)
        print("R: RMS Beam radius [mm]:", self.R_mm)

    def show_target(self):
        print("Target peak B [mT]:",self.target_Bpeak*1000)
        print("Target FWHM [mm]:", self.target_l*1000)
        print("Target f [cm]:",self.target_f*100)

    def exit(self):
        pass  # let the wrapper's del(self) handle it

    def describe(self, result):
        (B0, l, f, cs) = result
        print("Peak axial field:", B0*1000, "mT")
        print("Effective field length:", l*1000,"mm")
        print("Focal distance for given E:", f*100,"cm")
        print("Spherical aberration for given E:", cs)

    def scalc(self):  # placeholder output
        geometry = self.g
        scaling = self.s
        result = self.calc(scaling, geometry)
        return result

    def run_ctr(self, margin=5, maxiter=1000, ptol=6, verbose=2):
        constraints = self.define_ctr_constraints(t_Bpeak=self.target_Bpeak,
            t_l=self.target_l,
            t_f=self.target_f,
            t_p = [self.s,*self.g],
            margin=margin)
        out = self.ctr_minimize(self.s, self.g, constraints, max_iter=maxiter, ptol=ptol, verbose=verbose)
        self.s_opt = out.x[0]
        self.g_opt = out.x[1:]
        result = self.calc(self.s_opt, self.g_opt)
        self.describe(result)
        return out
