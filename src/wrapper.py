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

mm = 10**(-3)

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

    def exit(self):
        pass  # let the wrapper's del(self) handle it

    def scalc(self):  # placeholder output
        geometry = self.g
        scaling = self.s
        result = self.calc(scaling, geometry)
        return result

    def run_ctr(self, margin=5, maxiter=1000, ptol=6, verbose=2,
        target_Bpeak = "None",
        target_l = "None",
        target_f = "None",
        target_p = "None"
        ):
        # IS THIS NECESSARY? unconstrained should be allowed
        if target_Bpeak == "None": t_Bpeak = self.target_Bpeak
        else: t_Bpeak = target_Bpeak
        if target_l == "None": t_l = self.target_l
        else: t_l = target_l
        if target_f == "None": t_f = self.target_f
        else: t_f = target_f
        if target_p == "None": t_p = [self.s,*self.g]
        else: t_p = target_p

        constraints = self.define_ctr_constraints(
            t_Bpeak = t_Bpeak,
            t_l = t_l,
            t_f = t_f,
            t_p = t_p,
            margin=margin)
        out = self.ctr_minimize((self.s, *self.g), constraints, max_iter=maxiter, ptol=ptol, verbose=verbose)
        self.s_opt = out.x[0]
        self.g_opt = out.x[1:]
        self.last_message = out.message
        return (self.s_opt, self.g_opt)
