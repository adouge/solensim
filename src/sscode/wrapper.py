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

import sscode.backend.calc as calc
import sscode.backend.optim as optim
import sscode.backend.track as track
import sscode.handle as handle
from sscode.units import *

import numpy as np

def test_load_mcode_plugin():
    import plugins.mcode.wrapper as mwrapper
    o = mwrapper.OWrapper()
    return o

class Wrapper(handle.Handle):
    def __init__(self):
        handle.Handle.__init__(self)

class Old_Wrapper(handle.Legacy):
    """
    Backend I/O wrap methods
    """

    config_path = "WIP"

    def __init__(self):
        handle.Legacy.__init__(self)
        self.results = []

    def restart(self):
        """
        A restart method
        """
        print("Clearing result storage...")
        self.results = []  # flush result storage
        print("Reloading default settings...")
        self.load_default_settings()  # WIP
        print("Done.")

    def load_default_settings(self):
        """
            Load default settings from config file @ self.config_path (WIP)
        """
        wip()
        pass

### result storage
    def append_result(self):
        result = {
        "g" : self.g,
        "s" : self.s,
        "t_g" : self.target_g,
        "t_s" : self.target_s,
        "t_B" : self.target_Bpeak,
        "t_l" : self.target_l,
        "t_f" : self.target_f,
        "t_margin" : self.margin,
        "gopt": self.g_opt,
        "sopt": self.s_opt
        }
        self.results.append(result)

### main routine interlayer
    def run_ctr(self, maxiter=100, ptol=9, gtol=9, verbose=2,penalty=0):
        constraints = self.define_ctr_constraints()
        out = self.ctr_minimize(constraints,
            max_iter=maxiter,
            ptol=ptol, gtol=gtol,
            penalty=penalty,
            verbose=verbose)
        self.s_opt = out.x[0]
        self.g_opt = out.x[1:]
        self.last_message = out.message
        self.last_out = out
