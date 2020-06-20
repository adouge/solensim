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

import pycode.backend
import numpy as np

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

    def run_ctr(self, margin=5, maxiter=1000, ptol=6, gtol=6, verbose=2,penalty=0):
        constraints = self.define_ctr_constraints(margin=margin)
        out = self.ctr_minimize(constraints,
            max_iter=maxiter,
            ptol=ptol, gtol=gtol,
            penalty=penalty,
            verbose=verbose)
        self.s_opt = out.x[0]
        self.g_opt = out.x[1:]
        self.last_message = out.message
        self.last_out = out
