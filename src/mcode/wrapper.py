#########################################################################
#    Copyright 2020 Anton Douginets, Andrii Yanovets
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

import oct2py
import os
import matlab.engine

class OWrapper(oct2py.Oct2Py):
    """
    wrapper class for Oct2Py
    Octave instances start with mcode in PATH
    """
    mcode_path = work_dir = os.path.dirname(os.path.realpath(__file__))
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
    mcode_path = work_dir = os.path.dirname(os.path.realpath(__file__))
    engine.addpath(mcode_path)
    return engine
