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

import oct2py
import matlab.engine
import os

class Wrapper(oct2py.Oct2Py):
    def __init__(self):
        oct2py.Oct2Py.__init__(self)
        self.addpath("mcode")

    def restart(self):
        oct2py.Oct2Py.restart(self)
        self.addpath("mcode")

class MWrapper():
    mcode_path = os.path.dirname(os.path.realpath(__file__)) + "/mcode"
    def __init__(self):
        self.e = matlab.engine.start_matlab()
        self.e.cd(self.mcode_path, nargout=0)
