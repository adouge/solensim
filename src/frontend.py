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

def demo(params):
    """
    the future main() to be presented.
    """
    print("This should be the main script for presenting the project, with adjustable parameters: %s"%params)
