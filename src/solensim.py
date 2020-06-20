#########################################################################
#    Copyright 2020 Anton Douginets, Andrii Yanovets
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

# main program script

# config
#################
vstring = "0.1.0"
def_B = [75,125]
def_f = [50,50]
def_l = [45,55]
def_E = 3.5
def_Rbeam = 1
#################

import wrapper as wrap
import frontend as front
import numpy as np
import matplotlib.pyplot as plt


def tests():
    import mcode.wrapper as mwrap

    def test_matlab():
        print("Testing Matlab interface:")
        M = mwrap.mWrapper()
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
        O = mwrap.OWrapper()
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

    print("Here are some wrapper tests:")
    test_octave()
    print("----------------------------------------")
    test_matlab()

def main():
    print("solensim v%s Solenoid design tool"%vstring)
    print("========================================")
    e = front.API_iPython(def_E, def_Rbeam)
    e.target_Bpeak = def_B
    e.target_l = def_l
    e.target_f = def_f
    print("iPython interface handle initialized as \"e\".")

if __name__ == '__main__':
    main()
