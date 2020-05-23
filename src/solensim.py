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

import numpy as np
import wrapper


def matlab_test(engine):
    print("Testing Matlab interface:\nRunning Kamp's script:")
    engine.run("doNewCoil.m", nargout=0)
    print("done. \n Testing mwraptest.m:")
    A = engine.magic(5)
    B = engine.magic(5,5)
    C = engine.mwraptest(A,B)
    print(A,B)
    print("product:")
    print(C)

def wrapper_test():
    o = wrapper.OWrapper()
    print("Testing Oct2Py Wrapper:\nRunning Kamps's script:")
    o.run("doNewCoil")
    input("Press Enter to continue...")
    o.restart()
    print("done. \n Testing wraptest.m:")
    A = np.random.rand(5,5)
    B = np.random.rand(5,5)
    C = o.wraptest(A,B)
    print(A,B)
    print("product:")
    print(C)

def main():
    print("solensim v0.0.1 Solenoid simulation tool")
    print("========================================")
    print("Nothing to see here yet! \n Here's a wrapper test:")
    wrapper_test()
    print("----------------------------------------")
    m = wrapper.mWrapper()
    matlab_test(m)

if __name__ == '__main__':
    main()
