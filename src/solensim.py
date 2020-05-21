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
from wrappers import Wrapper

def wrapper_test():
    o = Wrapper()
    print("Testing Kamps's script:")
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
    print("solensim v0.0.0 Solenoid simulation tool")
    print("========================================")
    print("Nothing to see here yet! \n Here's a wrapper test:")
    wrapper_test()


if __name__ == '__main__':
    main()
