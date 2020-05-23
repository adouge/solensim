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
import frontend as fe

def main():
    [vstring] = fe.load_conf()
    print("solensim %s Solenoid simulation tool"%vstring)
    print("========================================")
    print("Nothing to see here yet! \n Here are some wrapper tests:")
    fe.test_octave()
    print("----------------------------------------")
    fe.test_matlab()
    print("========================================")

    fe.demo("Some Parameters")
if __name__ == '__main__':
    main()
