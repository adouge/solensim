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

import sscode.wrapper as wrapper
import sscode.frontend as front
import numpy as np
import matplotlib.pyplot as plt

# config
#################
vstring = "0.2.2"
def_B = [75,125]
def_f = [50,np.inf]
def_l = [45,55]
def_E = 3.5
def_Rbeam = 3
def_minRin = def_Rbeam*3
#################

print("solensim v%s Solenoid design & optimization tool"%vstring)
print("========================================")
e = front.legacy_opt_API()
e.E = def_E
e.R = def_Rbeam
e.minRin = def_minRin
e.target_Bpeak = def_B
e.target_l = def_l
e.target_f = def_f
print("Legacy interface handle initialized as \"e\".")
print("\nDefault settings:")
e.settings()
print("\nDefault targets:")
e.targets()
print("\nUse \"e.help()\" to view the help text.")

import demos
print("""
to run a few demos, try:
    demos.calc_REGAE(e) - describe the REGAE magnet, from T. Gehrke's thesis
    demos.opt_ok(e) - an OK approach to characteristics originally demanded
    demos.opt_ok2(e) - a second ok approach from different starting values
""")
