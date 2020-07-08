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
vstring = "0.2.0-RC1"
def_B = [75,125]
def_f = [50,50]
def_l = [45,55]
def_E = 3.5
def_Rbeam = 1
#################

import pycode.wrapper as wrap
import pycode.frontend as front
import numpy as np
import matplotlib.pyplot as plt

print("solensim v%s WIP Solenoid design tool"%vstring)
print("========================================")
e = front.API_iPython()
e.E = def_E
e.R = def_Rbeam
e.target_Bpeak = def_B
e.target_l = def_l
e.target_f = def_f
print("iPython interface handle initialized as \"e\".")
print("\nDefault settings:")
e.settings()
print("\nDefault targets:")
e.targets()
print("\nPlaceholder help text:")
e.help()
print("\nUse \"e.help()\" to view the startup help text again.")

import pycode.demos as demos
print("""
to run a few demos, try:
demos.opt_REGAE(e) - try to make the REGAE magnet even better (without considering the yoke);
    keyword arguments:
        margin - % bound on the deviation from source parameters, defaults to 5%
        maxiter - maximum iterations, defaults to 1000
        verbose - True/False - show optimizer's output to follow per-iteration progress, defaults to False
demos.calc_REGAE(e) - describe the REGAE magnet, from T. Gehrke's thesis

No "best" result to show yet, as we have problems with the general concept.
""")
