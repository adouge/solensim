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
import matplotlib.pyplot as plt
import pandas as pd

import solensim.wrapper as wrapper
import solensim.frontend as frontend
from solensim.units import *
import scipy.constants as const

# config
wrapper.load_ini()
print("(tried loading config)\n")
#################
vstring = "0.4.0"
#################

print("solensim v%s solenoid design & optimization tool"%vstring)
print("========================================")

print("Astra interface initialized as \"astra\"")
astra = frontend.Astra_Interface()

core = frontend.Core()
print("Core handle initialized as \"core\" (WIP).")

track = frontend.Tracker(astra)
track.linked_core = core  # bind core and track together
print("Tracker initialized as \"track\" (WIP).")

### BA section

from importlib import reload
import demos
import bachelor
