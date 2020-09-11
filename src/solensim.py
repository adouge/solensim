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

import solensim.wrapper as wrapper
import solensim.frontend as frontend

# config
wrapper.load_ini()
print("(tried loading config)\n")
#################
vstring = "0.3.2"
#################

print("solensim v%s solenoid design & optimization tool"%vstring)
print("========================================")

print("Astra interface initialized as \"astra\"")
astra = frontend.Astra_Interface()
print("astra.help():")
astra.help()
astra.presets()
print("\n")

track = frontend.Tracker(astra)
print("Tracker initialized as \"track\".")
core = frontend.Core()
print("API initialized as \"core\".")

import demos
demos.field_REGAE(core, astra)
astra.help()
astra.presets()
