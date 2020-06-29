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

# Demonstrations

import numpy as np
import matplotlib.pyplot as plt

def calc_REGAE(e):
    g_REGAE = [79.75, 41.8, 99.5]
    s_REGAE = 9*1000
    e.describe(s_REGAE, g_REGAE)
    e.illustrate(s_REGAE, g_REGAE)

def opt_REGAE(e, margin=10, maxiter=1000, verbose=False):
    g_REGAE = [79.75, 41.8, 99.5]
    s_REGAE = 9*1000
    e.g = g_REGAE
    e.s = s_REGAE
    e.target_g = (np.array(e.g)*(1-margin/100), np.array(e.g)*(1+margin/100))
    e.target_s = [7*1000, 11*1000]
    e.target_l = 102
    e.target_f = 27.7
    e.target_Bpeak = 80
    if verbose==True: v = 2
    else: v = 1
    e.settings()
    e.targets()
    e.run_ctr(maxiter=maxiter, margin=margin, verbose=v)

#def opt_best(e):
#     e.g = [30,40,50]
#     e.s = 6000
#     e.target_g = (np.array(e.g)*0.1, np.array(e.g)*2)
#     e.target_s = (2*1000,10*1000)
#     e.target_f = [50,np.inf]
#     e.target_l = [45,55]
#     e.target_Bpeak = [90,110]
#     e.run_ctr(maxiter=1000, verbose=1)
