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
    g_REGAE = [30, 99.5, 41.8]
    s_REGAE = 9*1000
    e.describe(s_REGAE, g_REGAE)
    e.illustrate(s_REGAE, g_REGAE)

def opt_REGAE(e, maxiter=100, verbose=False):
    g_REGAE = [30, 99.5, 41.8]
    s_REGAE = 9*1000
    e.g = g_REGAE
    e.s = s_REGAE
    e.target_g = e.g
    e.target_s = e.s
    e.target_l = 148
    e.target_f = 124.931
    e.target_Bpeak = 60.855
    if verbose==True: v = 2
    else: v = 1
    e.run_ctr(maxiter=maxiter, verbose=v)

def opt_ok(e):
    e.s = 8000
    e.g = [50,50,50]
    e.target_s = [4000,12000]
    e.target_g = [[10,10,10],[200,200,200]]
    e.target_f = [50, np.inf]

    e.target_Bpeak = [100,100]
    e.target_l = [45,55]
    e.run_ctr(maxiter=100, verbose=1)


def opt_ok2(e):
    e.s = 5000
    e.g = [100,100,100]
    e.target_s = [4000,12000]
    e.target_g = [[10,10,10],[200,200,200]]
    e.target_f = [50, np.inf]

    e.target_Bpeak = [100,100]
    e.target_l = [45,55]
    e.run_ctr(maxiter=100, verbose=1)
