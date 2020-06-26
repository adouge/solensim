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

    B = e.get_B(s_REGAE, g_REGAE)*1000
    z = np.linspace(-1,1,len(B))
    plt.figure()
    plt.plot(z, B)
    plt.axis([-0.2,0.2, 0, 100])
    plt.xlabel("Position on axis [m]")
    plt.ylabel("B_z [mT]")
    plt.show()

def opt_REGAE(e):
    g_REGAE = [79.75, 41.8, 99.5]
    s_REGAE = 9*1000
    e.g = g_REGAE
    e.s = s_REGAE
    e.target_g = (np.array(e.g)*0.5, np.array(e.g)*1.5)
    e.target_s = [7*1000, 11*1000]
    e.target_l = 102
    e.target_f = 27.7
    e.target_Bpeak = 80
    e.run_ctr(maxiter=1000, margin=10, verbose=0)

def opt_best(e):
     e.g = [30,40,50]
     e.s = 4000
     e.target_g = (np.array(e.g)*0.5, np.array(e.g)*1.5)
     e.target_s = (2*1000,6*1000)
     e.target_f = [50,50]
     e.target_l = [25,75]
     e.target_Bpeak = [100,100]
     e.run_ctr(maxiter=1000, verbose=0)
     B = e.get_B(e.s_opt, e.g_opt)*1000
     z = np.linspace(-1,1,len(B))
     plt.figure()
     plt.plot(z, B)
     plt.axis([-0.2,0.2, 0, max(B)])
     plt.xlabel("Position on axis [m]")
     plt.ylabel("B_z [mT]")
     plt.show()

def opt_free(e):
     e.g = [30,40,50]
     e.s = 4000
     e.target_g = ([10,10,5], [200,200,200])
     e.target_s = (1*1000,16*1000)
     e.target_f = [50,50]
     e.target_l = "None"
     e.target_Bpeak = "None"
     e.run_ctr(maxiter=1000, verbose=0)
     B = e.get_B(e.s_opt, e.g_opt)*1000
     z = np.linspace(-1,1,len(B))
     plt.figure()
     plt.plot(z, B)
     plt.axis([-0.2,0.2, 0, max(B)])
     plt.xlabel("Position on axis [m]")
     plt.ylabel("B_z [mT]")
     plt.show()