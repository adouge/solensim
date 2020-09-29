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
from solensim.units import *

def field_REGAE(handle, astra):
    print("Describing REGAE solenoid field (no yoke)\n via two-loop-approximation. (demos.field_REGAE(core, astra))")
    g_REGAE = [30, 99.5, 41.8]
    s_REGAE = 9*1000
    print("Rin, a, b [cm]:")
    print(g_REGAE)
    print("Scaling factor [A]: %d"%s_REGAE)

    p = (s_REGAE, *g_REGAE)
    z = np.linspace(-handle.zmax, handle.zmax, num = 2*10**handle.grain+1)
    B = handle.get_Bz(p)
    Bmax = handle.get_Bmax(p)
    fwhm = handle.get_fwhm(p)
    handle.E = 3.47
    focal = handle.get_f(p)
    print("Maximum field strength: %.3f mT, FWHM %.1f mm"%(Bmax/mm, fwhm/mm))
    print("Focal length for %.2f MeV energy: %.3f m"%(handle.E, focal))
    plt.figure(figsize=(7,4))
    plt.plot(z/cm, B/mm, "-k", label="Bz(z)")
    plt.xlabel("Axial position [cm]")
    plt.ylabel("On-axis field strength [mT]")
    plt.axis([-handle.zmax/cm, handle.zmax/cm, 0, Bmax*1.05/mm])
    plt.show()
    astra.write_field(z, B)
    print("Field saved to solenoid.dat (don't forget scaling in ASTRA runfile!)")
    #print("\n Running generator...")
    #astra.generate()
    #print("\n Running ASTRA...")
    #astra.run()
