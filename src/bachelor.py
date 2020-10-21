import numpy as np
import matplotlib.pyplot as plt
from solensim.aux import *
import scipy.constants as const

def line_rot(track, astra):
    astra.beam_preset = "line"
    astra.verbose = False
    astra.run()
    states = astra.read_states()
    s, zpos, parts, pref = track.process_states(states)

    beam = s.loc[2]
    beam0 = s.loc[0]
    rs = beam0["r"].values
    dphis = beam["dphi"].values
    turns = beam["turn"].values

    plt.figure(figsize=(8,5))
    plt.plot(rs*10**3, turns, ".k", label="%.2f MeV monochrome, N = 300, \nuniform horizontal line distribution"%(pref.get("pz").values[0]/10**6))
    plt.tick_params(axis="both",labelsize=12)
    plt.xlabel("Radial particle position at start [mm]", fontsize=20)
    plt.ylabel("Particle rotation [rad/pi]", fontsize=20)
    plt.legend(loc = "lower left", fontsize=16)
    plt.show()

    index_rmin = beam0["r"].idxmin()
    index_rmax = beam0["r"].idxmax()

    print("Rotation at %.3f: %.6f"%(rs[index_rmin], turns[index_rmin]))
    print("Rotation at %.3f: %.6f"%(rs[index_rmax], turns[index_rmax]))
