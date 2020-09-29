import numpy as np
import matplotlib.pyplot as plt
from solensim.units import *
import scipy.constants as const


def linefit(x, y, dy):
  # Lineare Regression (y = ax + b) mit Gewichtung nach dem Skript von U. Müller, S.41.
  # straightforward and easy, but has some drawbacks. (See Ref. 2)
  # takes:
  # x, y, dy vectors, with dy = s(y);
  # returns:
  # A, B - fit parameters,
  # dA, dB - their uncertainties.
  # Direct calculation by formulas in Mueller's Script (Ref. 1, p. 41)

  D = np.sum((x**2)/dy**2)*np.sum(1/dy**2) - (np.sum(x/dy**2))**2;
  A = (np.sum(x*y/dy**2)*np.sum(1/dy**2) - np.sum(x/dy**2)*np.sum(y/dy**2))/D;
  B = (np.sum((x**2)/dy**2)*np.sum(y/dy**2) - np.sum(x*y/dy**2)*np.sum(x/dy**2))/D;
  dA = np.sqrt(np.sum(1/dy**2)/D);
  dB = np.sqrt(np.sum((x**2)/dy**2)/D);

  return A, dA, B, dB


def line_rot(astra, track):
    print("Loaded line beam")
    astra.beam_preset = "line"
    astra.verbose = False
    print("Running ASTRA...")
    astra.run()
    states = astra.read_states()
    zpos, s, refs = track.process_states(states)

    beam = s.loc[zpos[-1]]
    beam0 = s.loc[zpos[0]]
    rs = beam0["r"].values
    dphis = beam["dphi"].values
    turns = beam["turn"].values
    plt.figure(figsize=(8,5))
    plt.plot(rs*10**3, turns, ".k", label="%.2f MeV monochrome, N = 300, \nuniform horizontal line distribution"%(refs.get("pz").values[0]/10**6))
    #plt.axis([0, np.max(rs)*1.05, np.min(turns)-np.max(dphis), np.max(turns)+np.max(dphis)])
    plt.tick_params(axis="both",labelsize=12)
    plt.xlabel("Radial particle position at start [mm]", fontsize=20)
    plt.ylabel("Particle rotation [rad/pi]", fontsize=20)
    plt.legend(loc = "upper left", fontsize=16)
    plt.show()

    index_rmin = beam0["r"].idxmin()
    index_rmax = beam0["r"].idxmax()

    print("Rotation at %.3f: %.6f"%(rs[index_rmin], turns[index_rmin]))
    print("Rotation at %.3f: %.6f"%(rs[index_rmax], turns[index_rmax]))


def ring_rot(astra, track):
    print("Loaded ring beam")
    astra.beam_preset = "ring"
    astra.verbose = False
    print("Running ASTRA...")
    astra.run()
    states = astra.read_states()
    zpos, s, refs = track.process_states(states)
    refpz = refs.get("pz").values[0]
    beam = s.loc[zpos[-1]]
    beam0 = s.loc[zpos[0]]
    dphis = track.calc_dphi_v(beam["phi"], beam0["phi"])
    pzs = beam0.query("particle>0").get("pz").values + refpz
    print("Passed output to track handle")
    plt.figure(figsize=(8,5))
    plt.plot(pzs/10**6, dphis[1:], ".k", label="%.2f MeV, sigmaE %.2f MeV, N = 100, \n uniform 1 mm ring distribution"%(refpz/10**6, refpz/10**6/10))
    plt.axis([2.9, 4.1, 1.13,1.19])
    plt.tick_params(axis="both",labelsize=12)
    plt.xlabel("Particle energy at start [MeV]", fontsize=20)
    plt.ylabel("Particle rotation [rad/pi]", fontsize=20)
    plt.legend(loc = "upper right", fontsize=16)
    plt.show()

def bmax_line_sweep(astra, track, sweep=[0.025,1], step = 0.025):
    print("Loaded line beam")
    astra.track_preset = "default"
    astra.beam_preset = "line"
    print("Sweeping between %d and %d mT, step: %d mT"%(sweep[0]*10**3, sweep[1]*10**3, step*10**3))
    astra.verbose = False
    print("Activating field scaling...")
    astra.runfile["solenoid"]["s_noscale"] = False
    astra.runfile["solenoid"]["s_smooth"] = 100

    astra.update_runfile()
    print("Running...")
    deltas = []
    Bs = []
    for Bmax in np.arange(sweep[0], sweep[1]+step, step=step):
        Bs.append(Bmax)
        astra.runfile["solenoid"]["maxb"] = Bmax
        astra.update_runfile()
        astra.run()
        states = astra.read_last()
        zpos, s, refs = track.process_states(states)
        beam0 = s.loc[zpos[0]]
        beam = s.loc[zpos[-1]]
        dphis = track.calc_dphi_v(beam["phi"], beam0["phi"])
        dphi = np.mean(dphis)
        deltas.append(dphi)
    deltas = np.array(deltas)
    Bs = np.array(Bs)
    plt.figure(figsize=(8,5))
    plt.plot(Bs*10**3, deltas, "--k")
    plt.plot(Bs*10**3, deltas, ".k", label="%.2f MeV monochrome, N = 100, \nuniform horizontal line distribution"%(refs.get("pz").values[0]/10**6))
    plt.axis([50,500,0,2])
    plt.tick_params(axis="both",labelsize=12)
    plt.xlabel("Peak field on axis [mT]", fontsize=20)
    plt.ylabel("Particle rotation [rad/pi]", fontsize=20)
    plt.legend(loc = "center left", fontsize=15)
    plt.title("Field scaling sweep form %d mT to %d mT, step size %d mT"%(sweep[0]*10**3, sweep[1]*10**3, step*10**3), fontsize=15)
    plt.show()
    return Bs, deltas

def pphi_visual(astra, track):
    print("Loaded line beam")
    astra.beam_preset = "line"
    astra.verbose = False
    print("Running ASTRA...")
    astra.run()
    states = astra.read_states()
    zpos, s, refs = track.process_states(states)
    plt.figure(figsize=(8,5))
    for z in zpos:
        plt.loglog(s.loc[z].get("r").values, s.loc[z].get("pphi")/10**3, ".", label="z=%.1f"%z)
    #plt.axis([0,1.75,1,1.25])
    plt.tick_params(axis="both",labelsize=12)
    plt.xlabel("Radial particle position at start [m]", fontsize=20)
    plt.ylabel("Particle rotation impulse [keV]", fontsize=20)
    plt.legend(loc = "lower right", fontsize=12)
    plt.show()
