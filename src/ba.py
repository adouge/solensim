"""BA scripts."""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from solensim.aux import *
import scipy.constants as const
from scipy.optimize import curve_fit
import pysnooper

def second_moment(a,b):
    n = len(a)
    return np.sum(a*b)/n - np.sum(a)*np.sum(b)/n**2

#@pysnooper.snoop()
def make_emits(track, core, label):
    m_e = const.m_e*const.c**2/const.e
    s = track.data[label]["s"].copy()
    zpos = track.data[label]["zpos"]
    emittances = pd.DataFrame(columns=["z", "eps_x", "eps_y", "eps_z", "eps_z_rel"])
    indices = ["x", "y", "z", "z_rel"]
    for z in zpos:
        b = s.loc[z].copy()
        #zref = track.data[label]["pref"].loc[(z, 0), "z"]
        #zref = np.mean(b["z"])
        #pzref = np.mean(b["pz"])
        for index in indices:
            i = b[index].values
            pi = b["p"+index].values
            m1 = second_moment(i, i)
            m2 = second_moment(pi, pi)
            m3 = second_moment(i, pi)
            eps2 = m1*m2 - m3**2
            if eps2 < 0:
                eps = 0
            else:
                eps = np.sqrt(eps2)#/m_e
            emittances.loc[z, "eps_"+index] = eps
            emittances.loc[z, index+"_ms"] = m1#np.sqrt(m1)
            emittances.loc[z, "p"+index+"_ms"] = m2#np.sqrt(m2)
            emittances.loc[z, "cov2_"+index] = m3**2
        emittances.loc[z, "z"] = z
        track.data[label]["eps"] = emittances

def plot_z_emits(track, core, labels):
    plt.figure(figsize=(9,9))
    color = {
        0: "r",
        1: "g",
        2: "b",
        3: "k"
    }
    center = 0
    shift = 0
    counter = 0
    for label in labels:
        eps = track.data[label]["eps"]
        maxe = (eps["eps_z_rel"].max())
        z_sol = track.runs.loc[label, "z_solenoid"]
        if z_sol == center:
            shift += 0.5
        center = z_sol
        z = track.data[label]["field_z"] + z_sol
        Bz = track.data[label]["field_Bz"]

        Bz = np.diff(Bz)/np.diff(z)
        z = z[1:] - np.diff(z)

        Bz *= maxe/np.max(Bz)/2
        plt.plot((eps["z"]+shift)*100, eps["eps_z_rel"].values, "-%s"%color[counter], label=label)
        plt.plot((z+shift)*100, -Bz, "--%s"%color[counter])
        counter += 1
    plt.xlabel("Arbitrary z position [cm]", fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel("RMS z-Emittance [pi*keV*mm]", fontsize=20)
    plt.legend(loc="upper right", fontsize=20)
    plt.axis([100,200,-0.0003,0.0006])
    plt.grid()
    plt.show()


def plot_all_emits(track, core, label):
    plt.figure(figsize=(9,9))
    color = {
        0: "r",
        1: "g",
        2: "b",
        3: "k"
    }
    eps = track.data[label]["eps"]
    maxe = (eps["eps_z"].max())
    z_sol = track.runs.loc[label, "z_solenoid"]
    z = track.data[label]["field_z"] + z_sol
    Bz = track.data[label]["field_Bz"]
    Bz *= maxe/np.max(Bz)/2
    for i in ["x", "y", "z"]:
        plt.plot((eps["z"])*100, eps["eps_z"].values, "-%s"%color[counter], label=label)

    plt.plot((z)*100, -Bz, "--k")
    plt.xlabel("Arbitrary z position [cm]", fontsize=20)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.ylabel("Norm. RMS Emittance [keV*mm]", fontsize=20)
    plt.legend(loc="upper right", fontsize=20)
    plt.show()




def focusing(track, core, label, compute=False):
    if compute:
        track.use_dat("plugins/astra/workspace/fields/"+lbl+".dat", normalize=True, label=lbl)
        track.overview_run()
        track.get_focal_region()
        try:
            track.focal_run()
        except:
            try:
                track.focal_run(step=0.15)
            except:
                track.focal_run(step=0.2)
        track.fit_focal_traj()
        track.fit_cs_expansion()
        track.check_felddurchgang(label=lbl, title="Smoothed, cut off wide field with soft edge")
        track.check_ray_fitting(label=lbl)



def larmor(track, core, compute=False):
    """kek."""
    labels_soft = [
        "thin_soft",
        "mid_soft",
        "wide_soft"
    ]
    labels_hard = [
        "thin_sharp",
        "mid_sharp",
        "wide_sharp"
    ]

    indices = {}

    fmts = [
        "or",
        "ob",
        "og",
        "dr",
        "db",
        "dg"
    ]
    fmt = {}
    disp_labels = [
        "thin/soft",
        "mid/soft",
        "wide/soft",
        "thin/sharp",
        "mid/sharp",
        "wide/sharp"
    ]
    disp_label = {}

    for i in range(3):
        fmt[labels_soft[i]] = fmts[i]
        fmt[labels_hard[i]] = fmts[i+3]
        disp_label[labels_soft[i]] = disp_labels[i]
        disp_label[labels_hard[i]] = disp_labels[i+3]
        indices[labels_soft[i]] = i
        indices[labels_hard[i]] = i+3


    zmax = {}
    s = {}
    idx = {}

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16, 8))
    fig.suptitle(
        "Larmor angle vs. initial radial position",
        fontsize=32)

    F1s = []
    phis = []
    As = {}
    Bs = {}
    pots = {}
    for label in [*labels_soft, *labels_hard]:
        if compute:
            track.use_dat(
                "plugins/astra/workspace/fields/"+label+".dat",
                normalize=True, label=label)
            track.overview_run()
            track.get_focal_region()
            core.sample_field(track.data[label]["field_z"], track.data[label]["field_Bz"])
            F1 = core.fint(1)
            track.runs.loc[label, "F1"] = F1
            F1s.append(F1)
        else:
            F1 = track.runs.loc[label, "F1"]
            F1s.append(F1)

        limit = track.runs.loc[label, "z_focal_left"]
        s = track.data[label]["s"].query("zpos<=@limit")
        idx = s.index[-1][0]
        phi0 = s.loc[idx, "turn"].min()
        phis.append(phi0)
        def parabola(r, A, B, pot=2):
            phi = A*r**pot + B
            return phi
        if label in labels_hard:
            axis = ax2
        else:
            axis = ax1
        y = s.loc[idx, "turn"].values
        x = s.loc[0, "r"].values
        bounds = [[0,0],[np.inf,phi0]]
        popt, pcov = curve_fit(parabola, xdata=x, ydata=y, p0=[0, phi0], bounds=bounds, sigma=x**4)
        As[label] = popt[0]
        Bs[label] = popt[1]
        #pots[label] = popt[2]
        #phis[indices[label]] = popt[1]
        rs = np.linspace(0, 30, 100)
        axis.plot(
            s.loc[0, "r"]*1000,
            (s.loc[idx, "turn"] - phi0)/phi0*100,
            fmt[label], markersize=8,
            label=disp_label[label]
        )
        axis.plot(
            rs,
            (parabola(rs/1000, *popt)-popt[1])/popt[1]*100,
            "-k"
        )

    for ax in (ax1, ax2):
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_xlabel("Initial radial position [mm]", fontsize=28)
        ax.plot(0,0,"-k",label="Expansion fit")
        ax.legend(loc="upper left", fontsize=24)
        ax.grid()
        #ax.plot([0,30], [0,0], "--k")
    #plt.axis([0,31,-0.025,1.025])
    ax1.set_ylabel("Deviation from min. observed PhiL [%]", fontsize=28)
    ax1.set_title("Soft edge", fontsize=24)
    ax2.set_title("Sharp edge", fontsize=24)
    plt.show()

    plt.figure(figsize=(16, 5))
    plt.title(
    "Residuals from the 2nd order expansion",
    fontsize = 32
    )
    for label in [*labels_soft, *labels_hard]:
        limit = track.runs.loc[label, "z_focal_left"]
        s = track.data[label]["s"].query("zpos<=@limit")
        idx = s.index[-1][0]
        y = s.loc[idx, "turn"].values
        x = s.loc[0, "r"].values
        plt.plot(
            x/mm,
            (y - parabola(x, As[label], Bs[label]))/y*100,
            fmt[label], label=disp_label[label],
            markersize=8
        )
    plt.legend(loc="upper left", fontsize=24)
    plt.xlabel("Initial radial position [mm]", fontsize=28)
    plt.xticks(fontsize=24)
    plt.ylabel("data - model [%]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.show()

    plt.figure(figsize=(16, 9))
    plt.title(
    "Larmor angle for paraxial particles vs. F1",
    fontsize=32)
    F1s = np.array(F1s)
    phis = np.array(phis)

    def slope(x, A):
        return A*x

    A0 = const.e/2/core.model.impuls_SI(track.E)/np.pi
    print("Guess for slope: %.2f"%A0)
    A, dA = curve_fit(slope, xdata=F1s, ydata=phis, p0=[A0])
    plt.plot(
        F1s/mm, slope(F1s, A)/mm,
        "-k", label="Slope fit: %.6f\npredicted: %.6f" % (A, A0),
        linewidth=2)
    for label in [*labels_soft, *labels_hard]:
        plt.plot(
            F1s[indices[label]]/mm, phis[indices[label]]/mm,
            fmt[label], markersize=24,
            label=disp_label[label])
    plt.plot(
        F1s/mm, slope(F1s, A)/mm,
        "-k",
        linewidth=2)
    plt.xlabel("F1 [mT*m]", fontsize=28)
    plt.xticks(fontsize=24)
    plt.ylabel("Minimum recorded PhiL [mrad/pi]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper left", fontsize=24)
    plt.show()
    print("dA: %.2f"%dA)

#def larmor_dphi_
