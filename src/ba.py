"""BA scripts."""

import numpy as np
import matplotlib.pyplot as plt
from solensim.aux import *
import scipy.constants as const
from scipy.optimize import curve_fit

def larmor_raw(track, core, compute=False):
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
        popt, pcov = curve_fit(parabola, xdata=x, ydata=y, p0=[0, phi0], bounds=bounds, sigma=x**3)
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

    #plt.axis([0,31,-0.025,1.025])
    ax1.set_ylabel("Deviation from PhiL @r=0 [%]", fontsize=28)
    ax1.set_title("Soft edge", fontsize=24)
    ax2.set_title("Sharp edge", fontsize=24)
    plt.show()
    print(As)
    print(Bs)
    print(phis)
    print(pots)

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
