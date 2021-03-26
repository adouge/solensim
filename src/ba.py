"""BA scripts."""

import numpy as np
import numpy.linalg as lina
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from solensim.aux import *
import scipy.constants as const
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import scipy.interpolate as interpolate

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

#import pysnooper

labels_soft = [
    "thin_soft",
    "mid_soft",
    "wide_soft"
]
labels_hard = [
    "thin_hard",
    "mid_hard",
    "wide_hard"
]

labels = [*labels_hard, *labels_soft]

indices = {}

fmts = [
    "or",
    "ob",
    "og",
    "sr",
    "sb",
    "sg"
]
fmt = {}
disp_labels = [
    "thin/soft",
    "mid/soft",
    "wide/soft",
    "thin/hard",
    "mid/hard",
    "wide/hard"
]
disp_label = {}

for i in range(3):
    fmt[labels_soft[i]] = fmts[i]
    fmt[labels_hard[i]] = fmts[i+3]
    disp_label[labels_soft[i]] = disp_labels[i]
    disp_label[labels_hard[i]] = disp_labels[i+3]
    indices[labels_soft[i]] = i
    indices[labels_hard[i]] = i+3

def generate_field(p, zmax, title, core, f=1.25, E=3.5):
    core.FM = "biswas"
    core.bcalc_zmax = zmax
    z = core.get_z()
    Bz = core.get_Bz(p)
    core.sample_field(z, Bz)
    Bz *= core.get_scale_factor(f, E)
    maxB = np.max(Bz)*1000
    maxZ = zmax*100
    core.sample_field(z, Bz)
    FWHM = core.get_fwhm()*1000
    F1 = core.fint(1)*10**3
    F2 = core.fint(2)*10**6
    plt.figure(figsize=(8, 5))
    plt.title(title, fontsize=24)
    plt.plot(z*100, Bz*1000, "-k")
    plt.axis([-maxZ, maxZ, 0, 1.1*maxB])
    plt.xlabel("On-axis position [cm]", fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("On-axis field [mT]", fontsize=22)
    lefttext = "f: %.2f @ %.1f MeV\nFWHM: %.0f mm" % (f, E, FWHM)
    righttext = "F1: %.2f m*mT  \nF2: %.1f m*mT2" % (F1, F2)
    plt.text(-0.95*maxZ, maxB*1.05, lefttext, fontsize=20, verticalalignment="top", horizontalalignment="left")
    plt.text(0.95*maxZ, maxB*1.05, righttext, fontsize=20, verticalalignment="top", horizontalalignment="right")
    plt.show()
    field = pd.DataFrame(np.array([z, Bz]).transpose(), columns=["z", "Bz"])
    return field

def show_field(z, Bz, title, core, f=1.25, E=3.5):
    core.sample_field(z, Bz)
    Bz2 = Bz*core.get_scale_factor(f, E)
    maxB = np.max(Bz2)*1000
    maxZ = np.max(z)*100
    core.sample_field(z, Bz2)
    FWHM = core.get_fwhm()*1000
    F1 = core.fint(1)*10**3
    F2 = core.fint(2)*10**6
    plt.figure(figsize=(8, 5))
    plt.title(title, fontsize=24)
    plt.plot(z*100, Bz2*1000, "-k")
    plt.axis([-maxZ, maxZ, 0, 1.1*maxB])
    plt.xlabel("On-axis position [cm]", fontsize=22)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel("On-axis field [mT]", fontsize=22)
    lefttext = "f: %.2f @ %.1f MeV\nFWHM: %.0f mm" % (f, E, FWHM)
    righttext = "F1: %.2f m*mT  \nF2: %.1f m*mT2" % (F1, F2)
    plt.text(-0.95*maxZ, maxB*1.05, lefttext, fontsize=20, verticalalignment="top", horizontalalignment="left")
    plt.text(0.95*maxZ, maxB*1.05, righttext, fontsize=20, verticalalignment="top", horizontalalignment="right")
    plt.show()
    field = pd.DataFrame(np.array([z, Bz2]).transpose(), columns=["z", "Bz"])
    return field

def second_moment(a,b):
    n = len(a)
    return np.sum(a*b)/n - np.sum(a)*np.sum(b)/n**2

#@pysnooper.snoop()
def make_emits(track, core, label):
    m_e = const.m_e*const.c**2/const.e
    s = track.data[label]["s"].copy()
    zpos = track.data[label]["zpos"]
    zmax = zpos[-1]
    tr = s.query("@zmax - 0.5 <= zpos")
    pz = tr.loc[zmax, "pz"].mean()
    emittances = pd.DataFrame(columns=["z", "eps_x", "eps_y"])
    indices = ["x", "y", "r"]
    for idx in indices:
        tr = tr.swaplevel()
        for part in tr.index.levels[0][1:]:
            tr.loc[(part, zmax), "%sprime" % idx] = np.mean(np.diff(tr.loc[part, idx])/np.diff(tr.loc[part, "z"]))
        tr = tr.swaplevel()
        i = np.array([*tr.loc[zmax, idx].values, *(-tr.loc[zmax, idx].values)])
        iprime = np.array([*tr.loc[zmax, "%sprime" % idx].values, *(-tr.loc[zmax, "%sprime" % idx].values)])
        cov = np.cov(i, iprime)
        print(cov)
        eps = pz/m_e*np.sqrt(lina.det(cov))
        track.runs.loc[label, "eps_%s_tr" % idx] = eps
    i = tr.loc[zmax, ["x", "y"]].values.transpose()
    iprime = tr.loc[zmax, ["xprime", "yprime"]].values.transpose()
    cov = np.cov(i, iprime)
    print(cov)
    eps = pz/m_e*np.sqrt(lina.det(cov))
    track.runs.loc[label, "eps_4d_tr"] = eps
    track.runs.loc[label, "eps_xy_tr"] = np.mean(track.runs.loc[label, ["eps_x_tr", "eps_y_tr"]])
    for z in zpos:
        b = s.loc[z].copy()
        #zref = track.data[label]["pref"].loc[(z, 0), "z"]
        #zref = np.mean(b["z"])
        #pzref = np.mean(b["pz"])
        for index in indices:
            i = b[index].values.astype(np.float64)
            i -= i.mean()
            pi = b["p"+index].values.astype(np.float64)
#            m1 = second_moment(i, i)
#            m2 = second_moment(pi, pi)
#            m3 = second_moment(i, pi)
#            eps2 = m1*m2 - m3**2
            cov = np.cov(i, pi)
#            eps2 = cov[0,0]*cov[1,1]-cov[0,1]*cov[1,0]
            eps2 = lina.det(cov)
            eps = np.sqrt(eps2)/m_e
            emittances.loc[z, "eps_"+index] = eps
#            emittances.loc[z, index+"_ms"] = m1#np.sqrt(m1)
#            emittances.loc[z, "p"+index+"_ms"] = m2#np.sqrt(m2)
#            emittances.loc[z, "cov2_"+index] = m3**2
        cov2 = np.cov(b.loc[:, ["x", "y"]].values.astype(np.float64).transpose(), b.loc[:, ["px", "py"]].values.astype(np.float64).transpose())
        emittances.loc[z, "eps_4d"] = np.sqrt(lina.det(cov2))/m_e**2
        emittances.loc[z, "z"] = z
    track.data[label]["eps"] = emittances

def emittances(track, core, label=None, compute=False):
    if label is not None:
        if type(label) == list:
            labels = label
        else:
            labels = [label]
    else:
        pass

    if compute:
        track.sig_r = 1
        track.N = 10000
        for lbl in labels:
            track.use_dat("plugins/astra/workspace/fields/"+lbl+".dat", normalize=True, label=lbl)
            track.overview_run()
            track.field_width=track.runs.loc[lbl, "field_width"]
            z = track.data[lbl]["field_z"]
            Bz = track.data[lbl]["field_Bz"]
            core.sample_field(z, Bz)
            track.runs.loc[lbl, "F1"] = core.fint(1)
            dz = z[:-1] - np.mean(np.diff(z))
            dBz = np.diff(Bz)/np.diff(z)
            ddz = dz[:-1] - np.mean(np.diff(dz))
            func = interpolate.interp1d(
                dz,
                dBz**2,
                fill_value=0,
                bounds_error=False
            )
            track.runs.loc[lbl, "F3"] = integrate.quad(func, -np.inf, np.inf)[0]
            track.data[lbl]["eps_astra"] = track.astra.read_tremit()
            make_emits(track, core, lbl)

    epsas = {}
    for lbl in labels:
        ea = track.data[lbl]["eps_astra"]
        ea["X"]["eps_y"] = ea["Y"]["eps_y_nrms"].copy()*1e-6
        ea = ea["X"].copy()
        ea["eps_x"] = ea["eps_x_nrms"]*1e-6
        ea["eps_4d"] = ea["eps_x"]*ea["eps_y"]
        ea["eps_xy"] = 0.5*(ea["eps_x"] + ea["eps_y"])
        track.data[lbl]["eps"]["eps_xy"] = 0.5*(track.data[lbl]["eps"]["eps_x"] + track.data[lbl]["eps"]["eps_y"])

        epsas[lbl] = ea.copy()

    fig, (a1, a2) = plt.subplots(1,2, figsize=(16, 6), sharey=True)
    fig.suptitle("Norm. RMS transverse emittance, own calculation vs. ASTRA", fontsize=32)
    for i in (0,1):
        lbl = ["wide_soft", "wide_hard"][i]
        eps = track.data[lbl]["eps"]
        epsa = epsas[lbl]
        z_sol = track.runs.loc[lbl, "z_solenoid"]
        (a1, a2)[i].plot((eps["z"]-z_sol)*1e2,eps["eps_xy"]*1e6, "--k", label="Own calculation")
        (a1, a2)[i].plot((epsa["z"]-z_sol)*1e2, epsa["eps_xy"]*1e6, "-k", label="ASTRA")
        (a1, a2)[i].set_xlabel("z [cm]", fontsize=24)
        (a1, a2)[i].tick_params(labelsize=24)
        (a1, a2)[i].set_xlim([-50,50])
        #(a1, a2)[i].set_title(["Wide field, soft edge", "Wide field, hard edge"][i], fontsize=28)
    a1.set_ylabel("[pi*mm*mrad]", fontsize=24)
    a1.set_ylim([-0.1,8])
    a1.text(-48, 7, "Not accounting for correlation", fontsize=24)
    a1.set_yticks([0, 2, 4, 8])
    a1.tick_params(size=10, axis="y")
    #a2.legend(fontsize=28, loc="upper right")
    plt.show()

    fig, (a1, a2) = plt.subplots(1,2, figsize=(16, 6), sharey=True)
    #fig.suptitle("Accounting for correlation", fontsize=32)
    for i in (0,1):
        lbl = ["wide_soft", "wide_hard"][i]
        eps = track.data[lbl]["eps"]
        epsa = epsas[lbl]
        z_sol = track.runs.loc[lbl, "z_solenoid"]
        (a1, a2)[i].plot((eps["z"]-z_sol)*1e2, eps["eps_4d"]**0.5*1e9, "--k", label="Own calculation")
        (a1, a2)[i].plot((epsa["z"]-z_sol)*1e2, epsa["eps_xy"]*1e9, "-k", label="ASTRA")
        (a1, a2)[i].set_xlabel("z [cm]", fontsize=24)
        (a1, a2)[i].tick_params(labelsize=24)
        (a1, a2)[i].set_xlim([-50,50])
        (a1, a2)[i].set_title(["Wide field, soft edge", "Wide field, hard edge"][i], fontsize=28)
    a1.set_ylabel("[pi*mm*mrad]*1e-3", fontsize=24)
    a1.set_ylim([0,0.6])
    a1.text(-48, 0.5, "Accounting for correlation", fontsize=24)
    a2.legend(fontsize=24, loc="lower right")
    plt.show()

    plt.figure(figsize=(16,6))
    plt.title("Norm. RMS x,y emittance by field [pi*mrad*mm]*1e-3", fontsize=32)
    plt.xlabel("Own calculation", fontsize=28)
    plt.ylabel("ASTRA", fontsize=28)

    own = [track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_4d"]**0.5*1e9 for lbl in labels]
    astra = [epsas[lbl].loc[epsas[lbl].index[-1], "eps_xy"]*1e9 for lbl in labels]
    slope = lambda x, A: A*x
    A, dA = curve_fit(slope, xdata=own, ydata=astra, p0=[1])
    print(A, dA)
    for label in labels:
        eps = track.data[label]["eps"]
        epsa = epsas[label]
        plt.plot(eps.loc[eps.index[-1], "eps_4d"]**0.5*1e9, epsa.loc[epsa.index[-1], "eps_xy"]*1e9, fmt[label], label=disp_label[label], markersize=24)
    x = np.linspace(0, 2.25)
    plt.plot(x, slope(x, A), "-k", label="Slope fit", linewidth=2)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(fontsize=28, loc="lower left", bbox_to_anchor=(0.4, 0))
    plt.axis([0,2.2,-0.1,2.2])
    plt.text(0.1, 1.9, "Slope: 1-%.2e" % (1-A), fontsize=28)
    plt.show()

    for lbl in labels:
        z = track.data[lbl]["field_z"]
        Bz = track.data[lbl]["field_Bz"]
        dz = z[:-1] - np.mean(np.diff(z))
        ddz = dz[:-1] - np.mean(np.diff(dz))
        dBz = np.diff(Bz)/np.diff(z)
        ddBz = np.diff(dBz)/np.diff(dz)
        func = interpolate.interp1d(
            dz,
            dBz**2,
            fill_value=0,
            bounds_error=False
        )
        func2 = interpolate.interp1d(
            ddz,
            ddBz**2,
            fill_value=0,
            bounds_error=False
        )
        track.runs.loc[lbl, ["F3", "F4"]] = [
            integrate.quad(func, -np.inf, np.inf)[0],
            integrate.quad(func2, -np.inf, np.inf)[0]
        ]
    c1s =  np.array([1/2*track.runs.loc[lbl, "F3"]/track.runs.loc[lbl, "F2"] for lbl in labels])
    c2s = np.array([5/64*track.runs.loc[lbl, "F4"]/track.runs.loc[lbl, "F2"] for lbl in labels])
    r = (track.astra.beam["x"]**2 + track.astra.beam["y"]**2)**0.5
    sig = np.mean((track.astra.beam["x"].std(), track.astra.beam["y"].std()))
    #sig = 1e-3
    print("Sigma: ", sig)
    epsapprox = np.sqrt(6)*sig**4*c1s/1.5
    pz = 3.4625e6
    m_e = const.m_e*const.c**2/const.e
    epstheory = pz/m_e*np.sqrt(6)*sig**4/1.5*np.sqrt(c1s**2 + 20*sig**2*c1s*c2s + 120*sig**4*c2s**2)
    #print("epstheory:    ", epstheory)
    #test = [track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_r"] for lbl in labels]
    #test = np.array(test)
    #print(np.array(test)/epstheory/2)
    #print(c1s)
    #print(c2s)
    #print(epsapprox*1e9)
    #print(epstheory*1e9)

    #print([track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_4d"]**0.5 for lbl in labels])

    plt.figure(figsize=(16,6))
    plt.title("Norm. RMS trace space emittance [pi*mrad*mm]*1e-3", fontsize=32)
    plt.ylabel("Own calculation", fontsize=28)
    plt.xlabel("Theoretical prediction", fontsize=28)

    own = [track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_4d"]**0.5 for lbl in labels]
    own2 = [track.runs.loc[lbl, "eps_xy_tr"] for lbl in labels]
    #print("=====================================")
    #print(own)
    #print(own2)
    #print(epstheory)
    #print("=====================================")
    #print(np.array(own)/test)
    slope = lambda x, A: A*x
    A, dA = curve_fit(slope, ydata=own2, xdata=epstheory)
    print(A, dA)
    for label in labels:
        eps = track.data[label]["eps"]
        epsa = epsas[label]
        plt.plot(epstheory[indices[label]]*1e9, track.runs.loc[label, "eps_xy_tr"]*1e9, fmt[label], label=disp_label[label], markersize=24)
    x = np.linspace(0, 8)
    plt.plot(x, slope(x, A), "-k", label="Slope fit", linewidth=2)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(fontsize=28, loc="lower left", bbox_to_anchor=(0.6, 0))
    plt.axis([0,7.5,-0.1,7.5])
    plt.text(0.1, 4.25, "Slope: %.2e" % A, fontsize=24)
    plt.show()

    plt.figure(figsize=(15,4))
    plt.title("Tracking results vs. theoretical predictions", fontsize=30)
    diffs = []
    diff0 = (epstheory[indices["thin_soft"]] - eps.loc[eps.index[-1], "eps_4d"]**0.5*A)/epstheory[indices[label]]*100
    for label in labels:
        diff = (-epstheory[indices[label]] + track.runs.loc[label, "eps_xy_tr"])/epstheory[indices[label]]*100
        eps = track.data[label]["eps"]
        plt.plot(
            indices[label] if label in labels_soft else indices[label] - 3,
            diff,
            fmt[label],
            markersize=20,
            label=({"mid_hard": "Hard edge", "mid_soft": "Soft edge"}[label] if label in ["mid_hard", "mid_soft"] else None)
        )
        diffs.append(diff)
    #plt.plot([-0.5, 4], [0, 0], "-k")
    diffs = np.array(diffs)
    #plt.axis([-1,3, -50, -56])
    plt.legend(loc="lower left", fontsize=22, bbox_to_anchor=(1.5/2.75, 0.4))
    plt.xlabel("Field width", fontsize=24)
    plt.xticks(fontsize=22)#, labels=["Thin", "Mid", "Wide"], ticks=[0, 1, 2])
    plt.ylabel("Deviation [%]", fontsize=24)
    plt.yticks(fontsize=22)#, ticks=[60,65,70,75])
    plt.grid()
    plt.show()

    for lbl in labels:
        own1 = track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_xy"]
        own2 = track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_4d"]**0.5
        astra = epsas[lbl].loc[epsas[lbl].index[-1], "eps_xy"]
        print("Diff: %s: %.3e%% XY, %.3e%% sqrt(4D)" % (lbl, (own1-astra)/astra*100, (own2-astra)/astra*100))

    for lbl in labels:
        print("%s: %.3e rad, %.3e transverse" % (lbl, track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_r"]*1e9, track.data[lbl]["eps"].loc[track.data[lbl]["eps"].index[-1], "eps_4d"]**0.5*1e9))

    print(sig*1e3)

def focusing(track, core, label=None, compute=False, expand=False, sigma=2, order=1):
    if label is not None:
        if type(label) == list:
            labels = label
        else:
            labels = [label]
    else:
        pass

    if compute:
        track.sig_r = 20/np.sqrt(2)
        track.N = 250
        for lbl in labels:
            track.use_dat("plugins/astra/workspace/fields/"+lbl+".dat", normalize=True, label=lbl)
            track.overview_run(beam="uniform")
            track.field_width=track.runs.loc[lbl, "field_width"]
            track.get_focal_region()
            try:
                track.focal_run()
            except:
                track.focal_run(step=0.2)
            track.fit_focal_traj(model="axial")
            z = track.data[lbl]["field_z"]
            Bz = track.data[lbl]["field_Bz"]
            core.sample_field(z, Bz)
            track.runs.loc[lbl, ["F1","F3","F4"]] = [core.fint(1), core.fint(3), core.fint(4)]
        track.runs["Delta_f"] = track.runs["f_max_observed"] - track.runs["f_min_observed"]

# DF vs. F1
    plt.figure(figsize=(16, 6))
    plt.title(
    "Min/max difference for observed foci vs. F1",
    fontsize=32)

    for label in labels:
        plt.plot(
            track.runs.loc[label, "F1"]/mm, track.runs.loc[label, "Delta_f"]*1e2,
            fmt[label], markersize=24,
            label=disp_label[label])
    plt.xlabel("F1 [mT$\cdot$m]", fontsize=28)
    plt.axis([4.5, 14.5, 0, 12])
    plt.xticks(fontsize=24)
    plt.ylabel("Focal region size [cm]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper right", fontsize=24)
    plt.show()

# True F vs. F1
    plt.figure(figsize=(16, 9))
    plt.title(
    "Max. observed focus vs. F1",
    fontsize=32)

    def slope(x, A):
        return A*x + 1.5

    A, dA = curve_fit(slope, xdata=track.runs["F1"], ydata=track.runs["f_max_observed"], p0=[0])
    for label in labels:
        plt.plot(
            track.runs.loc[label, "F1"]/mm, track.runs.loc[label, "f_max_observed"],
            fmt[label], markersize=24,
            label=disp_label[label])
    plt.plot(
        track.runs["F1"]/mm, slope(track.runs["F1"], A),
        "-k", label="Slope fit",
        linewidth=2)
    plt.text(10, 1.501, "Slope: %.6fe-3 [mT^-1]" % A, verticalalignment="bottom", fontsize=24)
    plt.xlabel("F1 [mT$\cdot$m]", fontsize=28)
    plt.axis([4.5, 14.5, 1.50, 1.55])
    plt.xticks(fontsize=24)
    plt.ylabel("Maximum observed focal length [m]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper left", fontsize=24)
    plt.show()
# Df expansion
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16, 8))
    fig.suptitle(
        "Focal length vs. initial radial position",
        fontsize=32)
    if expand or compute:
        for label in labels:
            track.run_label = label
            track.fit_cs_expansion(order=order, sigma_order=sigma)
    lim = []
    r = np.linspace(0,10, num=100)
    for label in labels:
        fits = track.data[label]["fits"]
        c = track.data[label]["exp_coeff"]
        f_real = track.runs.loc[label, "f"]
        if label in labels_soft:
            axis = ax1
        else:
            axis = ax2
        zsol = track.runs.loc[label, "z_solenoid"]
        axis.plot(fits["r0"]*1e3, ((fits["z_f"]-zsol)/f_real-1)*100, fmt[label], markersize=8, label=disp_label[label])
        y = (track.model.f_expansion(r/1e3, f_real, *c)/f_real-1)*100
        axis.plot(r, y, "-k")
        lim.append(np.min(y))

    for ax in (ax1, ax2):
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_xlabel("Initial radial position [mm]", fontsize=28)
        ax.plot(0,0,"-k",label="Expansion fit")
        ax.legend(loc="lower left", fontsize=24)
        ax.grid()
        ax.set_xlim([0, 10])
        #ax.plot([0,30], [0,0], "--k")
    #plt.axis([0,31,-0.025,1.025])
    ax1.set_ylabel("Deviation from max. observed f [%]", fontsize=28)
    ax1.set_ylim([np.min(lim), 0.1])
    ax1.set_title("Soft edge", fontsize=24)
    ax2.set_title("Hard edge", fontsize=24)
    plt.show()

    plt.figure(figsize=(16, 5))
    suffix = {
        1: "st",
        2: "nd",
        3: "rd"
    }
    plt.title(
        "Residuals from the %d%s order expansion, %s weighing" % (order, suffix[order], ("no" if sigma == 0 else "r%d" % sigma)),
        fontsize=32
    )
    diffs = []
    for label in labels:
        fits = track.data[label]["fits"]
        c = track.data[label]["exp_coeff"]
        f_real = track.runs.loc[label, "f"]
        modeled = track.model.f_expansion(fits["r0"], f_real, *c)
        y = fits["z_f"] - track.runs.loc[label, "z_solenoid"]
        diff = (y - modeled)/y*100
        plt.plot(
            fits["r0"]*1e3,
            diff*1e3,
            fmt[label], label=disp_label[label],
            markersize=8
        )
        diffs.append(diff)
    diffs = np.array(diffs)
    plt.axis([0,10,-15, 15])
    #plt.legend(loc="upper left", fontsize=24)
    plt.xlabel("Initial radial position [mm]", fontsize=28)
    plt.xticks(fontsize=24)
    plt.ylabel("data - model [0.001%]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.grid()
    plt.show()

    print("TABLE:")
    print("Field & $C_2\\unit{~[m^{-1}}$")
    for label in labels:
        print("%s & $%.2f$" % (disp_label[label], track.data[label]["exp_coeff"][0]))

def aberration(track, core, label=None, compute=False, expand=False, sigma=2, order=1):
    if label is not None:
        if type(label) == list:
            labels = label
        else:
            labels = [label]
    else:
        pass

    if compute:
        track.sig_r = 1
        track.N = 1000
        for lbl in labels:
            track.use_dat("plugins/astra/workspace/fields/"+lbl+".dat", normalize=True, label=lbl)
            track.overview_run()
            track.field_width=track.runs.loc[lbl, "field_width"]
            track.get_focal_region()
            try:
                track.focal_run()
            except:
                track.focal_run(step=0.2)
            track.fit_focal_traj(model="axial")
            z = track.data[lbl]["field_z"]
            Bz = track.data[lbl]["field_Bz"]
            core.sample_field(z, Bz)
            track.runs.loc[lbl, "F1"] = core.fint(1)
            dz = z[:-1] - np.mean(np.diff(z))
            ddz = dz[:-1] - np.mean(np.diff(dz))
            dBz = np.diff(Bz)/np.diff(z)
            ddBz = np.diff(dBz)/np.diff(dz)
            func = interpolate.interp1d(
                dz,
                dBz**2,
                fill_value=0,
                bounds_error=False
            )
            func2 = interpolate.interp1d(
                ddz,
                ddBz**2,
                fill_value=0,
                bounds_error=False
            )
            track.runs.loc[lbl, ["F3", "F4"]] = [
                integrate.quad(func, -np.inf, np.inf)[0],
                integrate.quad(func2, -np.inf, np.inf)[0]
            ]
        track.runs["Delta_f"] = track.runs["f_max_observed"] - track.runs["f_min_observed"]

    if expand or compute:
        for label in labels:
            track.run_label = label
            track.fit_cs_expansion(order=order, sigma_order=sigma)

# Df expansion
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16, 8))
    fig.suptitle(
        "Focal length expansion in the paraxial region by field",
        fontsize=32)

    lim = []
    r = np.linspace(0,3, num=100)
    for label in labels:
        fits = track.data[label]["fits"]
        c = track.data[label]["exp_coeff"]
        f_real = track.runs.loc[label, "f"]
        if label in labels_soft:
            axis = ax1
        else:
            axis = ax2
        zsol = track.runs.loc[label, "z_solenoid"]
        axis.plot(fits["r0"]*1e3, ((fits["z_f"]-zsol)/f_real-1)*100, fmt[label], markersize=8, label=disp_label[label])
        y = (track.model.f_expansion(r/1e3, f_real, *c)/f_real-1)*100
        axis.plot(r, y, "-k")
        lim.append(np.min(y))

    for ax in (ax1, ax2):
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_xlabel("Initial radial position [mm]", fontsize=28)
        ax.plot(0,0,"-k",label="Expansion fit")
        ax.legend(loc="lower left", fontsize=24)
        ax.grid()
        ax.set_xlim([0, 3])
        #ax.plot([0,30], [0,0], "--k")
    #plt.axis([0,31,-0.025,1.025])
    ax1.set_ylabel("Deviation from max. observed f [%]", fontsize=28)
    ax1.set_ylim([np.min(lim), 0.05])
    ax1.set_title("Soft edge", fontsize=24)
    ax2.set_title("Hard edge", fontsize=24)
    plt.show()

    plt.figure(figsize=(16, 5))
    suffix = {
        1: "st",
        2: "nd",
        3: "rd"
    }
    plt.title(
        "Residuals from the 1st order expansion, no weighing",
        fontsize=32
    )
    diffs = []
    for label in labels:
        fits = track.data[label]["fits"]
        c = track.data[label]["exp_coeff"]
        f_real = track.runs.loc[label, "f"]
        modeled = track.model.f_expansion(fits["r0"], f_real, *c)
        y = fits["z_f"] - track.runs.loc[label, "z_solenoid"]
        diff = (y - modeled)/y*100
        plt.plot(
            fits["r0"]*1e3,
            diff*1e3,
            fmt[label], label=disp_label[label],
            markersize=8
        )
        diffs.append(diff)
    diffs = np.array(diffs)
    plt.axis([0,3,-12, 12])
    #plt.legend(loc="upper left", fontsize=24)
    plt.xlabel("Initial radial position [mm]", fontsize=28)
    plt.xticks(fontsize=24)
    plt.ylabel("data - model [0.001%]", fontsize=28)
    plt.yticks(fontsize=24, ticks=[-8, -4, 0, 4, 8])
    plt.grid()
    plt.show()

# C2 vs. F1
    def cs_model(F3, F4):
        pz = (3.4625e6*const.e/const.c)
        R = 2e-3
        c_s = const.e**2 / 4 / pz**2 * R**4 * (F3 + const.e**2 / 3 / pz**2 * F4)
        return c_s*1.5

    plt.figure(figsize=(16, 6))
    plt.title(
    "First order expansion coefficient vs. F1",
    fontsize=32)

    for label in labels:
        plt.plot(
            track.runs.loc[label, "F1"]/mm, track.runs.loc[label, "c2"],
            fmt[label], markersize=24,
            label=disp_label[label])
    plt.xlabel("F1 [mT$\cdot$m]", fontsize=28)
    #plt.axis([4.5, 14.5, 0, 12])
    plt.xticks(fontsize=24)
    plt.ylabel("C2 [m^-2]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper right", fontsize=24)
    plt.show()

# C2 vs. F1
    def cs_model(F3, F2):
        c2 = 1 / 2 * F3 / F2
        return c2

    def cs_model2(F3, F4, F2):
        c2 = 1 / 2 * F3 / F2
        c4 = 5 / 64 * F4 / F2
        return c4

    cs_theory = [cs_model(*track.runs.loc[lbl, ["F3", "F2"]].values) for lbl in labels]

    def model(x, A):#, B):
        return A*x# + B

    #(A, B), (dA, dB) = curve_fit(model, ydata=cs_theory, xdata=[track.runs.loc[lbl, "c2"] for lbl in labels])
    A, dA = curve_fit(model, ydata=cs_theory, xdata=[track.runs.loc[lbl, "c2"] for lbl in labels])

    plt.figure(figsize=(16, 6))
    plt.title(
    "C2, theoretical vs. tracking results",
    fontsize=32)

    for label in labels:
        plt.plot(
            track.runs.loc[label, "c2"],
            cs_theory[indices[label]],
            fmt[label], markersize=24,
            label=disp_label[label])
    x = np.linspace(track.runs["c2"].min(), track.runs["c2"].max())
    plt.plot(
        x,
        model(x, A),
        "-k",
        label="Slope fit",
        linewidth=2
    )
    print([A, dA])
    plt.text(0, 700, "Slope: 1 + %.2e +- %.2e" % (A-1, dA), fontsize=24, verticalalignment="top")
    plt.xlabel("Measured C2 [m^-2]", fontsize=28)
    plt.axis([-25, 800, -25, 800])
    plt.xticks(fontsize=24)
    plt.ylabel("Theoretical C2 [m^-2]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="center left", bbox_to_anchor=(0.5,0.5), fontsize=24)
    plt.show()

    print("TABLE:")
    print("Field & $C_2\\unit{~[m^{-1}}$")
    for label in labels:
        print("%s & $%.2f$" % (disp_label[label], track.data[label]["exp_coeff"][0]))

    plt.figure(figsize=(16, 4))
    plt.title(
        "Theoretical predictions vs. tracking results",
        fontsize=32
    )
    diffs = []
    for label in labels:
        diff = (cs_theory[indices[label]] - track.runs.loc[label, "c2"])/cs_theory[indices[label]]*100
        plt.plot(
            indices[label] if label in labels_soft else indices[label] - 3,
            diff,
            fmt[label],
            markersize=20,
            label=({"mid_hard": "Hard edge", "mid_soft": "Soft edge"}[label] if label in ["mid_hard", "mid_soft"] else None)
        )
        diffs.append(diff)
    plt.plot([-0.5, 4], [0, 0], "-k")
    diffs = np.array(diffs)
    plt.axis([-0.5,2.5, -8, 8])
    plt.legend(loc="lower left", fontsize=24, bbox_to_anchor=(1.5/2.75, 0.5))
    plt.xlabel("Field width", fontsize=28)
    plt.xticks(fontsize=24, labels=["Thin", "Mid", "Wide"], ticks=[0, 1, 2])
    plt.ylabel("Model - Tracking [%]", fontsize=28)
    plt.yticks(fontsize=24, ticks=[-6, -4, -2,0,2,4,6])
    plt.grid()
    plt.show()


def larmor(track, core, compute=False):
    """kek."""

    zmax = {}
    s = {}
    idx = {}

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(16, 8))
    fig.suptitle(
        "Larmor angle vs. initial radial position",
        fontsize=32)

    F1s = []
    phis = []
    As = []
    A = {}
    B = {}
    Bs = []
    pots = {}

    def straighten_out(phi):
        if phi < 0:
            phi2 = phi + 1
        else:
            phi2 = phi - 1
        return phi2
    straighten_out_v = np.vectorize(straighten_out)

    for label in [*labels_soft, *labels_hard]:
        if compute:
            track.sig_r = 20 / np.sqrt(2)
            track.N = 250
            track.baseline_f = 1.5

            track.use_dat(
                "plugins/astra/workspace/fields/"+label+".dat",
                normalize=True, label=label)
            track.overview_run(beam_2d=True, beam="uniform")
            track.get_focal_region()
            core.sample_field(track.data[label]["field_z"], track.data[label]["field_Bz"])
            F1 = core.fint(1)
            track.runs.loc[label, "F1"] = F1
            F1s.append(F1)
        else:
            F1 = track.runs.loc[label, "F1"]
            F1s.append(F1)

        limit = track.runs.loc[label, "z_focal_left"] - 0.15
        s = track.data[label]["s"].query("zpos<=@limit")
        s = track.data[label]["s"].copy()
        idx = s.index[-1][0]
        closest = s.loc[0, "r"].idxmin()
        phi0 = straighten_out(s.loc[(idx, closest), "turn"])
        phis.append(phi0)
        def parabola(r, A, B, pot=2):
            phi = A*r**pot + B
            return phi
        if label in labels_hard:
            axis = ax2
        else:
            axis = ax1

        y = straighten_out_v(s.loc[idx, "turn"].values)
        x = s.loc[0, "r"].values
        maxc2 = (y.max() - phi0)/x.max()**2
        bounds = [[0.75*maxc2,0.75*phi0], [1.25*maxc2,1.25*phi0]]
        popt, pcov = curve_fit(parabola, xdata=x, ydata=y, p0=[maxc2, phi0], bounds=bounds, sigma=x**4)
        print(label, popt)
        As.append(popt[0])
        A[label] = popt[0]
        Bs.append(popt[1])
        B[label] = popt[1]
        #pots[label] = popt[2]
        #phis[indices[label]] = popt[1]
        maxr = 10
        rs = np.linspace(0, maxr, 100)
        axis.plot(
            s.loc[0, "r"]*1000,
            (straighten_out_v(s.loc[idx, "turn"]) - phi0)/phi0*100,
            fmt[label], markersize=8,
            label=disp_label[label]
        )
        axis.plot(
            rs,
            (parabola(rs/1000, *popt)-popt[1])/popt[1]*100,
            "-k"
        )

    Bs = np.array(Bs)
    As = np.array(As)
    plt.axis([0, maxr, 0, np.max((parabola(maxr/1e3, As, Bs)-Bs)/Bs*100)])
    ax1.set_xlim([0, maxr])
    for ax in (ax1, ax2):
        ax.tick_params(axis='both', which='major', labelsize=24)
        ax.set_xlabel("Initial radial position [mm]", fontsize=28)
        ax.plot(0,0,"-k",label="Expansion fit")
        ax.legend(loc="upper left", fontsize=24)
        ax.grid()
        #ax.plot([0,30], [0,0], "--k")
    #plt.axis([0,31,-0.025,1.025])
    ax1.set_ylabel("Deviation from prediction [\%]", fontsize=28)
    ax1.set_title("Soft edge", fontsize=24)
    ax2.set_title("Hard edge", fontsize=24)
    plt.show()

    plt.figure(figsize=(16, 9))
    plt.title(
    "Larmor angle for paraxial particles vs. $F_1$",
    fontsize=32)
    F1s = np.array(F1s)
    phis = np.array(phis)

    def slope(x, A):
        return A*x

    A0 = const.e/2/(3.4625e6*const.e/const.c)/np.pi*180*mm
    print("Guess for slope: %.2f"%A0)
    A, dA = curve_fit(slope, xdata=F1s/mm, ydata=np.array(Bs)*180, p0=[A0])
    plt.plot(
        F1s/mm, slope(F1s/mm, A),
        "-k", label="Slope fit:\ \ \ %.6f\npredicted: %.6f" % (A, A0),
        linewidth=2)
    for label in [*labels_soft, *labels_hard]:
        plt.plot(
            F1s[indices[label]]/mm, Bs[indices[label]]*180,
            fmt[label], markersize=24,
            label=disp_label[label])
    #plt.plot(
    #    F1s/mm, slope(F1s/mm, A),
    #    "-k",
    #    linewidth=2)
    plt.xlabel("$F_1$ [mT$\cdot$m]", fontsize=28)
#    plt.axis([4.5, 14.5, 60, 200])
    plt.xticks(fontsize=24)
    plt.ylabel("$\phi_L$ for paraxial particles [degrees]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper left", fontsize=24)
    plt.show()
    print("dA: %.2f"%dA)

    for lbl in labels:
        print("%s: %.2f" % (lbl, 180*phis[indices[lbl]]))

    plt.figure(figsize=(15, 6))
    plt.title(
    "$\phi_L$ after thin, hard edge field vs. $r_0$",
    fontsize=32)
    limit = track.runs.loc["thin_hard", "z_focal_left"] - 0.15
    s = track.data["thin_hard"]["s"].query("zpos<=@limit").copy()
    idx = s.index[-1][0]
    plt.plot(
        s.loc[0, "r"].values/mm, s.loc[idx, "turn"].values*180,
        ".k", label="Tracking data")
    predict = slope(F1s[indices["thin_hard"]]/mm, A0)
    plt.plot([0,10], [predict, predict], "--k", label="Paraxial approximation")
    plt.xlabel("Initial radial position $r_0$ [mm]", fontsize=28)
    plt.axis([0, 10, 12.729, 12.743])
    plt.xticks(fontsize=24)
    plt.ylabel("Larmor angle $\phi_L$ [degrees]", fontsize=28)
    plt.yticks(fontsize=24)
    plt.legend(loc="upper left", fontsize=24)
    plt.show()
    print("dA: %.2f"%dA)
