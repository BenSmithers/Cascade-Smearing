from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi, sin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt



f = open("cummulative_probs.dat",'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
raw_sorted = obj["raw_sorted"]
c_prob = obj["c_prob"]
chi2 = obj["chi2s"]

old = False

ps = [0.10, 0.01]
chis = [-2*log(p) for p in ps]

labels = ["90%", "99%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10)

evs = [4.47,]

def s2(theta):
    si = np.sin(2*theta)
    return si*si

print(chis)
for ev in evs:
    which_sliver = get_loc(ev, msqs)[0]
    chis = np.zeros(shape=(len(theta24s), len(theta34s)))
    for t24 in range(len(theta24s)):
        for t34 in range(len(theta34s)):
            chis[t24][t34] = chi2[t24][t34][which_sliver]

    plt.pcolormesh(theta24s*180/pi, theta34s*180/pi, chis.transpose(), vmin=0, vmax=10, cmap="bwr_r")
    #plt.pcolormesh(s2(theta24s), s2(theta34s), chis.transpose(), vmin=0, vmax=10, cmap="PuBu")
    cbar = plt.colorbar(extend='max')
    cbar.set_label(r"-2$\Delta$LLH", size=14)
    #ct = plt.contour(s2(theta24s), s2(theta34s), chis.transpose(), levels=[-2*log(0.10), -2*log(0.01)], cmap='cool')   
    ct = plt.contour(theta24s*180/pi, theta34s*180/pi, chis.transpose(), levels=[-2*log(0.10), -2*log(0.01)], cmap='cool')   
    set_lbls(ct)

    plt.title(r"90% CL Sensitivity with $\Delta m_{14}^{2}=$"+"{:.2f}".format(msqs[which_sliver]),size=16)
    plt.ylim([0,30])
    plt.xlim([0,30])
    plt.text(10,25, "Smithers Preliminary", color="r",size=14)
    plt.xlabel(r"$\theta_{24}$ [deg]",size=14)
    plt.ylabel(r"$\theta_{34}$ [deg]",size=14)
    plt.show()


t24s = [9*pi/180]
for t24 in t24s:
    which_sliver = get_loc(t24, theta24s)[0]
    chis = np.zeros(shape=(len(msqs), len(theta34s)))
    for t34 in range(len(theta34s)):
        for msq in range(len(msqs)):
            chis[msq][t34] = chi2[which_sliver][t34][msq]
    plt.pcolormesh(msqs, theta34s*180/pi, chis.transpose(), vmin=0, vmax=10, cmap="PuBu")
    cbar = plt.colorbar(extend='max')
    cbar.set_label(r"-2$\Delta$LLH", size=14)
    ct = plt.contour(msqs, theta34s*180/pi, chis.transpose(), levels=[-2*log(0.10), -2*log(0.01)], cmap='cool')   
    set_lbls(ct)
    plt.title(r"90% CL Sensitivity with $\theta_{24}=$"+"{:.2f}".format(theta24s[which_sliver]),size=16)
    plt.ylim([0,30])
    plt.text(2,25, "Smithers Preliminary", color="r",size=14)
    plt.xlabel(r"$\Delta m_{14}^{2}$ [eV$^{2}$]",size=14)
    plt.ylabel(r"$\theta_{34}$ [deg]",size=14)
    plt.show()



# these were written using the older likelihood setup 
if old:
    levels = [90.0, 68.0]

    # sum over the mass-squareds 
    c_prob_msq = np.zeros(shape=(len(theta24s), len(theta34s)))
    for entry in sorted_likelihoods:
        c_prob_msq[entry[0]][entry[1]] += entry[3]

    plt.contour( theta24s, theta34s, 100*c_prob_msq.transpose(), levels=levels)
    plt.xlabel(r"$\theta_{24}$",size=14)
    plt.ylabel(r"$\theta_{34}$",size=14)
    plt.show()


    c_prob_th34 = np.zeros(shape=(len(theta24s), len(msqs)))
    for entry in sorted_likelihoods:
        c_prob_th34[entry[0]][entry[2]] += entry[3]

    plt.contour( msqs, theta24s, 100*c_prob_th34, levels=levels)
    plt.ylabel(r"$\theta_{24}$",size=14)
    plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
    plt.show()

    c_prob_th24 = np.zeros(shape=(len(theta34s), len(msqs)))
    for entry in sorted_likelihoods:
        c_prob_th24[entry[1]][entry[2]] += entry[3]

    plt.contour( msqs, theta34s, 100*c_prob_th24, levels=levels)
    plt.ylabel(r"$\theta_{34}$",size=14)
    plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
    plt.show()


