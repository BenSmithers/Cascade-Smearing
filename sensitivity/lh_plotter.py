from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi
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

evs = [2, 4.47, 10]

print(chis)
for ev in evs:
    which_sliver = get_loc(ev, msqs)[0]
    chis = np.zeros(shape=(len(theta24s), len(theta34s)))
    for t24 in range(len(theta24s)):
        for t34 in range(len(theta34s)):
            chis[t24][t34] = chi2[t24][t34][which_sliver]


    ct = plt.contour(theta24s*180/pi, theta34s*180/pi, chis.transpose(), levels=[-2*log(0.10), -2*log(0.01)])
    set_lbls(ct)
    plt.title(r"90% CL Sensitivity with $\Delta m_{14}^{2}=$"+"{:.2f}".format(msqs[which_sliver]),size=16)
    plt.text(5,85, "Smithers Preliminary", color="r",size=14)
    plt.xlabel(r"$\theta_{24}$ [deg]",size=14)
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


