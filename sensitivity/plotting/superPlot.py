"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi, sin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")


#f = open("../cummulative_probs.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_nosys_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_0.0_0.0_0.0_0.0.dat",'rb')


f = open("/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_0.0_0.0_0.0_0.0.dat", 'rb')


#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1641e0/scaled_cummulative_probs_special_nosys_0.0_0.1641e0_0.2566e0_4.6416e0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probs_special_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_special_nosys_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_special_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1609e0/cummulative_probs_flat_smear_0.0_0.1609e0_0.2276e0_4.4700e0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_flat_smear_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1641e0/cummulative_probs_0.0_0.1641e0_0.2566e0_4.6410e0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_flat_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1609e0/cummulative_probs_0.0_0.1609e0_0.2276e0_4.4700e0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_compare_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1609e0/cummulative_probs_0.0_0.1609e0_0.0_4.4700e0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib/expectations/0.0/cummulative_probs_nosys_0.0_0.0_0.0_0.0.dat",'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
chi2 = np.array(obj["chi2s"]) 

deg = 180./3.1415926



ps = [0.10] #, 0.01]
chis_l = [4.605] # 2 degrees of freedom!! 

labels = ["90%"]#, "99%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((20,10),(30,5) )
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc)

evs = [0.5,3.0, 5.0, 10.0]
#evs = [1.1]

def s2(theta):
    si = np.sin(2*theta)
    return si*si


xs, ys= np.meshgrid( theta24s, theta34s)

scale_x = np.sin(xs)**2
scale_y = (np.cos(xs)**2)*(np.sin(ys)**2)

#minos_24 = 8.0/deg
#minos_34 = 37./deg

minos_24 = 7.0/deg
minos_34 = 26./deg

scaled_minos_24 = np.sin(minos_24)**2
t24s = np.logspace(-3, np.log10(scaled_minos_24), 1000)
scaled_minos_34 = [(np.sin(minos_34)**2)*(np.cos(theta)**2) for theta in t24s]

super_k = np.transpose(np.loadtxt("sk.dat", delimiter=","))
deepcore = np.transpose(np.loadtxt("deepcore.dat", delimiter=","))

for ev in evs:
    which_sliver = get_loc(ev, msqs)[0]
    chis = np.zeros(shape=(len(theta24s), len(theta34s)))
    

    for t24 in range(len(theta24s)):
        for t34 in range(len(theta34s)):
            chis[t24][t34] = chi2[t24][t34][which_sliver]

    plt.pcolormesh(scale_x,scale_y, chis.transpose(),vmin=0, vmax=10, cmap="PuBu")
    cbar = plt.colorbar(extend='max')
    ct = plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='cool')   
    plt.xscale('log')
    plt.title("At {:.2f} eV".format(msqs[which_sliver]) +r"$^{2}$")
#            plt.yscale('log')
    plt.ylim([0.00,0.30])
    plt.xlim([1e-3,0.5])
    plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{24}$",size=14)
    plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
#            plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)

    plt.plot(t24s, scaled_minos_34, color='k', label="Minos")
    plt.plot(super_k[0], super_k[1], color='blue', label="Super K")
    plt.plot(deepcore[0], deepcore[1], color='red', label="DeepCore")
    plt.vlines(x=scaled_minos_24, ymin=0, ymax=scaled_minos_34[-1], color='k')
    plt.legend(loc='upper left')
    plt.show()
