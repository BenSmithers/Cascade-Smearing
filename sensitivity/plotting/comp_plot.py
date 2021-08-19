"""
This puts the 90% contours for two different sensitivity results on the same plot
"""

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt 

plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

import numpy as np
import pickle
from math import pi, sin, log
deg = 180/pi

from cascade.utils import get_loc

logmode = False

def openfl(fn):
    f = open(fn, 'rb')
    data = pickle.load(f)
    f.close()
    return data

#open the data
#data1=openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probs_special_0.0_0.0_0.0_0.0.dat")
#data2=openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probs_special_nosys_0.0_0.0_0.0_0.0.dat")

data1 = openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_special_0.0_0.0_0.0_0.0.dat")
data2=  openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/cummulative_probs_special_nosys_0.0_0.0_0.0_0.0.dat")
theta24s = data1["theta24s"]
theta34s = data1["theta34s"]
msqs = data1["msqs"]
chi1 = np.array(data1["chi2s"])
chi2 = np.array(data2["chi2s"])


# relevant 90% chi squared for three DOF 
chis_l = [6.251]
labels = ["90%"]
assert(len(chis_l) == len(labels))

def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((20,10),(30,5) )
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc)

def s2(theta):
    si = np.sin(2*theta)
    return si*si

def get_slice(eV:float, chis:np.ndarray):
    """
    Get a slice of the chi-squared arrays at the requested mass squared difference 
    """
    which_sliver = get_loc(eV, msqs)[0]

    new_chis = np.zeros(shape=(len(theta24s), len(theta34s)))
    for i24 in range(len(theta24s)):
        for i34 in range(len(theta34s)):
            new_chis[i24][i34] = chis[i24][i34][which_sliver]
    return new_chis

at_ev = 5.0
which_sliver = get_loc(at_ev, msqs)[0]

slice1 = get_slice(at_ev, chi1)
slice2 = get_slice(at_ev, chi2)

if logmode:
    c1 = plt.contour(s2(theta24s), s2(theta34s), slice1.transpose(), levels=chis_l, cmap='summer')
    c2 = plt.contour(s2(theta24s), s2(theta34s), slice2.transpose(), levels=chis_l, cmap='autumn')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([1e-3,1])
    plt.xlim([1e-3,1])
    plt.xlabel(r"$\sin^{2}2\theta_{24}$",size=14)
    plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)
else:
    c1 = plt.contour(theta24s*deg, theta34s*deg, slice1.transpose(), levels=chis_l, cmap='summer')
    c2 = plt.contour(theta24s*deg, theta34s*deg, slice2.transpose(), levels=chis_l, cmap='autumn')
    plt.ylim([0,30])
    plt.xlim([0,30])
    plt.xlabel(r"$\theta_{24}$ [deg]",size=14)
    plt.ylabel(r"$\theta_{34}$ [deg]",size=14)


#set_lbls(c1)
#set_lbls(c2)
plt.title("At {:.2f} eV".format(msqs[which_sliver])+r"$^{2}$")
c1.collections[0].set_label("Yes Sys")
c2.collections[0].set_label("No Sys")
plt.legend()
plt.show()

