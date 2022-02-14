"""
This puts the 90% contours for two different sensitivity results on the same plot
"""

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt 

plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

import numpy as np
import pickle
from math import pi, sin, log, asin, sqrt
deg = 180/pi

from cascade.utils import get_loc, get_closest

logmode = True

def openfl(fn):
    f = open(fn, 'rb')
    data = pickle.load(f)
    f.close()
    return data

#open the data
#data1=openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probs_special_0.0_0.0_0.0_0.0.dat")
#data2=openfl("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probs_special_nosys_0.0_0.0_0.0_0.0.dat")

data1 = openfl("/home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_0.0_0.0_0.0_0.0.dat")
data2 = openfl("/home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_nosys_0.0_0.0_0.0_0.0.dat")
theta24s = data1["theta24s"]
theta34s = data1["theta34s"]
print(len(theta24s))
print(len(theta34s))
msqs = data1["msqs"]
chi1 = np.array(data1["chi2s"])
chi2 = np.array(data2["chi2s"])


th24_cut = asin(sqrt(0.03))
theta24s_cut = theta24s[:get_loc(th24_cut, theta24s)[0]]
chis1_cut = chi1[:get_loc(th24_cut, theta24s)[0]]
chis2_cut = chi2[:get_loc(th24_cut, theta24s)[0]]
xs_cut, ys_cut = np.meshgrid(theta24s_cut, theta34s)
scale_x_cut = np.sin(xs_cut)**2
scale_y_cut = (np.cos(xs_cut)**2)*(np.sin(ys_cut)**2)


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

interpolate = True
def get_slice(eV:float, chis:np.ndarray):
    """
    Get a slice of the chi-squared arrays at the requested mass squared difference 
    """
    which_sliver = get_loc(eV, msqs)
    print("slivers: {}".format(which_sliver))
    if interpolate:
        chis_final = np.zeros(shape=(2, len(chis), len(chis[0])))
        for t24 in range(len(chis)):
            for t34 in range(len(chis[t24])):
                #print("{} {}".format(t24, t34))
                chis_final[0][t24][t34] = chis[t24][t34][which_sliver[0]]
                chis_final[1][t24][t34] = chis[t24][t34][which_sliver[1]]
        return get_closest( eV, [msqs[which_sliver[0]], msqs[which_sliver[1]]], chis_final)
    else:
        which_sliver = which_sliver[0]
        new_chis = np.zeros(shape=(len(theta24s), len(theta34s)))
        for i24 in range(len(theta24s)):
            for i34 in range(len(theta34s)):
                new_chis[i24][i34] = chis[i24][i34][which_sliver]
        return new_chis

at_ev = 1.0
which_sliver = get_loc(at_ev, msqs)[0]

slice1 = get_slice(at_ev, chi1)
slice2 = get_slice(at_ev, chi2)

slice1_cut = get_slice(at_ev, chis1_cut)
slice2_cut = get_slice(at_ev, chis2_cut)

if logmode:
    c1 = plt.contour(scale_x_cut, scale_y_cut, slice1_cut.transpose(), levels=chis_l, cmap='summer')
    c2 = plt.contour(scale_x_cut, scale_y_cut, slice2_cut.transpose(), levels=chis_l, cmap='autumn')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim([3e-3,0.5])
    plt.xlim([1e-3,0.5])
    plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{24}$",size=14)
    plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
else:
    c1 = plt.contour(theta24s*deg, theta34s*deg, slice1.transpose(), levels=chis_l, cmap='summer')
    c2 = plt.contour(theta24s*deg, theta34s*deg, slice2.transpose(), levels=chis_l, cmap='autumn')
    plt.ylim([0,30])
    plt.xlim([0,30])
    plt.xlabel(r"$\theta_{24}$ [deg]",size=14)
    plt.ylabel(r"$\theta_{34}$ [deg]",size=14)


#set_lbls(c1)
#set_lbls(c2)
if interpolate:
    plt.title("At {:.2f} eV".format(at_ev)+r"$^{2}$")
else:
    plt.title("At {:.2f} eV".format(msqs[which_sliver])+r"$^{2}$")
c1.collections[0].set_label("Yes Sys")
c2.collections[0].set_label("No Sys")
plt.legend()
plt.show()

