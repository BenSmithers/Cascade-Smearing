"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from cascade.utils import gen_filename, config
from cascade.utils import get_loc, get_closest, get_color
from cascade.utils import shift_cmap

import pickle
from math import sqrt, asin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

interp = True

#f = open("../cummulative_probs.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_nosys_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_0.0_0.0_0.0_0.0.dat",'rb')

#/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_float_0.0_0.0_0.0_0.0.dat
# /home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_nosys_0.0_0.0_0.0_0.0.dat <- no sys
#f = open("/home/benito/software/data/cascade/hg_sib/0.1652e0/joint_likelihood_0.0_0.1652e0_0.2293e0_4.6416e0.dat", 'rb')


f2 = open("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_4_47eV_0.0_0.1609e0_0.0_4.4700e0.dat",'rb')
f = open("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_3_30eV_0.0_0.1609e0_0.0_3.3000e0.dat", 'rb')
# f = open("/home/benito/software/data/cascade/hg_sib/0.1609e0/joint_likelihood_0.0_0.1609e0_0.0_4.4700e0.dat", 'rb')

obj = pickle.load(f)
f.close()

theta14s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
chi2 = np.array(obj["chi2s"]) 

obj2 = pickle.load(f2)
chi22 = np.array(obj2["chi2s"])
f2.close()

deg = 180./3.1415926



ps = [0.10] #, 0.01]
chis_l = [4.605] # 2 degrees of freedom!! 
#chis_l = [6.251] # 3 degrees of freedom!! 

labels = ["90%"]#, "99%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((20,10),(30,5) )
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc)

evs = [1.0, 4.641, 10.0]
#evs = [1.1]


xs, ys= np.meshgrid( theta14s, theta34s)

from math import pi 
scale_x = np.sin(2*xs)**2
scale_y = np.sin(2*ys)**2
#(np.cos(xs)**2)*

chis = np.zeros(shape=(len(theta14s), len(theta34s)))
chis2= np.zeros(shape=(len(theta14s), len(theta34s)))

assert(len(chi2[0][0]) == 1)
for t24 in range(len(theta14s)):
    for t34 in range(len(theta34s)):
        chis[t24][t34] = chi2[t24][t34][0]
        chis2[t24][t34]=chi22[t24][t34][0]

#plt.pcolormesh(scale_x,scale_y, chis.transpose(),vmin=0, vmax=10, cmap="PuBu")
#cbar = plt.colorbar(extend='max')

#plt.xscale('log')
#plt.yscale('log')

plt.xlim([0,0.5])
plt.ylim([0,0.5])

#plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{14}$",size=14)
#plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)
plt.xlabel(r"$\sin^{2}2\theta_{14}$",size=14)


i_c = 0
def add_contour(filename, label):
    global i_c
    new_cmap = shift_cmap("cool", i_c*50./256)

    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    this_chi = np.array(data["chi2s"])
    these_chi = np.zeros(shape=(len(theta14s), len(theta34s)))

    assert(len(this_chi[0][0]) == 1)
    for t24 in range(len(theta14s)):
        for t34 in range(len(theta34s)):
            these_chi[t24][t34] = this_chi[t24][t34][0]

    ct = plt.contour(scale_x, scale_y, these_chi.transpose(), levels=chis_l, cmap=new_cmap,linestyles='--')
    ct.collections[0].set_label(label)
    i_c += 1



#ct = plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='cool',linestyles='-')
#ct2= plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='summer', linestyles='--')
plt.title("When "+r"$\theta_{24}=0.1609$")
add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_1_00eV_0.0_0.1609e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$")
add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_3_30eV_0.0_0.1609e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$")
add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_4_64eV_0.0_0.1609e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$")

print(np.nanmin(chis))
mindices = np.nanargmin(chis)
x,y = np.unravel_index(mindices, shape=np.shape(chis))
print(x,y)
# plt.scatter(np.sin(theta14s[x])**2, (np.cos(theta14s[x])**2)*(np.sin(theta34s[y])**2), s=80, marker=(5,1))
print(theta14s[x], theta34s[y])

#ct.collections[0].set_label(r"3.3eV$^{2}$")
#ct2.collections[0].set_label(r"4.5eV$^{2}$")

#plt.vlines(x=scaled_minos_24, ymin=0, ymax=scaled_minos_34[-1], color='k')
plt.legend(loc='upper right')

plt.savefig("bestPlot_1609.png", dpi=400)
plt.show()
