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

asimov = False


i_c = 0
def add_contour(filename, label, linestyle):
    global i_c
    new_cmap = shift_cmap("magma", (i_c*80.+50)/256)

    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    this_chi = np.array(data["chi2s"])
    these_chi = np.zeros(shape=(len(theta14s), len(theta34s)))

    
    def interpval(i, j):
        x = scale_x[i][i]
        y = scale_y[j][j]
        slope = (0.00698232-0.00817524)/(0.184218-0.162226)
        b = 0.00817524-slope*0.162226
        if y> slope*x+b:
            return 1e8
        else:
            return 1.0*chis_l[0]

    def other_int(i,j):
        x = scale_x[i][i]
        y = scale_y[j][j]
        p1 = (0.154637, 0.00721796)
        p2 = (0.163703, 0.00811777)
        slope = (p2[1]-p1[1])/(p2[0]-p1[0])
        b = p2[1]-slope*p2[0]
        if y>(slope*x+b):
            return 1e8
        else:
            return 1.0*chis_l[0]

    assert(len(this_chi[0][0]) == 1)
    for t24 in range(len(theta14s)):
        for t34 in range(len(theta34s)):
            use = this_chi[t24][t34][0]
            condition = this_chi[t24][t34][0] > 1e8 
            if condition:
                use= None # other_int(t24, t34)
            these_chi[t24][t34] = use

    ct = plt.contour(scale_x, scale_y, these_chi.transpose(), levels=chis_l, cmap=new_cmap,linestyles=linestyle)
    
    ct.collections[0].set_label(label)
    i_c += 1
    return ct



#ct = plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='cool',linestyles='-')
#ct2= plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='summer', linestyles='--')
small = False
if asimov:
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_1_00eV_smearing_0.0_0.1609e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.1609$", '-')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_3_30eV_smearing_0.0_0.1609e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.1609$", '-')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_4_64eV_smearing_0.0_0.1609e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.1609$", '-')
    i_c = 0

    ct = add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_1_00eV_smearing_0.0_0.3826e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.3826$", '--')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_3_30eV_smearing_0.0_0.3826e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.3826$", '--')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_4_64eV_smearing_0.0_0.3826e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.3826$", '--')
    new_color = list(ct.collections[0].get_color()[0])
    plt.plot([0.189364, 0.213617], [0.00808428, 0.00686016], ls='--', color=new_color)
else:
    #plt.vlines(x=[0.3,0.55], ymin=0, ymax=1.0, colors='black')
    plt.fill_between([0.3,0.55],0, 1.0,color=(0,0,0,0.1))
    plt.text(0.625, 0.09,"90\% CL, 2 DOF", fontsize='large')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_1_00eV_smearing_0.3555e0_0.1609e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.1609$", '-')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_3_30eV_smearing_0.3555e0_0.1609e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.1609$", '-')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_4_64eV_smearing_0.3555e0_0.1609e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.1609$", '-')
    i_c = 0

    ct = add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_1_00eV_smearing_0.3555e0_0.3826e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.3826$", '--')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_3_30eV_smearing_0.3555e0_0.3826e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.3826$", '--')
    add_contour("/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_4_64eV_smearing_0.3555e0_0.3826e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.3826$", '--')
    new_colr = list(ct.collections[0].get_color()[0])
    #plt.plot([0.154637, 0.161903], [0.00711, 0.00814], ls='--', color=(new_colr[0], new_colr[1], new_colr[2]))
"""
/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_1_00eV_smearing_0.3555e0_0.1609e0_0.0_1.0000e0.dat
/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_3_30eV_0.3555e0_0.3826e0_0.0_3.3000e0.dat
/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_1_00eV_0.3555e0_0.3826e0_0.0_1.0000e0.dat
/home/benito/software/data/cascade/hg_sib/0.3826e0/best_llh_4_64eV_0.3555e0_0.3826e0_0.0_4.6400e0.dat
/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_1_00eV_0.3555e0_0.1609e0_0.0_1.0000e0.dat
/home/benito/software/data/cascade/hg_sib/0.1609e0/best_llh_3_30eV_0.3555e0_0.1609e0_0.0_3.3000e0.dat
"""
if asimov:
    plt.xlim([0,0.5])
else:
    plt.xlim([0,1.0])
plt.ylim([0,0.1])
#plt.yscale('log')
#plt.ylim([1e-2, 1e0])

#plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{14}$",size=14)
#plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)
plt.xlabel(r"$\sin^{2}2\theta_{14}$",size=14)

print(np.nanmin(chis))
mindices = np.nanargmin(chis)
x,y = np.unravel_index(mindices, shape=np.shape(chis))
print(x,y)
# plt.scatter(np.sin(theta14s[x])**2, (np.cos(theta14s[x])**2)*(np.sin(theta34s[y])**2), s=80, marker=(5,1))
print(theta14s[x], theta34s[y])

#ct.collections[0].set_label(r"3.3eV$^{2}$")
#ct2.collections[0].set_label(r"4.5eV$^{2}$")

#plt.vlines(x=scaled_minos_24, ymin=0, ymax=scaled_minos_34[-1], color='k')
if asimov:
    plt.legend(loc='upper right', fancybox=True, frameon=True, framealpha=1.0,facecolor='white', edgecolor='black')
else:
    plt.legend(loc='upper left', fancybox=True, frameon=True,framealpha=1.0,facecolor='white', edgecolor='black')

plt.savefig("bestPlot_{}.png".format("asimov" if asimov else "not_asimov"), dpi=400)
plt.show()
