"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from os import set_inheritable
import os 
from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import get_loc, get_closest, get_color

import pickle
from math import sqrt, asin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))

interp = True

#f = open("../cummulative_probs.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_nosys_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.0/scaled_cummulative_probsnonorm_special_0.0_0.0_0.0_0.0.dat",'rb')

#/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_float_0.0_0.0_0.0_0.0.dat
# /home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_nosys_0.0_0.0_0.0_0.0.dat <- no sys
fname = gen_filename(config["datapath"], "joint_likelihood_smearing.dat" , SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.6416))

f = open(fname, 'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
chi2 = np.array(obj["chi2s"]) 

deg = 180./3.1415926



ps = [0.10] #, 0.01]
#chis_l = [4.605] # 2 degrees of freedom!! 
chis_l = [6.251] # 3 degrees of freedom!! 

#labels = ["90%"]#, "99%"]
move_it = 0
def set_lbls(ct_plot, labels, color):
    global move_it
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((1e-2 + (move_it*0.8e-2), 1e-1 - move_it*(1.7e-2)), )
    thingy = plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc, colors=(color,))
    move_it += 1

evs = [1.5,3.3,4.64,10]
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
antares = np.transpose(np.loadtxt("antares.dat", delimiter=","))

color_count = 0
for ev in evs:
    color_count+=1

    which_sliver = get_loc(ev, msqs)
    print(which_sliver)
    if interp:
        chis_f = np.zeros(shape=(2, len(theta24s), len(theta34s)))
        for t24 in range(len(theta24s)):
            for t34 in range(len(theta34s)):
                chis_f[0][t24][t34] = chi2[t24][t34][which_sliver[0]]
                chis_f[1][t24][t34] = chi2[t24][t34][which_sliver[1]]
        chis = get_closest( ev, [msqs[which_sliver[0]], msqs[which_sliver[1]]], chis_f)


    else:
        sliver = which_sliver[0]
        chis = np.zeros(shape=(len(theta24s), len(theta34s)))
        

        for t24 in range(len(theta24s)):
            for t34 in range(len(theta34s)):
                chis[t24][t34] = chi2[t24][t34][sliver]

    #plt.pcolormesh(scale_x,scale_y, chis.transpose(),vmin=0, vmax=10, cmap="PuBu")
    #cbar = plt.colorbar(extend='max')
    
    #plt.title("At {:.2f} eV".format(ev if interp else msqs[which_sliver[0]]) +r"$^{2}$")
#            plt.yscale('log')

    plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{24}$",size=14)
    plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
#            plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)

    #plt.plot(t24s, scaled_minos_34, color='k', label="Minos")
    #plt.plot(super_k[0], super_k[1], color=get_color(1,4,cmap="magma"), label="Super K")
    #plt.plot(deepcore[0], deepcore[1], color=get_color(2,4,"magma"), label="DeepCore")
    #plt.plot(antares[0], antares[1], color=get_color(3,4,"magma"), label="KM3NeT/ORCA NO")

    ct = plt.contour(scale_x, scale_y, chis.transpose(), levels=chis_l, cmap='cool',linestyles='dashed')   
    print(np.nanmin(chis))
    mindices = np.nanargmin(chis)
    x,y = np.unravel_index(mindices, shape=np.shape(chis))
    print(x,y)
    if ev==4.64:
        plt.scatter(np.sin(theta24s[x])**2, (np.cos(theta24s[x])**2)*(np.sin(theta34s[y])**2), s=80, marker=(5,1), color=get_color(color_count, len(evs)+1, "magma"))
    print(theta24s[x], theta34s[y])

    
    
    #set_lbls(ct, [str(ev)+r"eV$^{2}$"], get_color(color_count, len(evs), "viridis"))
    ct.collections[0].set_color(get_color(color_count, len(evs)+1, "magma"))
    ct.collections[0].set_label("{}".format(ev) + r" eV$^{2}$")

    #plt.vlines(x=scaled_minos_24, ymin=0, ymax=scaled_minos_34[-1], color='k')
    
plt.xlim([5e-3,0.2])
plt.xscale('log')
plt.yscale('log')
plt.ylim([9e-3, 5e-1])
plt.tight_layout()
plt.legend(loc='upper right')
plt.savefig("containment.png",dpi=400)
plt.show()
