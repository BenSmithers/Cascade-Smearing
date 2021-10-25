"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from os import set_inheritable
from cascade.utils import gen_filename, config
from cascade.utils import get_loc, get_closest, get_color

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
# f = open("/home/benito/software/data/cascade/hg_sib/0.1652e0/joint_likelihood_0.0_0.1652e0_0.2293e0_4.6416e0.dat", 'rb')
#
#f = open("/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_0.0_0.0_0.0_0.0.dat", 'rb')


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
#f = open("/home/benito/software/data/cascade/hg_sib/expectations/0.0/cummulative_probs_nosys_0.0_0.0_0.0_0.0.dat",'rb')\

# use this to just get the edges and stuff 
f = open("/home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_0.0_0.0_0.0_0.0.dat",'rb')
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

# evs = [1.5,3.3,4.64,10]
evs = [1.0]


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

ev = 1.0

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

def add_contour(filename, ls, cmap='cool', chop=False):
    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    this_chi = np.array(data["chi2s"]) 
    temp = np.zeros(shape=(2, len(theta24s), len(theta34s)))
    for t24 in range(len(theta24s)):
        for t34 in range(len(theta34s)):
            temp[0][t24][t34] = this_chi[t24][t34][which_sliver[0]]
            temp[1][t24][t34] = this_chi[t24][t34][which_sliver[1]]
    final_chi = get_closest( ev, [msqs[which_sliver[0]], msqs[which_sliver[1]]], temp)

    if chop:
        # need to chop off the right side :frown: 
        for i_x in range(len(scale_x)):
            for i_y in range(len(scale_y[i_x])):
                if scale_x[i_x][i_x]>6e-2:
                    final_chi[i_x][i_y]=None
    
    ct = plt.contour(scale_x, scale_y, final_chi.transpose(), levels=chis_l, cmap=cmap,linestyles=ls)   
    return ct

plt.xscale('log')
#plt.title("At {:.2f} eV".format(ev if interp else msqs[which_sliver[0]]) +r"$^{2}$")
#            plt.yscale('log')
plt.ylim([0.00,0.30])
plt.xlim([1e-3,0.5])
plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{24}$",size=14)
plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")
#            plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)

#plt.plot(t24s, scaled_minos_34, color='k', label="Minos")
plt.plot(super_k[0], super_k[1], color=get_color(1,4,cmap="magma"), label="Super K")
plt.plot(deepcore[0], deepcore[1], color=get_color(2,4,"magma"), label="DeepCore")
#plt.plot(antares[0], antares[1], color=get_color(3,4,"magma"), label="KM3NeT/ORCA NO")


ct = add_contour("/home/benito/software/data/cascade/hg_sib/0.0/joint_likelihood_0.0_0.0_0.0_0.0.dat", '--','cool',True)
ct.collections[0].set_label("Joint Fit")
ct = add_contour("/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_float_0.0_0.0_0.0_0.0.dat", '-','cool',False)
ct.collections[0].set_label("Cascades")


#set_lbls(ct, [str(ev)+r"eV$^{2}$"], get_color(color_count, len(evs), "viridis"))
#ct.collections[0].set_color(get_color(color_count, len(evs), "viridis"))

#plt.vlines(x=scaled_minos_24, ymin=0, ymax=scaled_minos_34[-1], color='k')
plt.yscale('log')
plt.ylim([3e-3, 5e-1])
plt.legend(loc='upper right')
root_name = "superPlot"
plt.savefig("{}_{:.2f}.png".format(root_name, ev if interp else which_sliver[0]),dpi=400)
plt.show()
