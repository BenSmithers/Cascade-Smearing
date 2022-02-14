"""
MANY PLOT

Draw our sensitivity contours at multiple values of \Delta m_{41}^{2}
"""

from os import set_inheritable
import os
from cascade.utils import gen_filename, config
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
root_dir = config["datapath"]
f = open(root_dir + "/0.0/joint_likelihood_smearing_0.0_0.0_0.0_0.0.dat",'rb')
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

xs, ys= np.meshgrid( theta24s, theta34s)

scale_x = np.sin(xs)**2
scale_y = (np.cos(xs)**2)*(np.sin(ys)**2)


def add_contour(filename, ls, cmap='cool', chop=False):
    f = open(filename, 'rb')
    data = pickle.load(f)
    f.close()
    this_chi = np.array(data["chi2s"]) 
    evs = [0.5, 1.0, 3.0, 10.0]

    count = 0
    for ev in evs:

        which_sliver = get_loc(ev, msqs)

        
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
                    if final_chi[i_x][i_y] < 0 or final_chi[i_x][i_y]>1000:
                        final_chi[i_x][i_y] = None
                    if scale_x[i_x][i_x]>6e-2:
                        final_chi[i_x][i_y]=None
                    
        
        ct = plt.contour(scale_x, scale_y, final_chi.transpose(), levels=chis_l, cmap=cmap,linestyles=ls) 
        ct.collections[0].set_color(get_color(count+1, len(evs)+1, "magma") )
        suffix = "Joint" if chop else "Cascades"
        ct.collections[0].set_label("{}, {:.1f}".format(suffix, ev) + r"eV$^{2}$")
        count += 1 
    return ct

plt.xscale('log')
plt.ylim([0.00,0.30])
plt.xlim([1e-3,0.7])
plt.xlabel(r"$\left|U_{\mu 4}\right|^{2}=\sin^{2}\theta_{24}$")
plt.ylabel(r"$\left| U_{\tau 4}\right|^{2}=\sin^{2}\theta_{34}\cdot\cos^{2}\theta_{24}$")

ct = add_contour(root_dir + "/0.0/joint_likelihood_0.0_0.0_0.0_0.0.dat", '--','cool',True)
ct = add_contour(root_dir + "/0.0/newSense_result_float_0.0_0.0_0.0_0.0.dat", '-','cool',False)


plt.yscale('log')
plt.ylim([3e-3, 5e-1])
plt.legend(loc='upper right')
plt.tight_layout()
plt.savefig("manyPlot.png",dpi=400)
plt.show()
