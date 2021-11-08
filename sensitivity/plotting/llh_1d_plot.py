"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from cascade.sensitivity.newSense import Scanner, JointLLH, doLLH
from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import get_loc, get_closest, get_color
from cascade.utils import shift_cmap

import pickle
from math import sqrt, asin
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))


"""
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_1_00eV_smearing_0.0_0.1609e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_3_30eV_smearing_0.0_0.1609e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_4_64eV_smearing_0.0_0.1609e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_1_00eV_smearing_0.0_0.3826e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.3826$", '--')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_3_30eV_smearing_0.0_0.3826e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.3826$", '--')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_4_64eV_smearing_0.0_0.3826e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.3826$", '--')
"""


msqs = [1.0, 3.3, 4.64]
th24s = [0.1609, 3.826]
th34s = [ 0.0 ]
thetas = np.concatenate(([0], np.arcsin(np.sqrt(np.logspace(-3,0,90)))/2))

filename = gen_filename(config["datapath"], "newllh_1d_scan.dat", params=SterileParams())
f = open(filename, 'rb')
data_dict = pickle.load(f)
f.close()

msqs = data_dict["msqs"]
th14s = data_dict["theta14s"]
th24s = data_dict["theta24s"]
th34s = data_dict["theta34s"]
chis = data_dict["chi2s"]
# chis = chis - np.nanmin(chis)

print("Chi Shape: {}".format(np.shape(chis)))

cl = 4 * 1.9
color_number = 0
for i24 in range(len(th24s)):
    if i24==2:
        continue
    for i34 in range(len(th34s)):
        if i34==2:
            continue
        if i24 == 0 :
            if i34==0:
                ls = 'solid'
            elif i34==1:
                ls = 'dotted'
            else:
                ls = ''
        elif i24== 1:
            if i34==0:
                ls = '-.'
            elif i34==1:
                ls = (0,(5,10))
            else:
                ls = 'dashdot'
        for ims in range(len(msqs)):

            color = get_color(ims, 4, "magma")
            llh_cut = [ chis[i14][i24][i34][ims] for i14 in range(len(th14s)) ]
            label = r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[i24])
            label+= r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[i34])
            label+= r"$\Delta m_{41}^{2}$: "+"{:.2f}".format(msqs[ims])

            plt.plot(np.sin(th14s)**2, llh_cut, ls=ls, color=color)


title_1 = r"\textbf{Mixing Angles: }"
title_2 = r"\textbf{Mass-Squared: }"
#title_2 = r"$ \Delta m_{41}^{2}:$"
all_labels=[title_1]
plt.gca().add_line(plt.Line2D([], [], color="none", label=title_1)) 


plt.plot([],[],color='black', ls='solid',label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[0]))
plt.plot([],[],color='black', ls='dotted',label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[1]))
#plt.plot([],[],color='black', ls='dashed',label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[2]))
plt.plot([],[],color='black', ls='-.',label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[0]))
plt.plot([],[],color='black', ls=(0,(5,10)),label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[1]))
#plt.plot([],[],color='black', ls='dashdot',label=r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[2]))
all_labels += [r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[0]), 
               r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[1]),
               #r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[0])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[2]),
               r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[0]), 
               r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[1])]
               #r"$\theta_{24}$: "+"{:.2f}, ".format(th24s[1])+ r"$\theta_{34}$: "+"{:.2f}, ".format(th34s[2])]

plt.gca().add_line(plt.Line2D([], [], color="none", label=title_2)) 

plt.plot([],[],color=get_color(0,4,"magma"), ls="solid",label=r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[0]))
plt.plot([],[],color=get_color(1,4,"magma"), ls="solid",label=r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[1]))
plt.plot([],[],color=get_color(2,4,"magma"), ls="solid",label=r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[2]))
all_labels += [r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[0]),
               r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[1]),
               r"$\Delta m_{41}^{2}$: "+"{:.2f}, ".format(msqs[2])]


def reorderLegend(ax=None, order=None):
    handles, labels = ax.get_legend_handles_labels()
    info = dict(zip(labels, handles))

    new_handles = [info[l] for l in order]
    return new_handles, order
handles, labels = reorderLegend(ax=plt.gca(), order=all_labels)

# plt.text(0.2, cl+0.5, "90\% CL, 4 DOF", color='red')
plt.hlines(cl, 0, 0.5, color='red')
plt.legend( bbox_to_anchor=(1.0, 1.1))
plt.xlabel(r"$\sin^{2}(\theta_{14})$")
plt.ylabel(r"-2(LLH$_{null}$ - LLH$_{injected})$")
plt.xlim([0, 0.5])
plt.ylim([1e0,1e3])
plt.yscale('log')
plt.tight_layout()
plt.savefig("legendary.png",dpi=400)
plt.show()

