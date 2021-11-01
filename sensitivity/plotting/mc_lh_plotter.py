from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi, sin
import numpy as np
import os

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))


#f = open("../cummulative_probs.dat_from_mc",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1609e0/cummulative_probs_from_mc_0.0_0.1609e0_0.0_4.4700e0.dat",'rb')
f = open("/home/bsmithers/software/data/hg_sib/expectations/0.0/cummulative_probs_from_mc_0.0_0.0_0.0_0.0.dat",'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
raw_sorted = obj["raw_sorted"]
c_prob = obj["c_prob"]
chi2 = obj["chi2s"]


ps = [ 0.10]
chis = [-2*log(p) for p in ps]

labels = ["90%"]#, "99%"]
#labels = ["95%", "68%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10, inline_spacing=-10)

def si2(values):
    sinsin = np.sin(2*values)
    return sinsin*sinsin

print(chis)


brazil_green = (0.0,156./255,59./255)
brazil_yellow = (1.0, 223./255, 0.0)
white = (1,1,1,0.0)

all_colors = [brazil_green for i in range(23)] + [brazil_yellow for i in range(46-23)] + [white for i in range(100-46)]

cmap = matplotlib.colors.ListedColormap(all_colors)
cmap.set_over('1.0')
cmap.set_under('0.75')
bounds = chis  
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)


chis = np.zeros(shape=(len(theta24s), len(msqs)))
for t24 in range(len(theta24s)):
    for msq in range(len(msqs)):
        chis[t24][msq] = chi2[t24][0][msq]
#plt.pcolormesh(si2(theta24s), msqs, chis.transpose(), vmin=0, vmax=10, cmap="PuBu")
#cbar = plt.colorbar(extend='max')
#cbar.set_label(r"-2$\Delta$LLH", size=14)
ct = plt.contour(si2(theta24s), msqs, chis.transpose(), levels=[-2*log(0.10)], cmap='Oranges_r')
plt.xlim([0.001, 1])
plt.ylim([0.01, 100])
#set_lbls(ct)
plt.xscale('log')
plt.yscale('log')
plt.text(2e-3, 3e1, "90\% CL, 2 DOF", color='black', size='x-large')
# plt.title(r"90% CL Sensitivity from Tracks", size=14)
#plt.text(2,25, "Smithers Preliminary", color="r",size=14)
plt.ylabel(r"$\Delta m_{14}^{2}$ [eV$^{2}$]",size=18)
plt.xlabel(r"$\sin^{2}(2\theta_{24})$",size=18)
plt.tight_layout()
plt.savefig("track_sensitivity.png", dpi=400)
plt.show()



