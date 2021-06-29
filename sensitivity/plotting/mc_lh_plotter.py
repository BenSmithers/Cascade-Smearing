from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi, sin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt



f = open("cummulative_probs.dat_from_mc",'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
raw_sorted = obj["raw_sorted"]
c_prob = obj["c_prob"]
chi2 = obj["chi2s"]


ps = [.10]
chis = [-2*log(p) for p in ps]

labels = ["90%"]#, "99%"]
#labels = ["95%", "68%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10)

def si2(values):
    sinsin = np.sin(2*values)
    return sinsin*sinsin

print(chis)

print(msqs)

brazil_green = (0.0,156./255,59./255)
brazil_yellow = (1.0, 223./255, 0.0)
white = (1,1,1)

cmap = matplotlib.colors.ListedColormap([brazil_green, brazil_yellow, white])
cmap.set_over('1.0')
cmap.set_under('0.75')
bounds = chis
norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)


chis = np.zeros(shape=(len(theta24s), len(msqs)))
for t24 in range(len(theta24s)):
    for msq in range(len(msqs)):
        chis[t24][msq] = chi2[t24][0][msq]
plt.pcolormesh(si2(theta24s), msqs, chis.transpose(), vmin=0, vmax=10, cmap=cmap)
cbar = plt.colorbar(extend='max')
cbar.set_label(r"-2$\Delta$LLH", size=14)
ct = plt.contour(si2(theta24s), msqs, chis.transpose(), levels=[-2*log(0.10)], cmap='cool')
plt.xlim([0.001, 1])
plt.ylim([0.01, 100])
set_lbls(ct)
plt.xscale('log')
plt.yscale('log')
plt.title(r"90% CL Sensitivity from Tracks", size=14)
#plt.text(2,25, "Smithers Preliminary", color="r",size=14)
plt.ylabel(r"$\Delta m_{14}^{2}$ [eV$^{2}$]",size=14)
plt.xlabel(r"$\sin^{2}(2\theta_{24})$",size=14)
plt.tight_layout()
plt.show()



