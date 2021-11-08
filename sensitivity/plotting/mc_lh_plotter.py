from cascade.utils import gen_filename, config
from cascade.utils import get_loc, get_closest

from cascade.utils import get_color

import pickle
from math import log, pi, sin, log10
import numpy as np
import os

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))


#f = open("../cummulative_probs.dat_from_mc",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib//expectations/0.1609e0/cummulative_probs_from_mc_0.0_0.1609e0_0.0_4.4700e0.dat",'rb')
f2 = open("/home/bsmithers/software/data/hg_sib/expectations/0.0/cummulative_probs_from_mc_0.0_0.0_0.0_0.0.dat",'rb')
f = open ("/home/bsmithers/software/data/hg_sib/0.0/newSense_mcfudgeresult_0.0_0.0_0.0_0.0.dat",'rb')
obj = pickle.load(f)
f.close()

eight_yr = pickle.load(f2)
f2.close()

yellow_left = np.transpose(np.loadtxt("brazil/yellow_left",delimiter=","))
yellow_right = np.transpose(np.loadtxt("brazil/yellow_right",delimiter=","))
green_left = np.transpose(np.loadtxt("brazil/green_left",delimiter=","))
green_right = np.transpose(np.loadtxt("brazil/green_right",delimiter=","))
print("min? {}".format(yellow_left[1][0]))
meow = np.transpose(np.loadtxt("meows_sense.dat", delimiter=","))

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
print(msqs)
chi2 = obj["chi2s"]

t24_8 = eight_yr["theta24s"]
msqs_8 = eight_yr["msqs"]
chi_8 = eight_yr["chi2s"]


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

yellow_left_interp = []
yellow_right_interp = []
green_left_interp = []
green_right_interp = []

for msq in msqs:
    if msq<yellow_left[1][0]:
        left = yellow_left[0][0]
    elif msq>yellow_left[1][-1]:
        left = yellow_left[0][-1]
    else:
        left = get_closest(msq, yellow_left[1], yellow_left[0])
    yellow_left_interp.append(left)


    if msq<yellow_right[1][0]:
        right = yellow_right[0][0]
    elif msq>yellow_right[1][-1]:
        left = yellow_right[0][-1]
    else:
        right = get_closest(msq, yellow_right[1], yellow_right[0])
    yellow_right_interp.append(right)


    if msq<green_left[1][0]:
        left = green_left[0][0]
    elif msq>green_left[1][-1]:
        left = green_left[0][-1]
    else:
        left = get_closest(msq, green_left[1], green_left[0])
    green_left_interp.append(left)

    if msq<green_right[1][0]:
        right = green_right[0][0]
    elif msq>green_right[1][-1]:
        left = green_right[0][-1]
    else:
        right = get_closest(msq, green_right[1], green_right[0])
    green_right_interp.append(right)
    

    
#yellow_left_interp = np.array([ get_closest(msq,yellow_left[1], yellow_left[0]) for msq in msqs ])
#yellow_right_interp = np.array([ get_closest(msq,yellow_right[1], yellow_right[0]) for msq in msqs ])

yellow_right_interp = np.array(yellow_right_interp)
yellow_left_interp = np.array(yellow_left_interp)

chis = np.zeros(shape=(len(theta24s), len(msqs)))
chis8 = np.zeros(shape=(len(t24_8), len(msqs_8)))
for t24 in range(len(theta24s)):
    for msq in range(len(msqs)):
        chis[t24][msq] = chi2[t24][0][msq]
for t24 in range(len(t24_8)):
    for msq in range(len(msqs_8)):
        chis8[t24][msq] = chi_8[t24][0][msq]

#plt.pcolormesh(si2(theta24s), msqs, chis.transpose(), vmin=0, vmax=10, cmap="PuBu")
#cbar = plt.colorbar(extend='max')
#cbar.set_label(r"-2$\Delta$LLH", size=14)

title_1 = r"\textbf{IC86 2016 1yr}"
title_2 = r"\textbf{This Work}"
all_labels=[title_1]

plt.gca().add_line(plt.Line2D([], [], color="none", label=title_1)) 

plt.fill_betweenx(msqs,yellow_right_interp, yellow_left_interp, color='yellow', label="95\%")
plt.fill_betweenx(msqs,green_right_interp, green_left_interp, color='green', label="68\%")
plt.plot(meow[0], meow[1], color='gainsboro' ,ls='--', label="Median")
all_labels += ["95\%", "68\%", "Median"]
#plt.plot(yellow_left_interp, msqs)
#plt.plot(yellow_right_interp, msqs)

plt.gca().add_line(plt.Line2D([], [], color="none", label=title_2)) 
ct = plt.contour(si2(theta24s), msqs, chis.transpose(), levels=[-2*log(0.10)], cmap='Oranges_r')
ct8 = plt.contour(si2(t24_8), msqs_8, chis8.transpose(), levels=[-2*log(0.10)], cmap='Oranges_r')
all_labels += [title_2, "1 year", "8 years"]
ct.collections[0].set_label("1 year")
ct.collections[0].set_color(get_color(2,4,'Oranges'))
ct8.collections[0].set_label("8 years")
ct8.collections[0].set_color(get_color(4,4,'Oranges'))
plt.xlim([0.001, 1])
plt.ylim([0.01, 10])
#set_lbls(ct)
plt.xscale('log')
plt.yscale('log')
# plt.title(r"90% CL Sensitivity from Tracks", size=14)
#plt.text(2,25, "Smithers Preliminary", color="r",size=14)
plt.ylabel(r"$\Delta m_{14}^{2}$ [eV$^{2}$]",size=18)
plt.xlabel(r"$\sin^{2}(2\theta_{24})$",size=18)


def reorderLegend(ax=None, order=None):
    handles, labels = ax.get_legend_handles_labels()
    info = dict(zip(labels, handles))

    new_handles = [info[l] for l in order]
    return new_handles, order
handles, labels = reorderLegend(ax=plt.gca(), order=all_labels)


plt.legend(handles=handles, labels=labels, loc='lower left', prop={'size':12}) #, fancybox=True, frameon=True, framealpha=0.25,facecolor='gray', edgecolor='black')
plt.tight_layout()
plt.savefig("track_sensitivity.png", dpi=400)
plt.show()



