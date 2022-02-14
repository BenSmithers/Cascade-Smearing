from cascade.utils import gen_filename, config
from cascade.utils import get_loc

import pickle
from math import log, pi, sin
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")


#f = open("../cummulative_probs.dat",'rb')
f = open("/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_fixed_0.0_0.0_0.0_0.0.dat",'rb')
#f = open("/home/benito/software/data/cascade/hg_sib/0.0/newSense_result_0.0_0.0_0.0_0.0.dat",'rb')
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
#f = open("/home/benito/software/data/cascade/hg_sib/expectations/0.0/cummulative_probs_nosys_0.0_0.0_0.0_0.0.dat",'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]

print(np.min(obj["chi2s"]))
chi2 = np.array(obj["chi2s"]) 
print("Max {} and min {}".format(np.max(chi2), np.min(chi2)))
print(np.argmin(chi2))

central = np.unravel_index(np.argmin(chi2), shape=np.shape(chi2))
print("It's at {}".format(central))
deg = 180./3.1415926
print("Best Fit: {} {} {}".format(deg*theta24s[central[0]], deg*theta34s[central[1]], msqs[central[2]]))


#chi2 = chi2 - np.min(chi2)

#flat = chi2.flatten()
#plt.hist(flat, bins=20)

old = False
log_mode = False


ps = [0.10] #, 0.01]
chis_l = [6.251]#, 11.345]

labels = ["90%"]#, "99%"]
def set_lbls(ct_plot):
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((20,10),(30,5) )
    plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc)

evs = [1.1, 5.0, 10.0]

def s2(theta):
    si = np.sin(2*theta)
    return si*si

for ev in evs:
    which_sliver = get_loc(ev, msqs)[0]
    chis = np.zeros(shape=(len(theta24s), len(theta34s)))
    for t24 in range(len(theta24s)):
        for t34 in range(len(theta34s)):
            chis[t24][t34] = chi2[t24][t34][which_sliver]

    if log_mode:
        plt.pcolormesh(s2(theta24s), s2(theta34s), chis.transpose(), vmin=0, vmax=20, cmap="PuBu")
        cbar = plt.colorbar(extend='max')
        ct = plt.contour(s2(theta24s), s2(theta34s), chis.transpose(), levels=chis_l, cmap='cool')   
        plt.xscale('log')
        plt.yscale('log')
        plt.ylim([1e-3,1])
        plt.xlim([1e-3,1])
        plt.xlabel(r"$\sin^{2}2\theta_{24}$",size=14)
        plt.ylabel(r"$\sin^{2}2\theta_{34}$",size=14)
    else:
        plt.pcolormesh(theta24s*180/pi, theta34s*180/pi, chis.transpose(), vmin=0, vmax=18, cmap="PuBu")
        cbar = plt.colorbar(extend='max')
        ct = plt.contour(theta24s*180/pi, theta34s*180/pi, chis.transpose(), levels=chis_l, cmap='Oranges_r')   

#        plt.ylim([0,30])
#        plt.xlim([0,30])
        plt.xlabel(r"$\theta_{24}$ [deg]",size=14)
        plt.ylabel(r"$\theta_{34}$ [deg]",size=14)


    cbar.set_label(r"-2$\Delta$LLH", size=14)


    plt.title(r"Sensitivity at $\Delta m_{14}^{2}=$"+"{:.2f}".format(msqs[which_sliver]),size=16)
    set_lbls(ct)
#    plt.text(10,25, "Smithers Preliminary", color="r",size=14)

    plt.tight_layout()
    plt.savefig("cascade_sens_{:.2f}.png".format(msqs[which_sliver]), dpi=400)
    plt.show()


t24s = [9*pi/180]
for t24 in t24s:
    which_sliver = get_loc(t24, theta24s)[0]
    chis = np.zeros(shape=(len(msqs), len(theta34s)))
    for t34 in range(len(theta34s)):
        for msq in range(len(msqs)):
            chis[msq][t34] = chi2[which_sliver][t34][msq]
    plt.pcolormesh(msqs, s2(theta34s), chis.transpose(), vmin=0, vmax=12, cmap="PuBu")
    #plt.pcolormesh(msqs, theta34s*180/pi, chis.transpose(), vmin=0, vmax=12, cmap="PuBu")
    cbar = plt.colorbar(extend='max')
    cbar.set_label(r"-2$\Delta$LLH", size=14)
    ct = plt.contour(msqs, s2(theta34s), chis.transpose(), levels=chis_l, cmap='cool')   
    #ct = plt.contour(msqs, theta34s*180/pi, chis.transpose(), levels=chis_l, cmap='cool')   
    set_lbls(ct)
    plt.title(r"90% CL Sensitivity with $\theta_{24}=$"+"{:.2f}".format(theta24s[which_sliver]),size=16)
    plt.yscale('log')
    plt.xscale('log')
    plt.xlim([1e-2, 1e2])
    plt.ylim([1e-3, 1])
    plt.xlabel(r"$\Delta m_{14}^{2}$ [eV$^{2}$]",size=14)
    plt.ylabel(r"$\sin^{2}2\theta_{34}$ [deg]",size=14)
    plt.tight_layout()
    plt.show()



# these were written using the older likelihood setup 
if old:
    levels = [90.0, 68.0]

    # sum over the mass-squareds 
    c_prob_msq = np.zeros(shape=(len(theta24s), len(theta34s)))
    for entry in sorted_likelihoods:
        c_prob_msq[entry[0]][entry[1]] += entry[3]

    plt.contour( theta24s, theta34s, 100*c_prob_msq.transpose(), levels=levels)
    plt.xlabel(r"$\theta_{24}$",size=14)
    plt.ylabel(r"$\theta_{34}$",size=14)
    plt.show()


    c_prob_th34 = np.zeros(shape=(len(theta24s), len(msqs)))
    for entry in sorted_likelihoods:
        c_prob_th34[entry[0]][entry[2]] += entry[3]

    plt.contour( msqs, theta24s, 100*c_prob_th34, levels=levels)
    plt.ylabel(r"$\theta_{24}$",size=14)
    plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
    plt.show()

    c_prob_th24 = np.zeros(shape=(len(theta34s), len(msqs)))
    for entry in sorted_likelihoods:
        c_prob_th24[entry[1]][entry[2]] += entry[3]

    plt.contour( msqs, theta34s, 100*c_prob_th24, levels=levels)
    plt.ylabel(r"$\theta_{34}$",size=14)
    plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
    plt.show()


