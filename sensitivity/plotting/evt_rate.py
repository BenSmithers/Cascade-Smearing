"""
Make a plot of our event rate at each bin
"""

f_name = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_smearedwell_0.0_0.0_0.0_0.0.dat"

import pickle
import numpy as np
from math import pi 

import matplotlib 
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

f = open(f_name, 'rb')
data = pickle.load(f)
f.close()

e_edges = data["e_edges"]
#a_edges = data["a_edges"]
a_edges = np.linspace(-1,1,11)
print(a_edges)
evt_r = data["event_rate"]

core_b = -0.98
mantle_b= -0.83

declination = True
def decd(obj):
    """
    declination convert. 
    """
    return np.sin(-(pi/2) + np.arccos(obj))

if declination:
    a_edges = decd(a_edges)
    core_b = decd(core_b)
    mantle_b = decd(mantle_b)

plt.pcolormesh(a_edges, e_edges, evt_r, cmap='viridis')

plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
cbar = plt.colorbar()
cbar.set_label("Events")
if declination:
    plt.xlabel(r"$\sin\delta^{reco}$")
else:
    plt.xlabel(r"$\cos\theta_{z}^{reco}$")
plt.ylabel(r"$E^{reco}$ [GeV]")
plt.yscale('log')
plt.ylim([1e2,1e6])
#plt.xlim([-1,0.2])
plt.show()
