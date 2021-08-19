"""
This makes two plots; flattened event rates! 
"""


f_name = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_smearedwell_0.0_0.0_0.0_0.0.dat"

import pickle
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

f = open(f_name, 'rb')
data = pickle.load(f)
f.close()

e_edges = data["e_edges"]
print(len(e_edges))
a_edges = np.linspace(-1,1,11)
print(a_edges)
evt_r = data["event_rate"]


fig, axes = plt.subplots(1,2,sharey=False,figsize=(11,5))

from math import pi
declination = np.sin(-(pi/2) + np.arccos(a_edges))


binned_in_e = np.sum(evt_r, axis=1)
binned_in_ang = np.sum(evt_r, axis=0)

axes[0].bar( x=e_edges[:-1], height=binned_in_e, width=e_edges[1:]-e_edges[:-1], align='edge')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlim([1e3,1e6])
axes[0].set_xlabel("$E^{reco}$ [GeV]")
axes[0].set_ylabel("Events")
axes[0].grid(which='both', alpha=0.5)
axes[0].set_ylim([5e-1,3e3])

axes[1].bar( x=declination[:-1], height=binned_in_ang, width=declination[1:]-declination[:-1], align='edge')
#axes[1].set_xlim([-1,0.2])
axes[1].set_xlabel(r"$\sin\delta^{reco}$")
axes[1].set_ylabel("Events")
axes[1].grid(which='both', alpha=0.5)
plt.show()
