"""
This makes two plots; flattened event rates! 
"""
from scipy.sparse import data
from cascade.utils import gen_filename, SterileParams, config
import pickle
import numpy as np

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")
from cascade.deporeco import DataReco

import os 

fn = "best_expected_flux.dat"
data_folder = os.path.join(config["datapath"], "expected_fluxes_reco")
f_name =  gen_filename(data_folder, fn, SterileParams())


f = open(f_name, 'rb')
data = pickle.load(f)
f.close()

e_edges = data["e_edges"]
print(len(e_edges))
a_edges = np.linspace(-1,1.0,11)
print(a_edges)
evt_r = data["event_rate"]

dataobj = DataReco(e_edges*(1e9), a_edges, e_edges*(1e9), a_edges)
tensor = [[[[ dataobj.get_energy_reco_odds(j,l)*dataobj.get_czenith_reco_odds(k,i,l) for i in range(10)] for j in range(20)] for k in range(10)] for l in range(20)]
evt_r = np.einsum('ij,klij', evt_r,tensor)

fig, axes = plt.subplots(1,2,sharey=False,figsize=(11,5))

from math import pi
declination = np.sin(-(pi/2) + np.arccos(a_edges))
declination = a_edges


binned_in_e = np.sum(evt_r, axis=1)
binned_in_ang = np.sum(evt_r, axis=0)

axes[0].bar( x=e_edges[:-1], height=binned_in_e, width=e_edges[1:]-e_edges[:-1], align='edge', color=(217/255,152/255,82/255 ))
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_xlim([1e3,1e6])
axes[0].set_xlabel(r"$E_{\nu}$"+" [GeV]")
axes[0].set_ylabel("Events")
#axes[0].grid(which='both', alpha=0.5)
axes[0].set_ylim([5e-1,3e3])

axes[1].bar( x=declination[:-1], height=binned_in_ang, width=declination[1:]-declination[:-1], align='edge', color=(141/255, 208/255, 214/255))
axes[1].set_xlim([-1.0,0.2])
#axes[1].set_xlabel(r"$\sin\delta^{reco}$")
axes[1].set_xlabel(r"$\cos\theta_{z}$")
axes[1].set_ylabel("Events")
#axes[1].grid(which='both', alpha=0.5)
plt.savefig("bothflat.png", dpi=400)
plt.show()
