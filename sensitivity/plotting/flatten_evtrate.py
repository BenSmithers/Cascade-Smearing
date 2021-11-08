"""
This makes two plots; flattened event rates! 
"""
from scipy.sparse import data
from cascade.utils import gen_filename, SterileParams, config
import pickle
import numpy as np
import os

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))

from cascade.deporeco import DataReco

import os 

fn = "best_expected_flux.dat"
data_folder = os.path.join(config["datapath"], "expected_fluxes_reco")
f_name =  gen_filename(data_folder, fn, SterileParams())
sterile_fname = gen_filename(data_folder, fn, SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.6416))

#sterile_string = r"$\sin^{2}(2\theta_{24})=0.1, \sin^{2}(2\theta_{34})=0.20, \Delta m_{41}^{2}=4.64$"
sterile_string = r"Sterile Neutrino"


def loadit(filename):
    f = open(filename, 'rb')
    tdata = pickle.load(f)
    f.close()
    return tdata

data = loadit(f_name)
sdata = loadit(sterile_fname)

e_edges = data["e_edges"]
print(len(e_edges))
a_edges = np.linspace(-1,1.0,11)
print(a_edges)
evt_r = data["event_rate"]

dataobj = DataReco(e_edges*(1e9), a_edges, e_edges*(1e9), a_edges)
tensor = [[[[ dataobj.get_energy_reco_odds(j,l)*dataobj.get_czenith_reco_odds(k,i,l) for i in range(10)] for j in range(20)] for k in range(10)] for l in range(20)]
evt_r = np.einsum('ij,klij', evt_r,tensor)
s_evt_r = np.einsum('ij,klij', sdata["event_rate"],tensor)

fig, axes = plt.subplots(2,1,figsize=(7,10))
fig.subplots_adjust(wspace=0.1)
fig.subplots_adjust(bottom=0.08, left=0.15, right=0.97, top=0.98)


from math import pi
declination = np.sin(-(pi/2) + np.arccos(a_edges))
declination = a_edges


binned_in_e = np.sum(evt_r, axis=1)
binned_in_ang = np.sum(evt_r, axis=0)

sbinned_in_e = np.sum(s_evt_r, axis=1)
sbinned_in_ang = np.sum(s_evt_r, axis=0)

axes[0].bar( x=e_edges[:-1], height=binned_in_e, width=e_edges[1:]-e_edges[:-1], align='edge', color='salmon',alpha=0.45, label=r"No Sterile")
axes[0].step( x=e_edges[1:], y=binned_in_e, color='salmon')
axes[0].step( x=e_edges[1:], y=sbinned_in_e, color='black',label=sterile_string)
axes[0].set_xscale('log')
axes[0].set_xlim([1e3,1e6])
axes[0].set_xlabel(r"$E_{\nu}$"+" [GeV]")
axes[0].set_ylabel("Events")
#axes[0].grid(which='both', alpha=0.5)
#axes[0].set_ylim([5e-1,1800])
axes[0].legend(prop={'size':18})

axes[1].bar( x=declination[:-1], height=binned_in_ang, width=declination[1:]-declination[:-1], align='edge', color='salmon', alpha=0.45, label=r"No Sterile")
axes[1].step( x=declination[1:], y=binned_in_ang, color='salmon')
plt.hlines(binned_in_ang[0], -1, declination[1], color='salmon')

axes[1].step( x=declination[1:], y=sbinned_in_ang, color='black', label=sterile_string)

plt.hlines(sbinned_in_ang[0], -1, declination[1], color='black')
axes[1].set_xlim([-1.0,0.2])
#axes[1].set_xlabel(r"$\sin\delta^{reco}$")
axes[1].set_xlabel(r"$\cos\theta_{z}$")
axes[1].set_ylabel("Events")
axes[1].set_ylim([0,1300])
axes[1].legend(prop={'size':18})
#axes[1].grid(which='both', alpha=0.5)
#plt.tight_layout()
plt.savefig("bothflat.png", dpi=400)
plt.show()
