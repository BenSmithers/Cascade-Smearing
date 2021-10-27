"""
Make a plot of our event rate at each bin. We will have two pannels! 
"""

#f_name_s = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_smearedwell_0.0_0.0_0.0_0.0.dat"
f_name_s = "/home/benito/Downloads/best_expected_flux_0.0_0.0_0.0_0.0.dat"
f_name_s = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_from_mc_0.0_0.0_0.0_0.0.dat"
#f_name_s = "/home/benito/software/data/cascade/hg_sib//expected_fluxes_reco/0.1609e0/best_expected_flux_0.0_0.1609e0_0.2247e0_4.4700e0.dat"
#f_name_s = "/home/benito/Downloads/best_expected_flux_0.0_0.1581e-1_0.1581e-1_0.0.dat"

f_name = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_0.0_0.0_0.0_0.0.dat"
#f_name = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_from_mc_0.0_0.0_0.0_0.0.dat"

import pickle
import numpy as np
from math import pi 

import matplotlib 
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

from cascade.utils import SterileParams, gen_filename
from cascade.deporeco import DataReco


def ldata(fname):
    f = open(fname, 'rb')
    data = pickle.load(f)
    f.close()
    return data

central_s = SterileParams()
sterile_s = SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.6416)
use_params = sterile_s

cascade_name_root = "best_expected_flux.dat"
track_name_root = "expected_flux_from_mc_smearedwell.dat"
datadir = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/"

cascade_fname = gen_filename(datadir, cascade_name_root, use_params)
track_fname = gen_filename(datadir, track_name_root, use_params)

data = ldata(f_name)
data2 = ldata(f_name_s)

cascade_datadict = ldata(cascade_fname)
track_datadict = ldata(track_fname)


e_edges = cascade_datadict["e_edges"]
a_edges = cascade_datadict["a_edges"]
n_a = len(a_edges)-1
n_e = len(e_edges)-1
_reco_obj = DataReco(e_edges*(1e9), a_edges, e_edges*(1e9), a_edges)
_reco_tensor = [[[[ _reco_obj.get_energy_reco_odds(j,l)*_reco_obj.get_czenith_reco_odds(k,i,l) for i in range(n_a)] for j in range(n_e)] for k in range(n_a)] for l in range(n_e)]
casc_e = np.einsum('ij,klij',cascade_datadict["event_rate"], _reco_tensor)

track_e = track_datadict["event_rate"]

core_b = -0.98
mantle_b= -0.83



declination = False
def decd(obj):
    """
    declination convert. 
    """
    return np.sin(-(pi/2) + np.arccos(obj))

if declination:
    a_edges = decd(a_edges)
    core_b = decd(core_b)
    mantle_b = decd(mantle_b)

fig, axes = plt.subplots(2,2,figsize=(10,6.5),sharey='row', gridspec_kw={'height_ratios':[0.05,1]})
fig.subplots_adjust(wspace=0.10)
# fig.subplots_adjust(bottom=0.15)
print(np.shape(axes))

pc = axes[1][0].pcolormesh(a_edges, e_edges, track_e, cmap='viridis', norm = matplotlib.colors.LogNorm())
axes[1][0].vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
axes[1][0].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[1][0].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
axes[1][0].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
fig.colorbar(pc,cax=axes[0][0], orientation='horizontal')

axes[1][0].set_xlabel(r"$\cos\theta_{z}^{reco}$")
axes[1][0].set_ylabel(r"$E^{reco}$ [GeV]")
axes[1][0].set_yscale('log')
axes[1][0].set_ylim([1e2,1e6])
axes[1][0].set_xlim([-1,0.2])

pc2 = axes[1][1].pcolormesh(a_edges, e_edges, casc_e, cmap='viridis', norm = matplotlib.colors.LogNorm(), vmin=1e-1)
axes[1][1].vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
axes[1][1].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[1][1].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
axes[1][1].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
fig.colorbar(pc2,cax=axes[0][1], orientation='horizontal')

axes[1][1].set_xlabel(r"$\cos\theta_{z}^{reco}$")
axes[1][1].set_yscale('log')
axes[1][1].set_ylim([1e2,1e6])
axes[1][1].set_xlim([-1,0.2])

plt.show()
plt.savefig("event_rate_sterile.png",dpi=400)
"""
plt.pcolormesh(a_edges, e_edges, evt_r, cmap='viridis', vmin=0, vmax=50000)

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
plt.xlim([-1,0.2])
plt.savefig("event_rate_plot.png", dpi=400)
plt.show()
"""