"""
Make a plot of our event rate at each bin. We will have two pannels! 
"""

from csv import excel_tab
import pickle
import numpy as np
from math import pi 
import os

import matplotlib 
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))

from cascade.utils import SterileParams, gen_filename, config
from cascade.deporeco import DataReco
from cascade.sensitivity.generate_all_integrated_fluxes import make_meta_flux

def ldata(fname):
    print("Loading {}".format(fname))
    f = open(fname, 'rb')
    data = pickle.load(f)
    f.close()
    return data

ratios = False

central_s = SterileParams()
#sterile_s = SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.6416)
sterile_s = SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.5)
use_params = central_s

cascade_name_root = "best_expected_flux.dat"
track_name_root = "expected_flux_from_mc_smearedwell.dat"
datadir = os.path.join(config["datapath"], "expected_fluxes_reco")
# datadir = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/"

sterile_casc_fname = gen_filename(datadir, cascade_name_root, sterile_s)
sterile_track_fname = gen_filename(datadir, track_name_root, sterile_s)


cascade_fname = gen_filename(datadir, cascade_name_root, central_s)
track_fname = gen_filename(datadir, track_name_root, central_s)

#data = ldata(f_name)
#data2 = ldata(f_name_s)
try:
    cascade_datadict = ldata(cascade_fname)
except IOError:
    cascade_datadict = make_meta_flux(central_s)
try:
    track_datadict = ldata(track_fname)
except IOError:
    track_datadict = make_meta_flux(central_s, True)
try:
    st_cascade_data = ldata(sterile_casc_fname)
except IOError:
    st_cascade_data = make_meta_flux(sterile_s)
try:
    st_track_data = ldata(sterile_track_fname)
except IOError:
    st_track_data = make_meta_flux(sterile_s, True)



e_edges = cascade_datadict["e_edges"]
a_edges = cascade_datadict["a_edges"]
n_a = len(a_edges)-1
n_e = len(e_edges)-1
_reco_obj = DataReco(e_edges*(1e9), a_edges, e_edges*(1e9), a_edges)
_reco_tensor = [[[[ _reco_obj.get_energy_reco_odds(j,l)*_reco_obj.get_czenith_reco_odds(k,i,l) for i in range(n_a)] for j in range(n_e)] for k in range(n_a)] for l in range(n_e)]


if ratios:
    ster_casc_e = np.einsum('ij,klij',st_cascade_data["event_rate"], _reco_tensor)
    cent_casc_e = np.einsum('ij,klij',cascade_datadict["event_rate"], _reco_tensor)
    casc_e = ster_casc_e/cent_casc_e
    track_e = st_track_data["event_rate"] / track_datadict["event_rate"]

else:
    casc_e = np.einsum('ij,klij',cascade_datadict["event_rate"], _reco_tensor)
    track_e = track_datadict["event_rate"]

core_b = -0.98
mantle_b= -0.83

fontsize=24

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

fig, axes = plt.subplots(2,1,figsize=(7,10), sharex=True) #, gridspec_kw={'height_ratios':[0.05,1], 'width_ratios':[1,1]})
fig.subplots_adjust(wspace=0.15, hspace=0.05, left=0.13, right=0.95,top=0.98, bottom=0.075)
# fig.subplots_adjust(bottom=0.15)

# norm = matplotlib.colors.LogNorm()
if ratios:
    pc = axes[0].pcolormesh(a_edges, e_edges, track_e, cmap='coolwarm', vmin=0.75, vmax=1.25 )
else:
    pc = axes[0].pcolormesh(a_edges, e_edges, track_e, cmap='PuBu', norm = matplotlib.colors.LogNorm())

    for i_e in range(len(e_edges) -1 ):
        for i_a in range(len(a_edges) - 1):
            axes[0].text( 0.5*a_edges[i_a]+0.5*a_edges[i_a+1], e_edges[i_e]*0.9+e_edges[i_e+1]*0.1, "{:.1f}".format(track_e[i_e][i_a]), color="white", fontsize="small", rotation=90)

#axes[0].vlines(core_b,ymin=1e2, ymax=10**6, colors="black", ls="-")
#axes[0].text(core_b+0.02, 5.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='black')
axes[0].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="black", ls="--")
axes[0].text(mantle_b+0.02, 5.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='black')

axes[0].set_ylabel(r"$E^{reco}$ [GeV]", size=fontsize)
axes[0].set_yscale('log')
axes[0].set_ylim([4e2,1e4])
axes[0].set_xlim([-1,0.2])
cbar = fig.colorbar(pc, ax=axes[0]) #, anchor=(0,0.9), panchor=(0,1))
cbar.set_label("Events")


#norm = matplotlib.colors.LogNorm()
if ratios:
    pc2 = axes[1].pcolormesh(a_edges, e_edges, casc_e, cmap='coolwarm', vmin=0.75, vmax=1.25)
else:
    pc2 = axes[1].pcolormesh(a_edges, e_edges, casc_e, cmap='PuBu', vmin=1e-1, norm = matplotlib.colors.LogNorm())
    for i_e in range(len(e_edges) -1 ):
        if e_edges[i_e]>=7e5:
            continue
        for i_a in range(len(a_edges) - 1):
            axes[1].text( 0.8*a_edges[i_a]+0.2*a_edges[i_a+1], e_edges[i_e]*0.8+e_edges[i_e+1]*0.2, "{:.1f}".format(casc_e[i_e][i_a]), color="white", fontsize="small", rotation=0)



#axes[1].vlines(core_b,ymin=1e2, ymax=10**6, colors="black", ls="-")
#axes[1].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='black')
axes[1].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="black", ls="--")
axes[1].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='black')

axes[1].set_xlabel(r"$\cos\theta_{z}^{reco}$", size=fontsize)
axes[1].set_ylabel(r"$E^{reco}$ [GeV]", size=fontsize)
axes[1].set_yscale('log')
axes[1].set_ylim([1e2,7e5])
axes[1].set_xlim([-1,0.2])
cbar=fig.colorbar(pc2,ax=axes[1] ) #, anchor=(0.5,0.9), panchor=(0.5,1))
cbar.set_label("Events")
#plt.tight_layout()

if ratios:
    plt.savefig("event_rate_ratios.png",dpi=400)
else:
    plt.savefig("event_rate.png",dpi=400)


plt.show()
