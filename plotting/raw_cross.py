"""
This script was written to take the propagated nusquids fluxes and convolve them with nusquids' total cross sections 

There was a little thing added in to handle some tau smearing - but that only really has an effect when we actually do this in a differential xs way
"""

import pickle
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")

from cascade.utils import SterileParams, config, gen_filename
from cascade.utils import bhist

from cascade.cross_section_test import get_total_flux as get_xs 
from cascade.nus_utils import get_flavor, get_neut, get_curr 

from cascade.raw_fluxes import raw_flux


width = 0.5

def _load_flux(params, tracks=False):
    name = gen_filename(config["datapath"], config["nu_flux"]+".dat", params)

    if False:
        f = open(name,'rb')
        all_data = pickle.load(f)
        f.close()
    else:
        kwargs ={"as_data":True}
        data = raw_flux(params, kwargs=kwargs) 
        flux = {}
        for key in data.get_keys(just_casc=(not tracks),just_tracks=tracks):
            flux[key] = np.array(data.fluxes[key])
        return (np.array(data.energies), np.array(data.angles), flux)

    return( all_data["e_true"], all_data["a_true"], all_data["flux"] )

null = SterileParams(0., 0., 0., 0.)
ster = SterileParams(0., 0.1609, 0.2205, 4.47)
#ster = SterileParams(0., 0.1609, 0.0, 4.47)
print("Loading files from {}".format(config['datapath']))
e_true, a_true, flux_null    = _load_flux(null, False)
e_true, a_true, flux_sterile = _load_flux(ster, False)

keys = flux_null.keys()
centers = bhist([e_true]).centers


for key in keys:
    flav = get_flavor(key)
    neut = get_neut(key)
    curr = get_curr(key) 

    for i_energy in range(len(centers)):
        xs = get_xs(centers[i_energy], flav, neut, curr)

        pcent = 1.0
        
        if False: #("Tau" in key) and ("CC" in key):
            pcent = 0.51

        flux_null[key][i_energy] *= xs*pcent
        flux_sterile[key][i_energy] *= xs*pcent


def get_ratio(nubar):

    just_nubar = True

    keep_key = "nuBar" if nubar else "nu_"

    
    ex = list(flux_null.keys())[0]
    null_total      = np.zeros(shape = np.shape(flux_null[ex]))
    sterile_total   = np.zeros(shape = np.shape(flux_null[ex]))

    for key in keys:

        if just_nubar and (keep_key not in key):
            print("Skip {}".format(key))
            continue

        null_total += flux_null[key]
        sterile_total += flux_sterile[key]

    ratio = 100*(sterile_total / null_total) -100

    return e_true, a_true, ratio

e_true, a_true, ratio_nubar = get_ratio(True)
e_true, a_true, ratio_nu = get_ratio(False)

#energies = np.array(bhist([e_true]).centers)
#czeniths = np.array(bhist([a_true]).centers)

fig, axes = plt.subplots(1,2,sharey=True,figsize=(11,5),
                        gridspec_kw={"width_ratios":[1,1.2]})

cf = axes[0].pcolormesh(a_true, e_true/(1e9), ratio_nu, cmap=cm.cividis, vmin=-100, vmax=100)
cf2 = axes[1].pcolormesh(a_true, e_true/(1e9), ratio_nubar, cmap=cm.cividis, vmin=-100, vmax=100)
core_b = -0.98
mantle_b= -0.83

axes[1].vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
axes[1].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[1].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
axes[1].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[1].text(0.05, 4e4, r"$\bar{\nu}$", fontsize='xx-large', color='white')


axes[0].vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
axes[0].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[0].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
axes[0].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
axes[0].text(0.05, 4e4, r"$\nu$", fontsize='xx-large', color='white')

axes[0].set_yscale('log')
axes[0].set_ylabel(r"E$_{\nu}^{true}$ [GeV]")
axes[0].set_xlabel(r"$\cos\theta_{z}^{true}$")
axes[1].set_xlabel(r"$\cos\theta_{z}^{true}$")

axes[0].set_ylim([10**2, 10**5])
axes[0].set_xlim([-1.0, 0.2])
axes[1].set_xlim([-1.0, 0.2])

from mpl_toolkits.axes_grid1.inset_locator import InsetPosition

cbar = fig.colorbar(cf2,ax=axes[1]) #cf,ticks=ticker.LogLocator())
cbar.set_label("Difference [%]")

fig.subplots_adjust(wspace=0,hspace=0)

plt.show()
