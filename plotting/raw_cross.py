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

from cascade.utils import SterileParams, config, gen_filename
from cascade.utils import bhist

from cascade.cross_section_test import get_total_flux as get_xs 
from cascade.nus_utils import get_flavor, get_neut, get_curr 

width = 0.2

def _load_flux(params):
    name = gen_filename(config["datapath"], config["nu_flux_downsize"]+".dat", params)

    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    return( all_data["e_true"], all_data["a_true"], all_data["flux"] )

just_nubar = True
keep_key = "_NC"


null = SterileParams(0., 0., 0., 0.)
ster = SterileParams(0., 0.1609, 0.2205, 4.47)
#ster = SterileParams(0., 0.1609, 0.0, 4.47)
print("Loading files from {}".format(config['datapath']))
e_true, a_true, flux_null    = _load_flux(null)
e_true, a_true, flux_sterile = _load_flux(ster)

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

ex = list(flux_null.keys())[0]
null_total      = np.zeros(shape = np.shape(flux_null[ex]))
sterile_total   = np.zeros(shape = np.shape(flux_null[ex]))

for key in keys:
    if just_nubar and (keep_key not in key):
        print("Skip {}".format(key))
        continue
    null_total += flux_null[key]
    sterile_total += flux_sterile[key]

ratio = sterile_total / null_total
energies = np.array(bhist([e_true]).centers)
czeniths = np.array(bhist([a_true]).centers)

cf = plt.pcolormesh(a_true, e_true/(1e9), ratio, cmap=cm.coolwarm, vmin=1.0-width, vmax=1.0+width)
plt.yscale('log')
plt.ylabel("True Energy [GeV]", size=14)
plt.xlabel(r"True $\cos\theta$",size=14)
plt.ylim([10**2, 10**6])
if just_nubar:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"Sterile Flux-XS / Null Flux-XS")
plt.show()
