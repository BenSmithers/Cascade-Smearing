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
import os

from cascade.utils import SterileParams, config, gen_filename
from cascade.utils import bhist

from cascade.cross_section_test import get_total_flux as get_xs 
from cascade.nus_utils import get_flavor, get_neut, get_curr 

width = 0.1

def _load_flux(name):

    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    return( all_data["e_true"], all_data["a_true"], all_data["flux"] )

just_nubar = True
keep_key = "Tau"


null = SterileParams(0., 0., 0., 0.)
ster = SterileParams(0., 0.1609, 0.2296, 4.47)

e_true, a_true, flux_null    = _load_flux("/home/benito/software/data/cascade/poly_sib/downsized_raw_det_flux_0.0_0.0_0.0_0.0.dat")
e_true, a_true, flux_tau = _load_flux("/home/benito/software/data/cascade/poly_sib/downsized_raw_det_flux_0.0_0.0_0.0_0.0Copy of .dat")

keys = flux_null.keys()
centers = bhist([e_true]).centers


for key in keys:
    flav = get_flavor(key)
    neut = get_neut(key)
    curr = get_curr(key) 

    for i_energy in range(len(centers)):
        xs = get_xs(centers[i_energy], flav, neut, curr)
    
        if "nubar" in key.lower():
            P = -1
        else:
            P = 1

        pcent = 1.0

        if False: #"tau" in key.lower():
            pcent = 0.5

        flux_null[key][i_energy] *= xs*pcent
        flux_tau[key][i_energy] *= xs*pcent

ex = list(flux_null.keys())[0]
null_total      = np.zeros(shape = np.shape(flux_null[ex]))
sterile_total   = np.zeros(shape = np.shape(flux_null[ex]))

for key in keys:
    if just_nubar and (keep_key not in key):
        print("Skip {}".format(key))
        continue
    null_total += flux_null[key]
    sterile_total += flux_tau[key]

ratio = sterile_total / null_total
energies = np.array(bhist([e_true]).centers)
czeniths = np.array(bhist([a_true]).centers)

cf = plt.pcolormesh(a_true, e_true/(1e9), ratio, cmap=cm.coolwarm, vmin=1.0-width, vmax=1.0+width)
plt.yscale('log')
plt.ylabel("Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_nubar:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"Tau Regen Flux-XS On/Off")
plt.savefig(os.path.join(config["datapath"], config["img_subdir"] + "tauregen.png"), dpi=400)
plt.show()
