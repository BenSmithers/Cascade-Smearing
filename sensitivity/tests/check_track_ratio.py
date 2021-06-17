from cascade.utils import Data, SterileParams, gen_filename, config
from cascade.sensitivity.make_from_mc import build_mc_flux
from cascade.sensitivity.eff_area_reader import quickload

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import os

# load the fluxes 
null_f = gen_filename(config["datapath"], "raw_det_flux.dat", SterileParams())
ster_f = gen_filename(config["datapath"], "raw_det_flux.dat", SterileParams(0.0, 0.1609, 0.0, 4.7))
null = Data(null_f)
ster = Data(ster_f)

# get the binning information 
filename = "effective_area.nu_mu.txt"
area_data = quickload(os.path.join(os.path.join(config["datapath"], "charm_search_supplemental/"), filename))
e_edges = np.array(area_data["e_reco"])
a_edges = np.array(area_data["cos_th"])

e_centers = 0.5*(e_edges[:-1] + e_edges[1:])
a_centers = 0.5*(a_edges[:-1] + a_edges[1:])

# get the flux in the centers of each of the bins
null_flux = np.zeros(shape=(len(e_centers), len(a_centers)))
ster_flux = np.zeros(shape=(len(e_centers), len(a_centers)))
keys = ["Mu_nu_CC", "Mu_nuBar_CC"]

for i_e in range(len(e_centers)):
    for i_a in range(len(a_centers)):
        for key in keys:
            null_flux[i_e][i_a] = null_flux[i_e][i_a] + null.get_flux(e_centers[i_e]*(1e9), key, angle=a_centers[i_a])
            ster_flux[i_e][i_a] = ster_flux[i_e][i_a] + ster.get_flux(e_centers[i_e]*(1e9), key, angle=a_centers[i_a])


plt.pcolormesh(a_edges, e_edges, ster_flux/null_flux)
plt.colorbar()
plt.xlabel("Cos theta")
plt.ylabel("Energy")
plt.yscale('log')
plt.show()
