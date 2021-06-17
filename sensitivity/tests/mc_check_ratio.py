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

# get the flux in the centers of each of the bins
null_flux_d = build_mc_flux(null)
ster_flux_d = build_mc_flux(ster)

null_flux = null_flux_d["event_rate"]
ster_flux = null_flux_d["event_rate"]

e_edges = null_flux_d["e_edges"]
a_edges = np.cos(null_flux_d["a_edges"])

plt.pcolormesh(a_edges, e_edges, ster_flux/null_flux)
plt.colorbar()
plt.xlabel("Cos theta")
plt.ylabel("Energy")
plt.yscale('log')
plt.show()
