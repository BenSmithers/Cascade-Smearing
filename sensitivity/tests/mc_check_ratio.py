from cascade.utils import Data, SterileParams, gen_filename, config
from cascade.sensitivity.make_from_mc import build_mc_flux

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import os
import pickle

# pass the fluxes to the MC reader
# it bins each of the events weight*flux 
# so we get an array of N_events in each bin
#   alternatively just load the thingy if this is already done 
def load(null_bool):
    filename = "null_exp.dat" if null_bool else "ster_exp.dat"

    if False: #os.path.exists(filename):
        f = open(filename,'rb')
        data = pickle.load(f)
        f.close()
        return data
    else:
        if null_bool:
            params = SterileParams()
        else:
            params = SterileParams(0.0, 0.1609, 0.0, 4.7)
        # load the fluxes 
        fu_filename = gen_filename(config["datapath"], "raw_det_flux.dat", params)
        flux = Data(fu_filename)

        data = build_mc_flux(flux)
        f = open(filename, 'wb')
        pickle.dump(data, f, -1)
        f.close()
        return data

null_flux_d = load(True)
ster_flux_d = load(False)


null_flux = null_flux_d["event_rate"]
ster_flux = ster_flux_d["event_rate"]

e_edges = null_flux_d["e_edges"]
a_edges = null_flux_d["a_edges"]

plt.pcolormesh(a_edges, e_edges, ster_flux/null_flux)
plt.colorbar()
plt.xlabel("Cos theta")
plt.ylabel("Energy")
plt.yscale('log')
plt.show()
