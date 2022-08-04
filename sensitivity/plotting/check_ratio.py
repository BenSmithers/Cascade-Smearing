"""
This file is used for plotting out the ratio between two expected fluxes
"""

from cascade.utils import SterileParams, gen_filename, config

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import os
import pickle
from math import pi

plt.style.use("/home/bsmithers/software/cascade/cascade.mplstyle")

# /home/benito/software/data/cascade/hg_sib//expected_fluxes_reco/0.1609e0/best_expected_flux_0.0_0.1609e0_0.2247e0_4.4700e0.dat

# pass the fluxes to the MC reader
# it bins each of the events weight*flux 
# so we get an array of N_events in each bin
#   alternatively just load the thingy if this is already done 
def load(param):
    """
    Use this to load in the fluxes. If they're there, load them. If not, make them

    pass on the sterileparams object so it knows what to load and what to pass to the data maker
    """
    #filename = "null_exp.dat" if null_bool else "ster_exp.dat"
    where = "/home/bsmithers/software/data/hg_sib/expected_fluxes_reco"

    #filename = gen_filename(where, "expected_flux_smeared.dat", param)
    filename = gen_filename(where, "best_expected_flux.dat", param)
    
    if os.path.exists(filename):
        print("Loading {}".format(filename))
        f = open(filename,'rb')
        data = pickle.load(f)
        f.close()
        return data
    else:
        # only import these if we need to 
        from cascade.sensitivity.make_from_mc import build_mc_flux
        from cascade.sensitivity.generate_all_integrated_fluxes import make_meta_flux

        data = make_meta_flux(param, do_mc=False)
        return data


sp= SterileParams(theta13=0.1609e0, theta23=0.2247e0, msq2=4.47e0)
null_flux_d = load(SterileParams())
ster_flux_d = load(sp)


null_flux = null_flux_d["event_rate"]
ster_flux = ster_flux_d["event_rate"]

print(np.shape(null_flux))
print(np.shape(ster_flux))

e_edges = null_flux_d["e_edges"]
a_edges = np.linspace(-1,1,11) # null_flux_d["a_edges"]

print("E bins {}".format(len(e_edges)))

print(e_edges)
print(a_edges)
amount = 100*(ster_flux/null_flux)

plt.pcolormesh(a_edges, e_edges, amount,vmin=70,vmax=130, cmap="coolwarm")
cf = plt.colorbar()
cf.set_label("Difference [%]")
#plt.title("Appearance with Bin-Breaking",size=14)
plt.xlabel(r"$\cos\theta_{z}^{reco}$")
plt.ylabel(r"$E_{\nu}^{reco}$ [GeV]")
plt.xlim([-1,0.2])
plt.ylim([1e2, 1e6])
plt.yscale('log')
plt.show()
