"""
This is a script to show the difference between two saved flux datafiles 

It shows a heatmap of their ratio as a function of the energy and cos(Zenith)
Note that you need to use the centers of the bins, but there's unfortunately a little disparity between how things are saved

downsized fluxes save the bin centers (TODO: fix this...)
reconstructed fluxes save the bin edges 
"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from cascade.utils import bhist, SterileParams, gen_filename
from cascade.utils import config

import pickle

livetime = 10*3600*24*365.

width = 0.1

def _load_flux(name):
    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_reco"]
    a_reco = all_data["a_reco"]
    flux = all_data["flux"]

    return( e_reco, a_reco, flux )


#null_pt = gen_filename(config["datapath"], config["nu_flux"], 0.,0.,0.)

null = SterileParams(0.,0.,0.,0.)
mud = SterileParams(0., 0.1609, 0.2205, 4.47)
eld = SterileParams(0.13, 0., 0.0, 1.3)

e_reco, a_reco, flux_null = _load_flux(gen_filename(config["datapath"], config["recon_flux"]+".dat", null))
e_reco, a_reco, flux_sterile = _load_flux(gen_filename(config["datapath"], config["recon_flux"]+".dat", mud))

ex = list(flux_null.keys())[0]

null_total = np.zeros(shape = np.shape(flux_null[ex]))
sterile_total = np.zeros(shape = np.shape(flux_null[ex]))

just_nubar = False

keep_key = "Tau"
for key in flux_null.keys():
    if just_nubar and (keep_key not in key):
        print("Skip {}".format(key))
        continue

    null_total += flux_null[key]
    sterile_total+=flux_sterile[key]

ratio = sterile_total / null_total

energies = np.array(bhist([e_reco]).centers)
czeniths = np.array(bhist([a_reco]).centers)

cf = plt.pcolormesh(czeniths, energies/(1e9), ratio, cmap=cm.coolwarm, vmin=1.0-width, vmax=1.0+width)
plt.yscale('log')
plt.ylabel("Reco Energy [GeV]", size=14)
plt.xlabel(r"Reco $\cos\theta$",size=14)
if just_nubar:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"Sterile Flux / Null Flux")
plt.savefig("flux_ratio" +("_{}".format(keep_key) if just_nubar else "")+"_muonD.png",dpi=400)
plt.show()

