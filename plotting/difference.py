import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from utils import bhist

import pickle

def _load_flux(name):
    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_reco"]
    a_reco = all_data["a_reco"]
    flux = all_data["flux"]

    return( e_reco, a_reco, flux )


e_reco, a_reco, flux_null = _load_flux(".flux_data_null.dat")
e_reco, a_reco, flux_sterile = _load_flux(".flux_data.dat")

ex = list(flux_null.keys())[0]

null_total = np.zeros(shape = np.shape(flux_null[ex]))
sterile_total = np.zeros(shape = np.shape(flux_null[ex]))

for key in flux_null.keys():
    null_total += flux_null[key]
    sterile_total+=flux_sterile[key]

ratio = sterile_total / null_total

energies = np.array(bhist([e_reco]).centers)
czeniths = np.array(bhist([a_reco]).centers)

cf = plt.pcolormesh(czeniths, energies/(1e9), ratio, cmap=cm.coolwarm)
plt.yscale('log')
plt.ylabel("Reco Energy [GeV]", size=14)
plt.xlabel(r"Reco $\cos\theta$",size=14)
plt.title("",size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"Sterile Flux / Null Flux")
plt.savefig("flux_ratio.png",dpi=400)
plt.show()

