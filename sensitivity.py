import numpy as np
import os
import pickle

from utils import sep_by_flavor, bhist

import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import ticker
import matplotlib.cm as cm

filename = ".flux_data.dat"
def _load_flux(filename):
    f = open(filename, 'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_reco"]
    a_reco = all_data["a_reco"]
    flux = all_data["flux"]

    return( e_reco, a_reco, flux )


e_reco, a_reco, flux = _load_flux(filename)

energies = np.array(bhist([e_reco]).centers)
czeniths = np.array(bhist([a_reco]).centers)

from_muon, from_not = sep_by_flavor(flux)

cf = plt.contourf( energies/(1e9), czeniths, np.transpose(from_not), cmap=cm.coolwarm, locator=ticker.LogLocator())
plt.xscale('log')
plt.xlabel("Event Energy [GeV]", size=14)
plt.ylabel(r"Reconstructed $\cos\theta$",size=14)
cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
cbar.set_label(r"The Flux")
plt.savefig("reco_flux_not_muon.png",dpi=400)
plt.show()


cf = plt.contourf( energies/(1e9), czeniths, np.transpose(from_muon), cmap=cm.coolwarm, locator=ticker.LogLocator())
plt.xscale('log')
plt.xlabel("Event Energy [GeV]", size=14)
plt.ylabel(r"Reconstructed $\cos\theta$",size=14)
cbar = plt.colorbar(cf,ticks=ticker.LogLocator())
cbar.set_label(r"The Flux")
plt.savefig("reco_flux_from_muon.png",dpi=400)
plt.show()

