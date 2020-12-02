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

if True:
    e_reco, a_reco, flux = _load_flux(filename)

    energies = np.array(bhist([e_reco]).centers)
    czeniths = np.array(bhist([a_reco]).centers)

    from_muon, from_not = sep_by_flavor(flux)

    cf = plt.pcolormesh(czeniths, energies/(1e9),  np.log10(from_not), cmap=cm.coolwarm) #, locator=ticker.LogLocator())
    plt.yscale('log')
    plt.ylabel("Event Energy [GeV]", size=14)
    plt.xlabel(r"Reconstructed $\cos\theta$",size=14)
    plt.title("Not Muons",size=14)
    cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$\log\Phi$ [GeV sr s cm$^{2}$]$^{-1}$")
    plt.savefig("reco_flux_not_muon.png",dpi=400)
    plt.show()


    cf = plt.pcolormesh(czeniths,energies/(1e9), np.log10(from_muon), cmap=cm.coolwarm) #, locator=ticker.LogLocator())
    plt.yscale('log')
    plt.ylabel("Event Energy [GeV]", size=14)
    plt.xlabel(r"Reconstructed $\cos\theta$",size=14)
    plt.title("Muons",size=14)
    cbar = plt.colorbar()# cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$\log\Phi$ [GeV sr s cm$^{2}$]$^{-1}$")
    plt.savefig("reco_flux_from_muon.png",dpi=400)
    plt.show()
else:
    e_reco, a_reco, flux = _load_flux(filename)

    flux = flux["Mu_nuBar_NC"]

    flux = np.ma.masked_where(flux<=0., flux)    

    energies = np.array(bhist([e_reco]).centers)
    czeniths = np.array(bhist([a_reco]).centers)

#    flux *= energies**2

#    from_muon, from_not = sep_by_flavor(flux)

    cf = plt.pcolormesh( czeniths, energies/(1e9), np.log10(flux), cmap=cm.coolwarm) #, locator=ticker.LogLocator())
    plt.yscale('log')
    plt.ylim([10**2, 10**5])
    plt.title("Muon-NuBar")
    plt.ylabel("Event Energy [GeV]", size=14)
    plt.xlabel(r"$\cos\theta_{z}$",size=14)
    cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
    cbar.set_label(r"$\log\Phi$ [GeV sr s cm$^{2}$]$^{-1}$")
    plt.savefig("reco_flux_not_muon.png",dpi=400)
    plt.show()


