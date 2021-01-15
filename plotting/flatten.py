import numpy as np
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

import pickle

from cascade.utils import config
from cascade.utils import gen_filename 
from cascade.utils import bhist 

def _load_flux(name):
    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_reco"]
    a_reco = all_data["a_reco"]
    flux = all_data["flux"]

    return( e_reco, a_reco, flux )

def _load_error(name):
    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_reco"]
    a_reco = all_data["a_reco"]
    error = all_data["error"]

    return( e_reco, a_reco, error )

e_reco, a_reco, kflux = _load_flux(gen_filename(config["datapath"], config["recon_flux"]+".dat", 0.1339, 0.0, 1.3))

e_reco, a_reco, kerror = _load_error(gen_filename(config["datapath"], config["flux_error"]+".dat", 0.1339, 0.0, 1.3))


flux = sum(kflux.values())
error = sum(kerror.values())

angle_widths = bhist([np.arccos(a_reco)]).widths
angle_centers = bhist([np.arccos(a_reco)]).centers
energy_centers = np.array(bhist([e_reco]).centers)
# first dim of flux is for energy 
# flux[energy][angle]

print(np.shape(flux))

summed_flux = np.array([ sum(flux[i_energy]*angle_widths) for i_energy in range(len(energy_centers)) ])
summed_error = np.array([ sum(error[i_energy]*angle_widths) for i_energy in range(len(energy_centers)) ])



plt.errorbar( x=energy_centers/(1e9), y=summed_flux, yerr=summed_error, capsize=5)
plt.xlabel("Energy [GeV]",size=14)
plt.ylabel("Flux", size=14)
plt.xscale('log')
plt.yscale('log')
plt.savefig("flattened.png",dpi=400)
plt.show()
