from cascade.utils import SterileParams, gen_filename, config, bhist
from cascade.cross_section_test import get_total_flux as get_xs
from cascade.nus_utils import get_flavor, get_neut, get_curr

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import pickle

from math import pi

def _load_flux(name):
    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    e_reco = all_data["e_true"]
    a_reco = all_data["a_true"]
    flux = all_data["flux"]

    return(e_reco, a_reco, flux)

e_reco, a_reco, flux = _load_flux(gen_filename(config["datapath"], config["nu_flux_downsize"]+".dat", SterileParams()))

energies = bhist([e_reco]).centers
a_widths = bhist([a_reco]).widths
angles = bhist([a_reco]).centers

keys = list(flux.keys())


def is_track(key):

    curr = key.split("_")[2].lower()
    if "nc"==curr:
        return(False)
    elif "cc"==curr:
        flavor = key.split("_")[0].lower()
        if flavor=="mu":
            return(True)
        elif flavor=="e":
            return(False)
        elif flavor=="tau":
            return(False)
        else:
            raise ValueError("Not sure what to do with {}".format(flavor))

    else:
        raise ValueError("Not sure what {} is".format(curr))
        
cascade_rate = np.zeros(shape=np.shape(flux[keys[0]])[0])
track_rate = np.zeros(shape=np.shape(flux[keys[0]])[0])


eff_width = (max(a_reco)-min(a_reco))/(2*len(a_reco))
widths_rad = [abs(np.arccos(ang+eff_width)-np.arccos(ang-eff_width)) for ang in np.linspace(min(a_reco), max(a_reco),len(a_reco))]

sterr = widths_rad*np.sin(np.arccos(a_reco))*2*pi

for key in keys:
    flav = get_flavor(key)
    neut = get_neut(key)
    curr = get_curr(key)

    for i_energy in range(len(energies)):
        xs = get_xs(energies[i_energy], flav, neut, curr)
        amount = sum(flux[key][i_energy]*sterr*xs)

        if is_track(key):
            track_rate[i_energy] += amount
        else:
            cascade_rate[i_energy] += amount

for i_energy in range(len(energies)):
    flav = get_flavor("Mu_nu_CC")
    curr = get_curr("Mu_nu_CC")

    neut_nu = get_neut("Mu_nu_CC")
    neut_nubar = get_neut("Mu_nuBar_CC")

    xs_nu = get_xs(energies[i_energy], flav, neut_nu, curr)
    xs_nubar = get_xs(energies[i_energy], flav, neut_nubar, curr)

    amount_nu = sum(flux["Mu_nu_NC"][i_energy]*sterr*xs_nu)
    amount_nubar=sum(flux["Mu_nuBar_NC"][i_energy]*sterr*xs_nubar)

    track_rate[i_energy]+= amount_nu+amount_nubar

volume = (1e3)**3
nucleon_density = 0.9168*pow(100,3)*(6.02e23)

plt.plot(e_reco/(1e9), track_rate*volume*nucleon_density, label="Track Rate")
plt.plot(e_reco/(1e9), cascade_rate*volume*nucleon_density, label="Cascade Rate")
plt.xscale('log')
plt.xlabel("Energy [GeV]", size=16)
plt.ylabel(r"Event Rate [GeV s]$^{-1}$",size=16)
plt.yscale('log')
plt.xlim([10**2, 10**6])
plt.legend()
plt.tight_layout()
plt.show()
