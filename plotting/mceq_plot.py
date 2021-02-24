"""
This plotting script makes two plots 
 - angle-integrated MCEq neutrino flux at atmosphere, separated by neutrino type
 - 2D heatmap of electron flux (cos(theta) vs energy)
"""
from cascade.utils import config 

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm

from math import acos, sin, pi

import numpy as np
import os
import pickle 

filename = os.path.join(config["datapath"], config["mceq_flux"])

if not os.path.exists(filename):
    raise IOError("It appears that the mceq data file has not been generated. Use the 'raw_fluxes.py' script to do that")

f = open(filename, 'rb')
inistate = pickle.load(f)
f.close()

n_e = 121
n_z = 50

energies = np.logspace( 1, 7, n_e)
zeniths = np.logspace(-0.999, 0., n_z)

z_width = sin(acos(0.999/n_z))*2*pi


e_flux = np.zeros(shape=(n_e, n_z))
mu_flux = np.zeros(shape=(n_e, n_z))
tau_flux = np.zeros(shape=(n_e, n_z))

for i_e in range(n_e):
    for i_z in range(n_z):
        for nutype in range(2):
            e_flux[i_e][i_z] += inistate[i_z][i_e][nutype][0]
            mu_flux[i_e][i_z] += inistate[i_z][i_e][nutype][1]
            tau_flux[i_e][i_z] += inistate[i_z][i_e][nutype][2]


# look at the energy distributions
e_sum = np.sum(e_flux, 1)*z_width
mu_sum = np.sum(mu_flux, 1)*z_width
tau_sum = np.sum(tau_flux, 1)*z_width

plt.plot(energies, e_sum, label=r"$\nu_{e}+\bar{\nu}_{e}$")
plt.plot(energies, mu_sum, label=r"$\nu_{\mu}+\bar{\nu}_{\mu}$")
plt.plot(energies, tau_sum, label=r"$\nu_{\tau}+\bar{\nu}_{\tau}$")
plt.xlabel("Energy [GeV]", size=14)
plt.ylabel(r"$\Phi_{\nu}$ [GeV cm$^{2}$ s]", size=14)
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.legend()
plt.show()
plt.clf()

cf = plt.pcolormesh(zeniths, energies, e_flux, cmap=cm.viridis, norm=colors.LogNorm())
cbar = plt.colorbar(cf)
#e_flux = all_data["nue_flux"] + all_data["nue_bar_flux"]
#mu_flux = all_data["numu_flux"] + all_data["numu_bar_flux"]
#tau_flux = all_data["nutau_flux"] + all_data["nutau_bar_flux"]

#plt.plot(energies, e_flux, label=r"$\nu_e$")
#plt.plot(energies, mu_flux, label=r"$\nu_{\mu}$")
#plt.plot(energies, tau_flux, label=r"$\nu_{\tau}$")
cbar.set_label(r"$\Phi_{e}$ [GeV cm$^{2}$ s sr]",size=14)

plt.ylabel(r"$E_{\nu}$ [GeV]",size=14)
plt.xlabel(r"$\cos\theta$", size=14)
plt.yscale('log')
plt.show()
