from cascade.utils import Data

import os
from math import log10

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from cascade.utils import config, gen_filename

null_pt = gen_filename(config["datapath"], config["nu_flux"], 0.,0.,0.)
sterile_pt = gen_filename(config["datapath"], config["nu_flux"],0.13388166, 0.0, 1.3)
null_dat = Data(null_pt, 3)
sterile_dat = Data(sterile_pt, 4)

n_bin = 100

_ang_range = (min(null_dat.angles), max(null_dat.angles))
_eng_range = (min(null_dat.energies), max(null_dat.energies))

angles = np.linspace(_ang_range[0], _ang_range[1], n_bin+1)
energies = np.logspace( log10(_eng_range[0]), log10(_eng_range[1]), n_bin)

#angles = null_dat._angles
#energies = np.array(null_dat._energies)

key_list = null_dat.get_keys()

#null_flux = np.zeros(shape=np.shape(null_dat._fluxes[key_list[0]]))
#sterile_flux = np.zeros(shape=np.shape(sterile_dat._fluxes[key_list[0]]))

null_flux=np.zeros(shape=(n_bin, n_bin+1))
sterile_flux=np.zeros(shape=(n_bin,n_bin+1))

for key in key_list:
    for i in range(len(energies)):
        for j in range(len(angles)):
        
            energy = energies[i]
            angle = angles[j]
            
#            null_flux[i][j] += null_dat._fluxes[key][i][j]
#            sterile_flux[i][j] += sterile_dat._fluxes[key][i][j]

            null_flux[i][j]     += null_dat.get_flux( energy, key, angle=angle)
            sterile_flux[i][j]  += sterile_dat.get_flux( energy, key, angle=angle)

cf = plt.pcolormesh(angles, energies/(1e9), np.log10(null_flux), cmap=cm.viridis)
plt.yscale('log')
plt.ylabel("True Energy [GeV]", size=14)
plt.xlabel(r"True $\cos\theta$",size=14)
cbar = plt.colorbar()
cbar.set_label("log( Null Flux )")
plt.savefig("null_flux_plot.png", dpi=400)
plt.show()
plt.clf()

cf = plt.pcolormesh(angles, energies/(1e9), np.log10(sterile_flux), cmap=cm.viridis)
plt.yscale('log')
plt.ylabel("True Energy [GeV]", size=14)
plt.xlabel(r"True $\cos\theta$",size=14)
cbar = plt.colorbar()
cbar.set_label("log( Sterile Flux )")
plt.savefig("sterile_flux_plot.png", dpi=400)
plt.show()
plt.clf()

cf = plt.pcolormesh(angles, energies/(1e9), sterile_flux / null_flux, cmap=cm.coolwarm, vmin=0.90, vmax=1.10)
plt.yscale('log')
plt.ylabel("True Energy [GeV]", size=14)
plt.xlabel(r"True $\cos\theta$",size=14)
cbar = plt.colorbar()
cbar.set_label("sterile / null")
plt.savefig("ratio_flux_plot.png", dpi=400)
plt.show()
plt.clf()
