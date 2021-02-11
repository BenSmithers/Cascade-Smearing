"""
This script here generates two fluxes at just the energy-deposition stage 

Then we sum them over true energy to get fluxes at each deposited energy/angle 
Then we plot each one individually
Then we plot their ratio
"""
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# these are like the MOST important things, damn
from cascade.utils import bhist, SterileParams, gen_filename, config
from cascade.deposit import generate_singly_diff_fluxes

null = SterileParams(0.,0.,0.,0.)
ster = SterileParams(0., 0.1609, 0.2296, 4.47)

raw_null = gen_filename(config["datapath"], config["nu_flux"]+".dat", null)
raw_ster = gen_filename(config["datapath"], config["nu_flux"]+".dat", ster)

n_bins = 40

true_e_edges, depo_e_edges, null_flux, czenith_edges, errs = generate_singly_diff_fluxes(n_bins, datafile=raw_null)
true_e_edges, depo_e_edges, ster_flux, czenith_edges, errs = generate_singly_diff_fluxes(n_bins, datafile=raw_ster)

true_e_widths = np.array(bhist([true_e_edges]).widths)

czeniths = np.array(bhist([czenith_edges]).centers)
energies = np.array(bhist([depo_e_edges]).centers)

keys = [str(key) for key in null_flux.keys()]

# the flux is saved like [deposited e][true e][angle]
# I want to sum over the event energies, so we're going to be writing this to a new shape 
# The new shape will involve the dimension of deposited energies, and czeniths 
base_shape = np.shape(null_flux[keys[0]])
new_shape = (base_shape[0], base_shape[2])

just_one = False
keep_key = "Mu"

null_new = np.zeros(shape=new_shape)
ster_new = np.zeros(shape=new_shape)

nulle = np.zeros(shape=(base_shape[0], base_shape[1]))
stere = np.zeros(shape=(base_shape[0], base_shape[1]))

nulle_cur = np.zeros(shape=new_shape)

for key in keys:
    if just_one and (keep_key not in key):
        print("Skip {}".format(key))
        continue

    # we subtract one from the lengths since the flux entries correspond 
    #       with bin WIDTHS and not bin EDGES
    null_new  += np.sum(null_flux[key], 1)
    ster_new += np.sum(ster_flux[key], 1)

    for i_e_depo in range(len(null_flux[key])-1):
        nulle_cur[i_e_depo] += null_flux[key][i_e_depo][10] 

        for i_e_true in range(len(null_flux[key][i_e_depo])-1):

            nulle[i_e_depo][i_e_true] += null_flux[key][i_e_depo][i_e_true][5]
            stere[i_e_depo][i_e_true] += ster_flux[key][i_e_depo][i_e_true][5]

# I still don't get why pcolormesh basically transposes things...


cf = plt.pcolormesh(czeniths, energies/(1e9), np.log10(nulle_cur))
plt.yscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"log$\Phi$")
plt.title("Single Slice!", size=14)
plt.show()
plt.clf()

cf = plt.pcolormesh(energies/(1e9), energies/(1e9), np.log10(nulle))
plt.yscale('log')
plt.xscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"Nu Energy [GeV]",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"log$\Phi$")
plt.show()
plt.clf()

cf = plt.pcolormesh(energies/(1e9), energies/(1e9), np.log10(stere))
plt.yscale('log')
plt.xscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"Nu Energy [GeV]",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"log$\Phi$")
plt.show()
plt.clf()

cf = plt.pcolormesh(czeniths, energies/(1e9), np.log10(null_new))
plt.yscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"log$\Phi$")
plt.show()
plt.clf()

cf = plt.pcolormesh(czeniths, energies/(1e9), np.log10(ster_new))
plt.yscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"log$\Phi$")
plt.show()
plt.clf()

cf = plt.pcolormesh(czeniths, energies/(1e9), ster_new/null_new, cmap=cm.coolwarm, vmin=0.9, vmax=1.1)
plt.yscale('log')
plt.ylabel("Depo Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_one:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"sterile / null")
plt.show()
plt.clf()
