"""
This plotting script makes two plots 
 - angle-integrated MCEq neutrino flux at atmosphere, separated by neutrino type
 - 2D heatmap of electron flux (cos(theta) vs energy)
"""
from cascade.utils import config 

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm

from math import acos, sin, pi

import numpy as np
import os
import pickle 
import sys

from cascade.utils import get_color

folders = ["/home/benito/software/data/cascade/hg_sib/", "/home/benito/software/data/cascade/poly_sib/"]

filename = os.path.join(folders[0], config["mceq_flux"])

if not os.path.exists(filename):
    raise IOError("It appears that the mceq data file has not been generated. Use the 'raw_fluxes.py' script to do that")

n_e = 121
n_z = 50

energies = np.logspace( 1, 7, n_e)
zeniths = np.logspace(-0.999, 0., n_z)

z_width = sin(acos(0.999/n_z))*2*pi


contour = False

def get_flux(filename):
    f = open(filename, 'rb')
    inistate = pickle.load(f)
    f.close()

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
    
    if contour:
        return( e_flux, mu_flux, tau_flux)
    else:
        return(e_sum, mu_sum, tau_sum)

def get_muon_flux(filename):
    f = open(filename, 'rb')
    inistate = pickle.load(f)
    f.close()

    mu_flux = np.zeros(shape=(n_e, n_z))
    mubar_flux = np.zeros(shape=(n_e, n_z))

    for i_e in range(n_e):
        for i_z in range(n_z):
            mu_flux[i_e][i_z] += inistate[i_z][i_e][0][1]
            mubar_flux[i_e][i_z] += inistate[i_z][i_e][1][1]
    
    mu_flux = np.sum(mu_flux,1)*z_width
    mubar_flux = np.sum(mubar_flux,1)*z_width

    return(mu_flux, mubar_flux)

e_flux, mu_flux, tau_flux = get_flux(filename)
pe_flux, pmu_flux, ptau_flux = get_flux(os.path.join(folders[1], config["mceq_flux"]))


if not contour:
    plt.plot(energies, e_flux, label=r"HG2012 $\nu_{e}$", ls='-',color=get_color(0,3,'plasma'))
    plt.plot(energies, mu_flux, label=r"HG2012 $\nu_{\mu}$", ls='-', color=get_color(1,3, 'plasma'))
    plt.plot(energies, tau_flux, label=r"HG2012 $\nu_{\tau}$", ls='-',color=get_color(2,3,'plasma'))

    plt.plot(energies, pe_flux, label=r"Poly $\nu_{e}$", ls='--', color=get_color(0,3,'plasma'))
    plt.plot(energies, pmu_flux, label=r"Poly $\nu_{\mu}$", ls='--', color=get_color(1,3,'plasma'))
    plt.plot(energies, ptau_flux, label=r"Poly $\nu_{\tau}$", ls='--', color=get_color(2,3,'plasma'))
    plt.xlabel("Energy [GeV]", size=16)
    #plt.ylabel(r"Flux/$\Phi_{\nu_{\mu}}$", size=14)
    plt.ylabel(r"$E^{2}\cdot \Phi$ [GeV cm$^{-2}$ s$^{-1}$] ", size=16)
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.legend(prop={'size':14})
    plt.show()
    plt.clf()

    sys.exit()

    print(filename)
    mu_flux, mubar_flux = get_muon_flux(filename)
    pmu_flux, pmubar_flux = get_muon_flux(os.path.join(folders[1], config["mceq_flux"]))

    plt.plot(energies, pow(energies,2)*mu_flux, label=r"Gaisser $\nu_{\mu}$", ls='-',color=get_color(0,2))
    plt.plot(energies, pow(energies,2)*mubar_flux, label=r"Gaisser $\bar{\nu}_{\mu}$", ls='-', color=get_color(2,2))
    plt.plot(energies, pow(energies,2)*pmu_flux, label=r"Poly $\nu_{\mu}$", ls='--', color=get_color(0,2))
    plt.plot(energies, pow(energies,2)*pmubar_flux, label=r"Poly $\bar{\nu}_{\mu}$", ls='--', color=get_color(2,2))
    plt.xlabel("Energy [GeV]", size=14)
    plt.ylabel(r"$E^{2}\cdot \Phi$ [GeV cm$^{-2}$ s$^{-1}$] ", size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.tight_layout()
    plt.legend()
    plt.show()
    plt.clf()

    plt.plot(energies, mu_flux/mubar_flux, label=r"gaisser", color=get_color(0,2))
    plt.plot(energies, pmu_flux/pmubar_flux, label=r"polygonato", color=get_color(2,2))
    plt.xlabel("energy [gev]", size=14)
    plt.ylabel(r"$\nu_{\mu}/\bar{\nu}_{\mu}$ ratio", size=14)
    plt.xscale('log')
    plt.tight_layout()
    plt.legend()
    plt.show()
    plt.clf()

    plt.plot(energies, mu_flux/pmu_flux, label=r"$\nu_{\mu}$", color=get_color(0,2))
    plt.plot(energies, mubar_flux/pmubar_flux, label=r"$\bar{\nu}_{\mu}$", color=get_color(2,2))
    plt.xlabel("energy [gev]", size=14)
    plt.ylabel(r"Gaisser/Poly ratio", size=14)
    plt.ylim([0,2])
    plt.xscale('log')
    plt.tight_layout()
    plt.legend()
    plt.show()
    plt.clf()


if contour:
    #cf = plt.pcolormesh(zeniths, energies, mu_flux/pmu_flux, cmap=cm.coolwarm,vmin=0.8,vmax=1.2)
    cf = plt.pcolormesh(zeniths, energies, (e_flux+mu_flux+tau_flux)/(pe_flux+pmu_flux+ptau_flux), cmap=cm.coolwarm,vmin=0.5,vmax=1.5)
    cbar = plt.colorbar(cf)
    #e_flux = all_data["nue_flux"] + all_data["nue_bar_flux"]
    #mu_flux = all_data["numu_flux"] + all_data["numu_bar_flux"]
    #tau_flux = all_data["nutau_flux"] + all_data["nutau_bar_flux"]

    #plt.plot(energies, e_flux, label=r"$\nu_e$")
    #plt.plot(energies, mu_flux, label=r"$\nu_{\mu}$")
    #plt.plot(energies, tau_flux, label=r"$\nu_{\tau}$")
    cbar.set_label(r"Gaisser/Polygonato",size=14)

    plt.ylabel(r"$E_{\nu}$ [GeV]",size=14)
    plt.xlabel(r"$\cos\theta$", size=14)
    plt.yscale('log')
    plt.show()
