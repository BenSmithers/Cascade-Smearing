from MCEq.core import MCEqRun
import crflux.models as crf

import numpy as np
import nuSQUIDSpy as nsq

from cascade.raw_fluxes import get_key
from cascade.utils import get_closest
from cascade.cross_section_test import xs_obj as xs
from cascade.cross_section_test import flavors, currents, neut_types, get_total_flux

un = nsq.Const()

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt


def make_plots():
    Emin = 1*un.GeV
    Emax = 10*un.PeV

    energies = nsq.logspace(Emin, Emax, 100)
    angles = nsq.linspace(-0.99, -0.98, 2)

    # just look at straight up/down
    # get MCEq flux there

    mag = 0.
    angle = 0. #degrees

    mceq = MCEqRun(
        interaction_model = 'SIBYLL23C',
	primary_model = (crf.HillasGaisser2012, 'H3a'),
	theta_deg = 0.
    ) 

    mceq.solve()

    n_nu = 3
    inistate = np.zeros(shape=(2, len(energies), 2, n_nu))
    flux = {}

    flux['e_grid'] = mceq.e_grid
    flux['nue_flux'] = mceq.get_solution('nue',mag)+mceq.get_solution('pr_nue',mag)
    flux['nue_bar_flux'] = mceq.get_solution('antinue',mag)+mceq.get_solution('pr_antinue',mag)
    flux['numu_flux'] = mceq.get_solution('numu',mag)+mceq.get_solution('pr_numu',mag)
    flux['numu_bar_flux'] = mceq.get_solution('antinumu',mag)+mceq.get_solution('pr_antinumu',mag)
    flux['nutau_flux'] = mceq.get_solution('nutau',mag)+mceq.get_solution('pr_nutau',mag)
    flux['nutau_bar_flux'] = mceq.get_solution('antinutau',mag)+mceq.get_solution('pr_antinutau',mag)

    # propagate it with nusquids
    for neut_type in range(2):
        for flavor in range(n_nu):
            flav_key = get_key(flavor, neut_type)
            if flav_key =="":
                continue
            for energy_bin in range(len(energies)):
                inistate[0][energy_bin][neut_type][flavor] = get_closest(energies[energy_bin]/un.GeV, flux['e_grid'], flux[flav_key])

                inistate[1][energy_bin][neut_type][flavor] = get_closest(energies[energy_bin]/un.GeV, flux['e_grid'], flux[flav_key])

    nus_atm = nsq.nuSQUIDSAtm(angles, energies, n_nu, nsq.NeutrinoType.both, True)
  	
    nus_atm.Set_MixingAngle(0,1,0.563942)
    nus_atm.Set_MixingAngle(0,2,0.154085)
    nus_atm.Set_MixingAngle(1,2,0.785398)

    nus_atm.Set_SquareMassDifference(1,7.65e-05)
    nus_atm.Set_SquareMassDifference(2,0.00247)

    nus_atm.SetNeutrinoCrossSections(xs)

    #settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    #nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_FUNCTIONS.GSL_STEP_RK4)

    # we load in the initial state. Generating or Loading from a file 
    nus_atm.Set_initial_state(inistate, nsq.Basis.flavor)
    print("Done setting initial state")

    # we turn off the progress bar for jobs run on the cobalts 
    nus_atm.Set_ProgressBar(False)
    nus_atm.Set_IncludeOscillations(True)

    print("Evolving State")
    nus_atm.EvolveState()
   

    # Evaluate the flux at the energies, flavors, and shit we care about
  
    eval_flux = {}
    for flavor in range(3):
        for nute in range(2):
            key = get_key(flavor, nute)
            if key=="":
                continue
            eval_flux[key] = np.zeros(len(energies))

            for energy in range(len(energies)):
                powed = energies[energy]
                
                eval_flux[key][energy] = nus_atm.EvalFlavor(flavor, -0.985, powed, nute)

    # multiply the fluxes with total cross sections at that energy 

    
    conv_flux = {}
    for flavor in flavors.keys():
        if flavor=="electron":
            i_flav = 0
        elif flavor=="muon":
            i_flav = 1
        else:
            i_flav = 2 
        for neut_type in neut_types.keys():
            if neut_type=="neutrino":
                i_neut = 0
            else:
                i_neut = 1

            eval_key = get_key(i_flav, i_neut)
            for current in currents.keys():
                conv_key = "_".join([flavor, neut_type, current])
                conv_flux[conv_key] = np.zeros(len(energies))
                for energy in range(len(energies)):
                    powed = energies[energy]

                    conv_flux[conv_key][energy] = eval_flux[eval_key][energy]*get_total_flux(powed, flavors[flavor], neut_types[neut_type], currents[current])


    # sum up the cascades and compare them with the cross sections 


    casc_fluxes = np.zeros(len(energies))
    track_fluxes = np.zeros(len(energies))
    for key in conv_flux:
        if is_cascade(key):
            casc_fluxes += conv_flux[key]
        else:
            track_fluxes+= conv_flux[key]


    plt.plot(energies/un.GeV, casc_fluxes, label="Cascades")
    plt.plot(energies/un.GeV, track_fluxes, label="Tracks")
    plt.xscale('log')
    plt.xlim([10**2, 10**7])
    plt.yscale('log')
    plt.legend()
    plt.xlabel("Energy [GeV]", size=14)
    plt.ylabel(r"Flux$\cdot\sigma_{total}$ [s GeV sr]$^{-1}$",size=14)
    plt.show()

def is_cascade(key):
    if "_NC_" in key:
        return(True)
    else:
        if "electron" in key:
            return(True)
        else:
            return(False) # CC muon and Taus 


if __name__=="__main__":
    make_plots()
