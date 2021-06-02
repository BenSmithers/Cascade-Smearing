"""
This script implements a function for setting up the astrophysical neutrino flux at the Earth's surface
"""
import os
import numpy as np
import pickle

from cascade.utils import SterileParams, gen_filename, config
from cascade.raw_fluxes import raw_flux
from cascade.oscillations.astro_mixing import get_flavor_ratio



def astro_initial_state(energies, zeniths, n_nu, kwargs):
    """
    This creates the initial state array that is sent to nuSQuIDS
    It assumes a 1-1-1 flavor ratio 
    """

    path = os.path.join(config["datapath"], config["astro_flux"])

    if os.path.exists(path):
        print("Loading Astro Flux")
        f = open(path, 'rb')
        inistate = pickle.load(f)
        f.close()
    else:
        print("Generating Astro Flux")
        e_bins = len(energies)
        a_bins = len(zeniths)

        inistate = np.zeros(shape=(a_bins, e_bins, 2, n_nu))
        # Using Analysis I fits
        # https://arxiv.org/pdf/2005.12943.pdf
        norm = 0.95 # I think this is [Gev sr cm2 s]^-1
        pivot = 100*(1e12) #TeV 
        gamma = -2.5 #unitless        
        if "norm" in kwargs:
            norm = (norm + kwargs["norm"])*(1e-18)
            print("Setting Norm to {}".format(norm))
        else:
            norm = norm*(1e-18)
        if "pivot" in kwargs:
            pivot = kwargs["pivot"]
            print("Setting Pivot to {}".format(pivot))
        if "dgamma" in kwargs:
            gamma = gamma + 0.11 + kwargs["dgamma"]
            print("Setting gamma to {}".format(gamma))
        else:
            gamma = gamma + 0.11
        def get_flux(energy):
            return norm*(energy/pivot)**(gamma)

        for i_e in range(e_bins):
            flux = get_flux(energies[i_e])
            for flavor in range(n_nu):
                for i_a in range(a_bins):
                    for neut_type in range(2):
                        inistate[i_a][i_e][neut_type][flavor] += flux*kwargs["flavor_ratio"][flavor]
        f = open(path, 'wb')
        pickle.dump(inistate, f, -1)
        f.close()
    return inistate 

def generate_astr_flux(params, **kwargs):
    """
    This script generates an astrophysical flux at the given SterileParameters

    It saves the flux to the datapath that the configuration file is configured to
    """
    print("Received {} point".format(params))
    flavor_ratio = get_flavor_ratio(params)
    forced_filename = gen_filename(config["datapath"], config["prop_astro_flux"], params)

    kwargs["flavor_ratio"] = flavor_ratio
    kwargs["state_setter"] = astro_initial_state
    kwargs["forced_filename"] =  forced_filename

    return raw_flux(params, kwargs)
