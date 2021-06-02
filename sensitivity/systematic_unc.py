"""
Several functions for calculating systematic sources of uncertainty


Each of these functions should return a tuple for the (-1 sigma, +1 sigma) shifts to the mean expectation 
"""

import numpy as np

from cascade.utils import get_loc, Data
from cascade.utils import SterileParams
from cascade.raw_fluxes import raw_flux

from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.sensitivity.eff_area_reader import build_flux

import pickle 
import os

perturb_folder = os.path.join(config["datapath"], "perturbed_rates")


def pickle_save(obj, filename):
    f = open(filename, 'wb')
    pickle.dump(obj, f, -1)
    f.close()
def pickle_load(filename):
    f = open(filename,'rb')
    obj = pickle.load(f)
    f.close()
    return obj

def ice_grad_0(data_dict):
    return _ice_grad(data_dict, 0)

def ice_grad_1(data_dict):
    return _ice_grad(data_dict, 1)

def _ice_grad(data_dict, grad_no):
    """
    Calculates the expected systematic uncertainty error in each of the bins from the ice gradients published in the MEOWS Paper
    """
    if not isinstance(data_dict, dict):
        raise TypeError("Expected {}, got {}".format(dict, type(data_dict)))

    new_uncertainty = np.zeros(shape=np.shape(data_dict["err"]))


    energies = data_dict["e_edges"]
    cos_th = data_dict["a_edges"]

    ice_grads = np.loadtxt(os.path.join(os.path.dirname(__file__), "resources","icegrad.dat"))
    ice_e = np.logspace(log10(500), 4, 13)

    for e_i in range(len(energies)-1):
        for a_i in range(len(cos_th)-1):
            ice_bin = get_loc( energies[e_i], ice_e )[0]
            
            new_uncertainty[e_i][a_i] = data_dict["event_rate"]*ice_grads[ice_bin][grad_no+1]/100.0
    return (-1*new_uncertainty, new_uncertainty)

def cr_perturb(dgamma = 0.0, dnorm=0.0):
    """
    Calculates the expected cosmic ray spectrum, then returns the expected gains/losses in each bin by perturbing the overall spectram index (dgamma) and normalization (dnorm) down/up 

    Returns a tuple containing np.ndarrays 
        The first is the event change from a negative perturbation
        The second is the event change from a positive perturbation
    """
    perturb_gamma = dgamma!=0.0
    perturb_norm = dnorm!=0.0

    if perturb_gamma == perturb_norm:
        raise ValueError("Cannot perturb both norm and gamma simultaneously. Also need to perturb at least one.")
    
    
    null = SterileParams()
    filename = os.path.join(perturb_folder, "cr_central.dat")
    if not os.path.exists(filename):
        flux_data = build_flux(Data(raw_flux(null)))
        pickle_save(flux_data, filename)
    else:
        flux_data = pickle_load(filename)


    p_plus = np.zeros(shape=np.shape(flux_data["stat_err"]))
    p_minus = np.zeros(shape=np.shape(flux_data["stat_err"]))

    if perturb_norm:
        mean_norm = 1.19
        return( (mean_norm-dnorm)*flux_data["event_rate"], (mean_norm+dnorm)*flux_data["event_rate"] )
    if perturb_gamma:
        # scale with 
        #phi(e) * (E/(2.2TeV))^deltagamma

        e_edges = flux_data["e_edges"]
        for i_e in range(len(e_edges)-1):
            for i_a in range(len(flux_data["a_edges"])-1):
                # we have the energies at the bin edges...
                # sooo
                mean_dgamma = 0.068
                effective_e = 0.5*(e_edges[i_e] + e_edges[i_e+1])
                p_plus[i_e][i_a] = flux_data["event_rate"][i_e][i_a]*(pow(effective_e/(2.2e3), mean_dgamma+dgamma) -1)
                p_minus[i_e][i_a] = flux_data["event_rate"][i_e][i_a]*(pow(effective_e/(2.2e3), mean_dgamma-dgamma) -1)
    
        return (p_minus - flux_data["event_rate"], p_plus-flux_data["event_rate"])



def astro_norm_unc(data_dict):
    """
    Calculates the expected astrophysical neutrino spectrum, then returns the expected gains/losses in each bin by perturbing the overall spectram index (dgamma) and normalization (dnorm) down/up 

    Returns a tuple containing np.ndarrays 
        The first is the event change from a negative perturbation
        The second is the event change from a positive perturbation
    """
    null = SterileParams()

    norm_p = 0.21
    shift_p = 0.19

#    unperturbed = generate_astr_flux(null)

    prefix = "astr_perturb"

    # for each of these... 
    #    generate the flux
    #    load it in as a data object 
    #    parse it through the effective area integrator
    # yay now we have fluxes in reconstruction space! 
    
    filename = os.path.join(perturb_folder, prefix+ "_central.dat")
    if os.path.exists(filename):
        central = pickle_load(filename)
    else:
        central = build_flux(Data(generate_astr_flux(null)))
        pickle_save(central, filename) 

    filename = os.path.join(perturb_folder, prefix+ "_norm_minus.dat")
    if os.path.exists(filename):
        norm_minus = pickle_load(filename)
    else:
        norm_minus = build_flux(Data(generate_astr_flux(null, norm=-1*norm_p)))
        pickle_save(norm_minus, filename)

    filename = os.path.join(perturb_folder, prefix+"_norm_plus.dat")
    if os.path.exists(filename):
        norm_plus = pickle_load(filename)
    else:
        norm_plus = build_flux(Data(generate_astr_flux(null, norm=norm_p)))
        pickle_save(norm_plus, filename)


    return(norm_minus["event_rate"]-central["event_rate"], norm_plus["event_rate"]-central["event_rate"])

def astro_shift_unc(data_dict):
    null = SterileParams()

    norm_p = 0.21
    shift_p = 0.19

#    unperturbed = generate_astr_flux(null)

    prefix = "astr_perturb"

    # for each of these... 
    #    generate the flux
    #    load it in as a data object 
    #    parse it through the effective area integrator
    # yay now we have fluxes in reconstruction space! 

    filename = os.path.join(perturb_folder, prefix+ "_central.dat")
    if os.path.exists(filename):
        central = pickle_load(filename)
    else:
        central = build_flux(Data(generate_astr_flux(null)))
        pickle_save(central, filename) 


    filename = os.path.join(perturb_folder, prefix+ "_gamma_minus.dat")
    if os.path.exists(filename):
        norm_minus = pickle_load(filename)
    else:
        norm_minus = build_flux(Data(generate_astr_flux(null, dgamma=-1*shift_p)))
        pickle_save(norm_minus, filename)

    filename = os.path.join(perturb_folder, prefix+"_gamma_plus.dat")
    if os.path.exists(filename):
        norm_plus = pickle_load(filename)
    else:
        norm_plus = build_flux(Data(generate_astr_flux(null, dgamma=shift_p)))
        pickle_save(norm_plus, filename)


    return(norm_minus["event_rate"], norm_plus["event_rate"])


