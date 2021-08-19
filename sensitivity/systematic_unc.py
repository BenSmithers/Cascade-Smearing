"""
Several functions for calculating systematic sources of uncertainty


Each of these functions should return a tuple for the (-1 sigma, +1 sigma) shifts to the mean expectation 
"""

import os 
from shutil import move
import numpy as np

from cascade.utils import get_loc, Data, config
from cascade.utils import SterileParams
from cascade.raw_fluxes import raw_flux

from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.sensitivity.eff_area_reader import build_flux, build_flux_sad
from cascade.sensitivity.make_from_mc import build_mc_flux

import pickle 
import os
from math import log10

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

def _flipper(central,perturbed):
    """
    Takes a tuple with two numpy arrays. 
    
    """

    u_minus = np.zeros(shape = np.shape(perturbed[0]))
    u_plus  = np.zeros(shape = np.shape(perturbed[1]))

    for i in range(len(u_minus)):
        for j in range(len(u_minus[i])):
            p1 = central[i][j] - perturbed[0][i][j]
            p2 = central[i][j] - perturbed[1][i][j]

            if p1>0:
                u_plus[i][j]  = abs(p1)
                u_minus[i][j] = abs(p2)
            else:
                u_plus[i][j]  = abs(p2)
                u_minus[i][j] = abs(p1)
    return(u_minus, u_plus)
    

def _ice_grad(data_dict, grad_no):
    """
    Calculates the expected systematic uncertainty error in each of the bins from the ice gradients published in the MEOWS Paper
    """
    if not isinstance(data_dict, dict):
        raise TypeError("Expected {}, got {}".format(dict, type(data_dict)))

    new_uncertainty = np.zeros(shape=np.shape(data_dict["stat_err"]))


    energies = data_dict["e_edges"]
    cos_th = data_dict["a_edges"]

    ice_grads = np.loadtxt(os.path.join(os.path.dirname(__file__), "resources","icegrad.dat"))
    ice_e = np.logspace(log10(500), 4, 13)

    for e_i in range(len(energies)-1):
        for a_i in range(len(cos_th)-1):
            if energies[e_i]<ice_e[0]:
                ice_bin = 0
            elif energies[e_i]>ice_e[-1]:
                ice_bin = len(ice_e)-1
            else:
                ice_bin = get_loc( energies[e_i], ice_e )[0]
            new_uncertainty[e_i][a_i] = ice_grads[ice_bin][grad_no+1]/100.0

    p_pert = (1+new_uncertainty)*data_dict["event_rate"]
    m_pert = (1-new_uncertainty)*data_dict["event_rate"]

    return _flipper( data_dict["event_rate"] , (m_pert, p_pert))

def cr_perturb(dgamma = 0.0, dnorm=0.0, use_mc = False, smearmode=False, special=False):
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
    
    
    filename = os.path.join(perturb_folder, "cr_central.dat")
    null = SterileParams()
    
    if special:
        flux_func = lambda *objs:build_flux_sad(*objs, good_angles=True)
    else:
        if smearmode:
            flux_func = build_flux_sad
        else:
            flux_func = build_mc_flux if use_mc else build_flux

    if True: #not os.path.exists(filename):
        kwargs = {'as_data':True}
        flux_data = flux_func(raw_flux(null, kwargs))
        pickle_save(flux_data, filename)
    else:
        flux_data = pickle_load(filename)

    p_plus = np.zeros(shape=np.shape(flux_data["stat_err"]))
    p_minus = np.zeros(shape=np.shape(flux_data["stat_err"]))

    if perturb_norm:
        mean_norm = 1.19
        return _flipper( flux_data["event_rate"], ((mean_norm-dnorm)*flux_data["event_rate"], (mean_norm+dnorm)*flux_data["event_rate"]))
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
                p_plus[i_e][i_a] = flux_data["event_rate"][i_e][i_a]*pow(effective_e/(2.2e3), -mean_dgamma-dgamma)
                p_minus[i_e][i_a] = flux_data["event_rate"][i_e][i_a]*pow(effective_e/(2.2e3),-mean_dgamma+dgamma)

        return _flipper(flux_data["event_rate"], (p_minus, p_plus)) 



def astro_norm_unc(use_mc = False, smearmode=False,special=False):
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

    prefix = "astr_perturb" + ("_frommc" if use_mc else "") 
    prefix+=("_smear" if smearmode else "")
    if special:
        prefix += "_smearedwell"

    # for each of these... 
    #    generate the flux
    #    load it in as a data object 
    #    parse it through the effective area integrator
    # yay now we have fluxes in reconstruction space! 
    
    if special:
        flux_func = lambda *objs:build_flux_sad(*objs, good_angles=True)
    else:
        if smearmode:
            flux_func = build_flux_sad
        else:
            flux_func = build_mc_flux if use_mc else build_flux

    filename = os.path.join(perturb_folder, prefix+ "_central.dat")
    if os.path.exists(filename):
        central = pickle_load(filename)
    else:
        central = flux_func(generate_astr_flux(null,as_data=True))
        pickle_save(central, filename) 

    filename = os.path.join(perturb_folder, prefix+ "_norm_minus.dat")
    if os.path.exists(filename):
        norm_minus = pickle_load(filename)
    else:
        norm_minus = flux_func(generate_astr_flux(null, norm=-1*norm_p,as_data=True))
        pickle_save(norm_minus, filename)

    filename = os.path.join(perturb_folder, prefix+"_norm_plus.dat")
    if os.path.exists(filename):
        norm_plus = pickle_load(filename)
    else:
        norm_plus = flux_func(generate_astr_flux(null, norm=norm_p,as_data=True))
        pickle_save(norm_plus, filename)

    return _flipper( central["event_rate"], (norm_minus["event_rate"], norm_plus["event_rate"]))
#    return(norm_minus["event_rate"]-central["event_rate"], norm_plus["event_rate"]-central["event_rate"])

def astro_shift_unc(use_mc=False, smearmode=False, special=False):
    null = SterileParams()

    norm_p = 0.21
    shift_p = 0.19

#    unperturbed = generate_astr_flux(null)

    prefix = "astr_perturb" + ("_frommc" if use_mc else "") 
    if smearmode:
        prefix += "_smear"
    if special:
        prefix += "_smearedwell"

    flux_func = build_mc_flux if use_mc else build_flux
    if special:
        flux_func = lambda *objs:build_flux_sad(*objs, good_angles=True)
    elif smearmode:
        flux_func = build_flux_sad

    # for each of these... 
    #    generate the flux
    #    load it in as a data object 
    #    parse it through the effective area integrator
    # now we have fluxes in reconstruction space! 

    filename = os.path.join(perturb_folder, prefix+ "_central.dat")
    if os.path.exists(filename):
        central = pickle_load(filename)
    else:
        central = flux_func(generate_astr_flux(null,as_data=True))
        pickle_save(central, filename) 


    filename = os.path.join(perturb_folder, prefix+ "_gamma_minus.dat")
    if os.path.exists(filename):
        norm_minus = pickle_load(filename)
    else:
        norm_minus = flux_func(generate_astr_flux(null, dgamma=-1*shift_p,as_data=True))
        pickle_save(norm_minus, filename)

    filename = os.path.join(perturb_folder, prefix+"_gamma_plus.dat")
    if os.path.exists(filename):
        norm_plus = pickle_load(filename)
    else:
        norm_plus = flux_func(generate_astr_flux(null, dgamma=shift_p,as_data=True))
        pickle_save(norm_plus, filename)


    return _flipper(central["event_rate"], (norm_minus["event_rate"], norm_plus["event_rate"]))


