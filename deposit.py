#!/usr/bin/python3.6
'''
This script plots the fluxes output by the raw_fluxes script. It also does some analysis and saves some arrays to disk 

But how does it work?
 1. Raw flux data from the included mceq+nuSQuIDS flux is loaded in by the Data Object. This object has functionality for sampling from the flux at arbitrary parent neutrino energy by interpolating between nsq points.
 2. This flux data is convolved with some cross sections to make a singly-differential flux array of deposited energy vs event energy 

The heavy lifting goes on in the generate_singly_diff_fluxes function 

Notes:
    + energies are generally in eV for consistency with nuSQuIDs
    + widths are in GeV for consistency with nuSQuIDS (why...)
    + the fluxes output by the Data object are in units of [GeV cm2 s]^-1 (incoming neutrino energy)
    + the cross sections are in units of cm2/GeV (outgoing neutrino energy)
'''

import sys

do_norm=False # deprecated 

# data analysis
import numpy as np

# file system, control
import os
import pickle
from warnings import warn
from math import log10

import nuSQUIDSpy as nsq

# specialty-made utility functions
from cascade.nus_utils import get_flavor, get_neut, get_curr
from cascade.utils import bhist, get_exp_std, get_width, get_nearest_entry_to
from cascade.utils import Data, get_index, get_loc, sci
from cascade.utils import config

from cascade.cross_section_test import get_diff_xs

# reconstruction data
from cascade.deporeco import DataReco

const = nsq.Const()
# load the data using the default filename, 'atmosphere.dat'

glob_angle = None

suffix=""
savefile = os.path.join( config["datapath"], config["all_fluxes"]+suffix+".dat")
def _load_data(glob_angle = None):
    """
    Loads the datafile. Returns tuple

    0 - parent energies
    1 - child energies
    2 - nuflux
    3 - cos(zenith) edges
    """
    print("Loading Data")
    f = open(savefile, 'rb')
    all_data = pickle.load(f)
    f.close()
    angle_edges = all_data["czenith"]
    nuflux = all_data["flux"]

    if glob_angle is not None:
        # figure out which slice to return
        # [:,:,N] would return the N'th angle bin
        if glob_angle<min(angle_edges) or glob_angle>max(angle_edges):
            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,0]*0.
        else:
            lower, upper = get_loc(glob_angle, angle_edges)

            print("Grabbed angle bin {}".format(lower))
            width = abs( np.arccos(angle_edges[upper]) - np.arccos(angle_edges[lower]))

            for key in nuflux.keys():
                nuflux[key] = nuflux[key][:,:,lower]*width
                

    return( all_data["e_true"], all_data["e_depo"], \
                nuflux, angle_edges )
    
def _save_data(filename, **kwargs):
    """
    Saves the generated data for use later. 
    """
    all_data = {}
    for kwarg in kwargs:
        all_data[kwarg] = kwargs[kwarg] # I hate this 

    f = open(savefile,'wb')
    pickle.dump( all_data, f, -1)
    f.close()   
    print("Data Saved!")



def do_for_key(event_edges,e_deposited_edges, key,data, angles):
    """
    This function takes the desired bin edges for the event energies and deposited energies along with the dictionary key corresponding to a specific combination of falvor, current, and neutrino type.

    It builds up the 2D flux array (singly differential), which it then returns along with a similar array but for flux uncertainties 
    """
    evt = bhist([event_edges])
    cas = bhist([e_deposited_edges])

    event_energies = evt.centers
    event_widths = evt.widths
    deposited_energies = cas.centers
    e_deposited_widths = cas.widths

    flav = key.split("_")[0]
    curr = key.split("_")[2]
    neut = key.split("_")[1]

    ang_list = bhist([angles]).centers

    #flux = bhist((e_deposited_edges, event_edges, angles))
    #err = bhist((e_deposited_edges, event_edges, angles))

    flux = np.zeros( shape=(len(e_deposited_edges)-1, len(event_edges)-1, len(angles)-1))
    err = np.zeros( shape=(len(e_deposited_edges)-1, len(event_edges)-1, len(angles)-1))

    for i_a in range(len(ang_list)):
        angle = ang_list[i_a]

        # in this case, knowing the cascade doesn't tell us anything about the event energy. 
        # so we loop over both, get the flux*differential_xs at each bin combination, and multiply by the widths of deposited-energy-bin to get the same units as in the CC case 
        for evt_bin in range(len(event_energies)):
            for dep_bin in range(len(deposited_energies)):
                deposited_energy = deposited_energies[dep_bin] #in hadronic shower
                lepton_energy = event_energies[evt_bin] - deposited_energies[dep_bin]

                if deposited_energy>event_energies[evt_bin]:
                    continue

                if curr=="CC":
                    if flav=="E":
                        scale = 1.0
                    elif flav=="Tau":
                        scale = 0.51 # this was manually calculated as the best value, regardless of the actual tau energy 
                    else:
                        continue

                    # in the charge current, some of the lepton's energy gets deposited too
                    # All of the electrons, and about half the tau's 
                    deposited_energy += scale*lepton_energy

                    try:
                        adj_dep_bin = get_loc( deposited_energy, e_deposited_edges )[0]
                    except ValueError:
                        continue
                else:
                    adj_dep_bin = dep_bin

                # we'll have nowhere to put these, so let's just skip this
                if deposited_energy < min(e_deposited_edges):
                    continue
                if deposited_energy > max(e_deposited_edges):
                    continue
        
                xs = get_diff_xs(event_energies[evt_bin], get_flavor(key), get_neut(key), get_curr(key), lepton_energy,0.0)*e_deposited_widths[dep_bin]*event_widths[evt_bin]

                amount =data.get_flux(event_energies[evt_bin],key, angle=angle)*xs*event_widths[evt_bin]
                amount_err = data.get_err(event_energies[evt_bin],key, angle=angle)*xs
                
                flux[adj_dep_bin][evt_bin][i_a] += amount/(e_deposited_widths[adj_dep_bin]*event_widths[evt_bin])
                err[adj_dep_bin][evt_bin][i_a] += amount_err/(e_deposited_widths[adj_dep_bin]*event_widths[evt_bin])

                #flux.register(amount, deposited_energy, event_energies[evt_bin], angle)
                #err.register(amount_err, deposited_energy, event_energies[evt_bin], angle)

    
    # build a new bhist in reconstruction space (Still with event energy too)
    # then scan through deposited-true angle space
    # and use the PDFS to smear the true values into reconstructed values, depositing them into the reconstruction bhist  

    return(flux, err)

def generate_singly_diff_fluxes(n_bins,datafile=config["nu_flux"]):
    """
    This is the newest, most streamlined function for generating the singly-differential flux arrays. 
    It has the same functionality as the old ones, but now it skips the step of filling in the regular flux array before swapping to the deposited energy flux array. 

    It goes over the different kind of fluxes and builds up a 2D or 3D array, convolving the entries with the relevant differential cross section
        the dimensionality depends on whether or we are integrating over the zenith angles
    """
    print("Depositing Energy")
    data = Data(datafile)

    e_min = 10*const.GeV
    e_max = 5*const.PeV
    extra = 0
    
    all_angles = data.angles
    
    event_edges       = np.logspace(np.log10(e_min), np.log10(e_max)+extra,n_bins+1)
    e_deposited_edges = np.logspace(np.log10(e_min)      , np.log10(e_max),n_bins+1)
    angle_edges = np.linspace(min(all_angles), max(all_angles), n_bins+1) # angle is in cos(zenith), so like -1->0
    
    nuflux = {}
    errs = {}

    for key in data.get_keys(): #flavor, current, interaction 
        fluxy, erry = do_for_key(event_edges,e_deposited_edges,key,data=data, angles=angle_edges)
        nuflux[key] = fluxy
        errs[key] = erry


    # if global variable "angle" isn't none, then we can separate out just a single angle

    return(event_edges,e_deposited_edges, nuflux, angle_edges, errs)

