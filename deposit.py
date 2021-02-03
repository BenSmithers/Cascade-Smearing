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



import nuSQUIDSpy as nsq

# specialty-made utility functions
from cascade.nus_utils import get_flavor, get_neut, get_curr
from cascade.utils import bhist, get_exp_std, get_width, get_nearest_entry_to
from cascade.utils import Data, get_index, get_loc, sci
from cascade.utils import config

from cascade.cross_section_test import get_diff_xs

# tau stuff
from cascade.tau_funcs import TauData

# reconstruction data
from cascade.deporeco import DataReco

const = nsq.Const()
# load the data using the default filename, 'atmosphere.dat'
tauData = TauData()

debug = False
glob_angle = None

# this block here is supposed to just plot all the raw fluxes 
if debug:
    #plotting imports
    import matplotlib
    # Need to use agg since Tk isn't on the cobalts??? 
    matplotlib.use('TkAgg', force=True)
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    from matplotlib import ticker #used for log-scale contourplots 

    # colormap thing
    cmap = plt.get_cmap('coolwarm')
    n_colors = 6
    def get_color(which, how_many=n_colors):
        return( cmap( float(which)/how_many ) )

    #data = Data("atmosphere_sterile.dat")
    data = Data("atmosphere_null.dat", 3)
    scale_e = np.array(data.energies)

    print("Making Debug plot")
    binny = len(scale_e)
    taus = np.zeros(binny)
    muon = np.zeros(binny)
    ele = np.zeros(binny)

    angle_bin = 0
    print("Angle: {}".format(data.angles[angle_bin]))

    for key in data.fluxes:
        flav=str(key).split('_')[0]
        if 'Tau'==flav:
            taus+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        elif "E"==flav:
            ele+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        elif "Mu"==flav:
            muon+=[data.fluxes[key][i][angle_bin] for i in range(len(data.fluxes[key]))]
        else:
            raise Exception("You might have done steriles? {} is unrecognized".format(flav))


    plt.plot( scale_e/const.GeV, muon, color=get_color(0,3),label="Muons")
    plt.plot(scale_e/const.GeV, taus, color=get_color(1,3),label="Taus")
    plt.plot(scale_e/const.GeV, ele, color=get_color(2,3), label="Electrons")
    plt.legend()
    plt.xlabel("Event Energy [GeV]",size=14)
    plt.ylabel("Flux [GeV sec cm$^{2}$ sr]$^{-1}$",size=14)
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig('raw_fluxes_{:.2f}.png'.format(glob_angle), dpi=400)
    print("saved raw fluxes")
    plt.close()


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
    
def _save_data(e_true, e_depo, flux, czenith):
    """
    Saves the generated data for use later. 
    """
    all_data = {"e_true": e_true,
                "e_depo": e_depo,
                "czenith": czenith,
                "flux": flux}
    f = open(savefile,'wb')
    pickle.dump( all_data, f, -1)
    f.close()   
    print("Data Saved!")



def do_for_key(event_edges,cascade_edges, key,data, angles=None):
    """
    This function takes the desired bin edges for the event energies and deposited energies along with the dictionary key corresponding to a specific combination of falvor, current, and neutrino type.

    It builds up the 2D flux array (singly differential), which it then returns along with a similar array but for flux uncertainties 
    """
    evt = bhist([event_edges])
    cas = bhist([cascade_edges])
    reco = bhist([cascade_edges])

    event_energies = evt.centers
    event_widths = evt.widths
    cascade_energies = cas.centers
    cascade_widths = cas.widths

    flav = key.split("_")[0]
    curr = key.split("_")[2]
    neut = key.split("_")[1]

    if angles is None:
        ang_list = [None]
    else:
        ang_list = bhist([angles]).centers

    if angles is None:
        flux = bhist((cascade_edges, event_edges))
        err = bhist((cascade_edges, event_edges))
    else:
        flux = bhist((cascade_edges, event_edges, angles))
        err = bhist((cascade_edges, event_edges, angles))

    for angle in ang_list:
        # need special Tau treatment 
        if curr=="CC":
            # deposit all the energy. Hadronic and Leptonic (event) contribute 
            # Since the energy always all gets deposited, we only need to do one loop!
            # So, for a given "deposited energy" (cascade_energy), we already know the total energy. 
            # Therefore we just get the total cross section * flux there... the units end up as [s GeV in sr]^-1 
            for cas_bin in range(len(cascade_energies)): 
                deposited_energy = cascade_energies[cas_bin] # energy going to the leptonic component! 
                xs = get_diff_xs(deposited_energy, get_flavor(key), get_neut(key), get_curr(key))*cascade_widths[cas_bin] # deposited (cascade) and event are the same here 

                amount =data.get_flux(deposited_energy,key, angle=angle)*xs                
                amount_err = data.get_err(energy=deposited_energy, key=key, angle=angle)*xs


                if flav.lower()=='tau':
                    # Etau is cascade_energies[cas_bin]
                    # How much energy is visible in the various tau decays though? 
                    # going from zero->deposited energy
                    if neut.lower()=="nubar":
                        deposited_energy = tauData(deposited_energy/const.GeV,-1)
                    else: # nu
                        deposited_energy = tauData(deposited_energy/const.GeV, 1)
                
                try:
                    if angle is None:
                        flux.register( amount, cascade_energies[cas_bin], deposited_energy) # add it in! 
                        err.register( amount_err, cascade_energies[cas_bin], deposited_energy)
                    else: 
                        flux.register( amount, cascade_energies[cas_bin], deposited_energy, angle)
                        err.register( amount_err, cascade_energies[cas_bin], deposited_energy, angle)
                except ValueError:
                    if flav.lower()!='tau':
                        raise Exception("It wasn't tau. Something's wrong")

        else:
            # in this case, knowing the cascade doesn't tell us anything about the event energy. 
            # so we loop over both, get the flux*differential_xs at each bin combination, and multiply by the widths of deposited-energy-bin to get the same units as in the CC case 
            for evt_bin in range(len(event_energies)):
                for cas_bin in range(len(cascade_energies)):
                    lepton_energy = event_energies[evt_bin] - cascade_energies[cas_bin]

                    # we'll have nowhere to put these, so let's just skip this
                    if lepton_energy < min(cascade_energies):
                        continue
                    if lepton_energy > max(cascade_energies):
                        continue
            
                    xs = get_diff_xs(event_energies[evt_bin], get_flavor(key), get_neut(key), get_curr(key), lepton_energy,0.0)*cascade_widths[cas_bin]*event_widths[evt_bin]

                    amount =data.get_flux(event_energies[evt_bin],key, angle=angle)*xs
                    amount_err = data.get_err(event_energies[evt_bin],key, angle=angle)*xs
                    if angle is None:
                        flux.register(amount, cascade_energies[cas_bin], event_energies[evt_bin])
                        err.register(amount_err, cascade_energies[cas_bin], event_energies[evt_bin])

                    else:
                        flux.register(amount, cascade_energies[cas_bin], event_energies[evt_bin], angle)
                        err.register(amount_err, cascade_energies[cas_bin], event_energies[evt_bin], angle)

    
    # build a new bhist in reconstruction space (Still with event energy too)
    # then scan through deposited-true angle space
    # and use the PDFS to smear the true values into reconstructed values, depositing them into the reconstruction bhist  

    return(flux.fill, err.fill)

def generate_singly_diff_fluxes(n_bins,debug=False, datafile=config["nu_flux"]):
    """
    This is the newest, most streamlined function for generating the singly-differential flux arrays. 
    It has the same functionality as the old ones, but now it skips the step of filling in the regular flux array before swapping to the deposited energy flux array. 

    It goes over the different kind of fluxes and builds up a 2D or 3D array, convolving the entries with the relevant differential cross section
        the dimensionality depends on whether or we are integrating over the zenith angles
    """
    data = Data(datafile)

    e_min = 10*const.GeV
    e_max = 1*const.PeV
    extra = 2
    
    all_angles = data.angles
    
    event_edges = np.logspace(np.log10(e_min), np.log10(e_max)+extra,n_bins+1)
    cascade_edges = np.logspace(np.log10(e_min), np.log10(e_max),n_bins+1)
    angle_edges = np.linspace(min(all_angles), max(all_angles), n_bins+1) # angle is in cos(zenith), so like -1->0
    
    nuflux = {}
    errs = {}

    for key in data.get_keys(): #flavor, current, interaction 
        fluxy, erry = do_for_key(event_edges,cascade_edges,key,data=data, angles=angle_edges)
        nuflux[key] = fluxy
        errs[key] = erry


    # if global variable "angle" isn't none, then we can separate out just a single angle

    return(event_edges,cascade_edges, nuflux, angle_edges, errs)

