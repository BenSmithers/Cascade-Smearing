#!/usr/bin/python3.6
'''
This takes the energy deposited/true array and converts it into a reconstructed/true one

This requires a lot of summing 
'''

from optparse import OptionParser
import sys
mode_str = "\nModes\n\
Debug - Raw fluxes for the three flavors. No cross sections\n \
4 - PLOT 2d HIST OF parent vs cascade, Ratio plot\n \
5 - The plot showing most probable energy for muon events and non-muon events\n \
6 - The plot just showing must probable energy for all of them\n \
8 - Makes two 2D histograms of parent vs cascade fluxes\n \
9 - Plot of Median event energy for all three flavors\n \
"

parser = OptionParser()
#mode_str = "0 - plot muon vs not muon\n1 - plot all the keys\n 2 - plot 2D hist of parent vs cascade"
parser.add_option("-d", "--debug", 
                dest="debug",
                default=False,
                action="store_true",
                help="activates debug options")
parser.add_option("-n", "--nbins",
                dest="n_bins",
                default=200,
                type=int,
                help="Number of bins to use for each axis")
options, args = parser.parse_args()
debug = options.debug
n_bins = options.n_bins

print("Configuration...")
print("    Using {} bins".format(n_bins))

# data analysis
import numpy as np

# file system, control
import os
import pickle
from warnings import warn

from cascade.cross_section_test import get_diff_xs

# specialty-made utility functions
from nus_utils import get_flavor, get_neut, get_curr
from utils import bhist, get_exp_std, get_width, get_nearest_entry_to
from utils import Data, get_index, get_loc, sci
from utils import config 

from deposit import generate_singly_diff_fluxes

# reconstruction data
from deporeco import DataReco
suffix= "_null"

savefile = os.path.join(config["datapath"], config["all_fluxes"]+suffix+".dat")
fluxfile = os.path.join(config["datapath"], config["recon_flux"]+suffix+".dat")

def _save_data(e_reco, e_true, a_reco, a_true, flux):
    """
    Saves the generated data for use later. 
    """
    all_data = {"e_reco": e_reco,
                "e_true": e_true,
                "a_reco": a_reco,
                "a_true": a_true,
                "flux": flux}
    f = open(savefile,'wb')
    pickle.dump( all_data, f, -1)
    f.close()   
    print("Saved {}".format(savefile))

def _save_flux(e_reco, a_reco, flux):
    all_data = {"e_reco": e_reco,
                "a_reco": a_reco,
                "flux": flux}
    f = open(fluxfile, 'wb')
    pickle.dump( all_data, f, -1)
    f.close()
    print("Saved {}".format(fluxfile))

def incorporate_recon(event_edges, cascade_edges, nuflux, angle_edges, just_flux=True):
    """
    This takes in the results from `generate_singly_diff_fluxes` and incorporates reconstruction uncertainties

    Should take in a list or array of energies (true, deposited), in units of eV
    And also take in a list of true cos(zenith) edges 
    """
    e_min = min(cascade_edges)
    e_max = max(cascade_edges)

    z_min = min(angle_edges)
    z_max = max(angle_edges)

    # we need to get all the centers for these bins with the given edges. 
    # these will be associeed with each of the bins in nuflux 
    
    cascade_centers = bhist([cascade_edges]).centers
    depo_widths = bhist([cascade_edges]).widths
    true_e_centers = bhist([event_edges]).centers
    true_ang_centers = bhist([angle_edges]).centers
   
    true_e_widths = bhist([event_edges]).widths
    true_ang_widths = bhist([angle_edges]).widths

    # these are reconstruction objects 
    r_energy = bhist([ np.logspace(np.log10(e_min), np.log10(e_max), int(len(cascade_edges)/2)) ])
    r_angle  = bhist([ np.linspace( z_min, z_max, int(len(angle_edges)/2))])
    
    print("Reconstruction Parameters")
    print("    Energy: {} to {} GeV".format(sci(e_min), sci(e_max)))
    print("    cos(t): {} to {} ".format(z_min, z_max))

    r_energy_centers = r_energy.centers
    r_energy_widths = r_energy.widths
    r_angle_centers = r_angle.centers
    r_angle_widths = r_angle.widths

    #build the data object
    # this thing take in those edges and centers and correctly builds normalized probabilities for the given bins 
    dataobj = DataReco(r_energy.edges, r_angle.edges, cascade_edges, angle_edges)

    # may god have mercy on our souls 
    recoflux = {}
    total_flux = {}
    for key in nuflux.keys():
        print("Reconstructing {} Flux".format(key))
        # energy x, angle y
        if not just_flux:
            recoflux[key] = np.zeros(shape=(len(r_energy_centers),len(true_e_centers), len(r_angle_centers),len(true_ang_centers)))
        total_flux[key] = np.zeros(shape=(len(r_energy_centers), len(r_angle_centers)))
        for i_e_reco in range(len(r_energy_centers)):
            for i_e_depo in range(len(cascade_centers)):
                depo_odds = dataobj.get_energy_reco_odds(i_e_depo, i_e_reco) # unitless
                if depo_odds<=0.:
                    continue
                for i_a_true in range(len(true_ang_centers)):
                    for i_a_reco in range(len(r_angle_centers)):
                        ang_odds = dataobj.get_czenith_reco_odds(i_a_true, i_a_reco,i_e_depo) # unitless 
                        if ang_odds<0.:
                            continue
                        for i_e_true in range(len(true_e_centers)):
                            amt = nuflux[key][i_e_depo][i_e_true][i_a_true]*depo_odds*ang_odds 
                            amt *= true_ang_widths[i_a_true]*true_e_widths[i_e_true]*depo_widths[i_e_depo]
                            if amt>=0:
                                if not just_flux:
                                    recoflux[key][i_e_reco][i_e_true][i_a_reco][i_a_true] += amt 
                                total_flux[key][i_e_reco][i_a_reco] += amt
    if not just_flux:
        _save_data(r_energy.edges, event_edges, r_angle.edges, angle_edges, recoflux)
    _save_flux(r_energy.edges, r_angle.edges, total_flux)


if __name__=="__main__":
    print("Running...")
    a,b,c,d = generate_singly_diff_fluxes(n_bins, debug=False, datafile="")
    incorporate_recon(a,b,c,d)

