"""
Ben Smithers

So this is kinda the master file. 
The gen_flux function handles generating the actual final state recon-recon flux files. 
"""

from cascade.utils import config
from cascade.utils import gen_filename, SterileParams

from cascade.raw_fluxes import raw_flux
from cascade.deposit import generate_singly_diff_fluxes
from cascade.reconstruct import incorporate_recon

import os
import sys

n_bins = config["n_bins"]

def gen_flux(params):
    """
    First we make a neutrino flux file at the detector (true flux binning)
    Then we get the fluxes, binned in energy deposited
    Then we incorporate the detector response 
    """
    if not isinstance(params, SterileParams):
        raise TypeError("Expected {} for params, not {}".format(SterileParams, type(params)))
    

    raw_flux_name = gen_filename(config["datapath"], config["nu_flux"]+".dat", params)
    if os.path.exists(raw_flux_name) and config["use_pregen_mceq_flux"]:
        pass
    else:
        raw_flux_name = raw_flux(params)
 
    a,b,c,d,err = generate_singly_diff_fluxes( config["n_bins"], debug=False, datafile=raw_flux_name)
    incorporate_recon(a,b,c,d,errors=err,params=params, just_flux=True)
    

if __name__=="__main__":
    #This just lets us run this from CLI or as a module 
    
    if len(sys.argv)!=5:
        print("Received incorrect number of args: {}".format(len(sys.argv)-1))

    theta03 = float(sys.argv[1])
    theta13 = float(sys.argv[2])
    theta23 = float(sys.argv[3])
    msq2    = float(sys.argv[4])
    
    print("Running with:")
    print("    theta03: {}".format(theta03))
    print("    theta13: {}".format(theta13))
    print("    theta23: {}".format(theta23))
    print("    msq2:    {}".format(msq2))

    params = SterileParams(theta03=theta03, theta13=theta13, theta23=theta23, msq2=msq2)

    gen_flux(params)
