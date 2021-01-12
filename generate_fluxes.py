
from utils import config
from utils import gen_filename

from raw_fluxes import raw_flux
from deposit import generate_singly_diff_fluxes
from reconstruct import incorporate_recon

import os
import sys

n_bins = config["n_bins"]

def gen_flux(theta13, theta23, msq3):    
    raw_flux_name = gen_filename(config["datapath"], config["nu_flux"], theta13, theta23, msq3)
    if os.path.exists(raw_flux_name) and config["use_pregen_mceq_flux"]:
        pass
    else:
        raw_flux_name = raw_flux(theta13, theta23, msq3)
 
    a,b,c,d = generate_singly_diff_fluxes( config["n_bins"], debug=False, datafile=raw_flux_name)
    incorporate_recon(a,b,c,d, just_flux=True, theta13=theta13, theta23=theta23, msq3=msq3)


if __name__=="__main__":
    
    if len(sys.argv)!=4:
        print("Received incorrect number of args: {}".format(len(sys.argv)-1))

    theta13 = float(sys.argv[1])
    theta23 = float(sys.argv[2])
    msq3    = float(sys.argv[3])
    
    print("Running with:")
    print("    theta13: {}".format(theta13))
    print("    theta23: {}".format(theta23))
    print("    msq3:    {}".format(msq3))

    gen_flux(theta13, theta23, msq3)
