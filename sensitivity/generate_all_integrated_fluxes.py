from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.sensitivity.eff_area_reader import build_flux 

from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import Data
from cascade.raw_fluxes import raw_flux

import numpy as np

import pickle
from time import time, localtime

def make_meta_flux(params):
    # look for the atmospheric fluxes. These should all be pre-generated 
    print("Loading Fluxes at {}".format(params))
    atmo_file = raw_flux(params)
    atmo_data = Data(atmo_file)

    # we need to load in the data ojects! 
    astr_filepath = generate_astr_flux(params) # will load file if it exists 
    astr_data = Data(astr_filepath)

    print("Calculating Expected Binned Flux at {}".format(params))
    # now we use these two to build the full expected flux
    full_flux = build_flux(atmo_data, astr_data) #dict 
    
    # save the object
    new_filename = gen_filename(config["datapath"]+ "/expected_fluxes_reco/", "expected_flux.dat", params)
    f = open(new_filename ,'wb')
    pickle.dump(full_flux, f, -1)
    f.close()


theta_24s = np.linspace(0,90, 90)
theta_34s = np.linspace(0,90, 90)
msqs = np.linspace(0,20,40)


if __name__=="__main__":
    import sys

    th24 = float(sys.argv[1])
    th34 = float(sys.argv[2])

    for msq in msqs:
        pm = SterileParams(theta13=th24, theta23=th34, msq2=msq)
        make_meta_flux(pm)

