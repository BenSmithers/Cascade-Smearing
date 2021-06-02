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
msqs = np.linspace(0,10,20)

t_prime = time()
to_do = len(theta_24s)*len(theta_34s)*len(msqs)
done = 0
for i_24 in theta_24s:
    for i_34 in theta_34s:
        for msq in msqs:
            t_start = time()
            pm = SterileParams(theta13=i_24, theta23=i_34, msq2=msq)
            make_meta_flux(pm)

            t_end = time()
            taken = t_end - t_prime 
            done += 1 
            more = t_prime + taken*to_do/done
            print("Expected completion time {}".format(localtime(more)))

            
