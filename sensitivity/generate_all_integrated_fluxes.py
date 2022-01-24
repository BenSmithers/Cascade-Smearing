from cmath import pi
from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.sensitivity.eff_area_reader import build_flux 
from cascade.sensitivity.make_from_mc import build_mc_flux 
from cascade.sensitivity.interp_eff_area import get_expectation

from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import Data
from cascade.raw_fluxes import raw_flux



import numpy as np

import pickle
from time import time, localtime
import os

def make_meta_flux(params, do_mc = False):
    """
    This makes and saves 
    """
    # look for the atmospheric fluxes. These should all be pre-generated 
    start = time() 
    kwargs = {}
    kwargs["as_data"]=True
    atmo_data = raw_flux(params,kwargs=kwargs)
    astr_data = generate_astr_flux(params, as_data=True)

    print("Calculating Expected Binned Flux at {}".format(params))
    # now we use these two to build the full expected flux
    if do_mc:
        full_flux = build_mc_flux(atmo_data, astr_data)
    else:
        full_flux = get_expectation(atmo_data, astr_data) #dict 
    middle = time()
    # save the object
    filename = "expected_flux_from_mc_smearedwell.dat" if do_mc else "best_expected_flux.dat"

    new_filename = gen_filename(config["datapath"]+ "/expected_fluxes_reco/", filename, params)

    print("Saving to {}".format(new_filename))
    f = open(new_filename ,'wb')
    pickle.dump(full_flux, f, -1)
    f.close()
    end = time()

    print("Flux Sim took {:.1f} seconds".format(middle-start))
    print("Saving took {:.3f} seconds".format(end-middle))

    return full_flux
    

if __name__=="__main__":

    n_m = 1
    msqs = np.concatenate(( np.array([0]), np.logspace(-2,2,n_m) ))

    if n_m==1:
        msqs=msqs[:1]

    import sys

    th24 = float(sys.argv[1])
    th34 = float(sys.argv[2])
    if len(sys.argv)==4:
        mc_mode = int(sys.argv[3])==1
    else:
        mc_mode = False

    for msq in msqs:
        pm = SterileParams(theta13=th24, theta23=th34, msq2=msq)
        make_meta_flux(pm, do_mc=mc_mode)

