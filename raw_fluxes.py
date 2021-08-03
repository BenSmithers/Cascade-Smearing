"""
Ben Smithers

This is a reimplementation of the c++ version of this same thing
Uses MCEQ to generate neutrino fluxes at-atmosphere, 
then nuSQuIDS to propagate the fluxes to the detector. 

Super cool though, it only uses MCEQ to make the fluxes at-atmosphere if this hasn't already been done. So that saves a lot of time when you make a lot of fluxes 
"""

# nusquids...
import nuSQUIDSpy as nsq

angular_bins = 50
energy_bins = 121

import numpy as np # the zeros function

from math import pi, acos, log10 # simple math things 
from math import asin, sin
import os # filename and datapath stuff 

from cascade.utils import gen_filename, SterileParams
from cascade.utils import sci # used to generate save names
from cascade.utils import get_closest # interpolation
from cascade.utils import config 
from cascade.utils import backup
from cascade.utils import Data
from cascade.cross_section_test import xs_obj as xs

# simulates the cosmic ray showers 
from MCEq.core import MCEqRun
import crflux.models as crf

import pickle # saving the 4D MCEq Datafile 
from tempfile import NamedTemporaryFile

un = nsq.Const()

def get_key(flavor, neutrino):
    """
    Little function to get the key to access the mceq dictionary
    You just give it the flavor/neutrino index (starts at 0, duh)
    """
    # flavor 0, 1, 2, 3
    # neutrino 0, 1

    key = ""
    if flavor==0:
        key+="nue_"
    elif flavor==1:
        key+="numu_"
    elif flavor==2:
        key+="nutau_"
    elif flavor==3:
        return("") # sterile 
    else:
        raise ValueError("Invalid flavor {}".format(flavor))

    if neutrino==0:
        key+="flux"
    elif neutrino==1:
        key+="bar_flux"
    else:
        raise ValueError("Invalid Neutrino Type {}".format(neutrino))

    return(key)

def get_initial_state(energies, zeniths, n_nu, kwargs):
    """
    This either loads the initial state, or generates it.
    Loading it is waaaay quicker.

    Possible issue! If you run a bunch of jobs and don't already have this flux generated, bad stuff can happen. 
    I'm imagining issues where a bunch of jobs waste time making this, and then all try to write to the same file
    Very bad. Big crash. Very Fail 
    """
    path = os.path.join(config["datapath"], config["mceq_flux"])
    if os.path.exists(path): 
#        print("Loading MCEq Flux")
        f = open(path, 'rb')
        inistate = pickle.load(f)
        f.close()
    else: 
#        print("Generating MCEq Flux")
        inistate = np.zeros(shape=(angular_bins, energy_bins, 2, n_nu))
        mceq = MCEqRun(
                interaction_model = config["interaction_model"],
                primary_model = (crf.HillasGaisser2012, 'H3a'),
                theta_deg = 0.
                )

        r_e = 6.378e6 # meters
        ic_depth = 1.5e3 # meters 
        mag = 0. # power energy is raised to and then used to scale the flux
        for angle_bin in range(angular_bins):
            # get the MCEq angle from the icecube zenith angle 
            angle_deg = asin(sin(pi-acos(zeniths[angle_bin]))*(r_e-ic_depth)/r_e)
            angle_deg = angle_deg*180./pi
            if angle_deg > 180.:
                angle_deg = 180.

            print("Evaluating {} deg Flux".format(angle_deg))
            # for what it's worth, if you try just making a new MCEqRun for each angle, you get a memory leak. 
            # so you need to manually set the angle 
            mceq.set_theta_deg(angle_deg)
            mceq.solve()

            flux = {}
            flux['e_grid'] = mceq.e_grid

            flux['nue_flux'] = mceq.get_solution('nue',mag)+mceq.get_solution('pr_nue',mag)
            flux['nue_bar_flux'] = mceq.get_solution('antinue',mag)+mceq.get_solution('pr_antinue',mag)
            flux['numu_flux'] = mceq.get_solution('numu',mag)+mceq.get_solution('pr_numu',mag)
            flux['numu_bar_flux'] = mceq.get_solution('antinumu',mag)+mceq.get_solution('pr_antinumu',mag)
            flux['nutau_flux'] = mceq.get_solution('nutau',mag)+mceq.get_solution('pr_nutau',mag)
            flux['nutau_bar_flux'] = mceq.get_solution('antinutau',mag)+mceq.get_solution('pr_antinutau',mag) 
            
            for neut_type in range(2):
                for flavor in range(n_nu):
                    flav_key = get_key(flavor, neut_type)
                    if flav_key=="":
                        continue

                    for energy_bin in range(energy_bins):                    
                        # (account for the difference in units between mceq and nusquids! )
                        inistate[angle_bin][energy_bin][neut_type][flavor] = get_closest(energies[energy_bin]/un.GeV, flux['e_grid'], flux[flav_key])
        
        # save it now
        f = open(path, 'wb')
        pickle.dump(inistate, f, -1)
        f.close()

    return(inistate)

def raw_flux(params, kwargs={}):
    """
    This is the main function. It saves a data file for the flux with a unique name for the given physics 
    """
    if not isinstance(params,SterileParams):
        raise TypeError("Expected {} for params, not {}".format(SterileParams, type(params)))

   
    if "forced_filename" in kwargs:
        forced_filename = kwargs["forced_filename"]
    else:
        forced_filename = None
    if "state_setter" in kwargs:
        state_setter = kwargs["state_setter"]
    else:
        state_setter = get_initial_state

    if "osc" in kwargs:
        osc = kwargs["osc"]
    else:
        osc = True
    if not osc:
        print("NOT USING OSCILLATIONS")

    if "as_data" in kwargs:
        as_data = kwargs["as_data"]
    else:
        as_data = False

    if forced_filename is not None:
        if not isinstance(forced_filename, str):
            raise TypeError("Forced filename should be {}, or {}".format(str, None))
  
    print("Propagating Neutrinos at {}".format(params))
    n_nu = 4 
    Emin = 1.*un.GeV
    Emax = 10.*un.PeV
    cos_zenith_min = -0.999
    cos_zenith_max = 0.2

    use_earth_interactions = True

    zeniths = nsq.linspace(cos_zenith_min, cos_zenith_max, angular_bins)
    energies = nsq.logspace(Emin, Emax, energy_bins) # DIFFERENT FROM NUMPY LOGSPACE

    nus_atm = nsq.nuSQUIDSAtm(zeniths, energies, n_nu, nsq.NeutrinoType.both, use_earth_interactions)

    nus_atm.Set_MixingAngle(0,1,0.563942)
    nus_atm.Set_MixingAngle(0,2,0.154085)
    nus_atm.Set_MixingAngle(1,2,0.785398)

    #sterile parameters 
    nus_atm.Set_MixingAngle(0,3,params.theta03)
    nus_atm.Set_MixingAngle(1,3,params.theta13)
    nus_atm.Set_MixingAngle(2,3,params.theta23)
    nus_atm.Set_SquareMassDifference(3,params.msq2)


    nus_atm.Set_SquareMassDifference(1,7.65e-05)
    nus_atm.Set_SquareMassDifference(2,0.00247)

    nus_atm.SetNeutrinoCrossSections(xs)

    nus_atm.Set_TauRegeneration(True)

    #settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    #nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_FUNCTIONS.GSL_STEP_RK4)
    
    # we load in the initial state. Generating or Loading from a file 
    inistate = state_setter(energies, zeniths, n_nu, kwargs)
    nus_atm.Set_initial_state(inistate, nsq.Basis.flavor)
    
    # we turn off the progress bar for jobs run on the cobalts 
    nus_atm.Set_ProgressBar(False)
    nus_atm.Set_IncludeOscillations(osc)
    
    nus_atm.EvolveState()

    int_en = 700
    int_cos = 100
    int_min_e = log10(Emin)
    int_max_e = log10(Emax)
      
    filename=""
    if not as_data:
        
        if forced_filename is None: 
            filename = gen_filename(config["datapath"], config["nu_flux"]+".dat", params)
        else:
            filename = forced_filename
        print("Saving File to {}".format(filename))
        
        if not config["overwrite"]:
            backup(filename)

        obj = open(filename, 'wt')
    else:
        cobalt = os.environ.get("_CONDOR_SCRATCH_DIR")
        if cobalt==None or cobalt=="" or cobalt==".":
            this_dir = None
        else:
            this_dir = cobalt
        obj = NamedTemporaryFile(mode='wt',buffering=1, dir=this_dir)

    angle = cos_zenith_min
    energy = int_min_e
    obj.write("# log10(energy) cos(zenith) flux_nuE flux_nuMu flux_nuTau flux_nuEBar flux_nuMuBar flux_nuTauBar\n")
    while angle<cos_zenith_max:
        energy = int_min_e
        while energy<int_max_e:
            #obj.write( ~string~ )
            obj.write("{} {}".format(energy, angle))
            reg_energy = pow(10., energy)
            for flavor in range(n_nu):
                obj.write(" {}".format(nus_atm.EvalFlavor(flavor, angle, reg_energy, 0)))
            for flavor in range(n_nu):
                obj.write(" {}".format(nus_atm.EvalFlavor(flavor, angle, reg_energy, 1)))

            energy += (int_max_e-int_min_e)/int_en
            obj.write("\n")

        angle += (cos_zenith_max-cos_zenith_min)/int_cos
    
    if as_data:
        data_obj = Data(obj.name)
        obj.close() # deletes tempfile 
        return data_obj
    else:
        obj.close()
        return(filename)
