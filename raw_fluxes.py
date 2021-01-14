"""
Ben Smithers

This is a reimplementation of the c++ version of this same thing
Uses MCEQ to generate neutrino fluxes at-atmosphere, 
then nuSQuIDS to propagate the fluxes to the detector 
"""

# nusquids...
import nuSQUIDSpy as nsq

angular_bins = 50
energy_bins = 121

import numpy as np # the zeros function

from math import pi, acos, log10 # simple math things 
import os # filename and datapath stuff 

from cascade.utils import gen_filename
from cascade.utils import sci # used to generate save names
from cascade.utils import get_closest # interpolation
from cascade.utils import config 
from cascade.utils import backup
from cascade.cross_section_test import xs_obj as xs

# simulates the cosmic ray showers 
from MCEq.core import MCEqRun
import crflux.models as crf

import pickle # saving the 4D MCEq Datafile 

un = nsq.Const()

def get_key(flavor, neutrino):
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
        return("")
    else:
        raise ValueError("Invalid flavor {}".format(flavor))

    if neutrino==0:
        key+="flux"
    elif neutrino==1:
        key+="bar_flux"
    else:
        raise ValueError("Invalid Neutrino Type {}".format(neutrino))

    return(key)

def get_initial_state(energies, zeniths, n_nu):
    """
    This either loads the initial state, or generates it.
    """
    path = os.path.join(config["datapath"], config["mceq_flux"])
    if os.path.exists(path):
        print("Loading MCEq Flux")
        f = open(path, 'rb')
        inistate = pickle.load(f)
        f.close()
    else: 
        print("Generating MCEq Flux")
        inistate = np.zeros(shape=(angular_bins, energy_bins, 2, n_nu))
        mceq = MCEqRun(
                interaction_model = 'SIBYLL23C',
                primary_model = (crf.HillasGaisser2012, 'H3a'),
                theta_deg = 0.
                )
        mag = 0.
        for angle_bin in range(angular_bins):
            # now we do the thingy for the 
            angle_deg = 180. - (acos(zeniths[angle_bin])*180./pi)
            if angle_deg > 180.:
                angle_deg = 180.

            print("Evaluating {} deg Flux".format(angle_deg))
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

def raw_flux(theta13=0.13388166, theta23=0.0, msq3 = 1.3):
    """
    This is the main function. It saves a data file for the flux
    """
    if not isinstance(theta13, (int, float)):
        raise TypeError("Expected {}, not {}, for theta13".format(theta13))
    if not isinstance(theta23, (int, float)):
        raise TypeError("Expected {}, not {}, for theta23".format(theta13))
    if not isinstance(msq3, (int, float)):
        raise TypeError("Expected {}, not {}, for msq3".format(theta13))


    n_nu = 4 
    Emin = 1.*un.GeV
    Emax = 10.*un.PeV
    cos_zenith_min = -0.999
    cos_zenith_max = 0.

    use_earth_interactions = True

    zeniths = nsq.linspace(cos_zenith_min, cos_zenith_max, angular_bins)
    energies = nsq.logspace(Emin, Emax, energy_bins) # DIFFERENT FROM NUMPY LOGSPACE

    nus_atm = nsq.nuSQUIDSAtm(zeniths, energies, n_nu, nsq.NeutrinoType.both, use_earth_interactions)

    nus_atm.Set_MixingAngle(0,1,0.563942)
    nus_atm.Set_MixingAngle(0,2,0.154085)
    nus_atm.Set_MixingAngle(1,2,0.785398)

    #sterile parameters 
    nus_atm.Set_MixingAngle(1,3,theta13)
    nus_atm.Set_MixingAngle(2,3,theta23)
    nus_atm.Set_SquareMassDifference(3,msq3)


    nus_atm.Set_SquareMassDifference(1,7.65e-05)
    nus_atm.Set_SquareMassDifference(2,0.00247)

    nus_atm.SetNeutrinoCrossSections(xs)

    #settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    #nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_FUNCTIONS.GSL_STEP_RK4)

    # LOAD INISTATE
    inistate = get_initial_state(energies, zeniths, n_nu)

    # set the initial flavor
    nus_atm.Set_initial_state(inistate, nsq.Basis.flavor)
    print("Done setting initial state")
    
    nus_atm.Set_ProgressBar(True)
    nus_atm.Set_IncludeOscillations(True)
    
    print("Evolving State")
    nus_atm.EvolveState()

    int_en = 700
    int_cos = 100
    int_min_e = log10(Emin)
    int_max_e = log10(Emax)
    
    filename = gen_filename(config["datapath"], config["nu_flux"], theta13, theta23, msq3)
    print("Saving File to {}".format(filename))

    if not config["overwrite"]:
        backup(filename)

    obj = open(filename, 'wt')
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
    obj.close()
    return(filename)
