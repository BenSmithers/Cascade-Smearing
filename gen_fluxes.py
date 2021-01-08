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

from cascade.utils import sci # used to generate save names
from cascade.utils import get_closest # interpolation
from cascade.utils import config 
from cascade.utils import backup

# simulates the cosmic ray showers 
from MCEq.core import config, MCEqRun
import crflux.models as crf

def mc_flux(angle_deg):
    mag = 1.
    flux = {}
    mceq = MCEqRun(
            interaction_model = 'SIBYLL23C',
            primary_model = (crf.HillasGaisser2012, 'H3a'),
            theta_deg = angle_deg
            )
    mceq.solve()
    
    flux['e_grid'] = mceq.e_grid

    flux['nue_flux'] = mceq.get_solution('nue',mag)+mceq.get_solution('pr_nue',mag)
    flux['nue_bar_flux'] = mceq.get_solution('antinue',mag)+mceq.get_solution('pr_antinue',mag)
    flux['numu_flux'] = mceq.get_solution('numu',mag)+mceq.get_solution('pr_numu',mag)
    flux['numu_bar_flux'] = mceq.get_solution('antinumu',mag)+mceq.get_solution('pr_antinumu',mag)
    flux['nutau_flux'] = mceq.get_solution('nutau',mag)+mceq.get_solution('pr_nutau',mag)
    flux['nutau_bar_flux'] = mceq.get_solution('antinutau',mag)+mceq.get_solution('pr_antinutau',mag)
  
    return(flux)

def main(theta13=0.13388166, theta23=0.0, msq3 = 1.3):
    """
    This is the main function. It saves a data file for the flux
    """
    if not isinstance(theta13, (int, float)):
        raise TypeError("Expected {}, not {}, for theta13".format(theta13))
    if not isinstance(theta23, (int, float)):
        raise TypeError("Expected {}, not {}, for theta23".format(theta13))
    if not isinstance(msq3, (int, float)):
        raise TypeError("Expected {}, not {}, for msq3".format(theta13))


    un = nsq.Const()

    n_nu = 4 
    Emin = 1.*un.GeV
    Emax = 10.*un.PeV
    cos_zenith_min = -0.999
    cos_zenith_max = 0.

    use_earth_interactions = True

    zeniths = nsq.linspace(cos_zenith_min, cos_zenith_max, angular_bins)
    energies = nsq.logspace(Emin, Emax, energy_bins)

    nus_atm = nsq.nuSQUIDSAtm(zeniths, energies, n_nu, nsq.NeutrinoType.both, use_earth_interactions)

    nus_atm.Set_MixingAngle(0,1,0.563942)
    nus_atm.Set_MixingAngle(0,2,0.154085)
    nus_atm.Set_MixingAngle(1,2,0.785398)

    #sterile parameters 
    nus_atm.Set_MixingAngle(1,3,0.13388166)
    nus_atm.Set_MixingAngle(2,3,0.0)
    nus_atm.Set_SquareMassDifference(3,1.3)


    nus_atm.Set_SquareMassDifference(1,7.65e-05)
    nus_atm.Set_SquareMassDifference(2,0.00247)


    #settting some zenith angle stuff 
    nus_atm.Set_rel_error(1.0e-6)
    nus_atm.Set_abs_error(1.0e-6)
    #nus_atm.Set_GSL_step(gsl_odeiv2_step_rk4)
    nus_atm.Set_GSL_step(nsq.GSL_STEP_FUNCTIONS.GSL_STEP_RK4)

    inistate = np.zeros(shape=(angular_bins, energy_bins, 2, n_nu))
    for angle_bin in range(angular_bins):
        # now we do the thingy for the 
        angle_deg = 180. - (acos(zeniths[angle_bin])*180./pi)
        if angle_deg > 180.:
            angle_deg = 180.

        print("Evaluating {} deg Flux".format(angle_deg))
        flux = mc_flux(angle_deg)

        for neut_type in range(2):
            for flavor in range(n_nu):
                flav_key = ""
                if flavor==0:
                    flav_key = "nue"
                elif flavor==1:
                    flav_key = "numu"
                elif flavor==2:
                    flav_key = "nutau"
                else: # sterile
                    continue
                if neut_type==1: # bar
                    flav_key += "_bar"
                
                flav_key += "_flux"

                flav_key = "nue_bar_flux"
                # now we get the flux at this zenith/energy combo

                for energy_bin in range(energy_bins):                    
                    # (account for the difference in units between mceq and nusquids! )
                    inistate[angle_bin][energy_bin][neut_type][flavor] = get_closest(energies[energy_bin]/un.GeV, flux['e_grid'], flux[flav_key])
        print()

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

    # determine output file 

    split  = config["nu_flux"].split(".")
    assert(len(split)==2) # should be like a filename 
    prefix = split[0]
    suffix = split[1] 
    filename_partial = "_".join((prefix,sci(theta13), sci(theta23),sci(msq3)))
    filename = os.path.join(config["datapath"], filename_partial)
    print("Saving File to {}".format(filename))
    
    if not config["overwrite"]:
        backup(filename)

    f = open(filename, 'wt')
    angle = cos_zenith_min
    energy = int_min_e
    while angle<cos_enith_max:
        while energy<int_max_e:
            #obj.write( ~string~ )
            obj.write("{} {}".format(energy, angle))
            reg_energy = pow(10., energy)
            for flavor in range(n_nu):
                obj.write(" {}".format(nus_atm.EvalFlavor(flavor, angle, reg_energy, 0)))
            for flavor in range(n_nu):
                obj.write(" {}".format(nus_atm.EvalFlavor(flavor, angle, reg_energy, 1)))
            obj.write("\n")

            energy += (int_max_e-int_min_e)/int_en
        angle += (cos_zenith_max-cos_zenith_min)/int_cos

main()
