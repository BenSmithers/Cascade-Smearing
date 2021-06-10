"""
Defines some stuff to do the likelihood calculations 
"""
import pickle
import numpy as np

from cascade.utils import get_closest, SterileParams, gen_filename, config, Data
from cascade.sensitivity.eff_area_reader import build_flux

null =  SterileParams()
null_flux = Data(gen_filename(config["datapath"], config["nu_flux"]+".dat", null))

from cascade.sensitivity.systematic_unc import astro_norm_unc
from cascade.sensitivity.systematic_unc import astro_shift_unc, cr_perturb
from cascade.sensitivity.systematic_unc import ice_grad_0, ice_grad_1

class LLHMachine:
    def __init__(self):

        # load null flux
        f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", null)
        f = open(f_name, 'rb')
        self.expectation = pickle.load(f)
        f.close()

        self.astr_norm_shift = astro_norm_unc(null_flux)
        self.astr_gamma_shift = astro_shift_unc(null_flux)

        self.cr_norm_shift = cr_perturb(dnorm=0.05)
        self.cr_gamma_shift = cr_perturb(dgamma=0.012)

        self.ice_grad_0 = ice_grad_0(null)
        self.ice_grad_1 = ice_graD_1(null)

        self.net_error_m = self.expectation["stat_err"][0]**2 + self.astr_norm_shift[0]**2 + self.astr_gamma_shift[0]**2
        self.net_error_m+= self.cr_norm_shift[0]**2 + self.cr_gamma_shift[0]**2 + self.ice_grad_0[0]**2 + self.ice_grad_1[0]**2
        self.net_error_m = np.sqrt(self.net_error_m)
        
        self.net_error_p = self.expectation["stat_err"][1]**2 + self.astr_norm_shift[1]**2 + self.astr_gamma_shift[1]**2
        self.net_error_p+= self.cr_norm_shift[1]**2 + self.cr_gamma_shift[1]**2 + self.ice_grad_0[1]**2 + self.ice_grad_1[1]**2
        self.net_error_p = np.sqrt(self.net_error_p)


    def get_likelihood(self, flux):
        """
        Calculates the odds that this flux was measured given the LLHMachine's expectation 
        """
        pass

    def _eval_probability(i_cth, i_e, bflux):
        """
        For a specific
        """



