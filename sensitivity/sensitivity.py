"""
Defines some stuff to do the likelihood calculations 
"""
import pickle
import numpy as np
import os

from cascade.utils import get_closest, SterileParams, gen_filename, config, Data
from cascade.sensitivity.eff_area_reader import build_flux

null =  SterileParams()
null_flux = Data(gen_filename(config["datapath"], config["nu_flux"]+".dat", null))

from cascade.sensitivity.systematic_unc import astro_norm_unc
from cascade.sensitivity.systematic_unc import astro_shift_unc, cr_perturb
from cascade.sensitivity.systematic_unc import ice_grad_0, ice_grad_1

from math import sqrt, pi, exp, log
from math import lgamma
twopi = sqrt(2/pi)

def llh_func(exp, cts, sigma):
    kappa = exp/sigma
    total = (kappa*kappa+1)*(log(exp)-2*log(sigma))
    total = total + lgamma(kappa*kappa + cts + 1)
    total = total - (cts + kappa*kappa+1)*log(cts*cts*(1+kappa/sigma))
    total = total - lgamma(kappa*kappa + 1)
    return total

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

        self.ice_grad_0 = ice_grad_0(self.expectation)
        self.ice_grad_1 = ice_grad_1(self.expectation)


        # distinctly add up the plus/minus errors in quadrature 
        self._net_error_m = self.expectation["stat_err"]**2 + self.astr_norm_shift[0]**2 + self.astr_gamma_shift[0]**2
        self._net_error_m = self._net_error_m +self.cr_norm_shift[0]**2 + self.cr_gamma_shift[0]**2 + self.ice_grad_0[0]**2 + self.ice_grad_1[0]**2
        self._net_error_m = np.sqrt(self._net_error_m)
        
        self._net_error_p = self.expectation["stat_err"]**2 + self.astr_norm_shift[1]**2 + self.astr_gamma_shift[1]**2
        self._net_error_p = self._net_error_p + self.cr_norm_shift[1]**2 + self.cr_gamma_shift[1]**2 + self.ice_grad_0[1]**2 + self.ice_grad_1[1]**2
        self._net_error_p = np.sqrt(self._net_error_p)
 
    def get_chi2(self, flux):
        return -2*self.get_llh(flux)
    
    def get_likelihood(self,flux):
        return exp(self.get_llh(flux))

    def get_llh(self, flux):
        """
        Calculates the odds that this flux was measured given the LLHMachine's expectation 

        So we scan over the bins and calculate the odds of measuring "flux" given the expected one
        """
        llh = 0.0
        for i_e in range(len(self.expectation["e_edges"])-1):
            for i_cth in range(len(self.expectation["a_edges"])-1):
                llh += self._eval_llh_bin(i_cth, i_e, flux["event_rate"][i_e][i_cth])
        return llh


    def _eval_llh_bin(self, i_cth, i_e, bflux):

        """
        We assume an asymmetric gaussian distribution
        """

        expected = self.expectation["event_rate"][i_e][i_cth]

        #cutting this out until I find a better asymmetric likelihood function
        if bflux < expected:
            sigma = self._net_error_m[i_e][i_cth]
        else:
            sigma = self._net_error_p[i_e][i_cth]
        
#        sigma = self._mean_error[i_e][i_cth]

        if sigma<0.0:
            raise ValueError("Negative signa??? {}".format(sigma))

        if expected==0.0: #  and bflux==0:
            return 0.0
        
        # likelihood! 
        # https://arxiv.org/pdf/1901.04645.pdf

        elif sigma==0 or (np.isnan(sigma)):
            print("Bad Bin! {}, {}. Expected {}. Got {}".format(i_cth, i_e, expected, bflux))
            return 0.0
        
#        return llh_func(expected, bflux, sigma)

        thingy = -(0.5*((bflux - expected)*(bflux-expected))/(sigma*sigma))
        return log(twopi/(self._net_error_m[i_e][i_cth] + self._net_error_p[i_e][i_cth])) + thingy 

if __name__=="__main__":
    # we will now build up an llh machine, and build up all the likelihood values 
    likelihooder = LLHMachine()

    # We can either scan over all the files found, or just hard-code in the ones...
    hardcode = True
    if hardcode:
        msqs = np.linspace(0,20,40)
        theta24s = np.linspace(0,0.5*pi, 90)
        theta34s = np.linspace(0,0.5*pi, 90)
    else:
        # although more reliable, this is slower. 
        # it's an O(n) process (or worse depending on the glob) 
        from cascade.utils import parse_folder
        theta14s, theta24s, theta34s, msqs = parse_folder(config["datapath"] + "/expected_fluxes_reco/", "exp*.dat")

    likelihoods = [[0,0,0,0.0] for i in range(len(msqs)*len(theta24s)*len(theta34s))]
    chi2 = np.zeros(shape=(len(theta24s), len(theta34s), len(msqs)))

    f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/", "expected_flux.dat", SterileParams()),'rb')
    central_chi = -2*likelihooder.get_llh(pickle.load(f))
    f.close()

    skipped = 0
    found = 0
    for th24 in range(len(theta24s)):
        for th34 in range(len(theta34s)):
            for msq in range(len(msqs)):
                pam = SterileParams(theta13=theta24s[th24], theta23=theta34s[th34], msq2=msqs[msq])

                f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", pam)
                index = th24 + th34*len(theta24s) + msq*len(theta24s)*len(theta34s)
                if not os.path.exists(f_name):
                    skipped+=1 
                    if skipped%1000==0:
                        print("Skipped {}, Found {}".format(skipped, found))
                    likelihoods[index][3] = 0.0
                    chi2[th24][th34][msq] = 1e8 #arbitrarily big...
                else:
                    f = open(f_name, 'rb')
                    try:
                        flux = pickle.load(f)
                    except EOFError:
                        likelihoods[index][3] = 0.0
                        chi2[th24][th34][msq] = 1e8 #arbitrarily big...
                        continue
                    found += 1

                    f.close()

                    llh = likelihooder.get_llh(flux)
                    likelihoods[index][3] = exp(llh)
                    chi2[th24][th34][msq] = -2*llh-central_chi
                   


                likelihoods[index][0] = th24 
                likelihoods[index][1] = th34
                likelihoods[index][2] = msq 


    # old code :frowning: 
    total_likelihood = np.sum(likelihoods, axis=0)[3]
    sorted_likelihoods = sorted(likelihoods, key=lambda entry:entry[3])
    for i in range(len(sorted_likelihoods)):
        sorted_likelihoods[i][3] = sorted_likelihoods[i][3]/total_likelihood

    # make a bunch of cummulative probabilities 
    c_prob = np.zeros(shape=(len(theta24s), len(theta34s), len(msqs)))
    for entry in sorted_likelihoods:
        c_prob[entry[0]][entry[1]][entry[2]] = entry[3]
     
    

    likelihood_dict = {
                "theta24s":theta24s,
                "theta34s":theta34s,
                "msqs":msqs,
                "raw_sorted":sorted_likelihoods,
                "c_prob":c_prob,
                "chi2s":chi2
            }
    f = open("cummulative_probs.dat",'wb')
    pickle.dump(likelihood_dict, f, -1)
    f.close()
