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

from math import sqrt, pi, exp

twopi = 1./sqrt(2*pi)

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


        # distinctly add up the plus/minus errors in quadrature 
        self.net_error_m = self.expectation["stat_err"][0]**2 + self.astr_norm_shift[0]**2 + self.astr_gamma_shift[0]**2
        self.net_error_m+= self.cr_norm_shift[0]**2 + self.cr_gamma_shift[0]**2 + self.ice_grad_0[0]**2 + self.ice_grad_1[0]**2
        self.net_error_m = np.sqrt(self.net_error_m)
        
        self.net_error_p = self.expectation["stat_err"][1]**2 + self.astr_norm_shift[1]**2 + self.astr_gamma_shift[1]**2
        self.net_error_p+= self.cr_norm_shift[1]**2 + self.cr_gamma_shift[1]**2 + self.ice_grad_0[1]**2 + self.ice_grad_1[1]**2
        self.net_error_p = np.sqrt(self.net_error_p)

    
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
                llh += self._eval_llh_bin(i_cth, i_e, flux[i_e][e_cth])
        return llh


    def _eval_llh_bin(i_cth, i_e, bflux):
        """
        We assume an asymmetric gaussian distribution
        """

        expected = self.expectation["event_rate"][i_e][i_cth]

        if bflux < expected:
            sigma = self._net_error_m[i_e][i_cth]
        else:
            sigma = self._net_error_p[i_e][i_cth]
        return log(twopi/sigma) - (0.5*((bflux - expected)**2)/(sigma*sigma))

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

    f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/", "expected_flux.dat", SterileParams()))
    central_chi = -2*likelihooder.get_llh(pickle.load(f))
    f.close()

    for th24 in range(len(theta24s)):
        for th34 in range(len(theta34s)):
            for msq in range(len(msqs)):
                pam = SterileParams(theta13=theta24s[th24], theta32=theta34s[th34], msq2=msqs[msq])

                f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", pam)
                if not os.path.exists(f_name):
                    print("File not found, skipping {}".format(f_name))
                    likelihoods[index][3] = 0.0
                    chi2[th24][th34][msq] = 1e8 #arbitrarily big...
                else:
                    f = open(f_name, 'rb')
                    flux = pickle.load(f)
                    f.close()

                    llh = likelihooder.get_llh(flux)
                    likelihoods[index][3] = exp(llh)
                    chi2[th24][th34][msq] = -2*llh-central_chi
                   
                index = th24 + th34*len(theta24s) + msq*len(theta24s)*len(theta34s)

                likelihoods[index][0] = th24 
                likelihoods[index][1] = th34
                likelihoods[index][2] = msq 


    # old code :frowning: 
    total_likelihood = np.sum(likelihoods, axis=0)[3]
    sorted_likelihoods = sorted(likelihoods, key=lambda entry:entry[3])
    for i in range(len(sorted_likelihoods))
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
    f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/","cummulative_probs.dat"),'wb')
    pickle.dump(likelihood_dict, f, -1)
    f.close()
