"""
Defines some stuff to do the likelihood calculations 
"""
import pickle
import numpy as np
import os

from scipy.optimize import minimize 

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


import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
def showplot(exp, thing, title):
    plt.pcolormesh(exp["a_edges"], exp["e_edges"], thing[0])
    plt.title("{}: Minus".format(title))
    plt.xlabel("Cos th")
    plt.ylabel("Energy")
    plt.yscale('log')
    plt.colorbar()
    plt.show()

    plt.pcolormesh(exp["a_edges"], exp["e_edges"], thing[1])
    plt.title("{}: Plus".format(title))
    plt.xlabel("Cos th")
    plt.ylabel("Energy")
    plt.yscale('log')
    plt.colorbar()
    plt.show()

class _generic_LLHMachine:
    def __init__(self):
        raise NotImplementedError("Use Inherited Class")

        self.f_name = ""

    def _configure(self):
        f = open(self.f_name, 'rb')
        self.expectation = pickle.load(f)
        f.close()


        self.astr_norm_shift = astro_norm_unc(null_flux, use_mc=self.use_mc)
        #showplot(self.expectation, self.astr_norm_shift, "Astro Norm")
        self.astr_gamma_shift = astro_shift_unc(null_flux, use_mc=self.use_mc)
        #showplot(self.expectation, self.astr_gamma_shift, "Astro Gamma")

        self.cr_norm_shift = cr_perturb(dnorm=0.05, use_mc = self.use_mc)
        self.cr_gamma_shift = cr_perturb(dgamma=0.012, use_mc = self.use_mc)
        #showplot(self.expectation, self.cr_norm_shift, "CR Norm")
        #showplot(self.expectation, self.cr_gamma_shift, "CR Gamma")

        self.ice_grad_0 = ice_grad_0(self.expectation)
        self.ice_grad_1 = ice_grad_1(self.expectation)
        #showplot(self.expectation, self.ice_grad_0, "Icegrad 0")
        #showplot(self.expectation, self.ice_grad_1, "Icegrad 1")

        # distinctly add up the plus/minus errors in quadrature 
        self._net_error_m = self.expectation["stat_err"]**2 + self.astr_norm_shift[0]**2 + self.astr_gamma_shift[0]**2
        self._net_error_m = self._net_error_m  + self.cr_gamma_shift[0]**2 + self.ice_grad_0[0]**2 + self.ice_grad_1[0]**2
        self._net_error_m = np.sqrt(self._net_error_m)
        
        self._net_error_p = self.expectation["stat_err"]**2 + self.astr_norm_shift[1]**2 + self.astr_gamma_shift[1]**2
        self._net_error_p = self._net_error_p  + self.cr_gamma_shift[1]**2 + self.ice_grad_0[1]**2 + self.ice_grad_1[1]**2
        self._net_error_p = np.sqrt(self._net_error_p)
 
    def get_chi2(self, flux):
        return -2*self.get_llh(flux)
    
    def get_likelihood(self,flux):
        return exp(self.get_llh(flux))

    def get_llh(self, flux):
        """
        Like below, but we allow the normalization to flutuate 
        """
        norm_central=1.0

        opt = minimize(lambda norm : -2*self.get_llh_fixed(flux, norm[0]), [norm_central])
        return self.get_llh_fixed(flux, opt.x[0])

    def get_llh_fixed(self, flux, norm=1.0):
        """
        Calculates the odds that this flux was measured given the LLHMachine's expectation 

        So we scan over the bins and calculate the odds of measuring "flux" given the expected one

        Overall normalization is fixed in this one
        """
        llh = 0.0
        for i_e in range(len(self.expectation["e_edges"])-1):
            for i_cth in range(len(self.expectation["a_edges"])-1):
                llh += self._eval_llh_bin(i_cth, i_e, norm*flux["event_rate"][i_e][i_cth])
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
        return thingy #log(twopi/(self._net_error_m[i_e][i_cth] + self._net_error_p[i_e][i_cth])) + thingy 

class LLHMachine_from_mc(_generic_LLHMachine):
    def __init__(self):
        # load null flux
        self.f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux_from_mc.dat", null)
        self.use_mc = True
        self._configure()
 
class LLHMachine(_generic_LLHMachine):
    def __init__(self):
        # load null flux
        self.f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", null)
        self.use_mc = False
        self._configure()

if __name__=="__main__":
    import sys

    if len(sys.argv)>=2:
        use_mc = True if (sys.argv[1].lower() in ["mc", "1"]) else False
    else:
        use_mc = False
    print("{}sing MC fluxes".format("U" if use_mc else "Not u"))

    # we will now build up an llh machine, and build up all the likelihood values 
    likelihooder = LLHMachine_from_mc() if use_mc else LLHMachine()
    print("Built Likelihooder")

    # We can either scan over all the files found, or just hard-code in the ones...
    hardcode = True
    if hardcode:


        if use_mc:
            theta34s = [0]
            msqs = np.logspace(-2,2,40)
            theta24s = np.arcsin(np.sqrt(np.logspace(-3,0, 100)))/2
        else:
            theta24s = np.linspace(0,0.5*pi, 90)
            theta34s = np.linspace(0,0.5*pi, 90)
            msqs = np.linspace(0,20,40)

    likelihoods = [[0,0,0,0.0] for i in range(len(msqs)*len(theta24s)*len(theta34s))]
    chi2 = np.zeros(shape=(len(theta24s), len(theta34s), len(msqs)))

    f_name_init = "expected_flux_from_mc.dat" if use_mc else "expected_flux.dat"
    f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/", f_name_init, SterileParams()),'rb')
    central_chi = -2*likelihooder.get_llh(pickle.load(f))
    f.close()

    skipped = 0
    found = 0
    for th24 in range(len(theta24s)):
        for th34 in range(len(theta34s)):
            for msq in range(len(msqs)):
                pam = SterileParams(theta13=theta24s[th24], theta23=theta34s[th34], msq2=msqs[msq])

                f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", f_name_init, pam)
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
    suffix = "_from_mc" if use_mc else ""
    f = open("cummulative_probs.dat"+suffix,'wb')
    pickle.dump(likelihood_dict, f, -1)
    f.close()
