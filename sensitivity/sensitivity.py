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
#null_flux = Data(gen_filename(config["datapath"], config["nu_flux"]+".dat", null))
null_flux=""

from cascade.sensitivity.systematic_unc import astro_norm_unc
from cascade.sensitivity.systematic_unc import astro_shift_unc, cr_perturb
from cascade.sensitivity.systematic_unc import ice_grad_0, ice_grad_1

from cascade.sensitivity.generate_all_integrated_fluxes import make_meta_flux

from cascade.sensitivity.load_expectation import load as loade

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
    def __init__(self, expectation=SterileParams(), use_syst=True, flatten=False):
        if not isinstance(expectation, SterileParams):
            raise TypeError("Expected {}, got {}".format(SterileParams, type(expectation)))
        if not isinstance(use_syst, bool):
            raise TypeError("Expected {}, got {}".format(bool, type(use_syst)))
        if not isinstance(flatten, bool):
            raise TypeError("Expected {}, got {}".format(bool, type(flatten)))
        self._use_systematics=use_syst
        self._flatten = flatten
        if not os.path.exists(self.f_name):
            make_meta_flux(expectation, self.use_mc)
        self._configure()

    def set_expectation(self, f_expect, scale_to_n_years=None):
        """
        Changes all the preconfigured expectations, but doesn't change the systematic sources of error 
        That'd be a bit too computationally expensive 
        """
        if isinstance(f_expect, str):
            self.f_name = f_expect

            f = open(self.f_name, 'rb')
            self.expectation = pickle.load(f)
            f.close()
        elif isinstance(f_expect, dict):
            self.expectation = f_expect
        else:
            raise TypeError("Don't know what to do with a {}".format(type(f_expect)))
        
        err = self.expectation["event_rate"]
        # distinctly add up the plus/minus errors in quadrature 
        if scale_to_n_years is not None:
            if not isinstance(scale_to_n_years, (int, float)):
                raise TypeError("Expected {} or {}, got {}".format(int, float, type(scale_to_n_years)))

            if self.use_mc:
                err = err*scale_to_n_years/8.
                self.expectation["event_rate"]*=scale_to_n_years/8.
            else:
                err = err*scale_to_n_years/10. 
                self.expectation["event_rate"]*=scale_to_n_years/10.

        err = np.sqrt(err)
        self._net_error_m_stat = err
        self._net_error_p_stat = err


    def _configure(self, scale_to_n_years=None):
        self.set_expectation(self.f_name, scale_to_n_years)
        
        if self._use_systematics:
            #showplot(self.expectation, self.astr_norm_shift, "Astro Norm")
            self.astr_gamma_shift = astro_shift_unc(use_mc=self.use_mc)
            #showplot(self.expectation, self.astr_gamma_shift, "Astro Gamma")

            self.cr_gamma_shift = cr_perturb(dgamma=0.012, use_mc = self.use_mc)
            #showplot(self.expectation, self.cr_gamma_shift, "CR Gamma")

            self.ice_grad_0 = ice_grad_0(self.expectation)
            self.ice_grad_1 = ice_grad_1(self.expectation)
            #showplot(self.expectation, self.ice_grad_0, "Icegrad 0")
            #showplot(self.expectation, self.ice_grad_1, "Icegrad 1")

            self._net_error_m_sys = self.astr_gamma_shift[0]**2
            self._net_error_m_sys = self._net_error_m_sys  + self.cr_gamma_shift[0]**2 + self.ice_grad_0[0]**2 + self.ice_grad_1[0]**2
            self._net_error_m_sys = np.sqrt(self._net_error_m_sys)

            self._net_error_p_sys = self.astr_gamma_shift[1]**2
            self._net_error_p_sys = self._net_error_p_sys  + self.cr_gamma_shift[1]**2 + self.ice_grad_0[1]**2 + self.ice_grad_1[1]**2
            self._net_error_p_sys = np.sqrt(self._net_error_p_sys)
        else:
            self._net_error_m_sys = np.zeros(shape=np.shape(self._net_error_m_stat))
            self._net_error_p_sys = np.zeros(shape=np.shape(self._net_error_p_stat))
 
    @property
    def _net_error_p(self):
        return np.sqrt(self._net_error_p_stat**2 + self._net_error_p_sys**2)

    @property
    def _net_error_m(self):
        return np.sqrt(self._net_error_m_stat**2 + self._net_error_m_sys**2)


    def get_chi2(self, flux):
        return -2*self.get_llh(flux)
    
    def get_likelihood(self,flux):
        return exp(self.get_llh(flux))

    def get_llh(self, flux):
        """
        Like below, but we allow the normalization to fluctuate 
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
            if self._flatten:
                llh += self._eval_llh_bin_flat(i_e, norm*sum(flux["event_rate"][i_e]))
            else:
                for i_cth in range(len(self.expectation["a_edges"])-1):
                    llh += self._eval_llh_bin(i_cth, i_e, norm*flux["event_rate"][i_e][i_cth])
        return llh

    def _eval_llh_bin_flat(self, i_e, bflux):
        """
        Same as before, but we're flattening along the angular bins
        """
        our_expected = sum(self.expectation["event_rate"][i_e])

        if our_expected==0.0:
            return 0.0

        if bflux<our_expected:
            sigma = np.sqrt(our_expected + sum(self._net_error_m_sys[i_e]**2))
        else:
            sigma = np.sqrt(our_expected + sum(self._net_error_p_sys[i_e]**2))

        if sigma<0:
            raise ValueError("Found negative sigma {}".format(sigma))
        elif sigma==0 or np.isnan(sigma):
            print("Bad Bin: sigma {}".format(sigma))

        return -0.5*(bflux - our_expected)*(bflux - our_expected)/(sigma*sigma)

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
            raise ValueError("Negative sigma??? {}".format(sigma))

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
    def __init__(self, expectation=SterileParams(), use_syst=True, flatten=False):
        self.f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux_from_mc.dat", expectation)
        self.use_mc = True
        _generic_LLHMachine.__init__(self, expectation, use_syst, flatten)

class LLHMachine(_generic_LLHMachine):
    def __init__(self, expectation=SterileParams(), use_syst=True, flatten=False):
        self.f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", expectation)
        self.use_mc = False
        _generic_LLHMachine.__init__(self, expectation, use_syst, flatten)

class LLHMachine_Data(_generic_LLHMachine):
    def __init__(self):
        self.f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", "expected_flux.dat", expectation)
        self.use_mc = False
        _generic_LLHMachine.__init__(self, SterileParams(), True, flatten=True)



if __name__=="__main__":
    import sys
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument("-m","--is_mc", dest="is_mc",
                        type=str, required=False, default="False",
                        help="Should we use the MC data (only for tracks)?")
    parser.add_argument("-s","--systematics", dest="use_syst",
                        type=str, required=False, default="True",
                        help="Should we apply systematic uncertainties?")
    parser.add_argument("--th14", dest="th14",
                        type=float, required=False,
                        default=0.0,
                        help="Value to use for theta-14 for the expectation")
    parser.add_argument("--th24", dest="th24",
                        type=float, required=False,
                        default=0.0,
                        help="Value to use for theta-24 for the expectation")
    parser.add_argument("--th34", dest="th34",
                        type=float, required=False,
                        default=0.0,
                        help="Value to use for theta-34 for the expectation")
    parser.add_argument("--deltam", dest="deltam",
                        type=float, required=False,
                        default=0.0,
                        help="Value to use for delta-m squared for the expectation")
    parser.add_argument("-l","--is_linear", dest="is_linear",
                        type=str, required=False,
                        default="",
                        help="Is the point distribution linear? False means it's logarithmic with the sines")
    parser.add_argument("-f", "--flatten", required=False, dest="flatten",
                        type=str, default="False",
                        help="Flatten along the cos-theta axis")
    parser.add_argument("--compare", dest="compare",
                        type=str, default="False",
                        help="Load the data file and get actual results?")
    args = parser.parse_args()
    use_mc = str(args.is_mc).lower()=='true'
    th14 = args.th14
    th24 = args.th24
    th34 = args.th34
    deltam = args.deltam
    use_syst = str(args.use_syst).lower()=="true"
    if args.is_linear!="":
        is_linear = str(args.is_linear).lower()=='true'
    else:
        is_linear = not use_mc
    flatten = str(args.flatten).lower()=="true" 
    compare = str(args.compare).lower()=="true"
    
    true_data = {}
    if compare:
        flatten = True
        use_mc = False
        true_data["e_edges"]= np.logspace(2,7,29)
        true_data["event_rate"] = loade()

    expectation = SterileParams(theta03=th14, theta13=th24, theta23=th34, msq2=deltam)

    print("{}sing MC fluxes".format("U" if use_mc else "Not u"))
    print("Using expectation {}".format(expectation))
    print("Using {} distribution of param points".format("linear" if is_linear else "logarithmic"))
    print("{}sing systematics".format("U" if use_syst else "Not u"))
    print("{}lattening".format("F" if flatten else "Not f"))

    # we will now build up an llh machine, and build up all the likelihood values 
    if compare:
        likelihooder = LLHMachine_Data()
    else:
        likelihooder = LLHMachine_from_mc(expectation, use_syst) if use_mc else LLHMachine(expectation,use_syst)
    print("Built Likelihooder")

    # We can either scan over all the files found, or just hard-code in the ones...
    hardcode = True
    if hardcode:
        if is_linear:
            theta24s = np.linspace(0,0.5*pi, 90)
            theta34s = np.linspace(0,0.5*pi, 90)
            msqs = np.linspace(0,20,40)
        else:
            n_p = 100 if use_mc else 90 # I don't know why I did this...
            min_s = -3 if use_mc else -2.5
            if use_mc:
                theta34s = [0]
            else:
                theta34s = np.arcsin(np.sqrt(np.logspace(min_s,0, n_p)))/2
            msqs = np.logspace(-2,2,40)
            theta24s = np.arcsin(np.sqrt(np.logspace(min_s,0, n_p)))/2


    likelihoods = [[0,0,0,0.0] for i in range(len(msqs)*len(theta24s)*len(theta34s))]
    chi2 = np.zeros(shape=(len(theta24s), len(theta34s), len(msqs)))

    central_chi = 0.0

    f_name_init = "expected_flux_from_mc.dat" if use_mc else "expected_flux.dat"
    if not compare:
        f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/", f_name_init, expectation),'rb')
        central_chi = -2*likelihooder.get_llh(pickle.load(f))
        f.close()
    else:
        likelihooder.set_expectation(true_data)
        central_chi = -2*likelihooder.get_llh(true_data)

    skipped = 0
    found = 0
    for th24 in range(len(theta24s)):
        for th34 in range(len(theta34s)):
            for msq in range(len(msqs)):

                pam = SterileParams(theta13=theta24s[th24], theta23=theta34s[th34], msq2=msqs[msq])
                f_name = gen_filename(config["datapath"] +"/expected_fluxes_reco/", f_name_init, pam)
                
                index = th24 + th34*len(theta24s) + msq*len(theta24s)*len(theta34s)
                if (skipped+found)%2000==0:
                    print("Skipped {}, Found {}".format(skipped, found))

                # this was implemented when only half of the fluxes were done running.
                if not os.path.exists(f_name):
                    skipped+=1 
                    likelihoods[index][3] = 0.0
                    chi2[th24][th34][msq] = 1e8 #arbitrarily big...
                else:
                    found += 1
                    if compare: #here, we have a data file with the measurement. So we use the available fluxes as an expectation 
                        likelihooder.set_expectation(f_name, 2.0)
                        llh = likelihooder.get_llh(true_data)
                    else:
                        f = open(f_name, 'rb')
                        try:
                            flux = pickle.load(f)
                        except EOFError:
                            likelihoods[index][3] = 0.0
                            chi2[th24][th34][msq] = 1e8 #arbitrarily big...
                            continue


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
    running_total = 0.0
    c_prob = np.zeros(shape=(len(theta24s), len(theta34s), len(msqs)))
    for entry in sorted_likelihoods:
        running_total = running_total + entry[3]
        c_prob[entry[0]][entry[1]][entry[2]] = running_total
     
    

    likelihood_dict = {
                "theta24s":theta24s,
                "theta34s":theta34s,
                "msqs":msqs,
                "raw_sorted":sorted_likelihoods,
                "c_prob":c_prob,
                "chi2s":chi2
            }

    if not compare:
        suffix = "_from_mc" if use_mc else ""
        suffix += "_flat" if flatten else ""
    else:
        suffix = "_compare"
    suffix += "" if use_syst else "_nosys"


    filename = gen_filename(config["datapath"]+"/expectations/", "cummulative_probs"+suffix+".dat", expectation)
    print("Writing Sensitivity File to {}".format(filename))
    f = open(filename,'wb')
    pickle.dump(likelihood_dict, f, -1)
    f.close()
