"""
Re-write of the old sensitivity code

This script is used to do likelihood scans 
"""

import numpy as np
import pickle 
import os
from numbers import Number
from math import sqrt

from scipy.optimize import minimize

from cascade.utils import SterileParams, gen_filename, config


from cascade.sensitivity.systematic_unc import unc_wrapper, Syst

data_folder = os.path.join(config["datapath"], "expected_fluxes_reco")

class generalLLH:
    def __init__(self):
        pass

    def get_llh(self, params):
        return 0.0


class doLLH(generalLLH):
    """
    Calculates likelihoods.We pass it a central likelihood we use to develop the uncertainties 

    There are a couple things you'll find useful 
        minimize - pass it a flux and it'll tell you the normalization that minimizes the delta-LLH
        load_file - pass it some params and it'll give you the actual flux dictionary
        get_llh - pass it some params, and it'll return the delta-llh based off the central expectation 
        eval_llh_norm - pass it a flux and a normalization, and it'll tell you the delta-llh 
    """
    def __init__(self,filenames, central_exp=SterileParams(), options={}):
        self._is_mc = False
        self._flatten = False
        self._upgoing = True
        self._skip_missing = False
        self._use_syst = True
        self._fix_norm = -1
        self._fixing_norm = False
        if not isinstance(filenames, str):
            raise TypeError("Filename should be {}, not {}".format(str, type(filenames)))
        self._filenames = filenames

        # do the initial configuration 
        self.set_central(central_exp)

        self._parse_options(options)
        self._options = options
        
        
        
        self._net_error_m_sys = np.zeros(shape=np.shape(self._net_error_m_stat))
        self._net_error_p_sys = np.zeros(shape=np.shape(self._net_error_p_stat))

        if self._use_syst:
            fn = gen_filename(data_folder, self._filenames, central_exp)
            f = open(fn,'rb')
            expectation = pickle.load(f)
            f.close()

            ice_grad_0 = unc_wrapper(expectation, Syst.icegrad0, options)
            ice_grad_1 = unc_wrapper(expectation, Syst.icegrad1, options)
            astro_grad = unc_wrapper(expectation, Syst.astrogam, options)
            cr_grad    = unc_wrapper(expectation, Syst.crgam   , options)
             
            self._net_error_m_sys = astro_grad[0]**2
            self._net_error_m_sys = self._net_error_m_sys  + cr_grad[0]**2 + ice_grad_0[0]**2 + ice_grad_1[0]**2
            self._net_error_m_sys = np.sqrt(self._net_error_m_sys)

            self._net_error_p_sys = astro_grad[1]**2
            self._net_error_p_sys = self._net_error_p_sys  + cr_grad[1]**2 + ice_grad_0[1]**2 + ice_grad_1[1]**2
            self._net_error_p_sys = np.sqrt(self._net_error_p_sys)

    def set_central(self, params):
        """
        Changes the expectation from which LLHs are calculated!

        Loads the expectation in, sets it 
        """
        fn =gen_filename(data_folder, self._filenames, params)
        f = open(fn, 'rb')
        data = pickle.load(f)
        f.close()
        
        self._expectation = data['event_rate']
        self.set_central_from_expectation(self._expectation)

    def set_central_from_expectation(self, exp):
        """
        Actually does the setting, recalculates the statistical uncertainty 
        """
        if not isinstance(exp, (np.ndarray, tuple, list)):
            raise TypeError()
        if not (np.shape(self._expectation)==np.shape(exp)):
            raise ValueError("Expected shape {}, got {}".format(self._expectation, exp))

        self._expectation = exp
        err = np.sqrt(exp)
        self._net_error_m_stat = err
        self._net_error_p_stat = err

    def _parse_options(self, options):
        """
        Parses the options passed to the doLLH object. 
        Sets the "is_mc" "flatten" and "upgoing" along with a lot of others 

          - upgoing, only look at those with cos(theta)<0.2 if True
          - fudge_factor: scales the statistics by a fudge factor: either a number-like or a np.ndarray with the same shape as the expectation 
        """
        known = ["is_mc", "flatten", "upgoing", "skip_missing", "use_syst", "fix_norm", "fudge"]

        for key in options.keys():
            if key not in known:
                raise KeyError("Not sure what to do with the key {}".format(key))
        
        def read(index, dtype):
            if known[index] in options:
                if not isinstance(options[known[index]], dtype):
                    raise TypeError(known[index])
                return options[known[index]]

        self._is_mc        = self._is_mc   if read(0, bool) is None else read(0, bool) 
        self._flatten      = self._flatten if read(1, bool) is None else read(1, bool)
        self._upgoing      = self._upgoing if read(2, bool) is None else read(2, bool)
        self._skip_missing = self._skip_missing if read(3,bool) is None else read(3,bool)
        self._use_syst     =self._use_syst if read(4, bool) is None else read(4, bool)
        self._fix_norm     =self._fix_norm if read(5,float) is None else read(5,float)
        
        # the fudge needs special handling
        self._fudge_factor = 1.0
        if  known[6] in options:
            if isinstance(options[known[6]], Number):
                self._fudge_factor = float(options[known[6]])
                self._fudge_unc = sqrt(self._fudge_factor)
            elif isinstance(options[known[6]], (np.ndarray, list, tuple)):
                if not np.shape(options[known[6]])==np.shape(self._expectation):
                    raise ValueError("Can't work with shapes {} and {}".format(np.shape(options[known[6]]), np.shape(self._expectation)))
                if isinstance(options[known[6]],np.ndarray):
                    self._fudge_factor = options[known[6]]
                else:
                    self._fudge_factor = np.array(options[known[6]])
                
                self._fudge_unc = np.sqrt(self._fudge_factor)
            else:
                raise TypeError("Unknown fudge type: {}".format(type(options[known[6]])))


        if self._fix_norm==-1:
            self._fixing_norm = False
        else:
            if self._fix_norm <=0.0:
                raise ValueError("Cannot fix norm to a negative value")
            self._fixing_norm = True

    def load_file(self, params):
        """
        Load the corresponding file in! return the event rate
        """
        fn = gen_filename(data_folder, self._filenames, params) 
        if not os.path.exists(fn):
            raise IOError("Couldn't find file {}".format(fn))
        f = open(fn, 'rb')
        try:
            data = pickle.load(f)
        except IOError:
            print("Error Loading File {}".format(fn))
            raise IOError("Error loading file, but not skipping these")
        f.close()
        return data

    def get_llh(self, params):
        """
        Calculates the LLH for the provided set of parameters

        We allow the normalization of the flux we're loading to float, and minimize the LLH over that normalization
        """

        # Load the file in. If there's an issue, check if we care
        # If we're skipping problematic files 
        try:
            data = self.load_file(params)
        except IOError as e:
            if self._skip_missing:
                return 1e8
            else:
                raise IOError(e)


        if self._fixing_norm:
            return self.eval_llh_norm(data["event_rate"], self._fix_norm)
        else:
            return self.eval_llh_norm(data["event_rate"], self.minimize(data["event_rate"]))

    def minimize(self, evt_rate):
        """
        return normalization that minimizes likelihood 
        """
        norm_guess = 1.0
        optimized = minimize(lambda norm: -2*self.eval_llh_norm(evt_rate, norm[0]), [norm_guess])
        return optimized.x[0]


    def eval_llh_norm(self, this_flux, norm:float):
        """
        Called by get_llh. 

        We evaluate the llh with the normalization fixed
        """
        llh = 0.0
        for i_e in range(len(this_flux)):
            if self._flatten:
                llh += self._eval_llh_flat_bin(i_e, norm*sum((self._fudge_factor*this_flux)[i_e]))
            else:
                for i_cth in range(len(this_flux[i_e])):
                    llh+= self._eval_llh_bin(i_e, i_cth, norm*(self._fudge_factor*this_flux)[i_e][i_cth])
        return llh

    def _eval_llh_flat_bin(self, i_e:int, bflux):
        raise NotImplementedError()

    def _eval_llh_bin(self, i_e:int, i_cth:int, bflux):
        """
        Calculates the LLH for a specific bin 
        """
        expected = (self._fudge_factor*self._expectation)[i_e][i_cth]

        if expected == 0.0:
            return 0.0

        if bflux < expected:
            sigma = self._net_error_m[i_e][i_cth]
        else:
            sigma = self._net_error_p[i_e][i_cth]

        if sigma<0.0:
            raise ValueError()

        if sigma==0.0 or np.isnan(sigma):
            raise ValueError("sigma: {}".format(sigma))

        return -0.5*((bflux-expected)*(bflux-expected))/(sigma*sigma)

    @property
    def _net_error_p(self):
        return np.sqrt((self._fudge_unc*self._net_error_p_stat)**2 + self._net_error_p_sys**2)

    @property
    def _net_error_m(self):
        return np.sqrt((self._fudge_unc*self._net_error_m_stat)**2 + self._net_error_m_sys**2)
    @property
    def options(self):
        return self._options

class JointLLH(generalLLH):
    """
    Use multiple different doLLH's to get a combined llh 

    So, if you have two different fluxes, make one of these and it'll do all the combining itself 

    It uses the first one to set an overall normalization too! 
    """
    def __init__(self, *args):
        for arg in args:
            if not isinstance(arg, doLLH):
                raise TypeError()

        self.doLLHs = list(args)

    def get_llh(self, params):
        """
        Gets the normalization from the first one, then use that to do the others 
        """
        data = self.doLLHs[0].load_file(params)
        norm = self.doLLHs[0].minimize(data["event_rate"])
        llh = self.doLLHs[0].eval_llh_norm(data["event_rate"], norm)

        for LLHobj in self.doLLHs[1:]:
            data = LLHobj.load_file(params)
            llh += LLHobj.eval_llh_norm(data["event_rate"], norm)

        return llh

class Scanner:
    """
    Pass a likelihood calculator and lists of sterile nu parameters
    This can then scan the expectations and return all the LLHs to you!
    """
    def __init__(self, llHooder : generalLLH, theta24s, theta34s, msqs):
        """
        Store the information, get ready to scan 
        """
        self.llh = llHooder

        if not isinstance(theta24s, (np.ndarray, list, tuple)):
            raise TypeError()
        if not isinstance(theta34s, (np.ndarray, list, tuple)):
            raise TypeError()
        if not isinstance(msqs, (np.ndarray, list, tuple)):
            raise TypeError()

        self.theta24s = theta24s
        self.theta34s = theta34s
        self.msqs = msqs

        self._compare = False
        self._central = None

    def set_central(self, central):
        """
        Signals to the scanner that we're planning on comparing the expectations to data 
        """
        self._compare = True
        self._central = central 

    def scan(self):
        """
        Scan over all the data files. This might take a while! 
        """
        chi2 = np.zeros(shape=(len(self.theta24s),len(self.theta34s), len(self.msqs)))
        
        n_todo = len(self.theta24s)*len(self.theta34s)*len(self.msqs)
        counter = 0
        pcent_i = 0
        pcents = np.concatenate((np.linspace(0,95,20), np.linspace(96,100,5)))


        how_long = "might take a while" if n_todo>20000 else "should be quick"
        print("Starting the scan! This "+how_long)

        for i24 in range(len(self.theta24s)):
            for i34 in range(len(self.theta34s)):
                for jm in range(len(self.msqs)):
                    counter += 1 
                    if 100*(counter/n_todo) > pcents[pcent_i]:
                        print("...{}%".format(pcents[pcent_i]))
                        pcent_i+=1 

                    pam = SterileParams(theta13 = self.theta24s[i24],
                                        theta23 = self.theta34s[i34],
                                        msq2 = self.msqs[jm])

                    if self._compare:
                        # if we _expect_ the sterile point, what's the llh of measuring what we measure?
                        self.llh._set_central(pam)
                        llh = self.llh.get_llh(self._central)
                    else:
                        # what are the odds of measuring PAM if we expect whatever the llh'er was configured to expect
                        llh = self.llh.get_llh(pam)

                    chi2[i24][i34][jm] = -2*llh
        if self._compare:
            chi2 = chi2 - np.min(chi2)

        return {
                "theta24s":self.theta24s,
                "theta34s":self.theta34s,
                "msqs":self.msqs,
                "chi2s":chi2,
                "options":self.llh.options
                }