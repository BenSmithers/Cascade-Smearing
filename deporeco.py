from math import exp, sqrt, pi, log, acos, log10, sinh, sin
import numpy as np
import os
from cascade.utils import bhist, get_loc, get_closest, config

import pickle
import scipy.optimize as opt
import sys 
"""
This script is here to approximate the uncertainties in going from "energy deposited" to "energy reconstructed"
Essentially this code here incorporates the detector reconstruction step
"""

# square roots are expensive. Let's only do this once 
rtwo = 1./sqrt(2*pi)

datafile = "reconstruction.txt"
data = np.loadtxt(os.path.join(os.path.dirname(__file__), datafile), dtype=float, delimiter=",")
data = data.transpose()
#data[1] = data[1]*data[0] 
# data[0] is energy
# data[1] is sigma



rtwolog = sqrt(2*log(2))

# Load in the datafile that has the angular error as a function of particle energy 
angdata = np.transpose(np.loadtxt(os.path.join(os.path.dirname(__file__),"angerror.txt"), delimiter=","))

# the fit is best when done in log-space of energy 
# so we do a little fitty-fit
logged = np.log10(angdata[0])
from scipy.optimize import curve_fit as fit
def func(x, a, gamma):
    return( a*(x**(-gamma)) )
popt, pcov = fit(func, logged, angdata[1])


def get_ang_error(energy):
    """
    Returns a best-guess for the angular error - the FIFTYth percentile

    Energy should be a float and in GeV
    """
    if not isinstance(energy, (int,float)):
        raise TypeError("Expected {}, got {}".format(float, type(energy)))
    
    try:
        # we really want to just use the data, so let's try interpolating first 
        error = get_closest( energy, angdata[0], angdata[1] )
    except ValueError: # failing that we use the fit (which isn't as good)
        error = func( log10(energy), popt[0], popt[1])
    return(min(max(error, 0), 180))

if False:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    plt.plot( angdata[0], angdata[1],'d', label="Data")
    xs = np.logspace(3.5,7.5,100)
    #ys = func(np.log10(xs), popt[0], popt[1])
    ys = [get_ang_error(en) for en in xs]
    plt.plot(xs, ys, label='Fit')
    plt.legend()
    plt.xscale('log')
    plt.xlabel("Energy Depo [GeV]",size=14)
    plt.ylabel("Angular Error [deg]",size=14)
    plt.show()
    sys.exit()

class DataReco:
    """
    This is the object that actually facilitates the energy reconstruction smearing 
    
    You pass it the edges of bins you have for energy/angle deposited/reconstructed, and then it builds up the normalized probabilities
    We keep this as a single object since building it is expensive, but accessing the probabilities is cheap! 
    """
    def __init__(self, reco_energy_edges, reco_czenith_edges, depo_energy_edges, true_czenith_edges):
        """
        Expects the energies in eV, but works in GeV
        """


        self._ereco = bhist([np.array(reco_energy_edges)*(1e-9)])
        self._edepo = bhist([np.array(depo_energy_edges)*(1e-9)])
        self._zreco = bhist([reco_czenith_edges]) # these two are in cos(Zenith)
        self._ztrue = bhist([true_czenith_edges])

        # these are now filled with the values of the probability DENSITY functions for each angle/energy combo 
        # TODO right now there is no assumed covariance ... this should be improved 
        self._energy_odds_array = np.array([[ get_odds_energy(deposited, reconstructed) for reconstructed in self.reco_energy_centers] for deposited in self.depo_energy_centers])
#        self._angle_odds_array = np.array([[ get_odds_angle(true, reconstructed) for reconstructed in self.reco_czenith_centers] for true in self.true_czenith_centers]) 
        self._angle_odds_array = np.array([[[ get_odds_angle(true, reconstructed, deposited) for reconstructed in self.reco_czenith_centers] for deposited in self.depo_energy_centers]for true in self.true_czenith_centers])

        # Avoid equating two floats. Look for a sufficiently small difference! 
        max_diff = 1e-12

        # normalize these things! 
        # so for each energy deposited... the sum of (PDF*width evaluated at each reco bin) should add up to 1. 
        for depo in range(len(self._energy_odds_array)):
            self._energy_odds_array[depo] *= 1./sum(self._energy_odds_array[depo]*self.reco_energy_widths)
            assert(abs(1-sum(self._energy_odds_array[depo]*self.reco_energy_widths)) <= max_diff)
    
        for true in range(len(self._angle_odds_array)):
            # for each of the possible true values
            for deposited in range(len(self._angle_odds_array[true])):
                # and each of the possible energies deposited
                self._angle_odds_array[true][deposited] *= 1./sum(self._angle_odds_array[true][deposited]*self.reco_czenith_widths)
                assert(abs(1-sum(self._angle_odds_array[true][deposited]*self.reco_czenith_widths)) <= max_diff)

    # these two are functions used to access those probabilities. 
    # they both take bins numbers
    # These use bin numbers _specifically_ to ensure that the user recognizes they need to use the exact same bins as they used to construct these 
    # these bin numbers should be for the bin centers/widths, NOT the edges
    def get_energy_reco_odds(self, i_depo, i_reco ):
        return(self._energy_odds_array[i_depo][i_reco]*self.reco_energy_widths[i_reco])

    def get_czenith_reco_odds(self, i_true, i_reco, i_e_true):
        return(self._angle_odds_array[i_true][i_e_true][i_reco]*self.reco_czenith_widths[i_reco])
   
    # here we have afew access functions to see what we used to build this (and what their associated bin centers/widths are)
    @property
    def reco_energy_centers(self):
        return(self._ereco.centers)
    @property
    def reco_czenith_centers(self):
        return(self._zreco.centers)
    @property
    def depo_energy_centers(self):
        return(self._edepo.centers)
    @property
    def true_czenith_centers(self):
        return(self._ztrue.centers)

    @property
    def reco_energy_widths(self):
        return(self._ereco.widths)
    @property
    def reco_czenith_widths(self):
        return(self._zreco.widths)
    @property
    def depo_energy_widths(self):
        return(self._edepo.widths)
    @property
    def true_czenith_widths(self):
        return(self._ztrue.widths)

# do some math, keep this in the global scope 
deg = np.pi/180.
mu = deg*6.1983e-1
fwhm = deg*(13.223 + 11.983)
sigma =  fwhm/rtwolog

def check_angle(value, cosmode=False):
    """
    A little functino to make sure that whatever is passed is both a number and in an acceptable range 
    """
    if not isinstance(value, (int, float)):
        raise TypeError("Expected {}, got {}".format(float, type(value)))
    if cosmode:
        if not (value<=1. and value>=-1.):
            raise ValueError("cos(theta) outside allowed bounds: {}".format(value))
    else:
        pass # value can be any number of radians 

class KappaGrabba:
    """
    This object calculates the smearing factor "kappa" for a bunch of possible angular errors, then saves them in a data file 
    """
    def __init__(self):
        self.numb = 100

        self.rad_range = np.linspace(5*pi/180., pi/2,self.numb)
        self.czen_range = np.cos(self.rad_range)

        self.datafile = os.path.join(config["datapath"], os.path.join(config["kappas"]))

        if os.path.exists(self.datafile):
            self.kappas = self._load()
        else:
            self.kappas = self.generate()
            self._save()

    def _load(self):
        f = open(self.datafile, 'rb')
        data = pickle.load(f)
        f.close()
        return(data)

    def _save(self):
        f = open(self.datafile, 'wb')
        pickle.dump(self.kappas, f, -1)
        f.close()
        print("Kappa values saved")
        
    def generate(self):
        """
        Should only be called once ever. This calculates all the kappa values! 
        """
        print("Generating Kappas for Angular uncertainty")
        kappas = np.zeros(shape=np.shape(self.rad_range))

        for i_czenith in range(len(self.czen_range)):
            czenith = self.czen_range[i_czenith]
            def funct(kappa):
                if kappa<=0:
                    return(100000)
                else:
                    return (exp(kappa) - exp(kappa*czenith))/(2*sinh(kappa)) - 0.5
            
            soln = opt.root(funct, np.array([10.]))
            kappas[i_czenith] = soln.x[0]

        return(kappas)

    def get_kappa(self, rad_error):
        """
        Accesses the stored calculated values for kappa
        """
        if not isinstance(rad_error, (int, float)):
            raise TypeError("Receied invalid datatype for rad error: {}".format(type(rad_error)))
        try:
            value = get_closest(rad_error, self.rad_range, self.kappas)
            return(value)
        except ValueError: # angle too far a way 
            if rad_error <= (5*pi/180.):
                return(max(self.kappas))
            elif rad_error>=(pi/2):
                return(min(self.kappas))
            else:
                raise Exception("Not sure how to work with {}".format(rad_error))

kappaCalc = KappaGrabba()

def get_spread(ang_unc):
    """
    This calculates "kappa" in the Kent distribution

    Value should be passed in RADIANS
    Returns unitless quantity
    """

    cos_unc = cos(ang_unc)
    inter = log(1./(exp(cos_unc)*sinh(1)/sinh(cos_unc) - 1))

def get_odds_angle( reco, true, energy_depo ):
    """
    Uses kent_pdf to get the zenith error. Uses a modified verion of Alex's work 

    MBSF: 
    The first three terms in the expansion of the Modifed Bessel Function of the First Kind where nu=0
    I_{0}(x)
    See: https://www.wolframalpha.com/input/?i=I_0%28x%29
    And: https://mathworld.wolfram.com/ModifiedBesselFunctionoftheFirstKind.html
    """
    check_angle(reco, True)
    check_angle(true, True)
    if not isinstance(energy_depo, (int, float)):
        raise TypeError("Energy of type {}, not {}".format(type(energy_depo), float))
    
    error = kappaCalc.get_kappa(get_ang_error( energy_depo )*deg)
    dreco = acos(reco)/deg
    dtrue = acos(true)/deg

    #return (rtwo/error)*exp(-0.5*pow((dreco-dtrue)/error,2))

    #this is the bessel part of the kent_pdf
    value = sqrt(1-reco*reco)*sqrt(1-true*true)
    mbsf = 1 + pow(value,2)/4 + pow(value,4)/64 # we use the polynomial verion since it should be faster
    
    try:
        prob = error/(2*sinh(error))
    except OverflowError:
        print(error)
        sys.exit()

    # Based on Alex's Kent Distribution work 
    return(prob*exp(error*reco*true)*mbsf)

def get_odds_angle_deprecated(true, reconstructed, energy_depo):
    """
    The returns the probability density of an RA/zenith angle of 'true' being reconstructed as 'reconstructed'

    It uses the 7-year cascade analysis' evaluation of angular reconstruction from Monopod 
    """
    sigma_fifty = get_ang_error( energy_depo )

    x = acos(reconstructed) - acos(true)
    
    value = (rtwo/sigma)*exp(-1*((x-mu)**2)/(2*sigma*sigma))
    return(value)


def get_odds_energy(deposited, reconstructed):
    """
    Takes an energy deposited and energy reconstructed.
    Loads the datafile and reads off the uncertainty from the second column. 
    The data is in %E_depo, so we scale this to be a ratio and then by that deposted energy
    """
    if not isinstance(deposited, (float,int)):
        raise Exception()
    if not isinstance(reconstructed, (float,int)):
        raise Exception()


    sigma = get_closest( deposited, data[0], data[1])*deposited*0.01
    
    s2 = np.log(1+ (sigma/deposited)**2)
    mu = np.log((deposited**2)/np.sqrt(deposited**2  + sigma**2))
    # now, we assume that the uncertainty follows a log normal distribution, and calculate the PDF here

    #prob = rtwo*(1./sigma)*exp(-0.5*((reconstructed - deposited)/sigma)**2)
    prob = rtwo*(1./np.sqrt(s2))*exp(-0.5*((np.log(reconstructed) - mu)**2)/s2)/reconstructed

    return(prob)

# this was here just for debugging purposes. Don't usually need it
doplot = False
if doplot:
    import matplotlib
    matplotlib.use('TkAgg')
    import matplotlib.pyplot as plt
    from matplotlib import ticker

    depos_e = np.logspace(0.1, 6.9, 100)*(1e9)
    recos_e = np.logspace(0.1, 6.9, 101)*(1e9)

    ang = np.linspace(-1,1,100)

    dataobj = DataReco( recos_e, ang, depos_e, ang)
    
    depos = dataobj.depo_energy_centers
    recos = dataobj.reco_energy_centers

    odds = np.array([[dataobj.get_energy_reco_odds(depo, reco) for depo in range(len(depos_e)-1)] for reco in range(len(recos_e)-1)])
    odds = np.log10(np.ma.masked_where(odds<=1e-10, odds))
    levels = np.logspace(-3,0,10)
    
    if False:
        print(np.min(odds))
        print(np.max(odds))

        figs, axes = plt.subplots(nrows=2, ncols=1, sharex=False, gridspec_kw={'height_ratios':[1,1]})

        ctf = axes[0].pcolormesh(depos,recos,odds, cmap='gist_yarg')
        #plt.pcolormesh(depos,recos,odds, cmap='gist_gray', locator=ticker.LogLocator())
        axes[0].set_xscale('log')
        axes[0].set_yscale('log')
        axes[0].set_xlim([10**0,10**7])
        axes[0].set_ylim([10**0,10**7])
        axes[0].set_xlabel("Deposited [GeV]")
        axes[0].set_ylabel("Reconstructed [GeV]")
        cbar = plt.colorbar(ctf) #, ticks=ticker.LogLocator())
        cbar.set_label("Probability Density")

        # slides !
        slices = 10**np.array([1,2,3,4,5])
        for cut in slices:
            


            left, right = get_loc(float(cut), depos)
            loc = left if abs(cut-depos[left])<abs(cut-depos[right]) else right
            width_factor = dataobj.depo_energy_widths[loc]
            odds = [ dataobj.get_energy_reco_odds(loc, val) for val in range(len(recos))]
            axes[1].plot(recos, odds, color=(.1,.1,.1))
        
        axes[1].set_xscale('log')
        axes[1].set_xlabel("Reconstructed [GeV]")       
        axes[1].set_ylabel("Prob. Density")
        
        figs.savefig("fig_deporeco.png",dpi=400)
        plt.show()
        
    choice_e, unused = get_loc( 1e4, depos)

    print("Error: {}".format(get_ang_error(1e4)))

    #ang_odds = np.array([[dataobj.get_czenith_reco_odds(true,reco, choice_e) for true in range(len(ang)-1)] for reco in range(len(ang)-1)])
    #ang_odds = np.array([ ang_odds[it]/sum(ang_odds[it]) for it in range(len(ang_odds))])

    #ang_odds=np.log10(np.ma.masked_where(ang_odds<=1e-10,ang_odds))
   
    choice_a, unused = get_loc(-0.99, ang)

    plt.figure(2)
    ang_odds = np.array([dataobj.get_czenith_reco_odds( choice_a, reco, choice_e) for reco in range(len(ang)-1)])
    print("Sum: {}".format(sum(ang_odds)))

    plt.plot( bhist([ang]).centers, ang_odds)
    #me = plt.pcolormesh(ang,ang,ang_odds,cmap='gist_yarg')
    plt.xlabel("Reconstructed Angle")
    plt.ylabel("Odds, normalized")
    #cbar = plt.colorbar(me)
    #cbar.set_label("Prob Density")
    plt.savefig("fig_angtruereco.png",dpi=400)
    plt.show()
    
    fluxes = np.array([(1.0 if true<-0.9 else 0.0) for true in bhist([ang]).centers])
    fluxes /= sum(fluxes)

    smeared = np.zeros(shape=np.shape(fluxes))
    for truth in range(len(fluxes)):
        for recon in range(len(fluxes)):
            smeared[recon] += fluxes[truth]*dataobj.get_czenith_reco_odds( truth, recon, choice_e)

    print("Smear sum: {}".format(sum(smeared)))
    plt.figure(3)
    plt.plot( bhist([ang]).centers, fluxes, label="True")
    plt.plot( bhist([ang]).centers, smeared, label="Smeared")
    plt.xlabel("Angle")
    plt.ylabel("Flux Density")
    plt.show()
