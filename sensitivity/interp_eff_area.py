"""
This file contains code for loading in the data representing the new DNN based event selection stuff 

So, we load the data in and do fits to them 
    the fits will be 4th order polynomials in log-log space 
    for the electron one, we also put in a Gaussian centered near the GR peak 

Then, we calculate the derivatives of these functions - this is easy since they're polynomials and Gaussians! 
So, we can then _call_ this object with some bin edges, and it'll return the new effective areas. Neato! 
"""

import os 
import numpy as np
from math import pi, sqrt, log10
from numpy.core.defchararray import not_equal

from scipy.optimize import curve_fit

from cascade.utils import SterileParams, get_loc
from cascade.raw_fluxes import raw_flux
from cascade.sensitivity.astro_flux_generator import generate_astr_flux

twop = 1/sqrt(2*pi)

def poly_3(x,*args):
    """
    This is for an order-3 polynomial
    """
    if not len(args)>=4:
        raise ValueError("Didn't get enough args... {}".format(len(args)))
    return args[0] + args[1]*x + args[2]*x*x  #+ args[3]*x*x*x 

def poly_3_prime(x, *args):
    """
    returns the derivative of the above poly_3 function
    """
    if not len(args)>=4:
        raise ValueError("Didn't get enough args: {}".format(args))
    return args[1] + 2*args[2]*x + 3*args[3]*x*x

mean = log10(6.3e6)
def poly_3_gauss(x,*args):
    """
    This is for an order-3 polynomial with a Gaussian superimposed 

    args
        0-3 are the first four coefficients in the 3d polyfit 
        4 is sigma, the width 
        5 is the scale 

        we fix the mean to the resonance peak
    """
    if not len(args)>=6:
        raise ValueError("Didn't get enough args... {}".format(len(args)))
    
    return args[5]*(twop/args[4])*np.exp(-0.5*(x-mean)*(x-mean)/(args[4]*args[4])) + poly_3(x,*args)

def poly_3_gauss_prime(x,*args):
    """
    Returns the derivative of the poly_3_gauss function
    """
    args[5]*(twop/args[4])*np.exp(-0.5*(x-mean)*(x-mean)/(args[4]*args[4]))*(x-mean)/(args[4]*args[4]) + poly_3_prime(x,*args)

file_dir = os.path.join(os.path.dirname(__file__), "new_areas")

class Fit:
    """
    This is a little class I use to make and hold the fits.
    Pass it the raw data, as loaded from the fit files, and pass it the flavor 

    It choose a fit type depending on the flavor, then fits the data. 
    Then you just call this object to get the fit-to value at that point! 
    """
    def __init__(self, data_raw, flavor):
        self.flavor = flavor

        self.energies = np.log10(data_raw[0])
        self.values = np.log10(data_raw[1])

        if self.flavor=="e":
            self._fit_f = poly_3_gauss
            params = (0, 0.5, -1, 0.0, 0.05, 1.5)
        else:
            self._fit_f = poly_3
            params = (0, 0.5, -1, 1)
        self.popt, self.pcov = curve_fit(self._fit_f, self.energies, self.values, params, maxfev=6000)

    def __call__(self, energy):
        """
        Energy should be in GeV! 
        Returns the best fit effective area there 
        """
        return 10**(self._fit_f(np.log10(energy), *self.popt))

    @property 
    def d(self):
        """
        Returns the derivative function
        """
        if self._fit_f==poly_3_gauss:
            return poly_3_gauss_prime
        else:
            return poly_3_prime

class AreaFit:
    """
    This is used to load in the published effective areas from that one thesis (https://ui.adsabs.harvard.edu/abs/2018PhDT........17N/abstract)
    We do fits to the effective areas, a fit for each angular bin and flavor, and store them

    Then you can get the effective area at arbitrary energy by making a call to this object! 
    """
    def __init__(self):
        """
        
        """
        self.a_edges = np.linspace(-1,1,11) # force this to match the files we're given 
        flavors = ["e", "mu", "tau"]

        self.fits = {}

        for flav in flavors:
            sub = os.path.join(file_dir, "nu"+flav)
            self.fits[flav] = []
            for i in range(len(self.a_edges)-1):
                angle_str = "{:.1f}_{:.1f}".format(self.a_edges[i], self.a_edges[i+1])
                filename = os.path.join(sub, "cz_"+angle_str+".csv")
                
                data_raw = np.transpose(np.loadtxt(filename, delimiter=","))
                
                self.fits[flav].append(Fit(data_raw, flav))

    def __call__(self, energy:float, cth:float, flavor:str):
        if flavor!="e" and flavor!="mu" and flavor!="tau":
            raise ValueError("Unrecognized flavor string: {}".format(flavor))
        
        i_ang = get_loc(cth, self.a_edges)[0]
        return self.fits[flavor][i_ang](energy)*1e4 #m2 to cm2

if __name__=="__main__":
    import matplotlib
    matplotlib.use('Qt5Agg')
    from matplotlib import pyplot as plt
    plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")
    from scipy.integrate import dblquad

    from cascade.utils import get_color
    import pickle

    flavors = ["e", "mu", "tau"]
    master_fit = AreaFit()

    do_fudge = True
    if do_fudge:
        core_b = -0.98
        mantle_b= -0.83


        # let's see how the new predicted event rate looks! 
        null = SterileParams()
        kwargs = {}
        kwargs["as_data"]=True
        atmo_data = raw_flux(null,kwargs=kwargs)
        astr_data = generate_astr_flux(null, as_data=True)
        n_e = 20
        n_a = 10
        e_bins = np.logspace(2,8,n_e+1)
        a_bins = np.linspace(-1,1,n_a+1)
        events = np.zeros(shape=(n_e, n_a))

        time = 3600.*24*365*10

        def funcy(energy, angle, flav):
            """
            This little function sums the contributions for each flux making up this combined flux 
            """
            nus = ["nu", "nuBar"]
            this_flav = flav[0].upper() + flav[1:].lower()
            keys = ["_".join([this_flav, nu, "NC"]) for nu in nus]
            net_f = 0.0
            # we take the average contribution from nu and nubar
            dobjs = (atmo_data, astr_data)
            for dobj in range(len(dobjs)):
                for key in keys:
                    ff = master_fit(energy, angle, flav)
                    if ff<0:
                        raise ValueError("Area fit {}".format(ff))
                    flux = 0.5*dobjs[dobj].get_flux((1e9)*energy, key, use_overflow=False, angle=angle)
                    if flux <0:
                        print("Flux is {} in {} flux".format(flux, "atmo" if dobj==0 else "astro"))
                    net_f = net_f + flux*ff
            return net_f*2*pi

        print("Doing those integrals now")
        for i_a in range(n_a):
            for i_e in range(n_e):
                for flav in flavors:
                    integral = dblquad(lambda energy, angle:funcy(energy, angle, flav), a=a_bins[i_a], b=a_bins[i_a+1], gfun=e_bins[i_e], hfun=e_bins[i_e+1])[0]*time
                    if integral<0:
                        print("{} {} {}".format(i_e, i_a, flav))
                        raise ValueError("Negative value in integral: {}".format(integral))

                    events[i_e][i_a] += integral

        plt.pcolormesh(a_bins, e_bins, events)
        #plt.xlim([-1,0.2])
        #plt.ylim([1e2,1e6])
        plt.yscale('log')
        plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
        plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
        plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
        plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
        for i_e in range(n_e):
            for i_a in range(n_a):
                plt.text(a_bins[i_a]+0.05,e_bins[i_e], "{:.1f}".format(events[i_e][i_a]), fontsize='x-small', color='white')

        cbar= plt.colorbar()
        cbar.set_label("Events")
        print("{} total events".format(np.sum(events)))
        plt.show()

        f_name_s = "/home/benito/software/data/cascade/hg_sib/expected_fluxes_reco/0.0/expected_flux_smearedwell_0.0_0.0_0.0_0.0.dat"
        f = open(f_name_s, 'rb')
        data = pickle.load(f)
        f.close()
        rate = data["event_rate"]


        np.save("full_fudge",events/rate)

        plt.pcolormesh(a_bins, e_bins, events/rate)
        #plt.xlim([-1,0.2])
       # plt.ylim([1e2,1e6])
        plt.yscale('log')
        plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
        plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
        plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
        plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
        for i_e in range(n_e):
            for i_a in range(n_a):
                plt.text(a_bins[i_a]+0.1,e_bins[i_e], "{:.1f}".format(events[i_e][i_a]/rate[i_e][i_a]), fontsize='x-small', color='white')

        cbar= plt.colorbar()
        cbar.set_label("New Fudge")
        plt.show()

    if not do_fudge:
        e_range = np.logspace(2,8,1000)


        for flav in flavors:
            fig = plt.figure()
            i_count = 0 
            for fit in master_fit.fits[flav]:
                ys = fit(e_range)
                plt.plot(e_range, ys, color=get_color(i_count, 9))
                i_count +=1 
            plt.title(flav,size=16)
            plt.xscale('log')
            plt.yscale('log')
            plt.xlim([1e2,1e8])
            plt.ylim([1e-2,1e3])
            plt.show()