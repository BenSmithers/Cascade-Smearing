"""
This file contains code for loading in the data representing the new DNN based event selection stuff 

So, we load the data in and do fits to them 
    the fits will be 4th order polynomials in log-log space 
    for the electron one, we also put in a Gaussian centered near the GR peak 

Then, we calculate the derivatives of these functions - this is easy since they're polynomials and Gaussians! 
So, we can then _call_ this object with some bin edges, and it'll return the new effective areas. Neato! 
"""

import os
from math import inf, log10, pi, sqrt

import numpy as np
from cascade.raw_fluxes import raw_flux
from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.utils import SterileParams, get_loc
from scipy.optimize import curve_fit

twop = 1/sqrt(2*pi)

def poly_3(x,*args):
    """
    This is for an order-3 polynomial
    """
    if not len(args)>=4:
        raise ValueError("Didn't get enough args... {}".format(len(args)))
    return args[0] + args[1]*x + args[2]*x*x*x*x  #+ args[3]*x*x*x

def sad_hat(x,*args):
    if not len(args)>=4:
        raise ValueError("Didn't get enough args... {}".format(len(args)))
    return args[0]*np.log(args[3]*(x-1))/(args[1]*(x-1))

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
    
    return args[5]*(twop/args[4])*np.exp(-0.5*(x-mean)*(x-mean)/(args[4]*args[4])) + sad_hat(x,*args)

def poly_3_gauss_prime(x,*args):
    """
    Returns the derivative of the poly_3_gauss function
    """
    return args[5]*(twop/args[4])*np.exp(-0.5*(x-mean)*(x-mean)/(args[4]*args[4]))*(x-mean)/(args[4]*args[4]) + poly_3_prime(x,*args)

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
        if False:
            if self.flavor=="e":
                self._fit_f = poly_3_gauss
                params = (0, 0.5, -1, 0.0, 0.05, 1.5)
            else:
                self._fit_f = poly_3
                params = (0, 0.5, -1, 1)
        else:
            if True: #self.flavor=="e":
                self._fit_f = poly_3_gauss
                params = (2, 1.5, 1.5, 2.0, 0.5, 1.5)
                bounds_l = (-inf, -inf,-inf,-inf,-inf,-inf)
                bounds_u = (inf, 3, 3, inf, inf, inf)
            else:
                self._fit_f = sad_hat
                params = (2, 1.5, 1.5, 2.0)
                bounds = [(-inf, inf) for i in params]
                bounds_l = (-inf, -inf,-inf,-inf)
                bounds_u = (inf, 1.7, 1.7, inf)
        
        self.popt, self.pcov = curve_fit(self._fit_f, self.energies, self.values, params, maxfev=7000)

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
    def __init__(self, debugMode=False):
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
                
                data_raw = sorted(np.loadtxt(filename, delimiter=","), key=lambda x: x[0])
                data_raw = np.transpose(data_raw)
                
                self.fits[flav].append(Fit(data_raw, flav))

                if debugMode:
                    dbx = np.logspace(2,8,2000)
                    dby = self.fits[flav][-1](dbx)
                    plt.plot(dbx,dby, label="fit")
                    plt.plot(data_raw[0], data_raw[1], label="data")
                    plt.xscale('log')
                    plt.legend()
                    plt.yscale('log')
                    plt.ylim([1e-4, 5e2])
                    plt.title("{}, {}".format(flav.upper(), angle_str))
                    plt.show()
                    

    def __call__(self, energy:float, cth:float, flavor:str):
        if flavor!="e" and flavor!="mu" and flavor!="tau":
            raise ValueError("Unrecognized flavor string: {}".format(flavor))
        
        i_ang = get_loc(cth, self.a_edges)[0]
        return self.fits[flavor][i_ang](energy)*1e4 #m2 to cm2

master_fit = AreaFit(False)

def get_expectation(*args):
    flavors = ["e", "mu", "tau"]

    n_e = 20
    n_a = 10
    e_bins = np.logspace(2,8,n_e+1)
    a_bins = np.linspace(-1,1,n_a+1)
    events = np.zeros(shape=(n_e, n_a))

    electron = np.zeros(shape=(n_e, n_a))
    muon = np.zeros(shape=(n_e, n_a))
    tau = np.zeros(shape=(n_e, n_a))


    def metaflux(energy, angle, flav):
        """
        This little function sums the contributions for each flux making up this combined flux 
        """
        nus = ["nu", "nuBar"]
        this_flav = flav[0].upper() + flav[1:].lower()
        keys = ["_".join([this_flav, nu, "NC"]) for nu in nus]
        net_f = 0.0
        # we take the average contribution from nu and nubar
        for dobj in range(len(args)):
            for key in keys:
                ff = master_fit(energy, angle, flav)
                if ff<0:
                    raise ValueError("Area fit {}".format(ff))
                flux = 0.5*args[dobj].get_flux((1e9)*energy, key, use_overflow=False, angle=angle)
                if flux <0:
                    print("Flux is {} in {} flux".format(flux, "atmo" if dobj==0 else "astro"))
                net_f = net_f + flux*ff
        return net_f*2*pi

    for flav in flavors:
        for i_a in range(n_a):
            for i_e in range(n_e):

                integral = dblquad(lambda energy, angle:metaflux(energy, angle, flav), a=a_bins[i_a], b=a_bins[i_a+1], gfun=e_bins[i_e], hfun=e_bins[i_e+1])[0]
                if integral<0:
                    print("{} {} {}".format(i_e, i_a, flav))
                    raise ValueError("Negative value in integral: {}".format(integral))

                events[i_e][i_a] += integral
                if flav=="e":
                    electron[i_e][i_a] += integral
                elif flav=="mu":
                    muon[i_e][i_a]+= integral
                else:
                    tau[i_e][i_a]+=integral

    livetime_ratio = 1.0
    seconds_up = livetime_ratio*10.*3600.*24.*365
    casc_flux = events*seconds_up
    casc_err = np.sqrt(casc_flux)

    return {
            "e_edges":e_bins,
            "a_edges":a_bins,
            "event_rate":casc_flux,
            "stat_err":casc_err,
            "electron":electron*seconds_up,
            "muon":muon*seconds_up,
            "tau":tau*seconds_up
           }

if __name__=="__main__":
    import matplotlib
    matplotlib.use('Qt5Agg')
    from matplotlib import pyplot as plt
    plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")
    

    import pickle

    from cascade.utils import get_color

    
    do_fudge = True
    flavors = ["e", "mu", "tau"]
    
    if do_fudge:
        core_b = -0.98
        mantle_b= -0.83


        # let's see how the new predicted event rate looks! 
        null = SterileParams()
        not_null= SterileParams(theta13=0.1609, theta23=0.2249, msq2=4.47)
        kwargs = {}
        kwargs["as_data"]=True
        atmo_data = raw_flux(null,kwargs=kwargs)
        astr_data = generate_astr_flux(null, as_data=True)
        n_e = 20
        n_a = 10

        raw_data = get_expectation(atmo_data, astr_data)

        e_bins = raw_data["e_edges"]
        a_bins = raw_data["a_edges"]
        events = raw_data["event_rate"]

        plt.pcolormesh(a_bins, e_bins, events)
        plt.xlim([-1,0.2])
        plt.ylim([1e2,1e6])
        plt.yscale('log')
        plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
        plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
        plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
        plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
        for i_e in range(n_e):
            if e_bins[i_e]>1e6:
                continue
            for i_a in range(n_a):
                if a_bins[i_a]>0.2:
                    continue
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
        
        e = raw_data["electron"]
        mu = raw_data["muon"]
        tau = raw_data["tau"]

        # do the actual rate stuff now 
        hertz = np.array([ sum(e_e) for e_e in events])/(10*365*24*3600)
        hertz_e = np.array([ sum(e_e) for e_e in e])/(10*365*24*3600)
        hertz_mu = np.array([ sum(e_e) for e_e in mu])/(10*365*24*3600)
        hertz_tau = np.array([ sum(e_e) for e_e in tau])/(10*365*24*3600)


        e_center = 0.5*(e_bins[1:] + e_bins[:-1])
        plt.plot(e_bins[:-1], hertz, drawstyle='steps',color="k", label="total")
        plt.plot(e_bins[:-1], hertz_e, drawstyle='steps', label=r"$\nu_e$")
        plt.plot(e_bins[:-1], hertz_mu, drawstyle='steps', label=r"$\nu_\mu$")
        plt.plot(e_bins[:-1], hertz_tau, drawstyle='steps', label=r"$\nu_\tau$")

        plt.xscale('log')
        plt.ylabel("Rate [Hz]")
        plt.xlabel(r"$E_{reco}$ [GeV]")
        plt.yscale('log')
        plt.xlim([1e2,1e7])
        plt.ylim([5e-10, 5e-4])
        plt.legend()
        plt.show()


        np.save("full_fudge",events/rate)

        plt.pcolormesh(a_bins, e_bins, events/rate)
        plt.xlim([-1,0.2])
        plt.ylim([1e2,1e6])
        plt.yscale('log')
        plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
        plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
        plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
        plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
        for i_e in range(n_e):
            if e_bins[i_e]>=1e6:
                continue
            for i_a in range(n_a):
                if a_bins[i_a]>=0.2:
                    continue
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
