"""
Get the expected flux at a given sterile neutrino point
"""

import numpy as np

from cascade.utils import get_closest, SterileParams, gen_filename, config, Data
from cascade.sensitivity.eff_area_reader import build_flux

null = SterileParams(0.,0.,0.,0.)
null_flux = Data(gen_filename(config["datapath"], config["nu_flux"]+".dat", null))

from cascade.sensitivity.systematic_unc import astro_norm_unc
from cascade.sensitivity.systematic_unc import astro_shift_unc, cr_perturb
print("doing norm")
astro_norm_unc(null_flux)
print("doing shift")
astro_shift_unc(null_flux)
print("COSMIC")
cr_perturb(dnorm=0.05)
cr_perturb(dgamma=0.012)

ster = SterileParams(0., 0.1609, 0.2205, 4.47)
ster_flux = Data(gen_filename(config["datapath"], config["nu_flux"]+".dat", ster))

test = build_flux(null_flux)
#test2 = build_flux(ster_flu)

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

if False:
    plt.pcolormesh(test["a_edges"], test["e_edges"],test2["event_rate"]/test["event_rate"], vmin=0.60, vmax=1.40, cmap="coolwarm")
    plt.xlim([-1.0,0.0])
    plt.ylim([10**2,10**6])
    plt.xlabel(r"$\cos\theta$",size=14)
    plt.ylabel("Energy [GeV]", size=14)
    plt.yscale('log')
    plt.colorbar()
    plt.show()


plt.pcolormesh(test["a_edges"], test["e_edges"],100*test["stat_err"]/test["event_rate"],vmin=0, vmax=50, cmap="inferno")
plt.xlim([-1.0,0.0])
plt.ylim([10**2,10**6])
plt.xlabel(r"$\cos\theta$",size=14)
plt.ylabel("Energy [GeV]", size=14)
plt.yscale('log')
plt.colorbar()
plt.show()
