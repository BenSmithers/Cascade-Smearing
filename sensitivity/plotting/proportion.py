from cascade.sensitivity.astro_flux_generator import generate_astr_flux
from cascade.sensitivity.eff_area_reader import build_flux, build_flux_sad

from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import Data
from cascade.raw_fluxes import raw_flux

from cascade.sensitivity.load_expectation import load as loade

import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as plt 
plt.style.use("/home/benito/software/cascade/cascade/cascade.mplstyle")
import numpy as np

"""
Comparing relative predicted fluxes for astrophysical and atmostpheric neutrinos 
"""

params = SterileParams()

kwargs = {}
kwargs["as_data"]=True

atmo_data = raw_flux(params,kwargs=kwargs)
astr_data = generate_astr_flux(params, as_data=True)

atmo_flux = build_flux_sad(atmo_data)
astr_flux = build_flux_sad(astr_data)

lt_scale = (641/365)/10.
atmo_summ = np.array([entry[0]*lt_scale for entry in atmo_flux["event_rate"]])
astr_summ = np.array([entry[0]*lt_scale for entry in astr_flux["event_rate"]])

e_edges = np.array(atmo_flux["e_edges"])
widths = e_edges[1:] - e_edges[:-1]

figs,axes = plt.subplots(nrows=2, ncols=1, sharex=True, gridspec_kw={'height_ratios':[3,1]})


axes[0].bar(x=e_edges[:-1], height=atmo_summ, width=widths, align='edge', label=r"$\nu_{atm}$ Exp",alpha=0.6)
axes[0].bar(x=e_edges[:-1], height=astr_summ, width=widths, align='edge',bottom=atmo_summ, label=r"$\nu_{astr}$ Exp",alpha=0.6)
axes[1].set_xlabel("Energy [GeV]")
axes[0].set_ylabel("N Events")
axes[0].set_xscale('log')
axes[0].set_xlim([5e2, 1e7])
axes[0].set_yscale('log')
axes[0].set_ylim([1e-1, 1e2])

occupation = np.array([entry[0] for entry in loade(True)])

axes[0].errorbar(x=0.5*(e_edges[1:]+e_edges[:-1]), y=occupation,xerr=0.5*(e_edges[1:]-e_edges[:-1]),yerr=np.sqrt(occupation),color='black',capsize=3, label="Observed",marker='d',ls='')
#plt.ylim([0.85,1.15])

axes[0].grid(which='major',alpha=0.4)
axes[1].set_ylabel("Obs/Exp")
axes[1].plot(e_edges[:-1], occupation/(atmo_summ+astr_summ), drawstyle='steps')
axes[1].set_ylim([0,2])
axes[0].legend()
plt.show()

