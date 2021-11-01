"""
This script will be used in making three-panel plots with P(nu_mu -> n_i), where we vary i over the tthree active neutrino flavors 

In each case, the initial state will be all muon. So finding the probabilities will be trivial! 

Then it also makes a plot showing the ->nu_tau and ->nu_tau_bar probabilities 
"""

from cascade.utils import SterileParams, config, gen_filename, Data
import nuSQuIDS as nsq

from cascade.raw_fluxes import raw_flux

import os 
import numpy as np
from matplotlib import pyplot as plt 
import pickle
plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))

params=SterileParams(theta13=0.1609, theta23=0.0, msq2=4.64)
#params=SterileParams(theta13=0.1609, theta23=0.2245, msq2=4.64)
#params=SterileParams(theta03=0.3555, theta13=0.1609, msq2=3.3)
#params= SterileParams()

colorscale = 'PuBu'

def state_setter(energies, zeniths, n_nu, kwargs):
    """
    This function prepares our initial flux. We set all the muon ones so we can easily extract 
        P(nu_mu -> nu_? )
        P(nu_mu_bar -> nu_?_bar )

    So everything is either a zero or a one! 
    """
    e_bins = len(energies)
    a_bins = len(zeniths)

    inistate = np.zeros(shape=(a_bins, e_bins, 2, n_nu))    

    for i_e in range(e_bins):
        for i_a in range(a_bins):
            for neut_type in range(2):
                inistate[i_a][i_e][neut_type][1] = 1.0
    return inistate 

# either generate the file or load it in! 
expected_fn = gen_filename(config["datapath"], "unitary_prob.dat", params)
if os.path.exists(expected_fn):
    final_probs = Data(expected_fn)
else:
    kwargs = {}
    kwargs["as_data"]=False
    kwargs["state_setter"] = state_setter
    kwargs["forced_filename"] = expected_fn
    final_probs = Data(raw_flux(params, kwargs=kwargs))

# now we need three keys 
curr = "CC" # doesn't matter 
neut = "nuBar"
root_flav = "Mu"
flavors = ['E', 'Mu', 'Tau']
# flavor neutrino current 
keys = ["_".join([flav, neut, curr]) for flav in flavors]

core_b = -0.98
mantle_b= -0.83

energies = np.logspace(2, 7, 200)
angles = np.linspace(-1,0.2, 100)

this_key = "Mu_nuBar_CC"
flux = [[ final_probs.get_flux(energy*(1e9),this_key, angle=angle) for angle in angles] for energy in energies]
plt.figure(figsize=(7,5))
mesh = plt.pcolormesh(angles, energies, flux, cmap=colorscale, vmin=0, vmax=1)
cbar = plt.colorbar(mesh)
cbar.set_label(r"$P(\bar{\nu}_{\mu}\to\bar{\nu}_{\mu})$")
plt.xlim([-1,0.2])
plt.ylim([1e2,1e6])
plt.yscale('log')
plt.ylabel(r"$E_{\nu}^{true}$")
plt.xlabel(r"$\cos\theta_{z}^{true}$")
plt.vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
plt.text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
plt.vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
plt.text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')
plt.tight_layout()
plt.savefig("nubar_survival.png", dpi=400)
plt.show()


fig, axes = plt.subplots(3,1,figsize=(5,10),sharex=True, gridspec_kw={'height_ratios':[1,1,1.45]})
fig.subplots_adjust(hspace=0.1)
fig.subplots_adjust(left=0.2)

for i in range(len(keys)):
    if neut=="nuBar":
        if i==0:
            stylized = r"$\bar{\nu}_{e}$"
        elif i==1:
            stylized = r"$\bar{\nu}_{\mu}$"
        elif i==2:
            stylized = r"$\bar{\nu}_{\tau}$" 
        elif i==3:
            stylized = r"$\bar{\nu}_{s}$"
        else:
            raise ValueError("What's goin on with this")
    else:
        if i==0:
            stylized = r"$\nu_{e}$"
        elif i==1:
            stylized = r"$\nu_{\mu}$"
        elif i==2:
            stylized = r"$\nu_{\tau}$" 
        elif i==3:
            stylized = r"$\nu_{s}$"
        else:
            raise ValueError("What's goin on with this")
    


    print(r"Doing $P(\nu_{\mu}\to$" + stylized + ")")

    axes[i].set_ylabel(r"$E_{\nu}^{true}$")
    axes[i].set_yscale('log')
    axes[i].set_ylim([1e2,1e6])
    axes[i].set_xlim([-1,0.2])
    flux = [[ final_probs.get_flux(energy*(1e9), keys[i], angle=angle) for angle in angles] for energy in energies]

    axes[i].text(0.0, 1e5, stylized, color='white' if i==1 else 'black', fontsize='x-large')
    pt = axes[i].pcolormesh(angles, energies, flux, cmap=colorscale, vmin=0.0, vmax=1.0)
    #axes[i].vlines(core_b,ymin=1e2, ymax=10**6, colors="white", ls="-")
    #axes[i].text(core_b+0.02, 1.5e2, "Inner/Outer Core Bdr",fontsize="x-small",rotation='vertical',color='white')
    axes[i].vlines(mantle_b,ymin=1e2, ymax=10**6, colors="white", ls="--")
    #axes[i].text(mantle_b+0.02, 1.5e2, "Core/Mantle Bdr",fontsize="x-small",rotation='vertical',color='white')

cbar = plt.colorbar(pt, orientation = 'horizontal', pad=0.22)
cbar.set_label(r"P($\bar{\nu}_{\mu}\to \bar{\nu}_{\alpha}$ )")
axes[-1].set_xlabel(r"$\cos\theta_{z}^{true}$")
#plt.tight_layout()
plt.savefig("survival_odds_{}.png".format(params),dpi=400)
plt.show()


# okay now we have the part that does the nu and nubar plot! 


fig, axes = plt.subplots(1,2,figsize=(12,5.5),sharey=True, gridspec_kw={'width_ratios':[1,1.125]})
fig.subplots_adjust(hspace=0.)

nus = ["nu", "nuBar"]
flavor = "Tau"
keys = ["_".join([flavor, nu, curr] ) for nu in nus]
for i in range(len(keys)):
    if i ==0 :
        stylized = r"$\nu_{\tau}$" 
    elif i==1:
        stylized = r"$\bar{\nu}_{\tau}$" 
    else:
        raise ValueError("How the hell...")

    axes[i].set_ylabel(r"$E_{\nu}^{true}$")
    axes[i].set_yscale('log')
    axes[i].set_ylim([1e2,1e6])
    axes[i].set_xlim([-1,0.2])
    axes[i].set_xlabel(r"$\cos\theta_{z}^{true}$")

    axes[i].text(0., 3e5, stylized, color='black', fontsize='xx-large')
    flux = [[ final_probs.get_flux(energy*(1e9), keys[i], angle=angle) for angle in angles] for energy in energies]

    pt = axes[i].pcolormesh(angles, energies, flux, cmap=colorscale, vmin=0.0, vmax=1.0)

cbar = plt.colorbar(pt)
cbar.set_label(r"P($\nu_{\mu}\to \nu_{tau}$ )")
plt.tight_layout()
plt.savefig("tau_odds_{}.png".format(params),dpi=400)
plt.show()