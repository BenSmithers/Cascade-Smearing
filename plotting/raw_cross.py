import pickle
import numpy as np
import matplotlib 
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from cascade.utils import SterileParams, config, gen_filename
from cascade.utils import bhist

from cascade.cross_section_test import get_total_flux as get_xs 
from cascade.nus_utils import get_flavor, get_neut, get_curr 

from cascade.tau_funcs import TauData
TD = TauData()

width = 0.1

def _load_flux(params):
    name = gen_filename(config["datapath"], config["nu_flux_downsize"]+".dat", params)

    f = open(name,'rb')
    all_data = pickle.load(f)
    f.close()

    return( all_data["e_true"], all_data["a_true"], all_data["flux"] )

just_nubar = False
keep_key = "Mu_nuBar_"


null = SterileParams(0., 0., 0., 0.)
ster = SterileParams(0., 0.1609, 0.2296, 4.47)

e_true, a_true, flux_null    = _load_flux(null)
e_true, a_true, flux_sterile = _load_flux(ster)

keys = flux_null.keys()
centers = bhist([e_true]).centers


for key in keys:
    flav = get_flavor(key)
    neut = get_neut(key)
    curr = get_curr(key) 

    for i_energy in range(len(centers)):
        xs = get_xs(centers[i_energy], flav, neut, curr)
    
        if "nubar" in key.lower():
            P = -1
        else:
            P = 1

        pcent = 1.0
        if "tau" in key.lower():
            pcent = (TD( centers[i_energy]/(1e9), P) - centers[i_energy])/centers[i_energy]

        flux_null[key][i_energy] *= xs*pcent
        flux_sterile[key][i_energy] *= xs*pcent

ex = list(flux_null.keys())[0]
null_total      = np.zeros(shape = np.shape(flux_null[ex]))
sterile_total   = np.zeros(shape = np.shape(flux_null[ex]))

for key in keys:
    if just_nubar and (keep_key not in key):
        print("Skip {}".format(key))
        continue
    null_total += flux_null[key]
    sterile_total += flux_sterile[key]

ratio = sterile_total / null_total
energies = np.array(bhist([e_true]).centers)
czeniths = np.array(bhist([a_true]).centers)

cf = plt.pcolormesh(a_true, e_true/(1e9), ratio, cmap=cm.coolwarm, vmin=1.0-width, vmax=1.0+width)
plt.yscale('log')
plt.ylabel("Energy [GeV]", size=14)
plt.xlabel(r"$\cos\theta$",size=14)
if just_nubar:
    plt.title("Only Looking at: {}".format(keep_key),size=14)
cbar = plt.colorbar() #cf,ticks=ticker.LogLocator())
cbar.set_label(r"Sterile Flux-XS / Null Flux-XS")
plt.show()