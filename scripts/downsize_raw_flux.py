"""
Discovers all the raw fluxes in the configuraed data folder, then downsizes them for faster plotting
"""

from glob import glob

import numpy as np
from math import log10
import pickle
import os 

from cascade.utils import config
from cascade.utils import bilinear_interp
from cascade.utils import Data

GeV = 1e9
PeV = (1e6)*GeV

Emin = 10*GeV
Emax = 10*PeV

zmin = -0.99
zmax = -0.01

Ebin = 100
zbin = 100

energies = np.logspace(log10(Emin), log10(Emax), Ebin)
zeniths = np.linspace(zmin, zmax, zbin)

def downsize(filename, destination):
    if os.path.exists(destination):
        print("Already done! Skipping {}".format(destination))
        return

    dataobj = Data(filename)

    flux = {}
    for key in dataobj.get_keys():
        flux[key] = np.zeros(shape=(Ebin, zbin))
        for i_e in range(Ebin):
            for i_z in range(zbin):
                flux[key][i_e][i_z] = dataobj.get_flux( energies[i_e], key=key, angle=zeniths[i_z])

    # okay now we pickle 
    all_data = {}
    all_data["e_true"] = energies
    all_data["a_true"] = zeniths
    all_data["flux"] = flux
    f = open(destination, 'wb')
    pickle.dump(all_data, f, -1)
    f.close()
    print("Saved {}".format(destination))

def main():
    # get the list of files to work with 
    which = os.path.join(config["datapath"], config["nu_flux"])
    all_files = glob(which+"*")
    print("Found {} files".format(len(all_files)))
    
    for each in all_files:
        dirname, filename = os.path.split(each)
        
        newname = ".".join(filename.split(".")[:-1])
        newname = "_".join([config["nu_flux_downsize"]] + newname.split("_")[3:])
        newname = os.path.join(config["datapath"], newname)
        try:
            downsize(filename, newname+".dat")
        except:
            print("Failed to read parse {}".format(filename))
            continue

if __name__=="__main__":
    main()
