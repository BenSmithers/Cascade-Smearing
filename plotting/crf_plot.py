import matplotlib
matplotlib.use("Qt5Agg")
import matplotlib.pyplot as plt

import numpy as np
from math import log10

import crflux.models as crf

from cascade.utils import get_color

hg = crf.HillasGaisser2012()
pg = crf.PolyGonato()

stuff ={14:"H",
        402:"He",
        1206:"CNO",
        2814:"MgAlSi",
        5426:"Fe"}

emin = 1e1
emax = 1e8
npoints = 100

e_grid = np.logspace(log10(emin), log10(emax), npoints)
evaled = {}
evaled_pg = {}
mag = 1

counter = 0
for cid in stuff.keys():
    cast = int(cid)


    n_nucl = hg.Z_A(cast)[1] #number of nucleons

    evaled[cast] = [hg.nucleus_flux(cast, e*n_nucl) for e in e_grid]
    evaled_pg[cast] = [pg.nucleus_flux(cast, e*n_nucl) for e in e_grid]

    plt.plot(e_grid, evaled[cast]*(e_grid**mag), label=stuff[cid], color=get_color(counter, len(list(stuff.keys()))-1), ls='-')
#    plt.plot(e_grid, evaled_pg[cast]*(e_grid**mag), color=get_color(counter, len(list(stuff.keys()))-1), ls='--')

    counter+=1

plt.xscale('log')
plt.yscale('log')
#plt.ylim([10**-1,10**4])
plt.xlabel(r"$E_{kinetic}$ [GeV/nucleon]", size=16)
plt.ylabel(r"E $dN/dE$", size=16)
plt.legend(prop={'size':14})
plt.tight_layout()
plt.show()

