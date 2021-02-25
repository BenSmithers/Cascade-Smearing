import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from math import pi, cos, sqrt

from cascade.oscillations.core import p_trans

# energy, baseline, flav_i, flav_f

d_earth = 2*(6.34e6) 
energy_range = np.logspace(3,8,100)
czeniths = np.linspace(-1, 0.01, 101)
thetas = np.arccos(czeniths)


#initial_flavor  = 1
#transmission    = [ [p_trans(en, d_earth, initial_flavor, flavor) for en in energy_range] for flavor in range(4) ] 

values = np.array([[ 1-p_trans(en, sqrt(2*d_earth*d_earth*(1-cos(pi-theta-theta))), 1, 1) for theta in thetas] for en in energy_range])

plt.pcolormesh(czeniths, energy_range, values,cmap=cm.inferno)
cbar = plt.colorbar()
cbar.set_label("test")
plt.title("Muon to not muon")
plt.ylabel("Energy [GeV]", size=14)
plt.xlabel(r"cos$\theta$", size=14)
plt.yscale('log')
plt.show()
