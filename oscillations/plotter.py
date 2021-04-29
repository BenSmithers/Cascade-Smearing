import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from math import pi, cos, sqrt

from cascade.oscillations.core import Oscillator
from cascade.utils import SterileParams

null_p = SterileParams()
ster_p = SterileParams(0., 0.1609, 0.2296, 4.47)
#ster_p = SterileParams(0.6597,0.160875,0.0, 4.47)
#ster_p = SterileParams(0.0,0.160875,0.0, 4.47)

null_osc = Oscillator(null_p)
ster_osc = Oscillator(ster_p)

# energy, baseline, flav_i, flav_f

d_earth = 2*(6.34e6) 
energy_range = np.logspace(2,6,100)
czeniths = np.linspace(-1, 0.01, 101)
thetas = np.arccos(czeniths)


#initial_flavor  = 1
#transmission    = [ [p_trans(en, d_earth, initial_flavor, flavor) for en in energy_range] for flavor in range(4) ] 

start_flav = 1
end_flav = 2

def calc_odds(which_osc):
    values = np.array([[ which_osc.p_trans(en, sqrt(2*d_earth*d_earth*(1-cos(pi-theta-theta))), start_flav, end_flav) for theta in thetas] for en in energy_range])
    
    return values

null_values = calc_odds(null_osc)
ster_values = calc_odds(ster_osc)

width = 0.1
plt.pcolormesh(czeniths, energy_range, ster_values-null_values,cmap=cm.coolwarm, vmin=-1*width, vmax=width)
cbar = plt.colorbar()
cbar.set_label("P(Sterile) - P(Null)")
plt.title(r"$\nu$"+str(start_flav)+r" to $\nu$"+str(end_flav))
plt.ylabel("Energy [GeV]", size=14)
plt.xlabel(r"cos$\theta$", size=14)
plt.yscale('log')
plt.show()
