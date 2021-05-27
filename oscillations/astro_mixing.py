from cascade.oscillations.core import Oscillator
from cascade.utils import SterileParams

import numpy as np

null = SterileParams()
#ster = SterileParams(0.0, 0.1609, 0.2205, 4.7)
ster = SterileParams(0.0, 0.1609, 0.0, 4.7)
osc = Oscillator(null)

initial = (1,2,0,0)


pivot = 100*(1e3) #GeV
norm = 0.787e-18  # [GeV sr s cm2]^-1
gamma = -2.5

def spectrum(E):
    """
    We assume this is the distribution of all the neutrino flavorinos 
    """
    return norm*(E/pivot)**gamma


# the energy reeeaaally doesn't matter... we're talking about stupendously large distances
# so consider that

arbitrarily_large_ratio = 1.0e25

final_ratio = np.array([0., 0., 0., 0.])

for i in range(4):
    for j in range(4):
        if initial[i]==0:
            continue

        for ratio in arbitrarily_large_ratio*np.logspace(-3, 3, 600):
            osc_prob = osc.p_trans(1.0, ratio, i, j)*initial[i]
            final_ratio[j] = final_ratio[j]+ osc_prob 
        print("{} to {} Mixing: {}".format(i, j, osc_prob))

final_ratio = final_ratio/np.sum(final_ratio)
print(final_ratio)

