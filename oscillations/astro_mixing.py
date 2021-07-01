from cascade.oscillations.core import Oscillator
from cascade.utils import SterileParams

import numpy as np

def get_flavor_ratio(params):
    """

    """
    if not isinstance(params, SterileParams):
        raise TypeError("Expected {}, not {}".format(SterileParams, type(params)))
    osc = Oscillator(params)

    initial = (1,2,0,0)

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

    # normalize these all to the muons!
    norm_flavor = 1
    norm_amt= final_ratio[norm_flavor]
    for i in range(4):
        final_ratio[i] = final_ratio[i]/norm_amt

    return final_ratio
    

