import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from math import pi, cos, sqrt

from cascade.oscillations.core import p_trans

energy = 1e3 # fixed 1TeV

positions = np.logspace(5,9, 300)
flavors = 4
values = [[p_trans(energy, pos, 1, flavor) for pos in positions] for flavor in range(flavors)]

for flavor in range(flavors):
    plt.plot(positions, values[flavor], label=r"$\nu$"+str(flavor))
plt.xlabel("Distance")
plt.ylabel("Probability")
plt.legend()
plt.xscale('log')
plt.show()
