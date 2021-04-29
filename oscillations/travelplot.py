import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm

import numpy as np
from math import pi, cos, sqrt

from cascade.oscillations.core import Oscillator
from cascade.utils import SterileParams, sci

ster_p = SterileParams(0.0,0.160875,0.0, 4.47)

osc = Oscillator(SterileParams())
#osc = Oscillator(ster_p)

energy = 5e1 # GeV



positions = np.logspace(5,8, 300)
flavors = 3
values = [[osc.p_trans(energy, pos, 1, flavor) for pos in positions] for flavor in range(flavors)]

names = ["e", r"$\mu$", r"$\tau$", "s"]
for flavor in range(flavors):
    plt.plot(positions, values[flavor], label=r"$\nu$"+names[flavor])

earth_r = 6.36e6

if min(positions)<(2*earth_r):
    plt.vlines(2*earth_r, ymin=0, ymax=1, colors="k", ls='--',label=r"$D_{earth}$")

plt.xlabel("Distance",size=14)
plt.ylabel("Probability",size=14)
plt.title("Muon Oscillations at {} GeV".format(energy),size=14)
plt.legend()
plt.xscale('log')
plt.show()
