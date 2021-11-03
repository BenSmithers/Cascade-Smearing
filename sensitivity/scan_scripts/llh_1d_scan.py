"""
SUPER PLOT

Draw our sensitivity contours and then overlay other experiments' sensitivities 
"""

from cascade.sensitivity.newSense import Scanner, JointLLH, doLLH
from cascade.utils import SterileParams, gen_filename, config
from cascade.utils import get_loc, get_closest, get_color
from cascade.utils import shift_cmap

import pickle
from math import sqrt, asin
import numpy as np
import os
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

plt.style.use(os.path.join(os.path.dirname(__file__), "..", ".." , "cascade.mplstyle"))


"""
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_1_00eV_smearing_0.0_0.1609e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_3_30eV_smearing_0.0_0.1609e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.1609e0/best_llh_4_64eV_smearing_0.0_0.1609e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.1609$", '-')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_1_00eV_smearing_0.0_0.3826e0_0.0_1.0000e0.dat", r"1.0eV$^{2}$, $\theta_{24}=0.3826$", '--')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_3_30eV_smearing_0.0_0.3826e0_0.0_3.3000e0.dat", r"3.3eV$^{2}$, $\theta_{24}=0.3826$", '--')
add_contour("/home/bsmithers/software/data/hg_sib/0.3826e0/best_llh_4_64eV_smearing_0.0_0.3826e0_0.0_4.6400e0.dat", r"4.64eV$^{2}$, $\theta_{24}=0.3826$", '--')
"""


msqs = [1.0, 3.3, 4.64]
th24s = [0.1609, 3.826]
th34s = [ 0.0 ]
thetas = np.concatenate(([0], np.arcsin(np.sqrt(np.logspace(-3,0,90)))/2))


systematics = True
smearing = True
central = SterileParams()

options ={
    "is_mc": False,
    "use_syst": systematics,
    "skip_missing":True,
    "smear":smearing
}
mc_options = {
    "is_mc" : True,
    "use_syst":systematics,
    "skip_missing":True,
    "smear":False
}

llhood = doLLH("best_expected_flux.dat",central_exp=central, options=options)
mcllhood = doLLH("expected_flux_from_mc_smearedwell.dat", central_exp=central, options=mc_options)
jointllh = JointLLH(mcllhood, llhood)

thetas = np.concatenate(([0], np.arcsin(np.sqrt(np.logspace(-3,0,90)))/2))

llh_dict = Scanner(jointllh, theta24s=th24s, theta34s=th34s, msqs=msqs, th14_mode=True, theta14s=thetas) 
results = llh_dict.scan()


outname = "llh_1d_scan.dat"
write_dir = gen_filename(config["datapath"], outname , central)
print("Wrote to {}".format(write_dir))
f = open(write_dir, 'wb')
pickle.dump(results, f, -1)
f.close()