"""
Use the scanner to actually do a likelihood scan! 
"""

from cascade.sensitivity.newSense import doLLH, Scanner, JointLLH
import numpy as np
import os 
import pickle
from math import asin, sqrt

from cascade.utils import get_loc, config, gen_filename, SterileParams
import sys
print("Running with : {}".format(sys.argv))
th14_mode = True
systematics = True
if len(sys.argv)<2 and th14_mode:
    raise ValueError()

#best_th14 = 0.3555
best_th14 = 0.0
best_th24 = 0.1609
#best_th24 = 0.3826
central = SterileParams(theta03=best_th14, theta13=best_th24, msq2=float(sys.argv[1]))

#central = SterileParams(theta13=0.1652, theta23=0.2293, msq2=4.6416)

#th03 = asin(sqrt(0.42))/2.0
#central = SterileParams(theta03 = th03, theta13=0.1609, msq2=4.47)
# central = SterileParams(theta13=0.1609, msq2=4.47)

if False:
    # load fudge
    fudge_file = os.path.join(os.path.dirname("__file__"), "fudge.dat")
    raw_fudge = np.transpose(np.loadtxt(fudge_file, delimiter=","))
    energies = raw_fudge[0]
    fudge = raw_fudge[1]

    bin_angles = np.linspace(-1,1,11)
    bin_energies = np.logspace(2,8,21)

    full_fudge = np.zeros(shape=(20,10))

    for e_i in range(20):
        if bin_energies[e_i] < energies[0]:
            j = 0
        elif bin_energies[e_i]>energies[-1]:
            j = -1
        else:
            j = get_loc(0.5*(bin_energies[e_i]+bin_energies[e_i+1]), energies)[0]
        for a_i in range(10):
            full_fudge[e_i][a_i] = fudge[j]

#true_fudge = np.load("full_fudge.npy")


smearing = True

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

x = np.concatenate(([0], np.arcsin(np.sqrt(np.logspace(-3,0,90)))/2))
if th14_mode:

    y = [float(sys.argv[1])]
else:
    y = np.concatenate(([0], np.logspace(-2,2,40)))
test = Scanner(jointllh, x, x, y, th14_mode) 

results = test.scan()
results["central"] = central

new_ev_string = "_".join("{:.2f}".format(y[0]).split("."))
smearstring = "_smearing" if smearing else ""
name = "best_llh_{}eV{}.dat".format(new_ev_string,smearstring) if th14_mode else "joint_likelihood{}.dat".format(smearstring)

write_dir = gen_filename(config["datapath"], name , central)
print("Wrote to {}".format(write_dir))
f = open(write_dir, 'wb')
pickle.dump(results, f, -1)
f.close()
