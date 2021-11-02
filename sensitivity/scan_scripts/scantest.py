"""
Use the scanner to actually do a likelihood scan! 
"""

from cascade.sensitivity.newSense import doLLH, Scanner  
import numpy as np
import os 
import pickle

from cascade.utils import get_loc, config, gen_filename, SterileParams


central = SterileParams()
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

#true_fudge = np.load("../full_fudge.npy")

# fix_norm
options ={
    "is_mc":False,
    "skip_missing":True,
    "use_syst":True,
    "smear":True

}


#llhood = doLLH("expected_flux_smearedwell.dat",central_exp=central, options=options)
llhood = doLLH("best_expected_flux.dat",central_exp=central, options=options)
x = np.concatenate(([0], np.arcsin(np.sqrt(np.logspace(-3,0,90)))/2))
y = np.concatenate(([0], np.logspace(-2,2,40))) 
test = Scanner(llhood, x, x, y) 

results = test.scan()

write_dir = gen_filename(config["datapath"], "newSense_result_smear.dat", central)
print("Wrote to {}".format(write_dir))
f = open(write_dir, 'wb')
pickle.dump(results, f, -1)
f.close()
