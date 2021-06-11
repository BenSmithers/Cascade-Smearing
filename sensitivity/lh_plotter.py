from cascade.utils import gen_filename, config

import pickle

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

levels = [90.0, 68.0]

f = open(gen_filename(config["datapath"]+"/expected_fluxes_reco/","cummulative_probs.dat"),'rb')
obj = pickle.load(f)
f.close()

theta24s = obj["theta24s"]
theta34s = obj["theta34s"]
msqs = obj["msqs"]
raw_sorted = obj["raw_sorted"]
c_prob = obj["c_prob"]

# sum over the mass-squareds 
c_prob_msq = np.zeros(shape=(len(theta24s), len(theta34s)))
for entry in sorted_likelihoods:
    c_prob_msq[entry[0]][entry[1]] += entry[3]

plt.contour( theta24s, theta34s, 100*c_prob_msq.transpose(), levels=levels)
plt.xlabel(r"$\theta_{24}$",size=14)
plt.ylabel(r"$\theta_{34}$",size=14)
plt.show()


c_prob_th34 = np.zeros(shape=(len(theta24s), len(msqs)))
for entry in sorted_likelihoods:
    c_prob_th34[entry[0]][entry[2]] += entry[3]

plt.contour( msqs, theta24s, 100*c_prob_th34, levels=levels)
plt.ylabel(r"$\theta_{24}$",size=14)
plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
plt.show()

c_prob_th24 = np.zeros(shape=(len(theta34s), len(msqs)))
for entry in sorted_likelihoods:
    c_prob_th24[entry[1]][entry[2]] += entry[3]

plt.contour( msqs, theta34s, 100*c_prob_th24, levels=levels)
plt.ylabel(r"$\theta_{34}$",size=14)
plt.xlabel(r"$\Delta m_{14}^{2}$",size=14)
plt.show()


