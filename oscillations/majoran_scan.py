from cascade.utils import SterileParams, get_color
from cascade.oscillations.core import Oscillator

import matplotlib
matplotlib.use('GTK3Agg')
from matplotlib import pyplot as plt 

import numpy as np
import os
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon
from scipy.optimize import minimize
from math import pi
from tqdm import tqdm



n_pt = 60
sinsq = np.linspace(0,1,n_pt)
radians =  np.arcsin(np.sqrt(sinsq))/2.0
dm2 = np.linspace(0,10, n_pt+1)

best_th14 = 0.0
best_dm2 = 1.3

M_LIGHTEST = 0.0

MAX_M_EFF = 0.155 # eV

"""
Scan over sterile neutrino space, get likelihood from best data
At each point, then get play with the phasese to minimize the effective mass. If it's beyond the 
"""



SKIP_SCAN = True
SAVE_NAME = "minimized_mass"
def get_minmas()->np.ndarray:
    fn = os.path.join(os.path.dirname(__file__), SAVE_NAME+"_"+"_".join([str(n_pt), str(n_pt+1) ])+".npy")

    if os.path.exists(fn) and SKIP_SCAN:
        return np.load(fn)
    else:
        minmass = np.zeros(shape=(len(sinsq), len(dm2)))
 
        for x,y in tqdm([(i,j) for j in range(len(dm2)) for i in range(len(sinsq))]): 
            def get_eff_mass(p):
                phase1, phase2, phase3 = p
                sp = SterileParams(theta03=radians[x], msq2=dm2[y],
                            maj_phase_1=phase1,
                            maj_phase_2=phase2,
                            maj_phase_3=phase3)
                osc = Oscillator(sp)
                mass = osc.get_eff_mass(M_LIGHTEST, 0)
                del osc
                del sp
                return mass
            
            bounds = ((0,2*pi),
                        (0,2*pi),
                        (0,2*pi))

            n_tries = 10
            true_min = 100000.
            for i in range(n_tries):
                p0 = tuple(np.random.rand(3)*2*pi) 
                minmin = minimize(get_eff_mass, p0,method='L-BFGS-B',  bounds=bounds, options={'eps':5e-3}).x

                mass = get_eff_mass(minmin)
                if mass<true_min:
                    true_min = mass
            minmass[x][y] = true_min
        np.save(fn, minmass)
        return minmass
inverted = True


_data = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "best_3s.dat"), delimiter=",")
three_sig = Polygon([Point(entry) for entry in _data])

two_s1 = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "best_2s-1.dat"), delimiter=",")
two_s2 = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "best_2s-2.dat"), delimiter=",")
two_sig = Polygon([Point(entry) for entry in two_s1]).union( Polygon([Point(entry) for entry in two_s2]))

one_s1 = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "best_1s-1.dat"), delimiter=",")
one_s2 = np.loadtxt(os.path.join(os.path.dirname(__file__), "data", "best_1s-2.dat"), delimiter=",")
one_sig = Polygon([Point(entry) for entry in one_s1]).union( Polygon([Point(entry) for entry in one_s2]))

lhs = np.zeros(shape=(len(sinsq), len(dm2)))
minmass = get_minmas()


for x,y in [(i,j) for j in range(len(dm2)) for i in range(len(sinsq))]:
    this_point = Point(sinsq[x], dm2[y])

    inside_1 = one_sig.contains(this_point)

    if inside_1:
        lhs[x][y] = 1.-0.682689492
        color = (0.0, 1.0, 0.0)
    else:
        inside_2 = two_sig.contains(this_point)
        if inside_2:
            lhs[x][y] = 1.- 0.954499736 # - 0.682689492
            color = (1.0, 1.0, 0.0)
        else:
            inside_3 = three_sig.contains(this_point)
            if inside_3:
                lhs[x][y] = 1.- 0.997300204 #- 0.954499736
                color = (1.0, 0.5, 0.0)
            else:
                color = (1.0,0.9,0.9)
                lhs[x][y] = 1.0-0.99993666


kamkam = np.zeros(shape=np.shape(minmass))
mask = minmass<MAX_M_EFF
kamkam[mask] = 0.90
kamkam[np.logical_not(mask)]=.10
plt.pcolormesh(sinsq, dm2, np.transpose(kamkam), vmin=0, vmax=1., cmap='magma')
plt.xlabel(r"$\sin^{2}2\theta_{14}$", size=14)
plt.ylabel(r"$\Delta m_{41}^{2}$", size=14)

for x,y in [(i,j) for j in range(len(dm2)) for i in range(len(sinsq))]:   
    break
    if minmass[x][y]>2*MAX_M_EFF:
        lhs[x][y] = lhs[x][y]*0.10
    elif minmass[x][y]>MAX_M_EFF:
        lhs[x][y] = lhs[x][y]*0.10
    else:
        lhs[x][y] = lhs[x][y]*0.90

llhs = np.log10(lhs)
chis =  -2*(llhs -np.max(llhs))



levels = [0.317, 0.045, 0.001][::-1]
#levels=[0.682689492, 
#        0.954499736 - 0.682689492 ,
#        0.997300204 - 0.954499736 ][::-1]

labels = [ "1s", "2s","3s"][::-1]
colors = [get_color(i+1,colormax=4) for i in range(3)]

def set_lbls(ct_plot, labels, color):
    global move_it
    fmt = {}
    for l,s in zip(ct_plot.levels, labels):
        fmt[l] = s
    loc = ((1e-2 + (move_it*0.8e-2), 1e-1 - move_it*(1.7e-2)), )
    thingy = plt.clabel(ct_plot, ct_plot.levels, inline=True, fmt=fmt, fontsize=10,manual=loc, colors=(color,))
    move_it += 1


ct = plt.contourf(sinsq, dm2, np.transpose(lhs), levels=levels, colors=colors, alpha=.5)
fmt={}
for l, s in zip(ct.levels, labels):
    fmt[l] = s

#plt.clabel(ct, ct.levels, fmt=fmt, inline=True)
plt.xlabel(r"$\sin^{2}2\theta_{14}$", size=14)
plt.ylabel(r"$\Delta m_{41}^{2}$", size=14)
plt.savefig("overlay.png",dpi=400)
plt.show()

