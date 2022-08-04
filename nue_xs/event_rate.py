"""
Here, we get an expected flux of neutrinos at the detector, and then use the cross sections to get an idea of what event rate we might expect 
"""

from cmath import pi
from cascade.raw_fluxes import raw_flux
from cascade.utils import SterileParams, get_loc
from cascade.nue_xs.cross_section import get_diff_xs, con
from math import acos, cos
import numpy as np

from cascade.cross_section_test import flavors, neut_types

from scipy.interpolate import interp2d
from scipy.integrate import quad

import matplotlib
matplotlib.use("GTK3Agg")
from matplotlib import pyplot as plt

import nuSQuIDS as nsq

n_bin = 100

_enus = np.logspace(2, 7, n_bin)*1e9 # eV
_zeniths = np.linspace(-1,1,int(n_bin/2)+1)

enus = 0.5*(_enus[1:] + _enus[:-1])
zeniths = 0.5*(_zeniths[1:] + _zeniths[:-1])

ewidths = _enus[1:] - _enus[:-1]
zwidths = np.abs(np.arccos(_zeniths[1:])-np.arccos( _zeniths[:-1]))
print(zwidths)
stangle = np.sin(np.arccos(zeniths))*zwidths*2*pi
print(stangle)

bjork_y = np.logspace(-6,0, n_bin, endpoint = False) # unitless 

print("... building interpolators")
all_xs_nodes = {}
# make a function (a scipy interpolator) that can get the differential xs at each point
for flavor in flavors.keys():
    all_xs_nodes[flavor] = {}
    for neut_type in neut_types.keys():
        nus_flav = flavors[flavor]
        nus_neut = neut_types[neut_type]

        all_xs_nodes[flavor][neut_type] = interp2d(enus, bjork_y, get_diff_xs(enus, bjork_y, nus_flav, nus_neut))

# now, we want to get the total cross sections... so we will integrate over Bjorken Y

print("... calculating total cross sections")
total_xs = {}
for flavor in flavors.keys():
    total_xs[flavor] = {}
    for neut_type in neut_types.keys():
        total_xs[flavor][neut_type] = np.array([ quad(lambda bjorky: all_xs_nodes[flavor][neut_type](energy, bjorky), a=1e-8, b=1.0-1e-8)[0] for energy in enus])

rate = np.zeros((len(zeniths), len(enus)))

nus_atm = raw_flux(SterileParams())

ice_rho = 917. # kg/m^3 
volume = (1000.)**3
m_ice= volume*ice_rho #kg
ice_molecule_mass = (15.999 + 2*1.008)*(1.66e-27) # kg/mol
n_molecules = m_ice/ice_molecule_mass # molecules 
e_density = float((2+8)/(2 + 16)) # electrons / nucleon 

n_targets = e_density*n_molecules

livetime = 10*365*24*3600
scale_factor = livetime*n_targets
print("{} target electrons".format(n_targets))

for i_z in range(len(zeniths)):
    for i_e in range(len(enus)):
        binwid = (1e-9)*ewidths[i_e]*stangle[i_z]

        for flavor in flavors.keys():
            for neut_type in neut_types.keys():
                nus_flav = flavors[flavor]
                nus_neut = neut_types[neut_type]

                rate[i_z][i_e] += total_xs[flavor][neut_type][i_e] * nus_atm.EvalFlavor(nus_flav, zeniths[i_z], enus[i_e], nus_neut) * scale_factor *binwid

#nus_atm.EvalFlavor(get_flavor(key), angle, energy, get_neut(key))

print("total rate : {}".format(np.sum(rate)))
plt.pcolormesh(_zeniths, _enus/(1e9), np.transpose(rate))
plt.yscale('log')
plt.ylabel("Energy [GeV]",size=14)
plt.xlabel(r"$\cos\theta$",size=14)
plt.colorbar()
plt.show()

