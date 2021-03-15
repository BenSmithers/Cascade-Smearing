"""
Builds the PMNS matrix for the 3+1 neutrino model and provides a survival probability function 
"""

import numpy as np
from math import sqrt, pi
from numpy import exp,cos, sin, arccos, arcsin

from cascade.utils import SterileParams, multimat

conj = np.conjugate

class Oscillator:
    """
    You provide this object a set of sterile neutrino parameters. It constructs a unitary 4-flavor mixing matrix according to those parameters, and you can ask it for oscillation probabilities 

    It just does vacuum oscillation probabilities right now 
    """
    def __init__(self, params):
        if not isinstance(params, SterileParams):
            raise TypeError("Expected {}, got {}".format(SterileParams))

        self.n_flavors = 4
        self._u_matrix = np.ones(shape=(self.n_flavors,self.n_flavors))*1j
        mixing_m = np.ones(shape=np.shape(self._u_matrix))

        mixing_m[0][1] = mixing_m[1][0] = 0.563942 
        mixing_m[0][2] = mixing_m[2][0] = 0.149575
        mixing_m[1][2] = mixing_m[2][1] = 0.855211

        mixing_m[0][3] = mixing_m[3][0] = params.theta03
        mixing_m[1][3] = mixing_m[3][1] = params.theta13
        mixing_m[2][3] = mixing_m[3][2] = params.theta23

        deltas = np.zeros(shape=np.shape(self._u_matrix))
        self.msq = np.zeros(shape=np.shape(self._u_matrix))

        self.msq[0][1] = 7.59e-5
        self.msq[1][0] = -self.msq[0][1]
        self.msq[1][2] = 2.32e-3
        self.msq[2][1] = -self.msq[1][2]
        self.msq[0][2] = 2.32e-3 # very approximate 
        self.msq[2][0] = -self.msq[0][2]

        self.msq[1][3] = params.msq2
        self.msq[3][1] = -self.msq[1][3]
        self.msq[0][3] = self.msq[0][1]+self.msq[1][3] 
        self.msq[2][3] = self.msq[1][3]-self.msq[1][2]
        self.msq[3][0] = -self.msq[0][3]
        self.msq[3][2] = -self.msq[2][3]

        # some macros! 
        def mcos(i,j):
            return( cos(mixing_m[i-1][j-1]))
        def msin(i,j):
            return( sin(mixing_m[i-1][j-1]))
        
        # we have these macros for creating 4x4 rotation matrices according to the mixing angles specified 
        def r_mat(i,j):
            matrix = np.identity(self.n_flavors)*1j
            matrix[i-1][i-1] = mcos(i,j)
            matrix[j-1][i-1] = msin(i,j)
            matrix[i-1][j-1] = -msin(i,j)
            matrix[j-1][j-1] = mcos(i,j)
            return(matrix)

        def r_hat(i,j):
            matrix = np.identity(self.n_flavors)*1j
            matrix[i-1][i-1] = mcos(i,j)
            matrix[j-1][i-1] = msin(i,j)*exp(-deltas[i-1][j-1]*1j)
            matrix[i-1][j-1] = -conj(msin(i,j)*exp(-deltas[i-1][j-1]*1j))
            matrix[j-1][j-1] = mcos(i,j)
            return(matrix)

        # we calculate the unitary rotation matrix via multiple sub-rotations!
        # guarantees that the matrix we build is unitary, regardless of what the specified mixing angles are 
        self._u_matrix = multimat( r_hat(2,4), r_mat(3,4), r_hat(1,4), r_mat(2,3), r_hat(1,3), r_mat(1,2))


        # hbar is 1
        self._hc = 0.1973269804*(1e-6) # eV * m
        self._c = 3e8 # m * s 

    @property
    def u_matrix(self):
        return self._u_matrix

    def p_trans(self, energy, baseline, flav_i, flav_f):
        """
        energy  - GeV
        baseline- meters
        flav_i  - initial flavor (0,3)
        flav_f  - final flavor (0,3)
        """

        if not isinstance(energy, (int, float)):
            raise TypeError("Expected {} for energy, not {}".format(float, type(energy)))
        if not isinstance(baseline, (int, float)):
            raise TypeError("Expected {} for baseline, not {}".format(float, type(baseline)))

        if not isinstance(flav_i, int):
            raise TypeError("Expected {} for initial flavor, not {}".format(int, type(flav_i)))
        if not isinstance(flav_f, int):
            raise TypeError("Expected {} for final flavor, not {}".format(int, type(flav_f)))

        if not (flav_i<=self.n_flavors and flav_i>=0):
            raise ValueError("Invalid initial flavor: {}".format(flav_i))
        if not (flav_f<=self.n_flavors and flav_f>=0):
            raise ValueError("Invalid final flavor: {}".format(flav_f))

        if not (energy>0):
            raise ValueError("Energy must be positive: {}".format(energy))
        if not (baseline>0):
            raise ValueError("Baseline must be positive: {}".format(baseline))

        prob = 1. if flav_i==flav_f else 0.

        energy      *= 1e9 # now in eV 
        baseline    *= (1/self._hc) # in eV^-1 

        for j in range(self.n_flavors):
            for i in range(self.n_flavors): 
                if i<j:
                    continue

                prob += 2*np.imag( conj(self.u_matrix[flav_f][j])*self.u_matrix[flav_f][i]*conj(self.u_matrix[flav_i][i])*self.u_matrix[flav_i][j]*sin(self.msq[i][j]*baseline/(2*energy)))

                prob -= 4*np.real( conj(self.u_matrix[flav_f][j])*self.u_matrix[flav_f][i]*conj(self.u_matrix[flav_i][i])*self.u_matrix[flav_i][j]*pow(sin(self.msq[i][j]*baseline/(4*energy)),2))

        return(prob)

