import numpy as np
from math import sqrt, pi
from numpy import exp,cos, sin, arccos, arcsin

n_flavors = 4
u_matrix = np.ones(shape=(n_flavors,n_flavors))*1j
mixing_m = np.ones(shape=np.shape(u_matrix))

mixing_m[0][1] = mixing_m[1][0] = 0.563942
mixing_m[0][2] = mixing_m[2][0] = 0.154085
mixing_m[1][2] = mixing_m[2][1] = 0.785398

mixing_m[0][3] = mixing_m[3][0] = 0.
mixing_m[1][3] = mixing_m[3][1] = 9.*pi/180
mixing_m[2][3] = mixing_m[3][2] = 16.*pi/180

deltas = np.zeros(shape=np.shape(u_matrix))

msq = np.ones(shape=np.shape(u_matrix))

def mcos(i,j):
    return( cos(mixing_m[i-1][j-1]))
def msin(i,j):
    return( sin(mixing_m[i-1][j-1]))

u_matrix[0][0] = mcos(1,2)*mcos(1,3)*mcos(1,4)
u_matrix[0][1] = mcos(1,3)*mcos(1,4)*msin(1,2)
u_matrix[0][2] = mcos(1,4)*msin(1,3)*exp(-1*deltas[0][2]*1j)
u_matrix[0][3] =           msin(1,4)*exp(-1*deltas[0][3]*1j)

u_matrix[1][0] = mcos(1,2)*(-1*mcos(2,4)*msin(1,3)*msin(2,3)*exp(-1*deltas[0][2]*1j) -mcos(1,3)*msin(1,4)*msin(2,4)*exp(-1*(deltas[1][3]-deltas[0][3])*1j))-mcos(2,3)*mcos(2,4)*msin(1,2)

u_matrix[1][1] = mcos(1,2)*mcos(2,3)*mcos(2,4)+msin(1,2)*(-1*exp(deltas[0][2]*1j)*mcos(2,4)*msin(1,3)*msin(2,3) - exp(-1*(deltas[1][3]-deltas[0][3])*1j)*mcos(1,3)*msin(1,4)*msin(2,4))


u_matrix[1][2] = mcos(1,3)*mcos(2,4)*msin(2,3) - msin(1,3)*msin(1,4)*msin(2,4)*exp(-1*(deltas[0][2]- deltas[0][3] + deltas[1][3])*1j)
u_matrix[1][3] = mcos(1,4)*msin(2,4)*exp(-1*deltas[1][3]*1j)

u_matrix[2][0] = -msin(1,2)*(-mcos(3,4)*msin(2,3)-exp(deltas[1][3]*1j)*mcos(2,3)*msin(2,4)*msin(3,4)) + mcos(1,2)*(-exp(deltas[0][3]*1j)*mcos(1,3)*mcos(2,4)*msin(1,4)*msin(3,4) - msin(1,3)*exp(deltas[0][2]*1j)*(mcos(2,3)*mcos(3,4) - exp(deltas[1][3]*1j)*msin(2,3)*msin(2,4)*msin(3,4)))
u_matrix[2][1] = -mcos(1,2)*(-mcos(3,4)*msin(2,3) - exp(deltas[1][3]*1j)*mcos(2,3)*msin(2,4)*msin(3,4)) + msin(1,2)*(-exp(deltas[0][3]*1j)*mcos(1,3)*mcos(2,4)*msin(1,4)*msin(3,4) - msin(1,3)*exp(deltas[0][2]*1j)*(mcos(2,3)*mcos(3,4) - exp(deltas[1][3]*1j)*msin(2,3)*msin(2,4)*msin(3,4)))
u_matrix[2][2] = exp((deltas[0][2]-deltas[0][3])*1j)*mcos(2,4)*msin(1,3)*msin(1,4)*msin(3,4) + mcos(1,3)*(mcos(2,3)*mcos(3,4) - exp(deltas[1][3]*1j)*msin(2,3)*msin(2,4)*msin(3,4))
u_matrix[2][3] = mcos(1,4)*mcos(2,4)*msin(3,4)

u_matrix[3][0] = -msin(1,2)*(-exp(deltas[1][3]*1j)*mcos(2,3)*mcos(3,4)*msin(2,4) + msin(2,3)*msin(3,4)) + mcos(1,2)*(-exp(deltas[1][3]*1j)*mcos(1,3)*mcos(2,4)*mcos(3,4)*msin(1,4) - exp(deltas[0][2]*1j)*msin(1,3)*(-exp(deltas[1][3]*1j)*mcos(3,4)*msin(2,3)*msin(2,4)-mcos(2,3)*msin(3,4)))
u_matrix[3][1] =  mcos(1,2)*(-exp(deltas[1][3]*1j)*mcos(2,3)*mcos(3,4)*msin(2,4) + msin(2,3)*msin(3,4)) + msin(1,2)*(-exp(deltas[1][3]*1j)*mcos(1,3)*mcos(2,4)*mcos(3,4)*msin(1,4) - exp(deltas[0][2]*1j)*msin(1,3)*(-exp(deltas[1][3]*1j)*mcos(3,4)*msin(2,3)*msin(2,4)-mcos(2,3)*msin(3,4)))
u_matrix[3][2] = -exp((deltas[0][2]-deltas[0][3])*1j)*mcos(2,4)*mcos(3,4)*msin(1,3)*msin(1,4) + mcos(1,3)*(-exp(deltas[1][3]*1j)*mcos(3,4)*msin(2,3)*msin(2,4)- mcos(2,3)*msin(3,4))
u_matrix[3][3] = mcos(1,4)*mcos(2,4)*mcos(3,4)

#print(u_matrix)

# numpy conjugate 
conj = np.conjugate

# hbar is 1
hc = 0.1973269804*(1e-6) # eV * m
c = 3e8 # m * s 

def p_trans(energy, baseline, flav_i, flav_f):
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

    if not (flav_i<=n_flavors and flav_i>=0):
        raise ValueError("Invalid initial flavor: {}".format(flav_i))
    if not (flav_f<=n_flavors and flav_f>=0):
        raise ValueError("Invalid final flavor: {}".format(flav_f))

    if not (energy>0):
        raise ValueError("Energy must be positive: {}".format(energy))
    if not (baseline>0):
        raise ValueError("Baseline must be positive: {}".format(baseline))

    prob = 1. if flav_i==flav_f else 0.

    energy      *= 1e9 # now in eV 
    baseline    *= (1/hc) # in eV^-1 

    for i in range(n_flavors):
        for j in range(n_flavors):
            if i<=j:
                continue

            prob+=2*np.imag( conj(u_matrix[flav_f][j])*u_matrix[flav_f][i]*conj(u_matrix[flav_i][i])*u_matrix[flav_i][j]*sin(msq[i][j]*baseline/(2*energy)))


            prob-=4*np.real( conj(u_matrix[flav_f][j])*u_matrix[flav_f][i]*conj(u_matrix[flav_i][i])*u_matrix[flav_i][j]*sin(msq[i][j]*baseline/(4*energy))**2)

    return(prob)

