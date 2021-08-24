import numpy as np
from math import pi

is_mc = True

n_grid = 1
theta03 =0

# sin^2(2theta) over 1e-3, 1e0

#theta13s = np.linspace(0, 0.5*pi, n_grid)
#theta23s = np.linspace(0, 0.5*pi, n_grid)

theta13s = np.arcsin(np.sqrt(np.logspace(-3, 0, n_grid)))/2.0
theta23s = np.arcsin(np.sqrt(np.logspace(-3, 0, n_grid)))/2.0

# use the zero point! 

theta13s = np.concatenate((np.array([0]) , theta13s))
theta23s = np.concatenate((np.array([0]) , theta23s))

filename = "susflux{}.dag".format("_from_mc" if is_mc else "")
obj = open(filename, "wt")

counter = 0
for theat23 in theta23s:
    for theat13 in theta13s:
        obj.write("JOB  susflux{:06d}  genflux{}.submit\n".format(counter, "_frommc" if is_mc else ""))
        counter+=1

counter = 0
for theta23 in theta23s:
    for theta13 in theta13s:
        obj.write('VARS  susflux{:06d}  arg=" {} {} {}"\n'.format(counter, theta13, theta23, 1 if is_mc else 0))
        counter+=1

obj.close()
