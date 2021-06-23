import numpy as np
from math import pi

n_grid = 90

theta13s = np.linspace(0, 0.5*pi, n_grid)
theta23 = 0.0

filename = "genflux_frommc.dag"
obj = open(filename, "wt")

counter = 0
for theta13 in theta13s:
    obj.write("JOB  genflux_frommc{:06d}  genflux_frommc.submit\n".format(counter))
    counter+=1

counter = 0
for theta13 in theta13s:
    obj.write('VARS  genflux_frommc{:06d}  arg=" {} {}"\n'.format(counter, theta13, theta23))
    counter+=1

obj.close()
