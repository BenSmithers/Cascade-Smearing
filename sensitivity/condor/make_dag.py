import numpy as np
from math import pi

n_grid = 20
theta03 =0

theta13s = np.linspace(0, 0.5*pi, 1)
theta23s = np.linspace(0, 0.5*pi, 1)

filename = "susflux.dag"
obj = open(filename, "wt")

counter = 0
for theat23 in theta23s:
    for theat13 in theta13s:
        obj.write("JOB  susflux{:06d}  genflux.submit\n".format(counter))
        counter+=1

counter = 0
for theta23 in theta23s:
    for theta13 in theta13s:
        obj.write('VARS  susflux{:06d}  arg=" {} {} {} {}"\n'.format(counter, theta13, theta23))
        counter+=1

obj.close()
