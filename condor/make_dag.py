import numpy as np
from math import pi

n_grid = 20
theta03s = np.linspace(0, pi, n_grid) #el
theta13c = 0.160875
span = 0.05
theta13s = [0.0,0.05, 0.1, 0.160875, 0.2]
#theta13s = np.linspace(theta13c-span , theta13c+span, 10)
theta23s = np.linspace(0, 40*pi/180., n_grid) #tau
msq = 4.47

filename = "genflux.dag"
obj = open(filename, "wt")

counter = 0
for theta03 in theta03s:
    for theat23 in theta23s:
        for theat13 in theta13s:
            obj.write("JOB  genflux{:06d}  genflux.submit\n".format(counter))
            counter+=1

counter = 0
for theta03 in theta03s:
    for theta23 in theta23s:
        for theta13 in theta13s:
            obj.write('VARS  genflux{:06d}  arg=" {} {} {} {}"\n'.format(counter,theta03, theta13, theta23, msq))
            counter+=1

obj.close()
