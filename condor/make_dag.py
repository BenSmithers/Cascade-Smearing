import numpy as np

theta13s = np.linspace(0, 2, 3)
theta23s = np.linspace(0, 2, 3)
msq3s = np.linspace(0,2.4,3)

filename = "genflux.dag"
obj = open(filename, "wt")

counter = 0
for theta13 in theta13s:
    for theat23 in theta23s:
        for msq3 in msq3s:
            obj.write("JOB  genflux{:06d}  genflux.submit\n".format(counter))
            counter+=1

counter = 0
for theta13 in theta13s:
    for theta23 in theta23s:
        for msq3 in msq3s:
            obj.write('VARS  genflux{:06d}  arg=" {} {} {}"\n'.format(counter, theta13, theta23, msq3))
            counter+=1

obj.close()
