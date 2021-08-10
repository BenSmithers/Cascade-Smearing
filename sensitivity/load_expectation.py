
import numpy as np
from math import pi,cos

import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

from cascade.utils import get_loc
def load(debug=False, e_edges = np.logspace(2,9,29)):

    filename = "/home/benito/software/data/cascade/hg_sib/charm_search_supplemental/event_properties.txt"
    data = np.loadtxt(filename, dtype=str)
    """
    Entries
        0 - event no
        1 - topology
        2 - photoelectrons
        3 - energy 
        4 - energy error, atmospheric
        5 - energy error, astro spectrum -2.47
        6 - energy error, astro spectrum -2
        7 - declination 
        8,9,10 - error on declination (same as before)
    """

    n_below = 0
    n_above = 0
    n_tracks = 0
    n_cascades = 0

    e_edges = np.logspace(2,7,29)
    #e_edges = np.logspace(2,8,21)

    occupation = np.zeros(shape=(len(e_edges)-1, 2))
    # one bin in cth 
    for event in data:
        if event[1].lower()=="track":
            n_tracks +=1
            continue
        else:
            n_cascades += 1
        # we have the direction the particles came from stored. 
        # Not the direction they are going! 
        cth = cos((pi/2) + float(event[7])*pi/180 )
        energy = float(event[3])*1e3 # TeV->GeV 
        if cth<0.2:
            n_below+=1
            i_bin = get_loc(energy, e_edges)[0]
            
            occupation[i_bin][0]+=1
        else:
            n_above+=1
#            occupation[i_bin][1]+=1
    print("Found {} tracks, {} cascades; {} below".format(n_tracks, n_cascades, n_below))
    if debug:
        return occupation
    else:
        return np.array(occupation)

if False:
    occupation = load(True)
    e_edges = np.logspace(2,8,21)
    print(e_edges)
    plt.bar(x=e_edges[:-1], height=occupation, width=(e_edges[1:]-e_edges[:-1]), align='edge')
    plt.xlabel("Energy [GeV]",size=14)
    plt.ylabel("N events")
    plt.xscale('log')
    plt.show()


