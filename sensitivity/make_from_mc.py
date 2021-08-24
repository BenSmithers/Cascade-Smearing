from math import acos
import os
import numpy as np

from cascade.utils import config
from cascade.utils import get_loc
from cascade.sensitivity.eff_area_reader import quickload_full

def parseline(line):
    """
    Reads one of the lines used in the sterile MC data release 
    """
    while line[0]==" ":
        line = line[1:]
        if len(line)==0:
            return

    if line[0]=="#":
        return 

    return [float(x) for x in line.split()]



def build_mc_flux(*dataobjs, livetime=10):
    """
    Takes some Data objects and uses the MC we have to build up expectation arrays 

    the optional livetime arg scales the expectation to to that number of years of uptime 
    """
    cobalt = os.environ.get("_CONDOR_SCRATCH_DIR")
    filename = "NuFSGenMC_nominal.dat"
    if cobalt==None or cobalt=="" or cobalt==".":
        file_dir = "/home/benito/Downloads/IC86SterileNeutrinoDataRelease/monte_carlo/" + filename
    else:
        file_dir = os.path.join(cobalt, "data", filename)
    
    def metaflux(energy, angle, key):
        net_f = 0.0
        for dobj in dataobjs:
            net_f = net_f + dobj.get_flux((1e9)*energy, key, angle=angle)
        return net_f

    filename = "effective_area.per_bin.nu_mu.cc.track.txt"
    area_data = quickload_full(os.path.join(os.path.join(config["datapath"], "charm_search_supplemental","effective areas/"), filename))

    e_edges = area_data["e_reco"]
    a_edges = area_data["cth_true"]
    net_flux = np.zeros(shape=(len(e_edges)-1,len(a_edges)-1))


    f = open(file_dir,'rt')
    print("Parsing MC")
    while True:
        line = f.readline()
        if line=="":
            break #EOF 
        
        parsed = parseline(line)
        if parsed is None:
            continue # comment line

        if parsed[0]<0:
            key = "Mu_nuBar_CC"
        else:
            key = "Mu_nu_CC"

        flux_here = metaflux(parsed[3], parsed[4], key)

        if flux_here==0.0:
            continue

        if parsed[1]<e_edges[0] or parsed[1]>e_edges[-1]:
            continue
        if parsed[2]<a_edges[0] or parsed[2]>a_edges[-1]:
            continue
        i_e = get_loc(parsed[1], e_edges)[0]
        i_a = get_loc(parsed[2], a_edges)[0]

        # let's let this event add tothe total in this bin!         
        net_flux[i_e][i_a] = net_flux[i_e][i_a] + livetime*(365./343.7)*parsed[5]*flux_here
        
    return {
            "e_edges": e_edges,
            "a_edges": a_edges,
            "event_rate": net_flux,
            "stat_err": np.sqrt(net_flux)
            }
