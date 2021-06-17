from math import acos
import os

from cascade.utils import get_loc
from cascade.sensitivity.eff_area_reader import quickload

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



def build_mc_flux(*dataobjs):
    """
    Takes some Data objects and uses the MC we have to build up expectation arrays 
    """
    file_dir = "/home/benito/Downloads/IC86SterileNeutrinoDataRelease/monte_carlo/NuFSGenMC_nominal.dat"
    
    bars = ["nu", "nuBar"]

    def metaflux(energy, angle, key):
        net_f = 0.0
        for dobj in dataobjs:
            net_f += dobj.get_flux((1e9)*energy, key, angle=np.cos(angle))
        return net_f

    filename = "effective_area.nu_{}.txt".format(flavor)
    area_data = quickload(os.path.join(os.path.join(config["datapath"], "charm_search_supplemental/"), filename))

    e_edges = area_data["e_reco"]
    a_edges = np.arccos(area_data["cos_th"])
    net_flux = np.zeros(shape=(28,10))


    f = open(file_dir,'rt')
    while True:
        line = f.readline()
        if line=="":
            break #EOF 
        
        parsed = parseline(line)
        if parsed is None:
            continue # comment line

        flux_here = metaflux(parsed[3], parsed[4], "Mu_nu_CC")
        flux_here+= metaflux(parsed[3], parsed[4], "Mu_nuBar_CC")

        i_e = get_loc(parsed[1], e_edges)[0]
        i_a = get_loc(parsed[2], a_edges)[0]

        # let's let this event add tothe total in this bin!         
        net_flux[i_e][i_a] = net_flux[i_e][i_a] + parsed[5]*flux_here
        
        
