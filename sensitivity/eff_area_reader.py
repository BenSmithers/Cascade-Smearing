import numpy as np
import os
from cascade.utils import Data, config
from scipy.integrate import dblquad

from math import pi 

def get_deets():
    filename = os.path.join(config["datapath"], "charm_search_supplemental/effective_area.nu_e.txt")

#    filename = "/home/benito/Downloads/IC86SterileNeutrinoDataRelease/monte_carlo/NuFSGenMC_nominal.dat"
    data = np.loadtxt(filename).transpose()
    print(np.shape(data))


    energies = np.concatenate((np.unique(data[0]), [data[1][-1]]))
    print("{} primary energies".format(len(energies)))
    print(energies)
    costh = np.concatenate((np.unique(data[2]), [data[3][-1]]))
    print("{} primary angles".format(len(costh)))
    print(costh)


def quickload(filename):
    data = np.loadtxt(filename)
    


    n_reco_e = 28
    n_cth = 10

    assert(n_reco_e*n_cth == len(data))

    energies = [data[i][0] for i in range(n_reco_e)] + [data[-1][1]]
    cos_thetas = [data[i*n_reco_e][2] for i in range(n_cth)] + [data[-1][3]]


    ldata = np.zeros(shape=(n_reco_e, n_cth))
    err = np.zeros(shape=np.shape(ldata))
    
    for i_e_reco in range(n_reco_e):
        for i_cth in range(n_cth):
            true_i = i_e_reco + n_reco_e*i_cth 
            
            ldata[i_e_reco][i_cth] = data[true_i][4]
            err[i_e_reco][i_cth] = data[true_i][5]

    # NOT THESE KEYS
    return {"e_reco": energies,
            "cos_th": cos_thetas,
            "data":ldata, # THIS IS THE QUICKLOAD BULL
            "err": err #STOP DON'T USE THESE
            }

m2_to_cm2 = 100*100

def build_flux(*dataobjs):
    file_dir = os.path.join(config["datapath"], "charm_search_supplemental/")
    flavors = ["e", "mu", "tau"]
    # flavors += ["e_bar", "mu_bar", "tau_bar"]
    interactions = ["cc", "nc"]

    # energies, cos_th 
    casc_flux = np.zeros(shape=(28,10))
    casc_err = np.zeros(shape=(28,10))
    
    def metaflux(energy, angle, key):
        net_f = 0.0
        for dobj in dataobjs:
            net_f += 2*pi*dobj.get_flux((1e9)*energy, key, angle=np.cos(angle))
        return net_f
    

    for flavor in flavors:
        for interaction in interactions:
            for barness in ["", "_bar"]:

                filename = "effective_area.nu_{}.txt".format(flavor)
                area_data = quickload(os.path.join(file_dir, filename))
                angles = np.arccos(area_data["cos_th"])
                energies = area_data["e_reco"]
                # we need the key in the dataobj 

                curr = interaction.upper()
                flav = flavor[0].upper() + flavor[1:] # E/Mu/Tau 
                nu = "nu" if barness=="" else "nuBar"
                key = flav + "_" + nu + "_" + curr

                if flav=="Mu" and curr.lower()=="cc":
                    continue
                print("Integrating {} Flux".format(key))
                for i_e in range(len(casc_flux)):
                    for i_cth in range(len(casc_flux[i_e])):
                        if np.cos(angles[i_cth])>0:
                            continue
                        

                        flux_integrated = dblquad( lambda energy, angle: metaflux(energy, angle, key), a=angles[i_cth+1], b=angles[i_cth],gfun=energies[i_e], hfun=energies[i_e+1] )[0]
                        
                        
                        casc_flux[i_e][i_cth] += area_data["data"][i_e][i_cth]*m2_to_cm2*flux_integrated
                        casc_err[i_e][i_cth] += (area_data["err"][i_e][i_cth]*m2_to_cm2*flux_integrated)**2
 #   casc_err  = np.sqrt(casc_err)

    seconds_up = 0.98*10*3600*24*365
    casc_flux=seconds_up*casc_flux
    casc_err = np.sqrt(casc_flux)
    return {
            "e_edges":energies,
            "a_edges":np.cos(angles),
            "event_rate": casc_flux,
            "stat_err": casc_err
            }
