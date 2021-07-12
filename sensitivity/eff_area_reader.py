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

def other_deets():
    filename = os.path.join(config["datapath"], "charm_search_supplemental/effective areas/effective_area.per_bin.nu_e_bar.cc.cascade.txt")

#    filename = "/home/benito/Downloads/IC86SterileNeutrinoDataRelease/monte_carlo/NuFSGenMC_nominal.dat"
    data = np.loadtxt(filename).transpose()
    print(np.shape(data))

    energies = np.concatenate((np.unique(data[0]), [data[1][-1]]))
    print("{} primary energies".format(len(energies)))
    print(energies)
    costh = np.concatenate((np.unique(data[2]), [data[3][-1]]))
    print("{} primary angles".format(len(costh)))
    print(costh)
    
    costh = np.concatenate((np.unique(data[4]), [data[5][-1]]))
    print("{} dep energy".format(len(costh)))
    print(costh)

    costh = np.concatenate((np.unique(data[6]), [data[7][-1]]))
    print("{} reco zenith angles".format(len(costh)))
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
            "err": err #STOP DON'T USE THESE KEYS
            }

def quickload_full(filename):
    data = np.loadtxt(filename)
    
    n_true_e = 100
    n_cth_true = 10
    n_dep_e = 20
    n_cth_r = 2

    assert(n_true_e*n_cth_true*n_dep_e*n_cth_r==len(data))

    e_true      = [data[i][0] for i in range(n_true_e)] + [data[-1][1]]
    cth_true    = [data[i*n_true_e][2] for i in range(n_cth_true)]+data[-1][3]
    e_depo      = [data[i*n_true_e*n_cth_true][4] for i in range(n_dep_e)] +data[-1][5]
    cth_reco    = [data[i*n_true_e*n_cth_true*n_dep_e][6] for i in range(n_cth_r)] + data[-1][7]

    ldata = np.zeros(shape=(n_true_e, n_dep_e,n_cth_true, n_cth_r))

    for i_e_t in range(n_true_e):
        for i_e_d in range(n_dep_e):
            for i_cth_t in range(n_cth_true):
                for i_cth_r in range(n_cth_r):
                    true_i = i_e_t + n_true_e*(i_cth_t + n_cth_true*(i_e_d + n_dep_e*i_cth_r))
                    
                    ldata[i_e_t][i_e_d][i_cth_t][i_cth_r] = data[true_i][8]
    return  { 
            "e_reco":e_depo,
            "e_true":e_true,
            "cth_true":cth_true,
            "cth_reco":cth_reco,
            "data":data
            }

m2_to_cm2 = 100*100

def build_flux_sad(*dataobjs):
    file_dir = os.path.join(config["datapath"], "charm_search_supplemental", "effective areas")
    flavors = ["e", "mu", "tau"]
    interactions=["cc", "nc"]
    
    casc_flux = np.zeros(shape=(20, 2)) 
    print("Integrating Flux over bins") 

    def metaflux(energy, angle, key):
        net_f = 0.0
        for dobj in dataobjs:
            net_f += 2*pi*dobj.get_flux((1e9)*energy, key, angle=np.cos(angle))
        return net_f
    
    for flavor in flavors:
        for interaction in interactions:
            for barness in ["","_bar"]:
                comb_flav = flavor+barness
                filename = "effective_area.per_bin.nu_{}.{}.track.txt".format(comb_flav, interaction)
                area_data= quickload_full(os.path.join(file_dir, filename))
                
                flux=area_data["data"]
                e_reco = area_data["e_reco"]
                e_true = area_data["e_true"]
                cth_true = area_data["cth_true"]
                cth_reco = area_data["cth_reco"]

                curr = interaction.upper()
                flav=flavor[0].upper() + flavor[1:].lower()
                nu="nu" if barness=="" else "nuBar"
                key = "_".join([flav, nu, curr])

                if flav=="Mu" and curr.lower()=="cc":
                    continue
                for i_e_d in range(len(e_reco)-1):
                    for i_e_t in range(len(e_true)-1):
                        for i_cth_t in range(len(cth_true)-1):
                            for i_cth_r in range(len(cth_reco)-1):
                                if np.cos(cth_reco[i_cth_r])>0:
                                    continue 
                                flux_integrated = dblquad(lambda energy, angle: metaflux(energy, angle, key), a=cth_true[i_cth_t+1], b=cth_true[i_cth_t], gfun=e_true[i_e_t], hfun=e_true[i_e_t+1])[0]

                                ldata[i_e_r][i_cth_r] = flux[i_e_t][i_e_d][i_cth_t][i_cth_r]*m2_to_cm2*flux_integrated

def build_flux(*dataobjs):
    file_dir = os.path.join(config["datapath"], "charm_search_supplemental/")
    flavors = ["e", "mu", "tau"]
    # flavors += ["e_bar", "mu_bar", "tau_bar"]
    interactions = ["cc", "nc"]

    # energies, cos_th 
    casc_flux = np.zeros(shape=(28,10))
    casc_err = np.zeros(shape=(28,10))
    print("Integrating Flux over bins") 
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
