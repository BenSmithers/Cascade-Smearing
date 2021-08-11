import numpy as np
import os
import sys
from cascade.utils import Data, config, get_closest, make_bin_probs
from scipy.integrate import dblquad
from scipy.special import i0 #modified bessel function of the first kind 
from math import pi, cos, sin, sinh, exp, acos
from cascade.deporeco import kappaCalc, get_ang_error

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
    dt = np.transpose(data)

    n_true_e = 100
    n_cth_true = 10
    n_dep_e = 20
    n_cth_r = 2

    assert(n_true_e*n_cth_true*n_dep_e*n_cth_r==len(data))

    e_true      = np.concatenate((np.unique(dt[0]), [dt[1][-1]]))
    cth_true    = np.concatenate((np.unique(dt[2]), [dt[3][-1]]))
    e_depo      = np.concatenate((np.unique(dt[4]), [dt[5][-1]]))
    cth_reco    = np.concatenate((np.unique(dt[6]), [dt[7][-1]]))

    ldata = np.zeros(shape=(n_true_e, n_dep_e,n_cth_true, n_cth_r))

    for i_e_t in range(n_true_e):
        for i_e_d in range(n_dep_e):
            for i_cth_t in range(n_cth_true):
                for i_cth_r in range(n_cth_r):
                    true_i = i_e_t + n_true_e*(i_cth_t + n_cth_true*(i_e_d + n_dep_e*i_cth_r))
                    
                    ldata[i_e_t][i_e_d][i_cth_t][i_cth_r] += data[true_i][8]
    return  { 
            "e_reco":e_depo,
            "e_true":e_true,
            "cth_true":cth_true,
            "cth_reco":cth_reco,
            "data":ldata
            }

m2_to_cm2 = 100*100

def get_cth_error(energy, true_angle, reco_angle):
    """
    Returns the probability for reconstructing one zenith angle, for a given true zenith angle, at a given energy
    """

    dtheta = get_ang_error(energy)*pi/180
    kappa = kappaCalc.eval(dtheta)
    value = i0(kappa*sin(true_angle)*sin(reco_angle)) #modified bessel first kind, order 0
    # i0 is faster than using the generic iv 

    return kappa/(2*sinh(kappa)) * exp(kappa*cos(true_angle)*cos(reco_angle))*value*sin(reco_angle)
       
    


def build_flux_sad(*dataobjs, good_angles = False):
    """
    Use the 'good_angles' argument to approximate the angular resolution! 
    """

    # we have a cobalt location and a regular location
    cobalt = os.environ.get("_CONDOR_SCRATCH_DIR")
    if cobalt==None or cobalt=="" or cobalt==".":
        file_dir = os.path.join(config["datapath"], "charm_search_supplemental", "effective areas")
    else:
        file_dir = os.path.join(cobalt, "data")

    flavors = ["e", "mu", "tau"]
    interactions=["cc", "nc"]
   
    n_true_e = 100
    n_cth_true = 10
    n_dep_e = 20
    n_cth_r = 2

    new_cth_edges = np.arccos(np.linspace(-1,1,n_cth_true+1))
    zero_point_two = acos(0.2) # precalculated so we don't do this every time in the loop

    if good_angles:
        casc_flux = np.zeros(shape=(n_dep_e,n_cth_true))
    else:
        casc_flux = np.zeros(shape=(n_dep_e,n_cth_r)) 
    print("Integrating Flux over bins") 

    def metaflux(energy, angle, key):
        """
        This little function sums the contributions for each flux making up this combined flux 
        """
        net_f = 0.0
        for dobj in dataobjs:
            net_f = net_f + dobj.get_flux((1e9)*energy, key, use_overflow=False, angle=cos(angle))
        return net_f*2*pi*sin(angle)
    
    # we need to iterate over all the files - this is how we do it 
    #   we could have just one loop to loop over the files themselves, but I found that might leave us a little blind
    #   so instead we pick files with conviction and drive 
    for flavor in flavors:
        for interaction in interactions:
            for barness in ["","_bar"]:
                for topo in ["cascade"]: # optionally add in the tracks!
                    comb_flav = flavor+barness

                    # load the file and make the data accessible 
                    filename = "effective_area.per_bin.nu_{}.{}.{}.txt".format(comb_flav, interaction,topo)
                    area_data= quickload_full(os.path.join(file_dir, filename))
                    eff_area=area_data["data"]
                    e_reco = area_data["e_reco"]
                    e_true = area_data["e_true"]
                    cth_true = area_data["cth_true"]
                    cth_reco = area_data["cth_reco"]
                    true_ang = np.arccos(cth_true)

                    # figure out the Data object key this corresponds to 
                    curr = interaction.upper()
                    flav=flavor[0].upper() + flavor[1:].lower()
                    nu="nu" if barness=="" else "nuBar"
                    key = "_".join([flav, nu, curr])
                    print("On key {}, fname {}".format(key, filename))

                    # now we iterate over the entries in the files 
                    for i_e_t in range(n_true_e):
                        e_center = 0.5*(e_true[i_e_t]+e_true[i_e_t+1])
                        for i_cth_t in range(n_cth_true):

                            # integrate the flux over incoming energy and angles 
                            flux_integrated = dblquad(lambda energy, angle: metaflux(energy, angle, key), a=true_ang[i_cth_t+1], b=true_ang[i_cth_t], gfun=e_true[i_e_t], hfun=e_true[i_e_t+1])[0]
                            if flux_integrated==0.0:
                                continue
                            # we'll scale this with the effective area later 

                            # and we need this bin center for splitting up the effective areas into multiple bins 
                            true_center = 0.5*(true_ang[i_cth_t]+true_ang[i_cth_t+1])
                            for i_e_d in range(n_dep_e):
                                for i_cth_r in range(n_cth_r):
                                    area = eff_area[i_e_t][i_e_d][i_cth_t][i_cth_r]
                                    if area==0.0:
                                        continue

                                    # for the "good angles" we split big bins into sub-bins
                                    #   up-going into six bins, down into four bins (each of width 0.2 costheta)
                                    if good_angles:
                                        if i_cth_r==0:
                                            shift = 0
                                            new_bins = new_cth_edges[0:7]
                                        elif i_cth_r==1:
                                            shift = 6
                                            new_bins = new_cth_edges[6:]
                                        else:
                                            raise ValueError("The reco bins changed in a way that'll require some restructuring :(")
                                       
                                        # we split the bins up according to our reconstruction efficiency 
                                        reco_bin_probs = make_bin_probs(lambda ang: get_cth_error(e_center, true_center, ang), new_bins, normalize=True)
                                        if not abs(1.0-sum(reco_bin_probs))<1e-10:
                                            raise ValueError("Something wrong {}".format(reco_bin_probs))
                                        # then split the effective areas, according to the above effiency, into multiple bins
                                        count_check = 0.0
                                        for i_cth_reee in range(len(reco_bin_probs)):
                                            index = shift + i_cth_reee
                                            count_check +=reco_bin_probs[i_cth_reee]
                                            casc_flux[i_e_d][index] = casc_flux[i_e_d][index] + area*m2_to_cm2*flux_integrated*reco_bin_probs[i_cth_reee]
                                        if abs(count_check-1)>(1e-10):
                                            raise ValueError("sum: {}, index is {}, i_cth_reee is {}".format(count_check, index, i_cth_reee))
                                    else:
                                        # don't split it up

                                        # should be in units of [s^-1] now
                                        casc_flux[i_e_d][i_cth_r] = casc_flux[i_e_d][i_cth_r] + area*m2_to_cm2*flux_integrated


    livetime_ratio = 1.0
    seconds_up = livetime_ratio*10.*3600.*24.*365
    casc_flux = casc_flux*seconds_up
    casc_err = np.sqrt(casc_flux)

    return {
            "e_edges":e_reco,
            "a_edges":cth_reco,
            "event_rate":casc_flux,
            "stat_err":casc_err
           }

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
            # the Data objects already 
            net_f = net_f + dobj.get_flux((1e9)*energy, key, angle=np.cos(angle))
        return net_f*2*pi #*np.sin(angle)
    


    for flavor in flavors:
        for interaction in interactions:
            for barness in ["", "_bar"]:

                filename = "effective_area.nu_{}.txt".format(flavor)
                area_data = quickload(os.path.join(file_dir, filename))
                angles = np.arccos(area_data["cos_th"])
                a_sins = np.sin(angles)
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
                        if cos(angles[i_cth])>0.0:
                            continue
                        
                        flux_integrated = dblquad( lambda energy, angle: metaflux(energy, angle, key), a=angles[i_cth+1], b=angles[i_cth],gfun=energies[i_e], hfun=energies[i_e+1] )[0] # double counting cc/nc
                         
                        
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
