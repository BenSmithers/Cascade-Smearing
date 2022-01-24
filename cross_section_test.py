"""
Ben Smithers 

Simple script to facilitate nuSQuIDS DIS cross section sampling
"""

import sys
legacy = False
try:
    import nuSQuIDS as nsq
except ImportError:
    import nuSQUIDSpy as nsq
    legacy = True
import numpy as np # useful for energy ranges
import time # so we can wait
import os 

constants = nsq.Const()

nBins = 100
eMin = 1.*constants.GeV
eMax = (10**12)*constants.GeV

energies = np.logspace( np.log10(eMin), np.log10(eMax), nBins)

# I use these dictionaries since the keys will be helpful in making plot labels 
flavors = {'electron' : nsq.NeutrinoCrossSections_NeutrinoFlavor.electron,
           'muon' : nsq.NeutrinoCrossSections_NeutrinoFlavor.muon,
           'tau': nsq.NeutrinoCrossSections_NeutrinoFlavor.tau
           }

currents = {'CC' : nsq.NeutrinoCrossSections_Current.CC, 
            'NC' : nsq.NeutrinoCrossSections_Current.NC
            }

neut_types = {'neutrino': nsq.NeutrinoCrossSections_NeutrinoType.neutrino, 
              'antineutrino': nsq.NeutrinoCrossSections_NeutrinoType.antineutrino
              }


# This uses some files that are stored on /data/user (or wherever)
# if we run a lot of jobs, the hdf5 files may be locked when we try accessing these 
# so we just wait a little bit and try again and again until it works or we run out of tries
tries_left = 100
sleep_time = 0.1 # seconds 
last_err = ""
cobalt = os.environ.get("_CONDOR_SCRATCH_DIR")
while tries_left > 0:
    try:
        if cobalt==None or cobalt=="" or cobalt==".":
            if legacy:
                xs_obj = nsq.NeutrinoDISCrossSectionsFromTables()
            else:
                xs_obj = nsq.loadDefaultCrossSections()
        else:
            xs_obj = nsq.NeutrinoDISCrossSectionsFromTables(os.path.join(cobalt, "data/csms_square.h5"))
        break
    except RuntimeError as err:
        last_err = str(err)
        time.sleep(sleep_time + np.random.rand()*sleep_time)

    tries_left -= 1

if tries_left<=0:
    raise RuntimeError("Ran out of tries while opening the DIS Cross Sections: {}".format(last_err))


def get_diff_xs( energy, flavor, neutrino, current, e_out=None, x=None):
    """
    This function is just a wrapper to the nusquids function to make sure I get my types right.... 
    """
    if not (isinstance( energy, int) or isinstance(energy, float)):
        raise TypeError("Expected {} for energy, got {}".format(float, type(energy)))
    if not isinstance( flavor, nsq.NeutrinoCrossSections_NeutrinoFlavor):
        raise TypeError("Expected {} for flavor, got {}".format( nsq.NeutrinoCrossSections_NeutrinoFlavor, type(flavor)))
    if not isinstance( current, nsq.NeutrinoCrossSections_Current):
        raise TypeError("Expected {} for current, got {}".format(nsq.NeutrinoCrossSections_Current, type(current)))
    if not isinstance( neutrino, nsq.NeutrinoCrossSections_NeutrinoType):
        raise TypeError("Expected {} for neutrino type, got {}".format(nsq.NeutrinoCrossSections_NeutrinoType, type(neutrino)))

    if (e_out is None) or (x is None):
        return( xs_obj.TotalCrossSection( energy, flavor, neutrino, current))
    else:
        if not (isinstance(e_out,float) or isinstance(e_out,int)):
            raise TypeError("Expected {} for e_out, got {}".format(float, type(e_out)))
        if not (isinstance(x,float) or isinstance(x,int)):
            raise TypeError("Expected {} for x, got {}".format(float, type(x)))
        return(xs_obj.SingleDifferentialCrossSection( energy, e_out, flavor, neutrino, current))

def get_total_flux( energy, flavor, neutrino, current):
    return(get_diff_xs(energy, flavor, neutrino, current))


# I don't want to totally remove this code, but I don't want it called as I import some of the funcitons above
# so I'm just commenting this all out for now. 
skip_plots = True
if not skip_plots:
    import matplotlib
    matplotlib.use('agg')
    import matplotlib.pyplot as plt

    total_xss = {}
    for current in currents:
        for nt in neut_types:
            for flavor in flavors:
                cross = np.array([ get_total_flux( en, flavors[flavor], neut_types[nt], currents[current]) for en in energies])
                total_xss[current+'_'+nt+'_'+flavor] = cross 


    muon_tracks = np.array([0. for i in energies])
    muon_cascades = np.array([0. for i in energies])
    for current in currents:
        for nt in neut_types:
            for flavor in flavors:
                if current=='CC':
                    if flavor=='muon':
                        muon_tracks += 0.5*total_xss[current+'_'+nt+'_'+flavor]
                else:
                    if flavor=='muon':
                        muon_cascades += 0.5*total_xss[current+'_'+nt+'_'+flavor]



    plt.plot(energies, muon_cascades,label="muon cascades")
    plt.plot(energies, muon_tracks,label='muon tracks')

    plt.yscale('log')
    plt.xscale('log')
    plt.xlabel('Neurino Energy [eV]',size=16)
    plt.ylabel(r'Total Cross Section [cm$^{2}$]',size=16)
    plt.legend()
    plt.savefig('test.png',dpi=400)
    plt.clf()

    # Let's try fixing the reconstructed energy 
    # and seeing how the differential cross section looks! 
    recon_energy = 100*constants.TeV
    by_min = -5
    by_max = 0
    in_energies = np.array([recon_energy/(1.-y) for y in (1-np.logspace(by_min, by_max, nBins))])
    diff_xs = [get_diff_xs(energy, flavors['muon'],neut_types['neutrino'],currents['NC'],recon_energy,0.) for energy in in_energies ]

    plt.plot(in_energies/constants.GeV, diff_xs)
    plt.xscale('log')
    plt.yscale('log')
    plt.title("Differential Cross Section from {:.2f} TeV Muon".format(recon_energy/constants.TeV))
    plt.xlabel("Muon Neutrino Energy [GeV]")
    plt.ylabel(r"differential xs [cm$^{2}$GeV$^{-1}$]")
    plt.savefig('test2.png',dpi=400)


