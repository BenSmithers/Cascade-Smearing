#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/combo/py3-v4.1.1/


##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.XLevel/build/
##METAPROJECT /data/user/bsmithers/metaprojects/snobo/py3-v4.0.1/ 



print(" ===== Welcome to SnowSuite =====")
print("")
print("                /\ ")
print("           __   \/   __")
print("           \_\_\/\/_/_/ ")
print("             _\_\/_/_ ")
print("            __/_/\_\__ ")
print("           /_/ /\/\ \_\ ")
print("                /\ ")
print("                \/ ")
print("")
print("=== Running Lepton Generation ===")

import time
from optparse import OptionParser

start_time = time.time()

#Create parser to accept arguments from user
parser = OptionParser()
parser.add_option("-s","--seed",dest="seed",default=1,type="int",   
                  help="Require a unique number under 4e9 for each job.")
parser.add_option("-o","--outfile",dest="outfile",type="string",    
                  help="Outfile name",
                  default = None)
parser.add_option("-i","--index",dest="index",type="string",        
                  help="Spectral index, positive", 
                  default = "2")
parser.add_option("--Emin", dest="Emin",type="string",              
                  help="Min energy to simulate [GeV]", 
                  default = "200")
parser.add_option("--Emax", dest="Emax",
                  type="string",
                  help="Maximum energy to smulate [GeV]", 
                  default = "1000000")
parser.add_option("-n","--nEvents", dest="nEvents",
                  type="string",
                  help="Number of events", 
                  default = "100000")
parser.add_option("--Zmin", dest="Zmin",
                  type="string",
                  help="Min zenith in degrees, 90 is horizon", 
                  default = "80")
parser.add_option("--Zmax", dest="Zmax",
                  type="string",
                  help="Max zenith in degrees, 180 is core", 
                  default = "180")
parser.add_option("--radius", dest="radius",
                  type="string",
                  help="The radius of the disk to radomly inject on.", 
                  default = "800")
parser.add_option("--length", dest="length",
                  type="string",          
                  help="The length of the cylinder to randomly inject in.", 
                  default = "1200")
parser.add_option("-t","--type", dest="type",
                  type="string",         
                  help="Primary neutrino type (SM only). Electron, Muon, Tau, or All", 
                  default = "Muon")
parser.add_option("-c","--interaction", dest="interaction",
                  type="string",
                  default="NC",
                  help="What kind of interaction do you want to simulate?")
options,args = parser.parse_args()

#pass the arguments to python variables
outfile     = options.outfile
seed        = int(options.seed)
Emin        = float(options.Emin)
Emax        = float(options.Emax)
index       = float(options.index)
nEvents     = int(options.nEvents)
Zmin        = float(options.Zmin)
Zmax        = float(options.Zmax)
radius      = float(options.radius)
length      = float(options.length)
nu_type     = options.type
interaction = options.interaction


print("=====> Importing necessary packages!")
# import python packages
import os
from os.path import expandvars
import numpy as np
import sys
import warnings

#import icecube ones
from I3Tray import *
from icecube import icetray, dataio, dataclasses, LeptonInjector
from icecube import phys_services
from icecube.icetray import I3Units

# verify an outfile is specified
if outfile is None:
    raise Exception("You must specify an outfile with the -o flag. Run `./1-process_gen.py --help' for help.")

print("=====> Generating Primaries")

print("")
print("Outfile: "+outfile)
print("Energy Range: "+str(Emin)+"GeV -> "+str(Emax)+"GeV")
print("Spectral Index: -"+str(index))
print("Number of Events: "+str(nEvents) + " "+nu_type+"s")
print("Z range: "+str(Zmin)+" -> "+str(Zmax))
print("Cylinder length: "+str(length)+", radius: "+str(radius))


tray = I3Tray()
tray.AddService("I3GSLRandomServiceFactory")(("Seed",seed))
tray.AddService("I3EarthModelServiceFactory","Earth")
tray.AddModule("I3InfiniteSource")(("Stream",icetray.I3Frame.DAQ))


#takes a string representing a particle, returns a pair of it and its anti-particle.
#   also accepts 'all', in which case it returns one of every charged lepton
def get_fs_particle_CC(nu_type):
    if nu_type=='Tau' or nu_type=='tau':
        return([ dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.TauPlus])
    elif nu_type=='Muon' or nu_type=='muon':
        return([ dataclasses.I3Particle.ParticleType.MuMinus , dataclasses.I3Particle.ParticleType.MuPlus])
    elif nu_type=='electron' or nu_type=='Electron':
        return([ dataclasses.I3Particle.ParticleType.EMinus  , dataclasses.I3Particle.ParticleType.EPlus])
    elif nu_type=='all' or nu_type =='All':
        return([ dataclasses.I3Particle.ParticleType.TauMinus, dataclasses.I3Particle.ParticleType.TauPlus,
                 dataclasses.I3Particle.ParticleType.MuMinus , dataclasses.I3Particle.ParticleType.MuPlus,
                 dataclasses.I3Particle.ParticleType.EMinus  , dataclasses.I3Particle.ParticleType.EPlus])
    else:
        raise Exception("I don't recogize what a '{}' is. Try Electron, Muon, Tau, or All.".format(nu_type))

def get_fs_particle_NC(nu_type):
    if nu_type=='Tau' or nu_type=='tau':
        return([ dataclasses.I3Particle.ParticleType.NuTau, dataclasses.I3Particle.ParticleType.NuTauBar])
    elif nu_type=='Muon' or nu_type=='muon':
        return([ dataclasses.I3Particle.ParticleType.NuMu , dataclasses.I3Particle.ParticleType.NuMuBar ])
    elif nu_type=='electron' or nu_type=='Electron':
        return([ dataclasses.I3Particle.ParticleType.NuE  , dataclasses.I3Particle.ParticleType.NuEBar  ])
    elif nu_type=='all' or nu_type =='All':
        return([ dataclasses.I3Particle.ParticleType.NuTau   , dataclasses.I3Particle.ParticleType.NuTauBar,
                 dataclasses.I3Particle.ParticleType.NuMu    , dataclasses.I3Particle.ParticleType.NuMuBar,
                 dataclasses.I3Particle.ParticleType.NuE     , dataclasses.I3Particle.ParticleType.NuEBar])
    else:
        raise Exception("I don't recogize what a '{}' is. Try Electron, Muon, Tau, or All.".format(nu_type))


def is_CC(pair):
    """
    Takes a pair of particles, tells if it's a CC interaction
    """
    if pair[1]==dataclasses.I3Particle.ParticleType.Hadrons:
        if (pair[0]==dataclasses.I3Particle.ParticleType.TauMinus\
        or pair[0]==dataclasses.I3Particle.ParticleType.TauPlus\
        or pair[0]==dataclasses.I3Particle.ParticleType.MuMinus\
        or pair[0]==dataclasses.I3Particle.ParticleType.MuPlus\
        or pair[0]==dataclasses.I3Particle.ParticleType.EMinus\
        or pair[0]==dataclasses.I3Particle.ParticleType.EPlus):
            return True
        else:
            return False
    else:
        return False 

def is_nubar(pair):
    if (pair[0]==dataclasses.I3Particle.ParticleType.TauMinus\
    or pair[0]==dataclasses.I3Particle.ParticleType.NuTau\
    or pair[0]==dataclasses.I3Particle.ParticleType.MuMinus\
    or pair[0]==dataclasses.I3Particle.ParticleType.NuMu\
    or pair[0]==dataclasses.I3Particle.ParticleType.EMinus\
    or pair[0]==dataclasses.I3Particle.ParticleType.NuE):
        return False
    else:
        return True


# need to prepare a list of inectors with which to instantiate the generator module
#   the cross sections and final state particles will all depend on which interaction we're doing. 

#dis_xs_folder = '/data/user/bsmithers/cross_sections/'
dis_xs_folder = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/xs_iso/'
injector_list = []

def get_xs(current, nuCharge):
    if not isinstance(current, str):
        raise TypeError("Expected {}, not {}".format(str, type(current)))
    if not isinstance(nuCharge, str):
        raise TypeError("Expected {}, not {}".format(str, type(nuCharge)))

    if not (nuCharge=="nu" or nuCharge=="nubar"):
        raise ValueError("Expected 'nu' or 'nubar', got {}".format(nuCharge))

    if current.lower() == "gr":
        return('{}single_xs_gr.fits'.format(dis_xs_folder) , '{}total_xs_gr.fits'.format(dis_xs_folder))
    elif current.lower()=="nc":
        return('{}dsdxdy_{}_CC_iso.fits'.format(dis_xs_folder,nuCharge), '{}sigma_{}_CC_iso.fits'.format(dis_xs_folder,nuCharge))
    elif current.lower()=="gr":
        return('{}dsdxdy_{}_NC_iso.fits'.format(dis_xs_folder,nuCharge), '{}sigma_{}_NC_iso.fits'.format(dis_xs_folder,nuCharge))

    else:
        raise ValueError("Unrecognized current: {}".format(currnet.lower()))

if interaction=='GR' or interaction=='gr':
    
    # in a GR interaciton, you have a NuEBar annihilate with an e-. Need to conserve charge + lepton no. 
    party_pairs = [[dataclasses.I3Particle.ParticleType.EMinus     , dataclasses.I3Particle.ParticleType.NuEBar    ],
                [dataclasses.I3Particle.ParticleType.MuMinus    , dataclasses.I3Particle.ParticleType.NuMuBar   ],
                [dataclasses.I3Particle.ParticleType.TauMinus   , dataclasses.I3Particle.ParticleType.NuTauBar  ],
                [dataclasses.I3Particle.ParticleType.Hadrons    , dataclasses.I3Particle.ParticleType.Hadrons   ]]
elif interaction=='CC' or interaction=='cc':
    
    party_pairs = [ [ entry , dataclasses.I3Particle.ParticleType.Hadrons ] for entry in get_fs_particle_CC(nu_type)]

elif interaction=='NC' or interaction=='nc':

    party_pairs = [ [ entry , dataclasses.I3Particle.ParticleType.Hadrons ] for entry in get_fs_particle_NC(nu_type)]

else:
    raise Exception("I don't know what kind of interaction '{}' is. Try GR, CC, or NC".format(interaction))

for pair in party_pairs:
    do_ranged = pair[0]==dataclasses.I3Particle.ParticleType.MuMinus or pair[0]==dataclasses.I3Particle.ParticleType.MuPlus

    xs_pair = get_xs(interaction.lower(), is_nubar(pair))

    injector_list.append( LeptonInjector.injector( 
            NEvents                             = int(nEvents / len(party_pairs)), # want this normalized to the number of events asked for! 
            FinalType1                          = pair[0], # the 
            FinalType2                          = pair[1],
            DoublyDifferentialCrossSectionFile  = xs_pair[0],
            TotalCrossSectionFile               = xs_pair[1],
            Ranged                              = do_ranged)
            )


tray.AddModule("MultiLeptonInjector",
    EarthModel = "Earth",
    Generators = injector_list,
    MinimumEnergy   =   Emin * I3Units.GeV,
    MaximumEnergy   =   Emax * I3Units.GeV,
    MaximumZenith   =   Zmax * I3Units.deg,
    MinimumZenith   =   Zmin * I3Units.deg,
    PowerLawIndex   =   index,
    InjectionRadius =   radius * I3Units.meter,
    EndcapLength    =   length * I3Units.meter
    )

# adds a header to the DAQ frames
event_id = 1;
def get_header(frame):
    global event_id
    header                  = dataclasses.I3EventHeader()
    header.event_id         = event_id
    header.run_id           = seed
    event_id                += 1
    frame["I3EventHeader"]  = header

tray.AddModule(get_header,streams=[icetray.I3Frame.DAQ]);
print("Creating lic file at " + os.path.abspath(os.path.join(outfile, os.pardir))  + "/Generation_data.lic")
outpath = ''

#this splits the outpath based around the slashes, and adds them all up until the second to last 
for path in outfile.split('/')[:-1]:
    outpath += path + '/'


tray.AddModule("InjectionConfigSerializer", OutputPath = outpath  + "/{}_{}_Generation_data.lic".format(interaction, nu_type))
tray.AddModule("I3Writer",Filename=outfile)
outputKeys = ['EventProperties','I3MCTree']

print("Preparations Complete, executing...")
tray.Execute()

end_time = time.time()
print("Done.")
print("That took "+str(end_time - start_time)+" seconds.")
