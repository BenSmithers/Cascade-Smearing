#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/simulation.trunk/build

##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/testing_meta/build/


#==============================
# Spencer N. Axani
# Questions? saxani@mit.edu
#==============================

import os
import sys,getopt,string
from os.path import expandvars
from icecube.simprod  import ipmodule
from icecube.simprod.util import CombineHits, DrivingTime, choose_max_efficiency
from icecube.simprod.util.fileutils import download,untar,isurl
from optparse import OptionParser
import logging
from icecube import icetray, dataclasses, dataio, phys_services, interfaces
from I3Tray import I3Tray,I3Units
from icecube.simprod.modules import ClSim
import subprocess
import time
from icecube.simprod import segments
from icecube import clsim

start_time = time.time()

usage = "usage: %prog [options] inputfile" 
parser = OptionParser(usage)
parser.add_option("--osg", type="string", default="False")

parser.add_option("-g","--gcdfile", 	type="string",default="/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz")

parser.add_option("-i","--infile", 		type="string",default="/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Nominal/Gen/00001-01000/Gen_00381.i3.zst") 

parser.add_option("-o","--outfile", 	type="string",default="/data/user/gparker/ALL_NSI_FLUXES/mu_tau/BFR_MC/process-Phot-out/test_Photon_level.i3.zst")

parser.add_option("-s","--seed", 		type="string",default="1")

parser.add_option("--usegpus",    		type="string",default="True")

parser.add_option("--maxparallelevents",      type="int",default="350")

parser.add_option("-e", "--generated_domeff",   	type="float",default="1.7")

parser.add_option("--ice_model",      	type="string",default = "spice_bfr-dv1_complete", help = 'spice_3.2, spice_mie')

parser.add_option("-m","--move",      	type="string",default = "True")

(options,args) = parser.parse_args() 


gcdfile 	= options.gcdfile
infile 		= options.infile
outfile 	= options.outfile
osg	 		= options.osg
move 		= options.osg
generated_domeff		= float(options.generated_domeff)
seed 		= options.seed
maxparallelevents = options.maxparallelevents

if options.usegpus == "True":
	usegpus = True
	usecpus = False
else:
	usegpus = False
	usecpus = True

if options.ice_model not in ['spice_3','spice_3_noAni','spice_3.2','spice_mie','SpiceHD','spice_lea','sys1_spice_3.2','sys2_spice_3.2','sys3_spice_3.2','sys4_spice_3.2','sys5_spice_3.2','sys6_spice_3.2','spice_bfr-dv1_complete']:
	print('Invalid Icemodel: '+options.ice_model)
	sys.exit()
	
icemodellocation = '/data/user/gparker/ALL_NSI_FLUXES/mu_tau/BFR_MC/'

ice_model = options.ice_model

print('GCD file: '+gcdfile)
print('Infile: '+infile)
print('Outfile: '+outfile)
print('DOM eff: ' + str(generated_domeff)) 
print('IceModel: '+ icemodellocation + ice_model)
print('Using GPUs: '+str(usegpus))

if osg == 'True':
	def copy_to_OSG(NPX_file):
		subprocess.check_call(['globus-url-copy','-nodcau','-rst', NPX_file,"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1])])
	def copy_to_NPX(NPX_file):
		subprocess.check_call(['globus-url-copy','-nodcau','-rst',"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1]), NPX_file])

	gcdfile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + gcdfile)
	gcdfile = str(os.getcwd()+'/'+gcdfile.split('/')[-1])

	infile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + infile)
	infile = str(os.getcwd()+'/'+infile.split('/')[-1])

	copy_to_OSG(gcdfile_NPX)
	copy_to_OSG(infile_NPX)

	infiles = [gcdfile,infile]
	outfile  = str('gsiftp://gridftp.icecube.wisc.edu' + outfile)
	outfile_temp = str(os.getcwd()+'/'+outfile.split('/')[-1])
else:
	infiles = [gcdfile, infile]
	outfile_temp = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/temp/'+outfile.split('/')[-1]

tray = I3Tray()

randomServicePhotons     = phys_services.I3GSLRandomService(seed = int(seed))
randomServicePropagators = phys_services.I3GSLRandomService(seed = int(seed))

tray.AddModule("I3Reader", "reader",filenamelist=infiles)

tray.AddModule("Rename",keys=["I3MCTree","I3MCTree_preMuonProp"])
tray.AddSegment(segments.PropagateMuons, "PropagateMuons",
	RandomService = randomServicePropagators)

tray.context['I3RandomService'] = randomServicePhotons
tray.AddSegment(clsim.I3CLSimMakePhotons, "MakingPhotons",
		UseCPUs								= usecpus,
		UseGPUs								= usegpus,
		UseOnlyDeviceNumber					= None,
		MCTreeName							= "I3MCTree",
		OutputMCTreeName					= None,
		FlasherInfoVectName					= None,
		FlasherPulseSeriesName				= None,
		MMCTrackListName					= "MMCTrackList",
		PhotonSeriesName					= "PhotonSeriesMap",
		ParallelEvents						= maxparallelevents,
		TotalEnergyToProcess				= 0,
		RandomService						= "I3RandomService",
		IceModelLocation					= icemodellocation + ice_model,
		DisableTilt							= False,
		UnWeightedPhotons					= False,
		UnWeightedPhotonsScalingFactor		= None,
		UseGeant4							= False,
		CrossoverEnergyEM					= None,
		CrossoverEnergyHadron				= None,
		UseCascadeExtension					= False,
		StopDetectedPhotons					= True,
		PhotonHistoryEntries				= 0,
		DoNotParallelize					= False,
		DOMOversizeFactor					= 1.0,
		UnshadowedFraction					= float(generated_domeff),
		HoleIceParameterization 			= expandvars('$I3_BUILD/ice-models/resources/models/angsens/as.nominal'),
		WavelengthAcceptance				= None,
		DOMRadius							= 0.16510*icetray.I3Units.m, # 13" diameter
		OverrideApproximateNumberOfWorkItems= None,
		IgnoreSubdetectors                  = ['IceTop'],
		ExtraArgumentsToI3CLSimModule		= dict(),
		If=lambda f: True )

i3streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')]
if move == "True":
	tray.AddModule("I3Writer","i3writer")(
		("Filename",outfile_temp),
		("Streams",i3streams)
		)
else:
	tray.AddModule("I3Writer","i3writer")(
		("Filename",outfile),
		("Streams",i3streams)
		)
tray.Execute()
tray.Finish()

if move == 'True':
	if osg == "True":
		copy_to_NPX(outfile)
	else:
		os.system(str("mv "+str(outfile_temp) + " " +str(outfile)))
else:
	if osg =="True":
		print('You need options. Move to be True if using the OSG')
		pass

print("Time: " + str(time.time() - start_time))
del tray

