#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/simulation.trunk/build

##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/simulation.trunk.Noise/build/

##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Metaprojects/simulation.V06-01-00/build/
##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Metaprojects/simulation.spe/build/
##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Metaprojects/simulation.V06-00-03/build/

import os, time, sys, subprocess
from I3Tray import *
from icecube import icetray
from icecube.simprod import segments
from optparse import OptionParser
from os.path import expandvars
from icecube import icetray, dataclasses, dataio, phys_services, simprod
from icecube.icetray import I3PacketModule
from icecube.sim_services import bad_dom_list_static
from icecube.simprod import util
from icecube.simprod.segments.MultiDomEffSample import *
from icecube.simprod import segments
from icecube import clsim
import string
from os.path import expandvars, exists, isdir, isfile

import numpy as np
from icecube.clsim import GetDefaultParameterizationList
from icecube.clsim import GetFlasherParameterizationList
from icecube.clsim import GetHybridParameterizationList

#from .I3CLSimMakePhotons import I3CLSimMakePhotons
#from .I3CLSimMakeHitsFromPhotons import I3CLSimMakeHitsFromPhotons

# This script runs the detector level
# generation -> photon -> DET -> level1 -> level2 -> xlevel -> h5

usage = "usage: %prog [options] inputfile"
parser = OptionParser(usage)

parser.add_option("-g", "--gcdfile",type="string", dest="gcdfile", help="Location of the GCD file.",
		default = "/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz")

parser.add_option("-i", "--infile",type="string", dest="infile", help="infile location",
		default = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Nominal/Phot/00001-01000/Phot_00327.i3.zst')

		#default = '/data/user/saxani/test_Photon_level.i3.zst')
parser.add_option("-o", "--outfile",type="string", dest="outfile", help="Outfile location.",
		default = '/data/user/saxani/test_Det.i3.zst')

parser.add_option("-s", "--seed",type="string",default='1',
		dest="seed", help="Seed for the random number generator. Must be below 2,000,000,000")

parser.add_option("-e", "--domeff",type="string", default="0.95",
		dest="domeff", help="Which DOM efficiency do you want to simulate?")

parser.add_option("--hole_ice",type="string", default="p1_0.30_p2_-1.0",#default="as.flasher_p1_0.30_p2_-1",
		dest="hole_ice", help="Which Hole Ice model do you want to use?")

parser.add_option("--generated_domeff",type="string", default="1.70",
		dest="generated_domeff", help="What was the maximum domeff that the photons were generated at?")

parser.add_option("-n", "--nFrames",type="int", default="0",
		dest="nFrames", help="The number of Frames to process")

parser.add_option("-m", "--move",type="string", default="False",
		dest="move", help="Do you want to write the file to a temp dir and then copy it to the correct location?")

parser.add_option("--osg",type="string", default="False",
		dest="osg", help="Do you want to run on the OSG?")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

#======================================
gcdfile 	= options.gcdfile
infile 		= options.infile
outfile 	= options.outfile
osg 		= options.osg
seed 		= int(options.seed)
domeff 		= float(options.domeff)
move 		= options.move
hole_ice 	= options.hole_ice
generated_domeff = float(options.generated_domeff)

# Randomly select one of the background files to inject events into simulation
if seed == '':
	print('Invalid Seed')
	sys.exit()
if domeff == '':
	print('Invalid domeff')
	sys.exit()

hole_ice_location = os.path.expandvars("$I3_BUILD/ice-models/resources/models/angsens_flasher/"+hole_ice)

if not os.path.isfile(hole_ice_location):
	hole_ice_location = os.path.expandvars("$I3_BUILD/ice-models/resources/models/angsens/"+hole_ice)
	if not os.path.isfile(hole_ice_location):
		print(hole_ice_location) 
		print('HoleIce file does not exist')
		sys.exit()
print('Using HoleIce: '+hole_ice_location)

if osg == 'True':
    def copy_to_OSG(NPX_file):
        subprocess.check_call(['globus-url-copy','-nodcau','-rst', NPX_file,"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1])])
        print('copy worked')

    def copy_to_NPX(NPX_file):
        subprocess.check_call(['globus-url-copy','-nodcau','-rst',"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1]), NPX_file])

    gcdfile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + gcdfile)
    gcdfile = str(os.getcwd()+'/'+gcdfile.split('/')[-1])
    infile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + infile)
    infile = str(os.getcwd()+'/'+infile.split('/')[-1])
    copy_to_OSG(gcdfile_NPX)
    copy_to_OSG(infile_NPX)
    outfile  = str('gsiftp://gridftp.icecube.wisc.edu' + outfile)
    outfile_temp = str(os.getcwd()+'/'+outfile.split('/')[-1])
    '''
    try:
        copy_to_OSG(outfile)
    except:
        print('Outfile does not exsist. Continuing...')
    if os.path.isfile(outfile_temp):
        print('Outfile already processed. Exiting.')
        sys.exit()
    '''
else:
	outfile_temp = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/temp/'+outfile.split('/')[-1]
	outfile = options.outfile

infiles = [gcdfile, infile]
#====================================

print('GCD file: ' 		+ infiles[0])
print('Infile: ' 		+ infiles[1])
print('Outfile: '		+ outfile)
print('Seed: '			+ str(seed))
print('DOM efficiency: '+ str(domeff))
print('Generated efficiency: '+ str(generated_domeff))

# set up a random number generator. Use GSL, random number has to be under 4,294,967,296, since it is 32 bit.
randomService = phys_services.I3GSLRandomService( seed = int(seed))

tray = I3Tray()
tray.context['I3RandomService'] = randomService
tray.Add("I3Reader", FilenameList=infiles)

# The only thing you might need to change is the name of the MC pe series map and the efficiency it was generated at.

if hasattr(dataclasses, "I3ModuleGeo"):
	tray.AddModule("I3GeometryDecomposer", "decomposeGeometry",
		If=lambda f: True)
		#If=lambda frame: If(frame) and ("I3OMGeoMap" not in frame) and ("I3ModuleGeoMap" not in frame))


tray.AddSegment(clsim.I3CLSimMakeHitsFromPhotons, "makeHitsFromPhotons",
	PhotonSeriesName	="PhotonSeriesMap",
	MCPESeriesName		="I3MCPESeriesMap",
	RandomService		='I3RandomService',
	DOMOversizeFactor	= 1.0,
    #GCDFile             = gcdfile,
    #IceModelLocation    =  expandvars("$I3_BUILD/ice-models/resources/models/spice_3.2.2/"),
	UnshadowedFraction	= generated_domeff,
	HoleIceParameterization	= hole_ice_location, #	expandvars("$I3_BUILD/ice-models/resources/models/angsens_flasher/"+hole_ice),
	If					= lambda f: True
	)

I3MCPESeriesMap = "I3MCPESeriesMap"

# Rename the selected MC pe series.
if I3MCPESeriesMap != 'I3MCPESeriesMap':
	def mover(frame):
		if 'I3MCPESeriesMap' in frame:
			del frame['I3MCPESeriesMap']
		frame['I3MCPESeriesMap'] = frame[mcpeseriesname]
	tray.Add(mover, Streams=[icetray.I3Frame.DAQ])

# This function will take the MC pe series map and downsample it. The output will be ''I3MCPESeriesMap_0.99''
tray.AddSegment(MultiDomEffSample, "multiDomEff",
	InputSeriesName 	= I3MCPESeriesMap,
	RandomService 		= randomService,
	GeneratedEfficiency = generated_domeff,
	SampleEfficiencies 	= [domeff],
    )

InputPESeriesMap="I3MCPESeriesMap_"+str(domeff)

def renameMCTree(frame):
	if frame.Has("I3MCTree_preMuonProp"):
		if frame.Has("I3MCTree"):
			del frame['I3MCTree']
		frame["I3MCTree"] = frame["I3MCTree_preMuonProp"]
		del frame["I3MCTree_preMuonProp"]
tray.AddModule(renameMCTree, "_renameMCTree", Streams=[icetray.I3Frame.DAQ])

MultisimAmplitudes = []
MultisimPhases = []
MultisimModes = []
IceXModelNumber = 0

'''
def getMultisim(frame):
	global MultisimAmplitudes
	global MultisimPhases
	global MultisimModes
	global IceXModelNumber
	MultisimAmplitudes = frame['MultisimAmplitudes']
	MultisimPhases = frame['MultisimPhases']
	MultisimModes = frame['MultisimModes']
	#IceXModelNumber = frame['IceXModelNumber']
	#del frame['MultisimAmplitudes']
	#del frame['MultisimPhases']
	#del frame['MultisimModes']
	#del frame['IceXModelNumber']
	#MultisimAmplitudes = frame['AmplitudeArray']
	#MultisimPhases = frame['PhaseArray']
	#MultisimModes = frame['ModeNumber']
	#del frame["AmplitudeArray"]
	#del frame["PhaseArray"]
	#del frame["ModeNumber"]
	#print(MultisimAmplitudes)
tray.AddModule(getMultisim, "_getMultisim", Streams=[icetray.I3Frame.Stream('M')])

def setMultisim(frame):
	global MultisimAmplitudes
	global MultisimPhases
	global MultisimModes
	del frame['MultisimAmplitudes']
	del frame['MultisimPhases']
	del frame['MultisimModes']
	frame.Put('MultisimAmplitudes'     ,dataclasses.I3VectorDouble(MultisimAmplitudes))
	frame.Put('MultisimPhases'     ,dataclasses.I3VectorDouble(MultisimPhases))
	frame.Put('MultisimModes'     ,dataclasses.I3VectorDouble(MultisimModes))
	#frame.Put('IceXModelNumber'     ,dataclasses.I3String(IceXModelNumber))
	#print(MultisimAmplitudes)
tray.AddModule(setMultisim, "_setMultisim", Streams=[icetray.I3Frame.Stream('Q')])
'''

# This first step does the noise.
from icecube import vuvuzela
from scipy.special import erfc
from numpy import log10
from copy import deepcopy

tray.AddModule("Inject", "AddNoiseParams",
               InputNoiseFile = expandvars("$I3_SRC/vuvuzela/resources/data/parameters.dat"),
                                      )
def VuvuzelaParamScaler(frame, Qd = 0.25 ):
    print "Modifying vuvuzela params according SPE efficiency" 
    calibration = deepcopy(frame['I3Calibration'])
    dom_cal = deepcopy(calibration.dom_cal)
    frame.Delete('I3Calibration')
    scalings = dataclasses.I3MapKeyDouble()   
    areaAbove = dataclasses.I3MapKeyDouble()   
    for oneOM in dom_cal.keys():
        if oneOM.om > 60 : continue # ignoring IceTop
        newSPE = dom_cal[oneOM].combined_spe_charge_distribution
        xi_val  = 0.0
        xi_val += (newSPE.exp1_width*newSPE.exp1_amp*np.exp( - Qd/newSPE.exp1_width)) 
        xi_val += (newSPE.exp2_width*newSPE.exp2_amp*np.exp( - Qd/newSPE.exp2_width))
        xi_val += (np.sqrt(2.*np.pi)*newSPE.gaus_amp*newSPE.gaus_width*(1.0 - 0.5*erfc((newSPE.gaus_mean - Qd)/(np.sqrt(2) * newSPE.gaus_width)))) 
        areaAbove[oneOM] = float(xi_val) # Area above discriminator
        #print oneOM , xi_val,
        xi_val = xi_val / 0.880149 
        #print xi_val
        scalings[oneOM] = xi_val # ratio of areas
        oldTR = float(dom_cal[oneOM].dom_noise_thermal_rate)
        oldDR = float(dom_cal[oneOM].dom_noise_decay_rate)
        oldSH = float(dom_cal[oneOM].dom_noise_scintillation_hits)
        oldSM = float(dom_cal[oneOM].dom_noise_scintillation_mean) # mean of the gaussian
        oldSS = float(dom_cal[oneOM].dom_noise_scintillation_sigma) # width of the gaussian
        # XXX: I am making conversion to float to avoid potential "pointer" problems
        # This might be not entirely necessary, bur rather a paranoid precaution
        newTR  = oldTR/xi_val
        newDR  = oldDR*(oldSH + 1.0)/(oldSH + xi_val)
        newSH  = oldSH/xi_val
        newSM  = oldSM + log10(xi_val) 
        newSM += log10(1.0 - (1.0 - xi_val)*oldSH/(oldSH+1.0))
        newSS  = oldSS

        # Writing parameters to the calibration
        calibration.dom_cal[oneOM].dom_noise_thermal_rate        = newTR
        calibration.dom_cal[oneOM].dom_noise_decay_rate          = newDR
        calibration.dom_cal[oneOM].dom_noise_scintillation_hits  = newSH
        calibration.dom_cal[oneOM].dom_noise_scintillation_mean  = newSM
        calibration.dom_cal[oneOM].dom_noise_scintillation_sigma = newSS
    frame['I3Calibration'] = calibration # saving new calibration
    frame['SPEScalingFactors'] = scalings # Putting parameters to the frame for book-keeping
    frame['SPEAbove'] = areaAbove
    return True

tray.AddModule(VuvuzelaParamScaler, "scaler", 
                Streams = [icetray.I3Frame.Calibration])

InputPESeriesMap_withoutNoise = InputPESeriesMap + "WithoutNoise"
tray.Add("Rename", "RenamePESeriesMap",
         Keys=[InputPESeriesMap, InputPESeriesMap_withoutNoise])

tray.AddModule("Vuvuzela", "vuvuzela_DetectorSim",
           InputHitSeriesMapName  = InputPESeriesMap_withoutNoise,
           OutputHitSeriesMapName = InputPESeriesMap,
           StartWindow            = -11.*I3Units.microsecond,
           EndWindow              = 11.*I3Units.microsecond,#Parameters from https://code.icecube.wisc.edu/projects/icecube/browser/IceCube/sandbox/hignight/LE_simulation_scripts/step_3_det_general_lowdt.py#L174
           IceTop                 = False,
           InIce                  = True,
           ScaleFactor            = 1.0,
           DeepCoreScaleFactor    = 1,
           DOMsToExclude          = [],
           RandomService          = randomService,
           SimulateNewDOMs        = True,
           DisableLowDTCutoff     = True,
           UseIndividual          = True
           )

# Run the detector simulation.
tray.AddSegment(segments.DetectorSim, "DetectorSim",
	RandomService = 'I3RandomService',
	RunID = seed,
	GCDFile = gcdfile,
	InputPESeriesMapName = InputPESeriesMap,
	KeepMCHits = True,
	SkipNoiseGenerator = True,# I'm skipping noise, and doing it after. simprod-scripts/python/segments/DetectorSim.py
	KeepPropagatedMCTree = True,
	TimeShiftSkipKeys=["LeptonInjectorProperties","EventProperties","GenerationSpec"])

remove_keys = ["PhotonSeriesMap","CalibratedWaveforms",InputPESeriesMap,InputPESeriesMap+'WithoutNoise']
tray.AddModule("Delete", "cleanup2", Keys = remove_keys)


#def RemoveMCPESeries(frame):
#	thekeys = frame.keys()
#	if frame.Has("I3MCTree_preMuonProp") and not frame.Has("I3MCTree"):
#		frame["I3MCTree"] = frame["I3MCTree_preMuonProp"]
#		del frame["I3MCTree_preMuonProp"]
#
#tray.AddModule(RemoveMCPESeries,Streams=[icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')],)

i3streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')]

print(outfile)
if move == "True":
	tray.AddModule("I3Writer","i3writer")(
		("Filename",outfile_temp),
		("Streams",i3streams)
		)
else:
	tray.AddModule("I3Writer","i3writer")(
			("Filename",outfile),
			("Streams",i3streams),
			#("compressionlevel", 16)
			)

	tray.AddModule("TrashCan", "the can")
icetray.logging.set_level('WARN')

if(options.nFrames==0):
	tray.Execute()
else:
	tray.Execute(int(options.nFrames))
tray.Finish()
del tray

end_time = time.time()
print("done ...")

#====================================
if options.move == 'True':
	if osg == "True":
		copy_to_NPX(outfile)
	else:
		os.system(str("mv "+str(outfile_temp) + " " +str(outfile)))
else:
	if osg =="True":
		print('You need options.move to be True if using the OSG')
		pass
#====================================
