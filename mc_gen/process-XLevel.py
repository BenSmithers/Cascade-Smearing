#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.XLevel/build/icerec.XLevel

##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.XLevel.Noise/build/

import time, sys, math, os.path,subprocess, os
from optparse import OptionParser
from shutil import copyfile
from os.path import expandvars
import numpy as np

from I3Tray import *
from icecube import icetray, dataio, dataclasses, STTools, MuonInjector
from icecube import phys_services, DomTools, simclasses, VHESelfVeto
from icecube import lilliput, gulliver, gulliver_modules, paraboloid
from icecube import linefit, MuonGun, WaveCalibrator,wavedeform
from icecube import photonics_service
from icecube import RandomStuff
from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
from icecube import TopologicalSplitter
import icecube.lilliput.segments
from icecube import tableio, hdfwriter
from icecube import MuonGun
from icecube.icetray import I3Units
from icecube.common_variables import direct_hits
from icecube.common_variables import hit_multiplicity
from icecube.common_variables import track_characteristics

load("bayesian-priors")
load("double-muon")
load("libmue")
load('MCParticleExtractor')

start_time = time.time()
print("Starting ...")

class option:
	def __init__(self,short,long,dest,description="",default=None,type="boolean"):
		self.short=short
		self.long=long
		self.dest=dest
		self.description=description
		self.default=default
		self.type=type

	def matches(self,flag):
		return flag==self.short or flag==self.long

	def parse_value(self,value):
		if(self.type=="boolean"):
			return self.parse_boolean(value)
		if(self.type=="int"):
			return self.parse_int(value)
		if(self.type=="float"):
			return self.parse_float(value)

	def parse_boolean(self,value):
		if(value==None or value=='' or value=='True' or value=='true' or value=='1' or value=='yes'):
			return True
		if(value=='False' or value=='false' or value=='0' or value=='no'):
			return False
		raise TypeError("option %s: invalid boolean value: %r" % (self.long, value))

	def parse_int(self,value):
		try:
			return int(value)
		except:
			raise TypeError("option %s: invalid integer value: %r" % (self.long, value))

	def parse_float(self,value):
		try:
			return float(value)
		except:
			raise TypeError("option %s: invalid floating point value: %r" % (self.long, value))

class OptionResult:
	pass

#======= Setup Program Arguments ========#
parser = OptionParser()
parser.add_option("-v","--verbose",dest="diagnoseCuts",default=0,type="int",help = 'Ignore this')
parser.add_option("-n","--nFrames", dest="nFrames", default=0, type="int", help = 'How many frames do you want to process')
parser.add_option("--ice_model", dest="ice_model", default='spice_bfr-dv1_complete', type="string", help = 'Ignore this')

parser.add_option("-g","--gcdfile",dest="gcdfile",type="string",
	default ="/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz")

parser.add_option("-i","--infile", dest="infile",type="string",
	default = '/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Ares/IC86.AVG/L2/domeff_0.97/00001-01000/L2_00_11_00991.i3.zst')

parser.add_option("-o","--outfile",dest="outfile",type="string",
	default = '/data/user/saxani/test_XLevel1.i3.zst')


parser.add_option("-m","--move",dest="move",default='False',type="string", help = 'Generate in tmp dir, then move to correct location?')

parser.add_option("-s","--osg",dest="osg",default='False',type="string", help = 'Is this running on the OSG?')
parser.add_option("-c","--cut",dest="cut",default='False',type="string", help = 'Do not touch this. Should be False.')

options,args = parser.parse_args()
#==================================================#

# Define a timewindow for getting the directhit and hit_multiplicity variables
dh_definitions = [direct_hits.I3DirectHitsDefinition("", -15*I3Units.ns, +200*I3Units.ns)]
ice_model = options.ice_model
if ice_model == 'spice_bfr-dv1_complete':
	ice_model_location = expandvars("$I3_BUILD/mue/resources/ice/spice_bfr-dv1_complete")
else:
	print('Invalid Ice model')
	sys.exit()

if options.cut not in ['True', 'False']:
	print("Cut options must be 'True' or 'False'")
	sys.exit()

# There is a difference in terminology between 2011 and all the other years.
if "NuFSGen2011" in options.infile or "IceCube/2011" in options.infile or "IC86NuFS" in options.infile:
    	print('The year is 2011')
	if 'pass2' in options.infile:
		print('And processing Pass2')
		InIcePulses 				= "InIcePulses"
		SRTInIcePulses 				= "SRTInIcePulses"
		SRTInIcePulses_NoDC_Qtot 	= "SRTInIcePulses_NoDC_Qtot"
		SRTInIcePulses_NoDC 		= "SRTInIcePulses_NoDC"
	else:
		InIcePulses 				= "OfflinePulses"
		SRTInIcePulses 				= "SRTInIcePulses"
		SRTInIcePulses_NoDC_Qtot 	= "SRTInIcePulses_NoDC_Qtot"
		SRTInIcePulses_NoDC 		= "SRTInIcePulses_NoDC"
else:
	InIcePulses 					= "InIcePulses"
	SRTInIcePulses 					= "SRTInIcePulses"
	SRTInIcePulses_NoDC_Qtot 		= "SRTInIcePulses_NoDC_Qtot"
	SRTInIcePulses_NoDC 			= "SRTInIcePulses_NoDC"
#==================================================#

def make_conj(f1,f2):
	def conj(frame):
		c1=f1(frame)
		if(not c1):
			return(False)
		return(f2(frame))
	return(conj)
def make_disj(f1,f2):
	def disj(frame):
		c1=f1(frame)
		if(c1):
			return(True)
		return(f2(frame))
	return(disj)

class selector(object):
	def __init__(self,func):
		self.func=func
	def __call__(self,frame):
		return(self.func(frame))
	def __and__(self,other):
		return(selector(make_conj(self.func,other.func)))
	def __or__(self,other):
		return(selector(make_disj(self.func,other.func)))
	def __invert__(self):
		return(selector(lambda frame: not self.func(frame)))

def fixWeightMap(frame):
	if(not frame.Has("CorsikaWeightMap")):
		return
	cwm=frame["CorsikaWeightMap"]
	expectedKeys=["AreaSum","Atmosphere","CylinderLength","CylinderRadius","DiplopiaWeight",
	"EnergyPrimaryMax","EnergyPrimaryMin","FluxSum","FluxSum0","Multiplicity","NEvents",
	"ParticleType","PrimarySpectralIndex","SpectralIndexChange","SpectrumType","TimeScale","Weight"]

	replace=False
	nwm=cwm
	for expectedKey in expectedKeys:
		if not expectedKey in cwm.keys():
			replace=True
			nwm[expectedKey]=NaN
	if(replace):
		frame.Delete("CorsikaWeightMap")
		frame["CorsikaWeightMap"]=nwm

#the existance of an I3MCTree is usually sufficient to distinguish simulated data
@selector
def is_simulation(frame):
    if frame.Has("I3MCTree") or frame.Has("I3MCTree_preMuonProp"):
	return True

def renameMCTree(frame):
    if frame.Has("I3MCTree_preMuonProp") and not frame.Has("I3MCTree"):
        frame["I3MCTree"] = frame["I3MCTree_preMuonProp"]
        del frame["I3MCTree_preMuonProp"]


#def CheckIfs(frame):

#	print("checking ifs")
#	print(is_simulation(frame),splitFrames(frame),afterpulses(frame))

#Select for data passing the muon filter in either year, and abstract away differences in
#processing by recording the name of that year's likelihood fit in a string which will be
#used in later modules

# Ensure that the frame passes the muon filter
@selector
def isMuonFilter(frame):
	if not frame.Has('FilterMask'):
		if(options.diagnoseCuts): print id,"--No filtermask"
                return(False)
	filterMask = filterMask = frame['FilterMask']

	if filterMask.has_key('MuonFilter_10'):
		muonFilter=filterMask['MuonFilter_10']
		frame["PoleMPEFitName"]=dataclasses.I3String("MPEFit_SLC")
		return(muonFilter.condition_passed and muonFilter.prescale_passed)
	elif filterMask.has_key('MuonFilter_11'):
		muonFilter=filterMask['MuonFilter_11']
		frame["PoleMPEFitName"]=dataclasses.I3String("MPEFit")
		return(muonFilter.condition_passed and muonFilter.prescale_passed)
	elif filterMask.has_key('MuonFilter_12'):
		muonFilter=filterMask['MuonFilter_12']
		frame["PoleMPEFitName"]=dataclasses.I3String("MPEFit")
		return(muonFilter.condition_passed and muonFilter.prescale_passed)
	elif filterMask.has_key('MuonFilter_13'):
		muonFilter=filterMask['MuonFilter_13']
		frame["PoleMPEFitName"]=dataclasses.I3String("MPEFit")
		return(muonFilter.condition_passed and muonFilter.prescale_passed)
	else:
		print("There is a problem! If none of the above exists, you are in the future. Hihao.")
		return(False)

# In case of weirdness, select only for events which have the necessary pulses
@selector
def hasTWSRTOfflinePulses(frame):
	return(frame.Has(SRTInIcePulses))

def getRecoPulses(frame,name):
	#print(sorted(frame.keys()))
	pulses = frame[name]
	if pulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
		cal = frame['I3Calibration']
		#print(cal)
		pulses = pulses.apply(frame)
		#pulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame,name)
	return pulses

deepCoreStrings = range(79,87)

def dumbOMSelection(frame, pulses, output, omittedStrings, IfCond):
	if not IfCond:
		return
	if not frame.Has(pulses):
		return
	ps = getRecoPulses(frame,pulses)
	nps = dataclasses.I3RecoPulseSeriesMap()
	for om in ps.keys():
		if not om.string in omittedStrings:
			nps[om] = ps[om]
	frame.Put(output,nps)

# Compute the 'Charge Weighted Distance' as defined by J. Feintzeig: The average distance
# of the pulses, weighted by charge, from the track hypothesis.

def ComputeChargeWeightedDist(frame, Pulses, Track):
	if(not frame.Stop==icetray.I3Frame.Physics):
		return
	if(not frame.Has(Pulses)):
		return
	if(not frame.Has(Track)):
		return
	pulses=getRecoPulses(frame,Pulses)
	track=frame[Track]
	if(track.__class__==dataclasses.I3String):
		Track=track.value
		if(not frame.Has(Track)):
			return
		track=frame[Track]
	geo=frame.Get('I3Geometry')
	omgeo=geo.omgeo
	Qtot=0
	AvgDistQ=0
	for dom in pulses:
		DomPosition=omgeo[dom[0]].position
		Dist=phys_services.I3Calculator.closest_approach_distance(track,DomPosition)
		Qdom=0
		for pulse in dom[1]:
			Qdom+=pulse.charge
		Qtot+=Qdom
		AvgDistQ+=Dist*Qdom
	if(Qtot==0):
		AvgDistQ=NaN
	else:
		AvgDistQ/=Qtot
	if not frame.Has(Track+'_AvgDistQ'):
		frame.Put(Track+"_AvgDistQ",dataclasses.I3Double(AvgDistQ))
	if not frame.Has(Pulses+"_Qtot"):
		frame.Put(Pulses+"_Qtot",dataclasses.I3Double(Qtot))

# Perform a set of very simple quality cuts to reduce data rate:
# For events whose cos(reconstructed zenith)>.2, make a qtot cut similar to the on used in the filter
# For all up-going events, make a cut demanding that the charge-weighted distance be smaller than
# 200 meters or the total charge be greater than 100 p.e.

@selector
def precut(frame):
	PoleMPEFitName=frame['PoleMPEFitName'].value
	if(not frame.Has(PoleMPEFitName) or not frame.Has(PoleMPEFitName+"_AvgDistQ") or not frame.Has(SRTInIcePulses_NoDC_Qtot)):
		if(options.diagnoseCuts): print id,"-- Missing PoleMPEFitName or SRTInIcePulses_NoDC_Qtot"
		return(False)
	
	track = frame.Get(PoleMPEFitName)
	qtot = frame.Get(SRTInIcePulses_NoDC_Qtot).value
	qwDist = frame.Get(PoleMPEFitName+"_AvgDistQ").value

	if(track.dir.zenith > (math.pi / 2.)):
		if(qtot < 100. and qwDist > 200.):
			if(options.diagnoseCuts): print id,"--Particle above horizon, total charge is too spall and AvgQdist is too large"
			return(False)
		return(True)
	else:
		c = math.cos(track.dir.zenith)
		lq = 0
		if(qtot > 0.0):
			lq = math.log10(qtot)
		if(c > .2):
			return(lq >= (0.6*(c-0.5) + 2.5))
		return(True)

@selector
def afterpulses(frame):
	return(frame.Has("IsAfterPulses"))

@selector
def splitFrames(frame):
	if(frame.Stop!=icetray.I3Frame.Physics):
		return(True)
	return (frame["I3EventHeader"].sub_event_stream=="TTrigger")

@selector
def basicRecosAlreadyDone(frame):
	if(frame.Stop!=icetray.I3Frame.Physics):
		return(True)
	return(frame.Has("MPEFit_TT"))

#In a fair number of events TopologicalSplitting does nothing, so we can save
#computation by just copying the results from the L2a reconstructions
class SamePulseChecker(icetray.I3PacketModule):
	copyKeys = {"MPEFit_SLC":"MPEFit_TT","MPEFit_SLCFitParams":"MPEFit_TTFitParams",
		"SPEFit4_SLC":"SPEFit4_TT","SPEFit4_SLCFitParams":"SPEFit4_TTFitParams",
		"SPEFitSingle_SLC":"SPEFitSingle_TT","SPEFitSingle_SLCFitParams":"SPEFitSingle_TTFitParams",
		"LineFit_SLC":"LineFit_TT"}

	def __init__(self, ctx):
		icetray.I3PacketModule.__init__(self, ctx, icetray.I3Frame.DAQ)
		self.oldPulseName = SRTInIcePulses
		self.newPulseName = "TTPulses"
		self.AddOutBox("OutBox")

	def Configure(self):
		pass

def add_basic_reconstructions(tray,suffix,pulses,condition):
	tray.AddSegment(linefit.simple,"LineFit"+suffix,
		If=condition,
		inputResponse=pulses,
		fitName="LineFit"+suffix
		)

	tray.AddSegment(icecube.lilliput.segments.I3SinglePandelFitter,"SPEFitSingle"+suffix,
		If=condition,
		domllh="SPE1st",
		pulses=pulses,
		seeds=["LineFit"+suffix],
		)

	tray.AddSegment(icecube.lilliput.segments.I3IterativePandelFitter,"SPEFit4"+suffix,
		If=condition,
		domllh="SPE1st",
		n_iterations=4,
		pulses=pulses,
		seeds=["SPEFitSingle"+suffix],
		)

	tray.AddSegment(icecube.lilliput.segments.I3SinglePandelFitter,"MPEFit"+suffix,
		If=condition,
		domllh="MPE",
		pulses=pulses,
		seeds=["SPEFit4"+suffix],
		)

	tray.AddModule(selectFinalFit,"FinalFit")

#Figure out which is the highest level track fit to succeed, and bless it as our track hypothesis
#Also record how many of the fits in the sequence failed
def selectFinalFit(frame):
	fits=["MPEFit_TT","SPEFit4_TT","SPEFitSingle_TT","LineFit_TT"]
	resultName="TrackFit"
	result=None
	params=None
	for fitName in fits:
		if(not frame.Has(fitName)):
			continue
		fit=frame.Get(fitName)
		if(fit.fit_status==dataclasses.I3Particle.OK):
			frame.Put(resultName,fit)
			frame.Put(resultName+"Fallback",icetray.I3Int(fits.index(fitName)))
			#print(fitName)
			if(frame.Has(fitName+"FitParams")):
				frame.Put(resultName+"FitParams",frame.Get(fitName+"FitParams"))
			break
	if(not frame.Has(resultName+"Fallback")):
		frame.Put(resultName+"Fallback",icetray.I3Int(len(fits)))

def firstHits(frame,inputPulses,outputPulses):
	# This function makes a mask of only the first pulses seen on every DOM.
	# see http://software.icecube.wisc.edu/offline_trunk/projects/dataclasses/masks.html
	if(frame.Has(inputPulses)):
		frame.Put(outputPulses,dataclasses.I3RecoPulseSeriesMapMask(frame,inputPulses,lambda OM,id,i3P: id==0))
	else:
		return True

def get_NChan(frame):
	# This is needed for Paraboloid sigma corrector.
	if frame.Has('TTPulses_hm'):
		pulse_hm = frame.Get('TTPulses_hm')
		#print(pulse_hm.n_hit_doms)
		frame.Put("NChanSource",icetray.I3Int(pulse_hm.n_hit_doms))
	else:
		return

def computeSimpleCutVars(tray,condition):
	tray.Add(firstHits,
		inputPulses = "TTPulses",
		outputPulses = "TTPulses_first")

	#Save the direct hit informaiton from the CommomVariables project
	tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'dh',
		DirectHitsDefinitionSeries       	= dh_definitions,
		PulseSeriesMapName              	= "TTPulses_first",
		ParticleName                     	= "TrackFit",
		OutputI3DirectHitsValuesBaseName 	= "TrackFit_dh",
		BookIt                           	= True,
		If 									= condition)

	# Save the hit multiplicity information as well.
	tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, '_hm',
		PulseSeriesMapName                	= "TTPulses_first",
		OutputI3HitMultiplicityValuesName 	= 'TTPulses_hm',
		BookIt                            	= True,
		If 									= condition)

	tray.Add(get_NChan)

	#=============Now do it without DeepCore================#

	# Select all the TTPulses that were not part of DeepCore
	tray.AddModule(dumbOMSelection,"NoDeepCore2",
		IfCond=condition,
		pulses="TTPulses",
		output="TTPulses_NoDC",
		omittedStrings=deepCoreStrings)

	tray.AddModule(firstHits,"FirstHits2",
		inputPulses = "TTPulses_NoDC",
		outputPulses = "TTPulses_NoDC_first")

	#Save the direct hit informaiton from teh CommomVariables project
	tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'dh_NoDC',
		DirectHitsDefinitionSeries       	= dh_definitions,
		PulseSeriesMapName               	= "TTPulses_NoDC_first",
		ParticleName                     	= "TrackFit",
		OutputI3DirectHitsValuesBaseName 	= "TrackFit_NoDC_dh",
		BookIt                           	= True,
		If 									= condition)

	# Save the hit multiplicity information as well.
	tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, 'hm_NoDC',
		PulseSeriesMapName                	= "TTPulses_NoDC_first",
		OutputI3HitMultiplicityValuesName 	= "TTPulses_NoDC_hm",
		BookIt                           	= True,
		If 									= condition)
	
	tray.AddModule(ComputeChargeWeightedDist,"CCWD2",Pulses="TTPulses_NoDC",Track="TrackFit")

#Only select events to calculate more complicated cut variables if they pass some very basic cuts
@selector
def doExpensiveRecos(frame):
	if(frame.Stop!=icetray.I3Frame.Physics):
		return(True)
	
	fallback=frame.Get("TrackFitFallback")
	if(fallback.value > 2):
		return(False)

	track = frame.Get("TrackFit")
	if(math.cos(track.dir.zenith) > 0.2):
		return(False)

	track_dh = frame.Get("TrackFit_dh")
	pulse_hm = frame.Get("TTPulses_hm")

	if(float(pulse_hm.n_hit_doms) < 15.0):
		#if(options.cut == 'True'):
		return(False)
	if(float(track_dh.n_dir_doms) <= 6.0):
		#if(options.cut == 'True'):
		return(False)
	if(float(track_dh.dir_track_length) < 200.0):
		#if(options.cut == 'True'):
		return(False)
	if(abs(float(track_dh.dir_track_hit_distribution_smoothness)) > 0.6):
		#if(options.cut == 'True'):
		return(False)

	return(True)

zenithPrior=None
def add_bayesian_reconstruction(tray,pulses,condition,seedTrack):
	seedLikelihood=seedTrack+"_BayesianLikelihoodSeed"
	tray.AddService('I3GulliverIPDFPandelFactory', seedLikelihood)(
		("InputReadout",pulses),
		("EventType","InfiniteMuon"),
		("Likelihood","SPE1st"),
		("PEProb","GaussConvolutedFastApproximation"),
		("JitterTime",15.*I3Units.ns),
		("NoiseProbability",10.*I3Units.hertz)
		)

	global zenithPrior
	if(zenithPrior==None):
		zenithPrior="ZenithPrior"
		tray.AddService("I3PowExpZenithWeightServiceFactory",zenithPrior)(
			("Amplitude",2.49655e-07),
			("CosZenithRange",[-1,1]),
			("DefaultWeight",1.383896526736738e-87),
			("ExponentFactor",0.778393),
			("FlipTrack",False),
			("PenaltySlope",-1000),
			("PenaltyValue",-200),
			("Power",1.67721)
			)

	tray.AddService("I3BasicSeedServiceFactory", seedTrack+"_BayesianSeed")(
		("InputReadout",pulses),
		("FirstGuesses",[seedTrack]),
		("TimeShiftType","TFirst")
		)

	tray.AddService("I3EventLogLikelihoodCombinerFactory",seedTrack+"_BayesianLikelihood")(
		("InputLogLikelihoods",[seedLikelihood,zenithPrior])
		)

	tray.AddModule("I3IterativeFitter",seedTrack+"Bayesian")(
		("OutputName",seedTrack+"Bayesian"),
		("RandomService","SOBOL"),
		("NIterations",8),
		("SeedService",seedTrack+"_BayesianSeed"),
		("Parametrization","default_simpletrack"),
		("LogLikelihood",seedTrack+"_BayesianLikelihood"),
		("CosZenithRange",[0, 1]),
		("Minimizer","default_simplex"),
		("If",condition)
		)

def add_paraboloid(tray,pulses,condition,seed):
	tray.AddService("I3BasicSeedServiceFactory",seed+"_ParaboloidSeed")(
		("InputReadout",pulses),
		("FirstGuesses",[seed]),
		("TimeShiftType","TFirst")
		)

	tray.AddService("I3GulliverIPDFPandelFactory",seed+"_ParaboloidLikelihood")(
		("InputReadout",pulses),
		("Likelihood","MPE"),
		("PEProb","GaussConvolutedFastApproximation"),
		("NoiseProbability",10 * I3Units.hertz),
		("JitterTime",4.0 * I3Units.ns),
		("EventType","InfiniteMuon")
		)

	tray.AddModule("I3ParaboloidFitter", seed+"Paraboloid")(
		("OutputName",seed+"Paraboloid"),
		("If",condition),
		("SeedService",seed+"_ParaboloidSeed"),
		("LogLikelihood",seed+"_ParaboloidLikelihood"),
		("MaxMissingGridPoints",1),
		("VertexStepSize",5.0 * I3Units.m),
		("ZenithReach",2.0 * I3Units.degree),
		("AzimuthReach",2.0 * I3Units.degree),
		("GridpointVertexCorrection",seed+"_ParaboloidSeed"),
		("Minimizer","default_simplex"),
		("NumberOfSamplingPoints",8),
		("NumberOfSteps",3),
		("MCTruthName","")
		)


def add_split_fits(tray,splitName,idx,pulses,condition):
	tray.AddSegment(linefit.simple,"LineFit"+splitName+idx,
		If=condition,
		inputResponse=pulses+splitName+idx,
		fitName="LineFit"+splitName+idx
		)

	tray.AddSegment(icecube.lilliput.segments.I3IterativePandelFitter,"SPEFit4"+splitName+idx,
		If=condition,
		domllh="SPE1st",
		n_iterations=4,
		pulses=pulses+splitName+idx,
		seeds=["LineFit"+splitName+idx]
		)

	add_bayesian_reconstruction(tray,pulses,condition,"SPEFit4"+splitName+idx)

	# The firstHit function is a mask that grabs only the first pulses on each DOM
	tray.Add(firstHits,
		inputPulses = pulses+splitName+idx,
		outputPulses = pulses+splitName+idx+"_first"
		)

	# Save the direct hit informaiton from the CommomVariables project
	tray.AddSegment(direct_hits.I3DirectHitsCalculatorSegment, 'dh_'+splitName+idx,
		DirectHitsDefinitionSeries       	= dh_definitions,
		PulseSeriesMapName               	= pulses+splitName+idx+"_first",
		ParticleName                     	= "SPEFit4"+splitName+idx,
		OutputI3DirectHitsValuesBaseName 	= "SPEFit4"+splitName+idx+'_dh',
		BookIt                           	= True,
		If 									= condition
		)
		
	# Save the hit multiplicity information as well.
	tray.AddSegment(hit_multiplicity.I3HitMultiplicityCalculatorSegment, 'hm_'+splitName+idx,
		PulseSeriesMapName                	= pulses+splitName+idx+"_first",
		OutputI3HitMultiplicityValuesName 	= pulses+splitName+idx+"_hm",
		BookIt                            	= True,
		If 									= condition
		)

	tray.AddModule("Delete",splitName+idx+"Cleanup")(
		("If",condition),
		("Keys",[pulses+splitName+idx+"_first"])
		)

#Split the pulses both geometrically and temporally, then reconstructions on all four resulting halves
def add_split_reconstructions(tray,pulses,condition,seed):
	tray.AddModule("I3ResponseMapSplitter",pulses+"SplitTime")(
		("If",condition),
		("InputPulseMap",pulses),
		("InputTrackName","") # deliberately blank
		)

	tray.AddModule("I3ResponseMapSplitter",pulses+"SplitGeo")(
		("If",condition),
		("InputPulseMap",pulses),
		("InputTrackName",seed)
		)

	add_split_fits(tray,"SplitTime","1",pulses,condition)
	add_split_fits(tray,"SplitTime","2",pulses,condition)
	add_split_fits(tray,"SplitGeo","1",pulses,condition)
	add_split_fits(tray,"SplitGeo","2",pulses,condition)

#======================================
gcdfile 	= options.gcdfile
print('CHANGIN GCD FILE!!!!!!!!!!!!!!!!!!!!!!!!!!')
print('I Had to use an older version of GCD file to be compatible with the metaproject.')
print('XLevel does not care about the GCD file containing the SPE templates.')
print('I cannot update the metraproject, since the version of MuEX is outdated.')
gcdfile = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/Noise/GeoCalibDetectorStatus_AVG_Fit_55697-57531_SPE_PASS2_Raw.i3.gz'

infile 		= options.infile
outfile 	= options.outfile
osg 		= options.osg

if osg == 'True':
    def copy_to_OSG(NPX_file):
        subprocess.check_call(['globus-url-copy','-nodcau','-rst', NPX_file,"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1])])

    def copy_to_NPX(NPX_file):
        subprocess.check_call(['globus-url-copy','-nodcau','-rst',"file:"+os.getcwd()+'/'+str(NPX_file.split('/')[-1]), NPX_file])

    gcdfile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + gcdfile)
    gcdfile = str(os.getcwd()+'/'+gcdfile.split('/')[-1])

    infile_NPX = str('gsiftp://gridftp.icecube.wisc.edu' + infile)
    infile = str(os.getcwd()+'/'+infile.split('/')[-1])

    spline_NPX = "gsiftp://gridftp.icecube.wisc.edu/data/ana/SterileNeutrino/IC86/HighEnergy/scripts/jobs/paraboloidCorrectionSpline.dat"

    copy_to_OSG(gcdfile_NPX)
    copy_to_OSG(infile_NPX)
    copy_to_OSG(spline_NPX)

    infiles = [gcdfile,infile]
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

    #outfile_temp = str(os.getcwd()+'/'+outfile.split('/')[-1])

    spline_path =  str(os.getcwd()+'/'+str(spline_NPX).split('/')[-1])
else:

	infiles = [gcdfile, infile]
	outfile_temp = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/temp/'+outfile.split('/')[-1]
	outfile = options.outfile
	spline_path = "/data/ana/SterileNeutrino/IC86/HighEnergy/scripts/jobs/paraboloidCorrectionSpline.dat"

#====================================
tray = I3Tray()
tray.AddService("I3GSLRandomServiceFactory","Random") # needed for bootstrapping
tray.AddModule("I3Reader","reader")(
		("FilenameList",infiles)
		)

exitStatus=RandomStuff.ExitStatus()

tray.AddModule(renameMCTree, "_renameMCTree", Streams=[icetray.I3Frame.DAQ])

tray.AddModule(isMuonFilter & hasTWSRTOfflinePulses,"selectValidData")
tray.AddModule(fixWeightMap,"patchCorsikaWeights")

truecondition=lambda frame: True

#======================================
# Compute precut variables and do precut
tray.AddModule(dumbOMSelection,"NoDeepCore",
	pulses			= SRTInIcePulses,
	output			= SRTInIcePulses_NoDC,
	omittedStrings	= deepCoreStrings,
	IfCond			= truecondition
	)

tray.AddModule(ComputeChargeWeightedDist,"CCWD",Pulses=SRTInIcePulses_NoDC,Track="PoleMPEFitName")

tray.AddModule(precut,"precut")
tray.AddModule("I3OrphanQDropper","OrphanQDropper")
#======================================
# Run TTrigger and do reconstructions

stConfigService = STTools.seededRT.configuration_services.I3DOMLinkSeededRTConfigurationService(
	allowSelfCoincidence    = False,            # Default: False.
	useDustlayerCorrection  = True,             # Default: True.
	dustlayerUpperZBoundary =    0*I3Units.m,   # Default: 0m.
	dustlayerLowerZBoundary = -150*I3Units.m,   # Default: -150m.
	ic_ic_RTTime            = 1000*I3Units.ns,  # Default: 1000m.
	ic_ic_RTRadius          =  150*I3Units.m    # Default: 150m.
	)

tray.AddModule("I3SeededRTCleaning_RecoPulse_Module", "SRTClean",
	InputHitSeriesMapName  = InIcePulses,
	OutputHitSeriesMapName = SRTInIcePulses,
	STConfigService        = stConfigService,
	#SeedProcedure         = "HLCCoreHits",
	NHitsThreshold         = 2,
	Streams                = [icetray.I3Frame.DAQ]
	)

tray.AddModule("I3TopologicalSplitter","TTrigger",
	SubEventStreamName 		= "TTrigger", #Spencer -- Is this ok?? I was getting warnigngs...
	InputName				= SRTInIcePulses,
	OutputName				= "TTPulses",
	Multiplicity			= 4,
	TimeWindow				= 4000*I3Units.ns,
	TimeCone				= 800*I3Units.ns,
	SaveSplitCount			= True
	)

tray.AddModule("AfterPulseSpotter","Afterpulses")(
	("StreamName","TTrigger"),
	("Pulses","TTPulses")
	)

tray.AddModule(SamePulseChecker,"SPC")

add_basic_reconstructions(tray,"_TT","TTPulses",splitFrames & ~afterpulses & ~basicRecosAlreadyDone)

#======================================
# Compute some simple cut variables
computeSimpleCutVars(tray,splitFrames & ~afterpulses)

#======================================
# Compute some expensive cut variables
add_bayesian_reconstruction(tray,"TTPulses",splitFrames & ~afterpulses & doExpensiveRecos,"TrackFit")
add_paraboloid(tray,"TTPulses",splitFrames & ~afterpulses & doExpensiveRecos,"TrackFit")
add_split_reconstructions(tray,"TTPulses",splitFrames & ~afterpulses & doExpensiveRecos,"TrackFit")

#======================================
# Write output

#def ParaboloidSigmaCorrector(frame):
#	NChan = frame.Get('TTPulses_hm').n_hit_doms
#	params = frame.Get('TrackFitParaboloidFitParams')
#	err1 = params.pbfErr1
#	err2 = params.pbfErr2
#	sigma = np.sqrt((err1*err1+err2*err2)/2.);
#	sigma *= correctionSpline(log10(*nChan))
#	frame.Put("CorrectedParaboloidSigma",dataclasses.I3Double(sigma))

tray.AddModule("ParaboloidSigmaCorrector","extractSigma")(
		("CorrectionSplinePath",spline_path),
		("ParaboloidParameters","TrackFitParaboloidFitParams"),
		("NChanSource",'NChanSource'), # with DeepCore!
		("Output","CorrectedParaboloidSigma"),)

@selector
def keepL3(frame):
	if(frame.Stop!=icetray.I3Frame.Physics):
		return(True)
	id = frame["I3EventHeader"].run_id,frame["I3EventHeader"].event_id,frame["I3EventHeader"].sub_event_id

	# We require that we had a good enough fit. The fallback.value represents what level of fit was sucessfull.
	fallback = frame.Get("TrackFitFallback")
	if(fallback.value > 2):
		if(options.diagnoseCuts): print id,"Fall back value too high"
		return(False)

	# Let's select some initial cuts.
	# dh: direct hit information
	# hm: hit multiplicity
	# Note: we typically cut while ignoring DeepCore.

	track 			= frame.Get("TrackFit")
	track_NoDC_dh 	= frame.Get("TrackFit_NoDC_dh")
	pulse_NoDC_hm 	= frame.Get("TTPulses_NoDC_hm")

	#frame.Put("trackfit_cuts_noDC_nchan",dataclasses.I3Double(frame.Get("TrackFitCuts_NoDC").nchan))
	#cuts = frame.Get("TrackFitCuts_NoDC")
	#frame.Put("trackfit_cuts_noDC_nchan",dataclasses.I3Double(frame.Get("TTPulses_NoDC_hm").n_hit_doms))
	#frame.Put("trackfit_cuts_noDC_nchan",dataclasses.I3Double(frame.Get("TTPulses_NoDC_hm").n_hit_doms))

	if(math.cos(track.dir.zenith) > 0.1):
		if(options.cut == 'True'):
			return(False)
	if(float(pulse_NoDC_hm.n_hit_doms) < 30.):
		if(options.cut == 'True'):
			return(False)
	if(math.cos(track.dir.zenith) > 0.0 and float(pulse_NoDC_hm.n_hit_doms) < 30.0 * math.exp(17.0*math.cos(track.dir.zenith))):
		if(options.cut == 'True'):
			return(False)

	#frame.Put("trackfit_cuts_nchan",dataclasses.I3Double(frame.Get("TrackFitCuts").nchan))

	if(not doExpensiveRecos(frame)):
		return(False)
	sigma = frame.Get("CorrectedParaboloidSigma")
	llhparams = frame.Get("TrackFitFitParams")
	frame.Put("RLogL",dataclasses.I3Double(llhparams.rlogl))
	frame.Put("Sigma",dataclasses.I3Double(sigma.value))
	if(sigma.value<2.5679769135e-2 and (llhparams.rlogl>(-42.967631051274708*sigma.value+8.6))):
		if(options.diagnoseCuts): print id,"fails due to paraboloid and rlogl (1):",sigma.value,llhparams.rlogl
		if(options.cut == 'True'):
			return(False)
	if(sigma.value>2.5679769135e-2 and (llhparams.rlogl>(-5.0*sigma.value+7.625))):
		if(options.diagnoseCuts): print id,"fails due to paraboloid and rlogl (2):",sigma.value,llhparams.rlogl
		if(options.cut == 'True'):
			return(False);
	bayes = frame.Get("TrackFitBayesian")
	bayesparams = frame.Get("TrackFitBayesianFitParams")
	speparams = frame.Get("SPEFit4_TTFitParams")

	if(bayes.fit_status==dataclasses.I3Particle.OK):
		frame.Put("bayes_fit_status",dataclasses.I3Double(1))
		llhdiff=bayesparams.logl-speparams.logl
		frame.Put("BayesLLHR",dataclasses.I3Double(llhdiff))
		if(math.cos(track.dir.zenith)<0.0 and llhdiff<33):
			if(options.diagnoseCuts): print id,"fails due to bayes diff (1):",llhdiff
			if(options.cut == 'True'):
				return(False)
		if(math.cos(track.dir.zenith)>=0.0):
			if(llhdiff<(33-86*math.cos(track.dir.zenith))):
				if(options.diagnoseCuts): print id,"fails due to bayes diff (2)",math.cos(track.dir.zenith),llhdiff
				if(options.cut == 'True'):
					return(False)
			if(llhdiff>(75.-45*math.sqrt(1-math.pow((math.cos(track.dir.zenith)-.1)/.1,2)))):
				if(options.diagnoseCuts): print id,"fails due to bayes diff (3)",math.cos(track.dir.zenith),llhdiff
				if(options.cut == 'True'):
					return(False)
	else:
		#print("bayes_fit_status FAILED")
		frame.Put("bayes_fit_status",dataclasses.I3Double(0))
		frame.Put("BayesLLHR",dataclasses.I3Double(float('NaN')))

	track_dh = frame.Get("TrackFit_dh")

	geosplit1 = frame.Get("SPEFit4SplitGeo1")
	geosplit2 = frame.Get("SPEFit4SplitGeo2")
	if(geosplit1.fit_status==dataclasses.I3Particle.OK and geosplit2.fit_status==dataclasses.I3Particle.OK):
		geocuts1 = frame.Get("SPEFit4SplitGeo1_dh")
		geocuts2 = frame.Get("SPEFit4SplitGeo2_dh")
		if((geocuts1.n_dir_doms + geocuts2.n_dir_doms) > 2 * track_dh.n_dir_doms):
			return(False)
		else:
			pass
	
	timesplit1 = frame.Get("SPEFit4SplitTime1")
	timesplit2 = frame.Get("SPEFit4SplitTime2")
	if(timesplit1.fit_status == dataclasses.I3Particle.OK and timesplit2.fit_status == dataclasses.I3Particle.OK):
		timecuts1 = frame.Get("SPEFit4SplitTime1_dh")
		timecuts2 = frame.Get("SPEFit4SplitTime2_dh")
		if((timecuts1.n_dir_doms + timecuts2.n_dir_doms) > 2 * track_dh.n_dir_doms):
			return(False)
		else:
			pass
	'''
	if geosplit1.fit_status==dataclasses.I3Particle.OK:
		frame.Put("geo_split1_fit_status",dataclasses.I3Double(1))
		frame.Put("geo_split1_dir_n",dataclasses.I3Double(frame.Get("SPEFit4SplitGeo1Cuts").ndir))
	else:
		#print("geosplit1_status FAILED")
		frame.Put("geo_split1_fit_status",dataclasses.I3Double(0))
		frame.Put("geo_split1_dir_n",dataclasses.I3Double(float('NaN')))

		if geosplit2.fit_status==dataclasses.I3Particle.OK:
				frame.Put("geo_split2_fit_status",dataclasses.I3Double(1))
		frame.Put("geo_split2_dir_n",dataclasses.I3Double(frame.Get("SPEFit4SplitGeo2Cuts").ndir))
		else:
		#print("geosplit2_status FAILED")
				frame.Put("geo_split2_fit_status",dataclasses.I3Double(0))
		frame.Put("geo_split2_dir_n",dataclasses.I3Double(float('NaN')))

	if timesplit1.fit_status==dataclasses.I3Particle.OK:
		frame.Put("time_split1_fit_status",dataclasses.I3Double(1))
		frame.Put("time_split1_dir_n",dataclasses.I3Double(frame.Get("SPEFit4SplitTime1Cuts").ndir))
	else:
		#print("timesplit1_status FAILED")
		frame.Put("time_split1_fit_status",dataclasses.I3Double(0))
		frame.Put("time_split1_dir_n",dataclasses.I3Double(float('NaN')))

		if timesplit2.fit_status==dataclasses.I3Particle.OK:
				frame.Put("time_split2_fit_status",dataclasses.I3Double(1))
			frame.Put("time_split2_dir_n",dataclasses.I3Double(frame.Get("SPEFit4SplitTime2Cuts").ndir))
	else:
		#print("timesplit2_status FAILED")
				frame.Put("time_split2_fit_status",dataclasses.I3Double(0))
		frame.Put("time_split2_dir_n",dataclasses.I3Double(float('NaN')))
	'''
	frame.Put("TTPulses_NChan",dataclasses.I3Double(frame.Get("TTPulses_hm").n_hit_doms))
	frame.Put("TTPulses_NoDC_NChan",dataclasses.I3Double(frame.Get("TTPulses_NoDC_hm").n_hit_doms))
	frame.Put("Dir_N",dataclasses.I3Double(frame.Get("TrackFit_dh").n_dir_doms))
	frame.Put("Dir_S",dataclasses.I3Double(frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness))
	frame.Put("Dir_L",dataclasses.I3Double(frame.Get("TrackFit_dh").dir_track_length))
	frame.Put("Dir_N_NoDC",dataclasses.I3Double(frame.Get("TrackFit_NoDC_dh").n_dir_doms))
	frame.Put("Dir_S_NoDC",dataclasses.I3Double(frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness))
	frame.Put("Dir_L_NoDC",dataclasses.I3Double(frame.Get("TrackFit_NoDC_dh").dir_track_length))
	return(True)

def planL3Cut(frame):
	if((splitFrames & (afterpulses | ~keepL3))(frame)):
		frame.Put("CutL3",icetray.I3Bool(True))

tray.AddModule(planL3Cut,"CutPlanner")

#=======================================================
# Inserting the Energy section
#=======================================================

@selector
def cutL3(frame):
	return(frame.Has("CutL3"))

@selector
def parabaloidCut(frame):
	if(frame["TrackFitParaboloidFitParams"].pbfStatus!=dataclasses.I3Particle.OK): return False
	return(True)

@selector
def afterpulses(frame):
	return(frame.Has("IsAfterPulses"))

@selector
def goodFit(frame):
	if(not frame.Has("TrackFit")):
		return(False)
	if(frame["TrackFit"].fit_status!=dataclasses.I3Particle.OK):
		return(False)
	if(frame.Has("TrackFitFallback") and frame["TrackFitFallback"].value > 2):
		return(False)
	return(True)

finalSample = splitFrames & ~cutL3 & ~afterpulses & goodFit

#=========================DO I NEED THIS???????????==============================
tray.AddModule("I3MCParticleExtractor","extractMCTruth")(
		("StoreLevel",2)
		)
#================================================================================

#tray.AddModule("I3WaveformTimeRangeCalculator","WaveformRange")
#tray.AddModule('I3WaveformTimeRangeCalculator', 'SuperDSTRange', If=lambda frame: not 'CalibratedWaveformRange' in frame)
#tray.AddModule(wavedeform.AddMissingTimeWindow,"PulseTimeRange")

tray.AddModule("I3WaveformTimeRangeCalculator","WaveformRange", If=lambda frame: not 'CalibratedWaveformRange' in frame)
tray.AddModule(wavedeform.AddMissingTimeWindow,"PulseTimeRange")

def findHighChargeDOMs(frame, pulses, outputList):
	pulsemap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, pulses)
	charges = np.array([sum([p.charge for p in pulses]) for pulses in pulsemap.itervalues()])
	qmean = charges.mean()
	baddies = [om for om, q in zip(pulsemap.keys(), charges) if q/qmean > 10.]
	frame[outputList] = dataclasses.I3VectorOMKey(baddies)

tray.AddModule(findHighChargeDOMs, "findHighChargeDOMs",
	pulses = InIcePulses,
	outputList = "BrightDOMs")

tray.AddModule("muex", "muex")(
	("pulses", "TTPulses"),#MuEx does not handle noise, so it should run on cleaned pulses
	("rectrk", "TrackFit"),
	("result", "MuEx"),
	("lcspan", 0),
	("repeat", 0),
	("rectyp", True),
	("usempe", False),
	("detail", False),
	("energy", True),
    #("compat", False),#New Addition bug fix
	("icedir", ice_model_location),
	("If",finalSample))

def add_energy(frame):
	if frame.Has('I3MCTree'):
		frame.Put("",dataclasses.get_most_energetic_muon(frame['I3MCTree']).energy)

def true_trackfit(frame):
	if frame.Has('I3MCTree'):
		# We want the true trajectory of the muon for this. But if there is no muon (ie NuGen Nue simulation), then take the cascade
		try:
			track = dataclasses.I3Particle(dataclasses.get_most_energetic_muon(frame['I3MCTree']))
		except:
			track = dataclasses.I3Particle(dataclasses.get_most_energetic_cascade(frame['I3MCTree']))
		frame.Put('true_trackfit',track)
	else:
		return
tray.Add(true_trackfit)

'''
tray.AddModule("muex", "muextrue")(
	("pulses", "TTPulses"),
	("rectrk", "true_trackfit"),
	("result", "MuExTrue"),
	("lcspan", 0),
	("repeat", 0),
	("rectyp", True),
	("usempe", False),
	("detail", False),
	("energy", True),
    #("compat",False),# New addition, bug fix
	#("icedir", expandvars("$I3_BUILD/mue/resources/ice/mie")),
	("icedir", ice_model_location),
	("If",finalSample & is_simulation))
'''
tray.AddModule("FiducialVolumeEntryPointFinder", "findEntryPoints",
	TopBoundaryWidth=0.,
	BottomBoundaryWidth=0.,
	SideBoundaryWidth=0.,
	AssumeMuonTracksAreInfinite=True,
	ParticleIntersectionsName="IntersectingParticles",
	If=is_simulation)

def _FindDetectorVolumeIntersections(recoParticle, geometry):
	intersectionPoints = VHESelfVeto.IntersectionsWithInstrumentedVolume(geometry, recoParticle)
	intersectionTimes = []
	for intersectionPoint in intersectionPoints:
		vecX = intersectionPoint.x - recoParticle.pos.x
		vecY = intersectionPoint.y - recoParticle.pos.y
		vecZ = intersectionPoint.z - recoParticle.pos.z

		prod = vecX*recoParticle.dir.x + vecY*recoParticle.dir.y + vecZ*recoParticle.dir.z
		dist = math.sqrt(vecX**2 + vecY**2 + vecZ**2)
		if prod < 0.: dist *= -1.

		if abs(prod-dist) > 1e-3*icetray.I3Units.m:
			raise RuntimeError("intersection points are not on track")

		intersectionTimes.append(dist/dataclasses.I3Constants.c + recoParticle.time)

	sortedTimes = sorted(intersectionTimes)
	return sortedTimes

def FindDetectorVolumeIntersections(frame, TrackName="", OutputTimeWindow=None, TimePadding=0.):
	if not (finalSample & is_simulation)(frame):
		return
	if OutputTimeWindow is not None:
		twName = OutputTimeWindow
	else:
		twName = TrackName + "TimeRange"

	theTrack = frame[TrackName]
	geometry = frame["I3Geometry"]

	times = _FindDetectorVolumeIntersections(theTrack, geometry)

	if len(times) == 0:
		#raise RuntimeError("track does not intersect the detector volume")
		frame[twName] = dataclasses.I3TimeWindow()
	elif len(times)	== 1:
		raise RuntimeError("tracks with only one intersection are not supported")
	else:
		tWindow = dataclasses.I3TimeWindow(times[0]-TimePadding, times[-1]+TimePadding)
		frame[twName] = tWindow

tray.AddModule(FindDetectorVolumeIntersections, "FindDetectorVolumeIntersections",
	TimePadding = 60.*I3Units.m/dataclasses.I3Constants.c,
	TrackName="TrackFit",
	OutputTimeWindow="ContainedTimeRange")

def EntryEnergy(frame):
	if not (finalSample & is_simulation)(frame):
		return
	mcTree = frame["I3MCTree"]
	#locate the highest energy muon
	muon = None
	for p in mcTree:
		if p.type in [dataclasses.I3Particle.ParticleType.MuPlus,
				dataclasses.I3Particle.ParticleType.MuMinus]:
			if(muon==None):
				muon = p
			elif(p.energy>muon.energy):
				muon = p
			#break

	#print "Muon energy is",muon.energy
	#sum all losses before the muon enters the detector
	timeWindow = frame["ContainedTimeRange"]
	#print "Contained time window is [",timeWindow.start,',',timeWindow.stop,']'
	containedEnergyMC = 0.
	earlyLosses = 0.
	for p in frame["I3MCTree"]:
		if p.shape == dataclasses.I3Particle.ParticleShape.Dark:
			continue
		if p.type in [dataclasses.I3Particle.ParticleType.MuPlus,
				dataclasses.I3Particle.ParticleType.MuMinus,
				dataclasses.I3Particle.ParticleType.TauPlus,
				dataclasses.I3Particle.ParticleType.TauMinus,
				dataclasses.I3Particle.ParticleType.NuE,
				dataclasses.I3Particle.ParticleType.NuEBar,
				dataclasses.I3Particle.ParticleType.NuMu,
				dataclasses.I3Particle.ParticleType.NuMuBar,
				dataclasses.I3Particle.ParticleType.NuTau,
				dataclasses.I3Particle.ParticleType.NuTauBar,
				]:
			continue # skip tracks for now
		if p.location_type != dataclasses.I3Particle.LocationType.InIce:
			continue
		if p.time > timeWindow.stop:
			continue
		if p.time < timeWindow.start:
			earlyLosses += p.energy
			continue
		containedEnergyMC += p.energy
		if muon is None:
			print('No muon, could be NC?')
			if not frame.Has("MuonEntryEnergy"):
				frame["MuonEntryEnergy"]=dataclasses.I3Double(NaN)
	'''
	else:
		frame["MuonEntryEnergy"]=dataclasses.I3Double(muon.energy-earlyLosses)
	frame["MuonEnergyLoss"]=dataclasses.I3Double(containedEnergyMC)
	lengthInDetector=dataclasses.I3Constants.c*(timeWindow.stop-timeWindow.start)/icetray.I3Units.m
	frame["MuonLengthInDetector"]=dataclasses.I3Double(lengthInDetector)
	frame["MuonAverage_dEdX"]=dataclasses.I3Double(containedEnergyMC/lengthInDetector)
	'''

#tray.AddModule("Dump","d")
tray.AddModule(EntryEnergy,"EntryEnergy")

def isNeutrinoType(pType):
	return(pType==dataclasses.I3Particle.ParticleType.NuE
		or pType==dataclasses.I3Particle.ParticleType.NuEBar
		or pType==dataclasses.I3Particle.ParticleType.NuMu
		or pType==dataclasses.I3Particle.ParticleType.NuMuBar
		or pType==dataclasses.I3Particle.ParticleType.NuTau
		or pType==dataclasses.I3Particle.ParticleType.NuTauBar)
def isChargedLepton(pType):
	return(pType==dataclasses.I3Particle.ParticleType.EMinus
		or pType==dataclasses.I3Particle.ParticleType.EPlus
		or pType==dataclasses.I3Particle.ParticleType.MuMinus
		or pType==dataclasses.I3Particle.ParticleType.MuPlus
		or pType==dataclasses.I3Particle.ParticleType.TauMinus
		or pType==dataclasses.I3Particle.ParticleType.TauPlus)


def finalNeutrino(frame):
	if not (finalSample & is_simulation)(frame):
		return
	mcTree = frame["I3MCTree"]
	try:
		primaryNeutrino = dataclasses.get_most_energetic_primary(mcTree)
	except:
		primaryNeutrino = mcTree.most_energetic_primary

	if(primaryNeutrino==None or not isNeutrinoType(primaryNeutrino.type)):
		return

	#walk down the tree to find the first daughter neutrino which is 'InIce'
	neutrino=primaryNeutrino
	while(neutrino.location_type!=dataclasses.I3Particle.LocationType.InIce):
		children=mcTree.get_daughters(neutrino)
		foundNext=False
		#take the first child which is a neutrino;
		#for in-Earth NC interactions it should be the only one anyway
		for child in children:
			if(isNeutrinoType(child.type)):
				neutrino=child
				foundNext=True
				break
		if(not foundNext):
			#print " did not find a daughter neutrino"
			return #bail out
	frame["InteractingNeutrino"]=neutrino
tray.AddModule(finalNeutrino,"finalNeutrino")

def interaction_type(frame):
	if not (finalSample & is_simulation)(frame):
		return
	mcTree = frame['I3MCTree']
	
	if 'NuFSGen' in options.infile:
		# for NuFSGen, the neutrino is the primary. There is no background.
		neutrino = dataclasses.get_most_energetic_primary(mcTree)
	else:
		# NuGEn is different. Let's find the highest energy neutrino and use that.
		neutrino =  dataclasses.get_most_energetic_neutrino(mcTree)

	# Let's go through the neutrino daughters to see what kind of interaction it was. If there are no children,
	# the interaction was probably just from a cosmic ray. It looks like the injected CORSIKA just propagates muons. so
	# I don't think we have any extra neutrinos.

	PrimaryMuonEnergy     = 0
	PrimaryMuonDir        = 0
	PrimaryMuonAzimuth    = 0
	PrimaryCascadeEnergy  = 0
	PrimaryCascadeDir     = 0
	PrimaryCascadeAzimuth = 0
	neutral_current       = 0
	neutral_current_type  = 0
	charged_current       = 0
	charged_current_type  = 0

	if not neutrino:
		CORSIKA = 1
		print('There is a problem. This better be CORSIKA, there are no neutrinos.')
		return

	children = mcTree.get_daughters(neutrino)
	if len(children) == 0:
		CORSIKA = 1
		#print('This neutrino had no children, CORSIKA. We had a neutrino, but it didn't interact. It could be a CORSIKA event with secondary neutrinos.')
		#print(frame['I3MCTree'])
		print('There is no children on the neutrino! Check this. Probably CORSIKA. Setting CORSIKA =1')
		return
	else:
		CORSIKA = 0

	neutrino_energy = neutrino.energy
	count = 0

	#let's start by talking the total energy of the mu+hadrons. This variable defines the NuFSGen energy range.
	total_energy = 0
	for child in children:
		#check for charged current. Neutral current won't have a Mu- or Mu+.
		if child.type == dataclasses.I3Particle.ParticleType.MuMinus or child.type == dataclasses.I3Particle.ParticleType.MuPlus:
			#if CC add the total energy
			for c in children:
				total_energy += c.energy
	frame.Put("CC_energy",dataclasses.I3Double(total_energy))

	for child in children:
		if(isNeutrinoType(child.type)):
			if child.type == dataclasses.I3Particle.ParticleType.NuE:
				neutral_current_type = 12
				neutral_current = 1
			elif child.type == dataclasses.I3Particle.ParticleType.NuEBar:
				neutral_current_type = -12
				neutral_current = 1
			elif child.type == dataclasses.I3Particle.ParticleType.NuMu:
				neutral_current_type = 14
				neutral_current = 1
			elif child.type == dataclasses.I3Particle.ParticleType.NuMuBar:
				neutral_current_type = -14
				neutral_current = 1
			elif child.type == dataclasses.I3Particle.ParticleType.NuTau:
				neutral_current_type = 16
				neutral_current = 1
			elif child.type == dataclasses.I3Particle.ParticleType.NuTauBar:
				neutral_current_type = -16
				neutral_current = 1
			else:
				neutral_current_type = 0
				neutral_current = 0
				charged_current_type = 0
				charged_current = 0

				PrimaryMuonEnergy = NaN
				PrimaryMuonZenith = NaN
				PrimaryMuonAzimuth = NaN
				PrimaryMuonType = NaN
				print('-- I think this is a Neutral Current interaction:')
				print(frame['I3MCTree'])
				break

		elif (isChargedLepton(child.type)):
			if child.type == dataclasses.I3Particle.ParticleType.EMinus:
				charged_current_type = 12
				charged_current = 1

			elif child.type == dataclasses.I3Particle.ParticleType.EPlus:
				charged_current_type = -12
				charged_current = 1

			elif child.type == dataclasses.I3Particle.ParticleType.MuMinus:
				charged_current_type = 14
				charged_current = 1
			
			PrimaryMuonEnergy = child.energy
			PrimaryMuonZenith = child.dir.zenith
			PrimaryMuonAzimuth = child.dir.azimuth
			PrimaryMuonType = 14

		elif child.type == dataclasses.I3Particle.ParticleType.MuPlus:
			charged_current_type = -14
			charged_current = 1
			PrimaryMuonEnergy = child.energy
			PrimaryMuonZenith = child.dir.zenith
			PrimaryMuonAzimuth = child.dir.azimuth
			PrimaryMuonType = -14

		elif child.type == dataclasses.I3Particle.ParticleType.TauMinus:
			charged_current_type = 16
			charged_current = 1

		elif child.type == dataclasses.I3Particle.ParticleType.TauPlus:
			charged_current_type = -16
			charged_current = 1
		else:
			print('There is a problem with how we are looking for CC events. This event had a lepton as a child, but unknown type.')
			print(frame['I3MCTree'])
	
	#NuGen seems to have a problem, some events don't have secondaries (except Hadrons). These will show up as 0000, and CORSIKA = 0.

	# both CC and NC have cascades.
	cascade_energy = dataclasses.get_most_energetic_cascade(frame['I3MCTree']).energy
	cascade_dir = dataclasses.get_most_energetic_cascade(frame['I3MCTree']).dir.zenith
	cascade_dir = dataclasses.get_most_energetic_cascade(frame['I3MCTree']).dir.azimuth
	'''
	frame.Put("PrimaryMuonEnergy",dataclasses.I3Double(PrimaryMuonEnergy))
	frame.Put("PrimaryMuonDir",dataclasses.I3Double(PrimaryMuonDir))
	frame.Put("PrimaryMuonAzimuth",dataclasses.I3Double(PrimaryMuonAzimuth))
	frame.Put("PrimaryMuonType",dataclasses.I3Double(PrimaryMuonType))
	'''
	frame.Put("injectedMuonEnergy",dataclasses.I3Double(PrimaryMuonEnergy))
	frame.Put("injectedMuonZenith",dataclasses.I3Double(PrimaryMuonDir))
	frame.Put("injectedMuonAzimuth",dataclasses.I3Double(PrimaryMuonAzimuth))
	frame.Put("primaryType",dataclasses.I3Double(PrimaryMuonType))

	frame.Put("PrimaryCascadeEnergy",dataclasses.I3Double(PrimaryCascadeEnergy))
	frame.Put("PrimaryCascadeDir",dataclasses.I3Double(PrimaryCascadeDir))
	frame.Put("PrimaryCascadeAzimuth",dataclasses.I3Double(PrimaryCascadeAzimuth))

	frame.Put("neutral_current",dataclasses.I3Double(neutral_current))
	frame.Put("neutral_current_type",dataclasses.I3Double(neutral_current_type))
	frame.Put("charged_current",dataclasses.I3Double(charged_current))
	frame.Put("charged_current_type",dataclasses.I3Double(charged_current_type))
	frame.Put("CORSIKA",dataclasses.I3Double(CORSIKA))

#tray.Add(interaction_type)


def overburden(frame):
	if frame.Has('MuEx'):
		muex = frame['MuEx']
		z = muex.pos.z
		zen = muex.dir.zenith
		a = 0.259
		b = 0.000363
		R = 12714000/2
		h = 1950
		d = np.sqrt(muex.pos.x**2+muex.pos.y**2)
		L = d/np.tan(zen)
		J = R - h + L + z
		p = np.sqrt(d**2 + (J-L)**2)
		beta = 3*np.pi/2 - zen - np.arctan(J/d)
		alpha = np.arcsin(p * np.sin(beta)/R)
		theta = np.pi - (beta + alpha)
		OB = R*np.sin(theta)/np.sin(beta)
		frame.Put("Overburden",dataclasses.I3Double(OB))

tray.Add(overburden)

def get_pulses(frame):
	if frame.Has('TTPulses') and frame.Has('TTPulses_NoDC'):
		TTpulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses')
		TTpulses_NoDC = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses_NoDC')
		TTnchan  = len(TTpulses)
		TTall_pulses = [p for i,j in TTpulses for p in j]
		TTtotal_charge = sum([p.charge for p in TTall_pulses])
		TTnchan_NoDC = len(TTpulses_NoDC)
		TTall_pulses_NoDC = [p for i,j in TTpulses_NoDC for p in j]
		TTtotal_charge_NoDC = sum([p.charge for p in TTall_pulses_NoDC])
		frame.Put("TTPulses_nchan"	,dataclasses.I3Double(TTnchan))
		frame.Put("TTPulses_qtot"	,dataclasses.I3Double(TTtotal_charge))
		frame.Put("TTPulses_NoDC_nchan"	,dataclasses.I3Double(TTnchan_NoDC))
		frame.Put("TTPulses_NoDC_qtot"	,dataclasses.I3Double(TTtotal_charge_NoDC))

outer_string_list = [31,41,51,60,68,75,76,77,78,72,73,74,67,59,50,40,30,21,13,6,5,4,3,2,1,7,14,22]
top_bottom_om_list = [1,60]

def outer_charge(frame):
	if frame.Has('TTPulses'):
		TTpulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses')
		Radial_charges = 0
		TB_charges = 0
		for omkey, pulse, in TTpulses:
			if omkey[0] in outer_string_list:
				for p in pulse:
					Radial_charges += p.charge
			if omkey[1] in top_bottom_om_list and not omkey[0] in outer_string_list:
				for p in pulse:
					TB_charges += p.charge
		TTPulses_total_outer_charge = Radial_charges+TB_charges
		all_TTpulses = [p for i,j in TTpulses for p in j]
		total_TTcharge = sum([p.charge for p in all_TTpulses])
		frame.Put("TTPulses_total_outer_charge",dataclasses.I3Double(TTPulses_total_outer_charge))
		frame.Put("TTPulses_outer_charge_ratio",dataclasses.I3Double(TTPulses_total_outer_charge/total_TTcharge))
tray.Add(outer_charge)

def get_cut_variables(frame):
	if frame.Has('TTPulses') and frame.Has('TrackFit'):
		geometry = frame["I3Geometry"]
		pulses_map = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses')
		particle = frame["TrackFit"]
		hit_multiplicity_values = hit_multiplicity.calculate_hit_multiplicity(geometry, pulses_map)
		direct_hits_map = direct_hits.calculate_direct_hits(geometry, pulses_map, particle)
		cylinder_radius = 150*I3Units.m
		track_characteristics_values = track_characteristics.calculate_track_characteristics(geometry,pulses_map,particle,cylinder_radius)

		frame.Put("dir_N_D_strings"		,dataclasses.I3Double(direct_hits_map['D'].n_dir_strings))
		frame.Put("dir_N_D_doms"		,dataclasses.I3Double(direct_hits_map['D'].n_dir_doms))
		frame.Put("dir_N_D_pulses"		,dataclasses.I3Double(direct_hits_map['D'].n_dir_pulses))
		frame.Put("dir_L_D"				,dataclasses.I3Double(direct_hits_map['D'].dir_track_length))
		frame.Put("dir_S_D"				,dataclasses.I3Double(direct_hits_map['D'].dir_track_hit_distribution_smoothness))
		frame.Put("dir_N_B_strings"		,dataclasses.I3Double(direct_hits_map['B'].n_dir_strings))
		frame.Put("dir_N_B_doms"		,dataclasses.I3Double(direct_hits_map['B'].n_dir_doms))
		frame.Put("dir_N_B_pulses"		,dataclasses.I3Double(direct_hits_map['B'].n_dir_pulses))
		frame.Put("dir_L_B"				,dataclasses.I3Double(direct_hits_map['B'].dir_track_length))
		frame.Put("dir_S_B"				,dataclasses.I3Double(direct_hits_map['B'].dir_track_hit_distribution_smoothness))
		frame.Put("trackfit_separation"	,dataclasses.I3Double(track_characteristics_values.track_hits_separation_length))
		return(True)
tray.Add(get_cut_variables)
tray.Add(get_pulses)

# Find the muon we will use for weighting
def getWeightingCandidate(frame):
	WeightingCandidatesThisRound=0
	if frame.Has('I3MCTree'):
		for i in frame['I3MCTree']:
			if((i.type==i.ParticleType.MuMinus) or (i.type==i.ParticleType.MuPlus)):
				WeightingCandidatesThisRound+=1
				Candidate=i
		if(WeightingCandidatesThisRound==1):
			frame.Put("WeightingMuon",dataclasses.I3Particle(Candidate))
		else:
			print "Warning! Multiple weighting candidates found! Storing none of them"
			print "This is fine for NuGen, NuFSGen shouldn't see this message Since there is always only one muon."



tray.AddModule(getWeightingCandidate,"weightingcand")

finalWithParab=finalSample & parabaloidCut

tray.AddModule(finalWithParab,"finalSample")

tray.AddModule('Delete', 'delkeys',
	Keys = ['InIceRawData',
		'I3MCPulseSeriesMap',
		'I3MCPulseSeriesMapParticleIDMap',
		'CleanIceTopRawData',
		'ClusterCleaningExcludedTanks',
		#'IceTopDSTPulses',
		'IceTopPulses',
		'IceTopRawData',
		'FilterMask_NullSplit0',
		'FilterMask_NullSplit1',
		'FilterMask_NullSplit2',
		'FilterMask_NullSplit3',
		'FilterMask_NullSplit4',
		'FilterMask_NullSplit5',
		'FilterMask_NullSplit6',
		'OfflineIceTopHLCTankPulses',
		'OfflineIceTopHLCVEMPulses',
		'OfflineIceTopSLCVEMPulses',
		'TankPulseMergerExcludedTanks'])

tray.AddModule("PacketCutter","Cutter")(
		("CutStream","TTrigger"),
		("CutObject","CutL3")
		)

i3streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics,icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation]

outfile = options.outfile
outputKeys = ["GenerationSpec",'MuEx']

if options.move == 'True':
	tray.AddModule("I3Writer","i3writer")(
		("Filename",outfile_temp),
		("Streams",i3streams)
		)
	if options.cut == 'True':
		tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
			("tableservice",hdfwriter.I3HDFTableService(outfile_temp.replace('.i3.bz2','_golden.h5'))),
			("SubEventStreams",["TTrigger"]),
			("keys",outputKeys)
			)
	else:
		pass
else:
	if osg =='True':
		print('You need options.move to be True if using the OSG')
	else:
		tray.AddModule("I3Writer","i3writer")(
			("Filename",outfile),
			("Streams",i3streams)
			)
		if options.cut == 'True':
			tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
				("tableservice",hdfwriter.I3HDFTableService(outfile_temp.replace('.i3.bz2','_golden.h5'))),
				("SubEventStreams",["TTrigger"]),
				("keys",outputKeys)
				)
		else:
			pass

if(options.nFrames==0):
	tray.Execute()
else:
	tray.Execute(options.nFrames)
tray.Finish()

# If the file finished processing, move it to the storage location.

if options.move == 'True':
	if osg == "True":
		copy_to_NPX('gsiftp://gridftp.icecube.wisc.edu' + outfile)
		if options.cut == 'True':
			copy_to_NPX('gsiftp://gridftp.icecube.wisc.edu' + outfile.replace('.i3.bz2','_golden.h5'))
	else:
		os.system(str("mv "+str(outfile_temp) + " " +str(outfile)))
		if options.cut == 'True':
			os.system(str("mv "+str(outfile_temp.replace('.i3.bz2','_golden.h5')) + " " +str(outfile.replace('.i3.bz2','_golden.h5'))))
else:
	if osg =="True":
		print('You need options.move to be True if using the OSG')
	pass

end_time = time.time()
print("Processing time: "+str(end_time-start_time))

exit(exitStatus.status)


