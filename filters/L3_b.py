#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v3/icetray-start
#METAPROJECT /data/user/nwandkowsky/tarballs/icerec.V05-02-00_legacy

# ./L3_b.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz -i L3a_nugen.i3.zst -o L3b_nugen.i3.zst


from I3Tray import *
from icecube import lilliput, gulliver, icetray, phys_services, dataio, dataclasses
import sys, math, numpy
import copy
import glob

import os

import time

start_time = time.asctime()
print 'Started:', start_time

# handling of command line arguments
from optparse import OptionParser
parser = OptionParser()
usage = """%prog [options]"""
parser.set_usage(usage)

parser.add_option("-i", "--input", action="store", type="string", default="", dest="INPUT", help="Output i3 file")
parser.add_option("-o", "--output", action="store", type="string", default="", dest="OUTPUT", help="Output i3 file")
parser.add_option("-g", "--gcd", action="store", type="string", default="", dest="GCD", help="GCD file for input i3 file")

# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

infile = options.INPUT
outfile = options.OUTPUT

tray = I3Tray()

tray.Add("I3Reader", 'reader', FilenameList=[options.GCD, infile])
################################
## ALL THE TRACK FIT GOODIES
################################
from icecube import DomTools
from icecube import level3_filter_muon
tray.AddModule("StaticDOMTimeWindowCleaning",
        InputPulses="TSInIcePulses",
        OutputPulses="TWTSInIcePulses",
        MaximumTimeDifference=3e3*I3Units.ns)

load("libmue")
tray.AddModule("muex", "muex_angular4",
        Pulses="TWTSInIcePulses",
        rectrk="",
        result="MuEXAngular4",
        lcspan=0,
        repeat=4,
        usempe=True,
        detail=False,
        energy=False,
        icedir=os.path.expandvars("$I3_BUILD/mue/resources/ice/mie"))

#from icecube import photonics_service_old
# spline MPE
from icecube.photonics_service import I3PhotoSplineService
infmuonprobsplinepath = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits"
infmuonampsplinepath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'
spline=I3PhotoSplineService(infmuonampsplinepath,infmuonprobsplinepath,4)
llh="MPE"

cascade_service_mie = I3PhotoSplineService(amplitudetable = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.abs.fits',\
        timingtable = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/ems_mie_z20_a10.prob.fits',\
        timingSigma = 0.)

tray.AddService("I3BasicSeedServiceFactory", "SplineSeed", FirstGuesses=['MuEXAngular4'])

from icecube import spline_reco
tray.AddService("I3SplineRecoLikelihoodFactory","LLHSpline",
        PhotonicsService=spline,
        Pulses="TWTSInIcePulses",
        Likelihood=llh,
        NoiseRate=10*I3Units.hertz)

tray.AddService("I3GulliverMinuitFactory", "Minuit",
        Algorithm="SIMPLEX",
        MaxIterations=1000,
        Tolerance=0.01)

tray.AddService("I3SimpleParametrizationFactory", "SimpleTrack",
        StepX = 20*I3Units.m,
        StepY = 20*I3Units.m,
        StepZ = 20*I3Units.m,
        StepZenith = 0.1*I3Units.radian,
        StepAzimuth= 0.2*I3Units.radian,
        BoundsX = [-2000*I3Units.m, 2000*I3Units.m],
        BoundsY = [-2000*I3Units.m, 2000*I3Units.m],
        BoundsZ = [-2000*I3Units.m, 2000*I3Units.m])

tray.AddModule("I3SimpleFitter", "SplineFit",
        SeedService="SplineSeed",
        Parametrization="SimpleTrack",
        LogLikelihood="LLHSpline",
        Minimizer="Minuit",)

# muex - differential energy
tray.AddModule("muex", "muex_differential",
        Pulses = "TWTSInIcePulses",
        rectrk = "SplineFit",
        result = "SplineMPEMuEXDifferential",
        detail = True, # differential
        energy = True,
        lcspan = 0,
        icedir=os.path.expandvars("$I3_BUILD/mue/resources/ice/mie"))

tray.Add('Delete',Keys=['SaturatedDOMs','BrightDOMs'])
from icecube.millipede import HighEnergyExclusions, MonopodFit
exclusions = tray.AddSegment(HighEnergyExclusions, Pulses="SplitInIcePulses")

# millipede
tray.AddModule("MuMillipede", "millipede",
        MuonPhotonicsService=None,
        CascadePhotonicsService=cascade_service_mie,
        PhotonsPerBin=5,
        ShowerRegularization=1e-9,
        MuonSpacing=0.*I3Units.m,
        ShowerSpacing=2.5*I3Units.m,
        SeedTrack="SplineFit",
        Output="Millipede_SplineMPE",
        ReadoutWindow="SplitInIcePulsesTimeRange",
        ExcludedDOMs=exclusions,
        DOMEfficiency=0.99,
        Pulses="SplitInIcePulses",
        PartialExclusion=True)

#tray.Add('Dump')

def determine_millipede_deposited_energy(frame):
    if "Millipede_SplineMPE" in frame:
        losses = frame["Millipede_SplineMPE"]
        first_loss_found=False
        deposited_energy = 0
        for i in range(0,len(losses)):
            if losses[i].pos.z<500 and losses[i].pos.z>-500 and math.sqrt(losses[i].pos.x*losses[i].pos.x+losses[i].pos.y*losses[i].pos.y)<550:
                deposited_energy = deposited_energy + losses[i].energy
                if first_loss_found==False and losses[i].energy>10:
                    frame["MillipedeFirstLoss"] = losses[i]
                    first_loss_found=True
    else:
        deposited_energy = 0.

    if "MillipedeFirstLoss" not in frame:
        frame["MillipedeFirstLoss"] = frame["MonopodFit4_lea_tilt"]
    frame["MillipedeDepositedEnergy"] = dataclasses.I3Double(deposited_energy)
    #print "deposited energy: ", frame["MillipedeDepositedEnergy"].value
tray.Add(determine_millipede_deposited_energy)

#####################################
## INCOMING TRACK VETO W/O CUTS
#####################################

from icecube.CascadeL3_IC79.level4.veto import ACausalHitSelector, TrackVeto, LayerVeto, VetoMarginCalculator
cleaned_pulses = tray.AddSegment(ACausalHitSelector, 'ACausalSplitInIcePulses_milli', pulses='SplitInIcePulses', vertex="MillipedeFirstLoss")
tray.AddSegment(TrackVeto, 'VetoTrack_milli', TimeWindow=(-15, 1000), Pulses=cleaned_pulses, Radius=100,
        Direction='down', TrackType=dataclasses.I3Particle.StoppingTrack, Vertex="MillipedeFirstLoss")

tray.AddSegment(TrackVeto, 'UpgoingTrack_milli', TimeWindow=(-30, 500), Pulses=cleaned_pulses, Radius=100,
        Direction='up', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MillipedeFirstLoss")

cleaned_pulses_hlc = tray.AddSegment(ACausalHitSelector, 'ACausalSplitInIcePulsesHLC_milli', pulses='SplitInIcePulsesHLC', vertex="MillipedeFirstLoss")
tray.AddSegment(TrackVeto, 'StartingTrackHLC_milli', TimeWindow=(-30, 500), Pulses=cleaned_pulses_hlc, Radius=100,
        Direction='all', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MillipedeFirstLoss")

tray.AddSegment(TrackVeto, 'UpgoingTrackHLC_milli', TimeWindow=(-30, 500), Pulses=cleaned_pulses_hlc, Radius=100,
        Direction='up', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MillipedeFirstLoss")

tray.AddModule(VetoMarginCalculator, 'VetoTrackMargin_milli', Vertex="MillipedeFirstLoss", Geometry="I3Geometry")

##################################################
## CASCADE/TRACK DISCRIMINATOR ESSENTIALS
##################################################
def getenergyconfinement(cascade,track_particles):
    key='reco'
    # due to a falsely configured muon spacing in millipede, need to split the cascade and track segments and then rebin
    if key == 'reco':
        particle_buffer = []
        track_dx_buffer, track_dy_buffer, track_dz_buffer, track_de_buffer = [], [], [], []
        for particle in track_particles:
            if particle.shape == dataclasses.I3Particle.Cascade:
                track_dx_buffer.append(particle.pos.x)
                track_dy_buffer.append(particle.pos.y)
                track_dz_buffer.append(particle.pos.z)
                track_de_buffer.append(particle.energy)
            else:
                particle_buffer.append(particle)
        track_dx_buffer = numpy.asarray(track_dx_buffer)
        track_dy_buffer = numpy.asarray(track_dy_buffer)
        track_dz_buffer = numpy.asarray(track_dz_buffer)
        track_de_buffer = numpy.asarray(track_de_buffer)
        for particle in particle_buffer:
            argmin = numpy.argmin(numpy.sqrt((track_dx_buffer - particle.pos.x)**2
                + (track_dy_buffer - particle.pos.y)**2
                + (track_dz_buffer - particle.pos.z)**2))
            track_de_buffer[argmin] += particle.energy

    distancethreshold = 40.
    mask = ( numpy.sqrt((track_dx_buffer-cascade.pos.x)**2
        + (track_dy_buffer-cascade.pos.y)**2
        + (track_dz_buffer-cascade.pos.z)**2) < distancethreshold )

    # calculate the energy confinement
    esum = numpy.sum(track_de_buffer)
    edep = numpy.sum(track_de_buffer[mask])
    econfinement = edep/esum
    return econfinement

def get_track_length(frame):
    if "Millipede_SplineMPE" in frame:
        losses = frame["Millipede_SplineMPE"]
        first_loss_found=False
        for i in range(0,len(losses)):
            if first_loss_found==False and losses[i].energy>10 and losses[i].pos.z<500 and losses[i].pos.z>-500 and math.sqrt(losses[i].pos.x*losses[i].pos.x+losses[i].pos.y*losses[i].pos.y)<550:
                first_loss = losses[i]
               first_loss_found=True
            if losses[i].energy>10 and losses[i].pos.z<500 and losses[i].pos.z>-500 and math.sqrt(losses[i].pos.x*losses[i].pos.x+losses[i].pos.y*losses[i].pos.y)<550:
                last_loss = losses[i]
        if first_loss_found==False:
            frame["TrackLength"] = dataclasses.I3Double(0.)
        else:
            dist = numpy.sqrt((first_loss.pos.x-last_loss.pos.x)*(first_loss.pos.x-last_loss.pos.x)+ (first_loss.pos.y-last_loss.pos.y)*(first_loss.pos.y-last_loss.pos.y)+ (first_loss.pos.z-last_loss.pos.z)*(first_loss.pos.z-last_loss.pos.z))
            frame["TrackLength"] = dataclasses.I3Double(dist)
    else:
        frame["TrackLength"] = dataclasses.I3Double(0.)
    #print "Track length: ", frame["TrackLength"].value
    frame["EnergyConfinement"] = dataclasses.I3Double(getenergyconfinement(frame["MonopodFit4_lea_tilt"],frame["Millipede_SplineMPE"]))
    #print "Energy confinement: ", frame["EnergyConfinement"].value
tray.Add(get_track_length)

tray.Add('I3LCPulseCleaning', 'cleaning_fortw', OutputHLC='TWTSInIcePulsesHLC', OutputSLC='', Input='TWTSInIcePulses')

def CleanDeepCore(frame):
    mask = dataclasses.I3RecoPulseSeriesMapMask(frame, "TWTSInIcePulsesHLC", lambda om, idx, pulse: om.string <= 78)
        frame["TWTSInIcePulsesHLC_NoDC"] = mask
tray.AddModule(CleanDeepCore, 'nodc_again')

def getRecoPulses(frame,name):
    pulses = frame[name]
        if pulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
            pulses = pulses.apply(frame)
        return pulses

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

tray.AddModule(ComputeChargeWeightedDist,"CCWD_Track",Pulses="TWTSInIcePulsesHLC_NoDC",Track="SplineMPE")
tray.AddModule(ComputeChargeWeightedDist,"CCWD_Cascade",Pulses="TWTSInIcePulsesHLC_NoDC",Track="MonopodFit4_lea_tilt")

def check_type(frame):
    if "MCPrimary" in frame:
        if frame["MCPrimary"].type==12 or frame["MCPrimary"].type==-12:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)

        if (frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14) and frame["I3MCWeightDict"]["InteractionType"]==2.0:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)
        elif (frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14) and frame["I3MCWeightDict"]["InteractionType"]==1.0:
            frame["IsCascade_true"]= icetray.I3Bool(False)
                frame["IsTrack_true"]= icetray.I3Bool(True)

        if (frame["MCPrimary"].type==16 or frame["MCPrimary"].type==-16) and frame["I3MCWeightDict"]["InteractionType"]==2.0:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)
        elif  (frame["MCPrimary"].type==16 or frame["MCPrimary"].type==-16) and frame["I3MCWeightDict"]["InteractionType"]==1.0:
            frame["IsCascade_true"]= icetray.I3Bool(False)
            frame["IsTrack_true"]= icetray.I3Bool(True)
tray.Add(check_type)

#########################################
## GATHER SOME MC TRUTH
#########################################

from icecube import VHESelfVeto
# find detector volume intersection for weighting
def _FindDetectorVolumeIntersections(frame, recoParticle, geometry):
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

    min_index = 0
    for i in range(len(intersectionTimes)):
        if i==0:
            min_time = intersectionTimes[i]
        else:
            if intersectionTimes[i]<min_time:
                min_index = i
                min_time = intersectionTimes[i]
    if len(intersectionTimes)==0:
        frame["IntersectionPoint"] = dataclasses.I3Position()
    else:
        frame["IntersectionPoint"] = dataclasses.I3Position(intersectionPoints[min_index])

    sortedTimes = sorted(intersectionTimes)
    return sortedTimes

def FindDetectorVolumeIntersections(frame, TrackName="", OutputTimeWindow=None, TimePadding=0.):
    if TrackName in frame:
        if OutputTimeWindow is not None:
            twName = OutputTimeWindow
        else:
            twName = TrackName + "TimeRange"

            theTrack = frame[TrackName]
            geometry = frame["I3Geometry"]

            times = _FindDetectorVolumeIntersections(frame, theTrack, geometry)

            if len(times) == 0:
                #raise RuntimeError("track does not intersect the detector volume")
                    frame[twName] = dataclasses.I3TimeWindow()
            elif len(times) == 1:
                raise RuntimeError("tracks with only one intersection are not supported")
            else:
                tWindow = dataclasses.I3TimeWindow(times[0]-TimePadding, times[-1]+TimePadding)
                    frame[twName] = tWindow
tray.AddModule(FindDetectorVolumeIntersections, "FindDetectorVolumeIntersections",
        TimePadding = 60.*I3Units.m/dataclasses.I3Constants.c,
        TrackName="MCTrack",
        OutputTimeWindow="ContainedTimeRange")

# Determine true deposited energy
def collectStats(frame):
    if "I3MCTree" in frame:
        timeWindow = frame["ContainedTimeRange"]
        validRows = []

        losses = 0
        emlosses = 0
        hadlosses = 0
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
            if p.time > timeWindow.start and p.time < timeWindow.stop:
                if p.type in [p.Hadrons, p.PiPlus, p.PiMinus, p.NuclInt]:
                    hadlosses += p.energy
                    if p.energy < 1*I3Units.GeV:
                        losses += 0.8*p.energy
                    else:
                        energyScalingFactor = 1.0 + ((p.energy/I3Units.GeV/0.399)**-0.130)*(0.467 - 1)
                                losses += energyScalingFactor*p.energy
                else:
                    emlosses += p.energy
                    losses += p.energy
        frame["DepositedEnergyMC"] = dataclasses.I3Double(losses)
        frame["DepositedEMEnergyMC"] = dataclasses.I3Double(emlosses)
        frame["DepositedHadEnergyMC"] = dataclasses.I3Double(hadlosses)
        #print "True deposited: ", frame["DepositedEnergyMC"].value
tray.AddModule(collectStats, "CollectStats")

#tray.Add('Delete',
tray.AddModule('I3Writer', 'writer',
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        filename=outfile)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

del tray

stop_time = time.asctime()

print 'Started:', start_time
print 'Ended:', stop_time
