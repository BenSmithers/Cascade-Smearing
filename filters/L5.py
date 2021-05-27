#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/combo/py3-v4.1.1/


##!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v1/icetray-start
##METAPROJECT /data/user/nwandkowsky/tarballs/icerec.V04-11-02

#./L5_nugen.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz -i /data/ana/Cscd/StartingEvents/NuGen/NuMu/low_energy/IC86_2011/l4/1/l4_00000002.i3.bz2 -o /data/user/nwandkowsky/simulations/2014/signal/condor/event_selection/data/nugen/IC86_2011/numu.i3
from I3Tray import *
from icecube import icetray, dataio, dataclasses
import sys, math, numpy
import copy
import glob

import os

import time

start_time = time.asctime()
print('Started:', start_time)

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

infiles=[options.GCD, options.INPUT]

tray = I3Tray()
tray.AddModule('I3Reader', 'reader', FilenameList=infiles)
print("Reading input file...", len(infiles))

def charge_cut(frame):
    if frame["HomogenizedQTot"].value<100.:
        return False
    else:
        return True
tray.Add(charge_cut)

def check_filters(frame):
    NPassedFilters = 0
    if frame["CascadeFilter"].value==True:
        NPassedFilters = NPassedFilters+1

    if frame["EHEFilter"].value==True:
        NPassedFilters = NPassedFilters+1

    if frame["MuonFilter"].value==True:
        NPassedFilters = NPassedFilters+1

    return NPassedFilters>1
tray.Add(check_filters)

def event_filter(frame):
    #print "Before filter: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
    del frame["IsUpgoingMuon"]
    del frame["IsCascade"]
    del frame["IsHese"]

    ###########
    ##  HESE
    ###########
    if frame["IsHESE_ck"].value==False and (frame["L4VetoLayer0"].value+frame["L4VetoLayer1"].value)==0 and frame["HomogenizedQTot"].value>6000. and frame["IsCascade_reco"].value==False:
        frame["IsHese"]=icetray.I3Bool(False)
        print("Jakob HESE event that's probably a muon removed...")
    elif (frame["L4VetoLayer0"].value+frame["L4VetoLayer1"].value==0 and frame["HomogenizedQTot"].value>6000.) or frame["IsHESE_ck"].value==True:
        #print "HESE event! ",frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
        frame["IsHese"]=icetray.I3Bool(True)
        frame["IsCascade"]=icetray.I3Bool(False)
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)
        return True
    else:
        frame["IsHese"]=icetray.I3Bool(False)
    if frame["IsHese"].value==False and frame["HomogenizedQTot"].value>6000.:
        return False


    ###########
    ##  CASCADES
    ###########

    # MONOPOD WITH OFFLINE PULSES
    if  frame['L4VetoTrackMarginOfflineTop'].value > 100 and frame['L4VetoTrackMarginOfflineSide'].value > 0:
        side = math.log10(frame['HomogenizedQTot'].value) > 3.40876826865 - ((frame['L4VetoTrackMarginOfflineSide'].value+5)**1.73562251277)/17266.4002818
        top  = math.log10(frame['HomogenizedQTot'].value) > 3.40214654218 -((frame['L4VetoTrackMarginOfflineTop'].value-100)**1.880562951)/23709.9897396
    else:
        side=False
        top=False
    inside_volume_mono_offline = (side and top)
    # MONOPOD WITH SPLIT PULSES
    if  frame['L4VetoTrackMarginSplitTop'].value > 100 and frame['L4VetoTrackMarginSplitSide'].value > 0:
        side = math.log10(frame['HomogenizedQTot'].value) > 3.40876826865 - ((frame['L4VetoTrackMarginSplitSide'].value+5)**1.73562251277)/17266.4002818
        top  = math.log10(frame['HomogenizedQTot'].value) > 3.40214654218 -((frame['L4VetoTrackMarginSplitTop'].value-100)**1.880562951)/23709.9897396
    else:
        side=False
        top=False
    inside_volume_mono_split = (side and top)

    # MILLIPEDE WITH OFFLINE PULSES
    if  frame['L4VetoTrackMarginMilliOfflineTop'].value > 100 and frame['L4VetoTrackMarginMilliOfflineSide'].value > 0:
        side = math.log10(frame['HomogenizedQTot'].value) > 3.40876826865 - ((frame['L4VetoTrackMarginMilliOfflineSide'].value)**1.73562251277)/17266.4002818
        top  = math.log10(frame['HomogenizedQTot'].value) > 3.40214654218 -((frame['L4VetoTrackMarginMilliOfflineTop'].value-100)**1.880562951)/23709.9897396
    else:
        side=False
        top=False
    inside_volume_milli_offline = (side and top)

    if frame["IsCascade_reco"].value==True:
        if frame["CascadeFilter"].value==True and frame["HomogenizedQTot_split"].value>100. and ((frame["L4VetoTrackOfflineVetoCharge"].value==0 and inside_volume_mono_offline) or (frame["L4VetoTrackOfflineVetoCharge"].value<6 and \
                inside_volume_mono_offline and inside_volume_milli_offline)) and \
                (numpy.log10(frame["L4MonopodFit"].energy)-numpy.cos(frame["L4MonopodFit"].dir.zenith)>2.2 or \
                numpy.cos(frame["L4MonopodFit"].dir.zenith)<0.65):
            frame["IsCascade"]=icetray.I3Bool(True)
            frame["IsUpgoingMuon"]=icetray.I3Bool(False)
            #print "CASCADE (cascade): ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,frame['HomogenizedQTot'].value#, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
            return True
        else:
            print("Cascade rejected: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,\
                    frame["L4VetoTrackOfflineVetoCharge"].value, inside_volume_mono_offline, \
                    inside_volume_milli_offline, (numpy.log10(frame["L4MonopodFit"].energy)-numpy.cos(frame["L4MonopodFit"].dir.zenith)), numpy.cos(frame["L4MonopodFit"].dir.zenith))
            frame["IsCascade"]=icetray.I3Bool(False)
    else:
        if (frame["MuonFilter"].value==True and frame["HomogenizedQTot_split"].value>100. and frame["L4VetoTrackMilliOfflineVetoCharge"].value<2 and \
            frame["L4VetoTrackOfflineVetoCharge"].value<2 and inside_volume_milli_offline and inside_volume_mono_offline):# or\
            #(frame["L4VetoTrackMilliSplitVetoCharge"].value<2 and inside_volume_milli_split):
            frame["IsCascade"]=icetray.I3Bool(True)
            frame["IsUpgoingMuon"]=icetray.I3Bool(False)
            #print "CASCADE (track): ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,frame['HomogenizedQTot'].value#, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
            return True
        else:
            print("Cascade (tracklike) rejected ",  frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,\
                    frame["L4VetoTrackOfflineVetoCharge"].value, inside_volume_mono_offline, \
                    frame["L4VetoTrackMilliOfflineVetoCharge"].value, inside_volume_milli_offline)
            frame["IsCascade"]=icetray.I3Bool(False)

    ###########
    ##  UPMU
    ###########
    if ( ( frame['L4UpgoingTrackOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackOfflineVetoCharge'].value > \
            frame['L4VetoTrackOfflineVetoCharge'].value and frame['L4UpgoingTrackOfflineVetoChannels'].value > 3 ) or \
            ( frame['L4UpgoingTrackSplitVetoCharge'].value > 10 and frame['L4UpgoingTrackSplitVetoCharge'].value > \
            frame['L4VetoTrackSplitVetoCharge'].value and frame['L4UpgoingTrackSplitVetoChannels'].value > 3 )  or \
            ( frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > \
            frame['L4VetoTrackMilliOfflineVetoCharge'].value and frame['L4UpgoingTrackMilliOfflineVetoChannels'].value > 3 and \
            frame['L4UpgoingTrackSplitVetoCharge'].value > 6 ) )\
            and frame["TrackFit"].dir.zenith>1.5 and frame["MuonFilter"].value==True:# and frame["OnlineL2Filter"].value==True:
                frame["IsUpgoingMuon"]=icetray.I3Bool(True)
        #print "Upgoing Muon event!", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
    else:
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)
    return frame["IsUpgoingMuon"].value==True or frame["IsHese"].value==True or frame["IsCascade"].value==True
tray.Add(event_filter)


def get_track_length(frame):
    losses = frame["MuonLosses"]
    first_loss_found=False
    for i in range(0,len(losses)):
        if first_loss_found==False and losses[i].energy>1 and losses[i].pos.z<500 and losses[i].pos.z>-500 and math.sqrt(losses[i].pos.x*losses[i].pos.x+losses[i].pos.y*losses[i].pos.y)<550:
            first_loss = losses[i]
            first_loss_found=True
        if losses[i].energy>1 and losses[i].pos.z<500 and losses[i].pos.z>-500 and math.sqrt(losses[i].pos.x*losses[i].pos.x+losses[i].pos.y*losses[i].pos.y)<550:
            last_loss = losses[i]
    if first_loss_found==False:
        frame["TrackLength"] = dataclasses.I3Double(0.)
    else:
        dist = numpy.sqrt((first_loss.pos.x-last_loss.pos.x)*(first_loss.pos.x-last_loss.pos.x)+ (first_loss.pos.y-last_loss.pos.y)*(first_loss.pos.y-last_loss.pos.y)+ (first_loss.pos.z-last_loss.pos.z)*(first_loss.pos.z-last_loss.pos.z))
        frame["TrackLength"] = dataclasses.I3Double(dist)
tray.Add(get_track_length)

def discriminate(frame):
    del frame["IsCascade_reco"]
    del frame["IsTrack_reco"]

    if frame["IsUpgoingMuon"].value==True:
        frame["IsCascade_reco"]=icetray.I3Bool(False)
        frame["IsTrack_reco"]=icetray.I3Bool(True)
        return True

    if ((frame["L4MonopodFit_AvgDistQ"].value-30)*0.9<frame["TrackFit_AvgDistQ"].value):
        is_track_reco=False
        is_cascade_reco=True
    else:
        is_track_reco=True
        is_cascade_reco=False

    if frame["IsCascade"].value==True:
        if is_track_reco==True:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
        elif frame["L4StartingTrackHLCOfflineVetoCharge"].value==0:
            frame["IsCascade_reco"]=icetray.I3Bool(True)
            frame["IsTrack_reco"]=icetray.I3Bool(False)
        elif frame["L4StartingTrackHLCOfflineVetoCharge"].value>1.5:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
        else:
            frame["IsCascade_reco"]=icetray.I3Bool(True)
            frame["IsTrack_reco"]=icetray.I3Bool(False)

    if frame["IsHese"].value==True:
        if is_track_reco==True and frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>50.:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
        elif frame["TrackLength"].value>550 and  frame["L4StartingTrackHLCOfflineVetoCharge"].value>30.:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
        elif frame["L4StartingTrackHLCOfflineVetoCharge"].value==0:
            frame["IsCascade_reco"]=icetray.I3Bool(True)
            frame["IsTrack_reco"]=icetray.I3Bool(False)
        elif ((frame["L4StartingTrackHLCOfflineVetoCharge"].value>1000. and frame["L4StartingTrackHLCOfflineVetoCharge"].value<10000.) or (frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>1500. and frame["L4StartingTrackHLCOfflineVetoCharge"].value>200.)) and  frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>frame["L4StartingTrackHLCOfflineVetoCharge"].value:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
        elif  frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>200 and frame["L4StartingTrackHLCOfflineVetoCharge"].value>60 and frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>frame["L4StartingTrackHLCOfflineVetoCharge"].value:
            frame["IsCascade_reco"]=icetray.I3Bool(False)
            frame["IsTrack_reco"]=icetray.I3Bool(True)
            #print frame["IsTrack_true"].value, frame["IsTrack_reco"].value, " HESE track (both charges>30)", frame["L4StartingTrackHLCOfflineVetoCharge"].value, frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value
        else:
            frame["IsCascade_reco"]=icetray.I3Bool(True)
            frame["IsTrack_reco"]=icetray.I3Bool(False)
tray.Add(discriminate)

def remove_coincident(frame):
    distance = numpy.sqrt( (frame["SPEFit4_split"].pos.x-frame["SplineMPE_split"].pos.x)**2+(frame["SPEFit4_split"].pos.y-frame["SplineMPE_split"].pos.y)**2 +(frame["SPEFit4_split"].pos.z-frame["SplineMPE_split"].pos.z)**2)
    frame["SPEFit4Distance"]=dataclasses.I3Double(distance)
    distance = numpy.sqrt( (frame["SPEFitSingle_split"].pos.x-frame["SplineMPE_split"].pos.x)**2+(frame["SPEFitSingle_split"].pos.y-frame["SplineMPE_split"].pos.y)**2 +(frame["SPEFitSingle_split"].pos.z-frame["SplineMPE_split"].pos.z)**2)
    frame["SPEFitSingleDistance"]=dataclasses.I3Double(distance)
    not_coinc = frame["SPEFit4Distance"].value<200 and frame["SPEFitSingleDistance"].value<200
    #print "Coincident Event: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
    if frame["IsUpgoingMuon"].value==True and not_coinc==False:
        return False

    if frame["IsHese"].value == False and frame["L4MonopodFit_AvgDistQ"].value>150. and frame["TrackFit_AvgDistQ"].value>110:
        return False
    else:
        return True
tray.Add(remove_coincident)


def coinc_cut(frame):
    if "IsUpgoingMuon" in frame:
        if frame["IsUpgoingMuon"].value==True and (frame["TrackFit"].pos.z>370. or frame["SPEFit4_split"].dir.zenith<1.5):
            #print "coinc event removed: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
            return False
        elif frame["IsCascade"].value==True and frame["TrackFit"].pos.z>370.:
            #print "Cascade removed, potentialli coincident: ",  frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
            return True
tray.Add(coinc_cut)

def track_cut(frame):
    def vector(zenith, azimuth):
        return [numpy.sin(zenith)*numpy.cos(azimuth), numpy.sin(zenith)*numpy.sin(azimuth), numpy.cos(zenith)]

    def unit_vector(vector):
        """ Returns the unit vector of the vector.  """
        return vector / numpy.linalg.norm(vector)

    def angle_between(v1, v2):
        v1_u = unit_vector(v1)
        v2_u = unit_vector(v2)
        return numpy.arccos(numpy.clip(numpy.dot(v1_u, v2_u), -1.0, 1.0))

    vec1 = vector(frame["SPEFit2_offline"].dir.zenith,frame["SPEFit2_offline"].dir.azimuth)
    vec2 = vector(frame["TrackFit"].dir.zenith,frame["TrackFit"].dir.azimuth)
    angle_diff = angle_between(vec1,vec2)

    if frame["IsUpgoingMuon"].value==True and angle_diff>0.5:
        #print "track event removed: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, numpy.cos(frame["TrackFit"].dir.zenith)
        return False

    vec1 = vector(frame["SPEFit4_split"].dir.zenith,frame["SPEFit4_split"].dir.azimuth)
    vec2 = vector(frame["TrackFit"].dir.zenith,frame["TrackFit"].dir.azimuth)
    angle_diff = angle_between(vec1,vec2)

    if frame["IsUpgoingMuon"].value==True and angle_diff>0.5:
        #print "track event removed: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, numpy.cos(frame["TrackFit"].dir.zenith)
        return False
    else:
        return True

tray.Add(track_cut)

def check_fits(frame):
    # check MonopodFit
    if frame["IsHese"].value==True or (frame["L4MonopodFit"].pos.z<600 and frame["L4MonopodFit"].pos.z>-600 and numpy.sqrt(frame["L4MonopodFit"].pos.x*frame["L4MonopodFit"].pos.x+frame["L4MonopodFit"].pos.y*frame["L4MonopodFit"].pos.y)<600):
        return True
    else:
        return False
    if "SplineMPEMuEXDifferential" in frame:
        if frame["IsHese"].value==False and frame["IsCascade_reco"].value == False and ( frame["SplineMPEMuEXDifferential"].energy == 0.0 or numpy.isnan(frame["SplineMPEMuEXDifferential"].energy) ):
            #print "Bad reco event..."
            return False
tray.Add(check_fits)

from icecube.CascadeL3_IC79.level3.reco import CascadeLlhVertexFit
tray.AddSegment(CascadeLlhVertexFit, 'L5CscdLlh', Pulses='TWOfflinePulsesHLC')

from icecube.photonics_service import I3PhotoSplineService
base = os.path.expandvars('$I3_DATA/photon-tables/splines/ems_mie_z20_a10.%s.fits')
pxs = I3PhotoSplineService(base % "abs", base % "prob", 0)

tray.Add('Delete', Keys=['BrightDOMs'])
from icecube import  millipede
exclusions = tray.AddSegment(millipede.HighEnergyExclusions, 'monopod_DOM_exclusions',
        Pulses='TWOfflinePulsesHLC',
        CalibrationErrata='CalibrationErrata',
        ExcludeSaturatedDOMs='OfflineInIceCalibrationErrata',
        ExcludeBrightDOMs='BrightDOMs',
        ExcludeDeepCore=False,
        BadDomsList='BadDomsList',
        SaturationWindows='SaturationTimes'
        )
millipede_config = dict(Pulses="OfflinePulses", PartialExclusion=False, CascadePhotonicsService=pxs, BadDOMs=exclusions, Parametrization='HalfSphere', MinTimeWidth=15)

"""
exclusions = tray.AddSegment(HighEnergyExclusions, Pulses='OfflinePulses')

millipede_config = dict(Pulses='OfflinePulses', CascadePhotonicsService=pxs, PartialExclusion=False, BadDOMs=exclusions,    Parametrization='HalfSphere')
"""
tray.Add(millipede.MonopodFit, 'L5MonopodAmpFit', Seed='L5CscdLlh',
        PhotonsPerBin=-1, **millipede_config)
tray.Add(millipede.MonopodFit, 'L5MonopodPreFit', Seed='L5MonopodAmpFit',
        PhotonsPerBin=5, **millipede_config)
tray.AddSegment(millipede.MonopodFit, 'L5MonopodFit4', Seed='L5MonopodPreFit',
        PhotonsPerBin=5, Iterations=4, **millipede_config)

from simsegments import TrackVetoes
tray.Add(TrackVetoes, 'L5Offline',pulses="OfflinePulses",Vertex='L5MonopodFit4')

def printy(frame):
    if frame["IsHese"].value==True and frame["IsHESE_ck"].value==False and frame["CascadeFilter"].value==False:
        del frame["IsHese"]
        frame["IsHese"]=icetray.I3Bool(False)
    if frame["IsCascade"].value==True and frame["IsCascade_reco"].value==True and (frame["L5MonopodFit4"].pos.z<-500 or frame["L5MonopodFit4"].pos.z>500 or numpy.sqrt(frame["L5MonopodFit4"].pos.x*frame["L5MonopodFit4"].pos.x+frame["L5MonopodFit4"].pos.y*frame["L5MonopodFit4"].pos.y)>550):
        print("Cascade event outside of detector... ",  frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        del frame["IsCascade"]
        frame["IsCascade"]=icetray.I3Bool(False)
    if frame["IsCascade"].value==True and frame["IsCascade_reco"].value==False and (frame["L4VetoTrackMilliOfflineVetoCharge"].value>2 or frame['L4VetoTrackMarginMilliOfflineSide'].value<125):
        del frame["IsCascade"]
        frame["IsCascade"]=icetray.I3Bool(False)
    if frame["IsUpgoingMuon"].value==True:
        print("UpMu event: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        return True
    elif frame["IsHese"].value==True:
        print("HESE event: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        return True
    elif frame["IsCascade"].value==True and frame["L4VetoTrackOfflineVetoCharge"].value<6 and frame["L4VetoTrackL5OfflineVetoCharge"].value<2:
        print("Cascade event: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        return True
    elif frame["IsCascade"].value==True and frame["L4VetoTrackOfflineVetoCharge"].value<2 and frame["L4VetoTrackL5OfflineVetoCharge"].value<3:
        print("Cascade event: ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
        return True
    else:
        if ( ( frame['L4UpgoingTrackOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackOfflineVetoCharge'].value > \
                frame['L4VetoTrackOfflineVetoCharge'].value and frame['L4UpgoingTrackOfflineVetoChannels'].value > 3 ) or \
                ( frame['L4UpgoingTrackSplitVetoCharge'].value > 10 and frame['L4UpgoingTrackSplitVetoCharge'].value > \
                frame['L4VetoTrackSplitVetoCharge'].value and frame['L4UpgoingTrackSplitVetoChannels'].value > 3 )  or \
                ( frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > \
                frame['L4VetoTrackMilliOfflineVetoCharge'].value and frame['L4UpgoingTrackMilliOfflineVetoChannels'].value > 3 and \
                frame['L4UpgoingTrackSplitVetoCharge'].value > 6 ) )\
                and frame["TrackFit"].dir.zenith>1.5 and frame["MuonFilter"].value==True:# and frame["OnlineL2Filter"].value==True:
            del frame["IsUpgoingMuon"]
            del frame["IsCascade"]
            frame["IsUpgoingMuon"]=icetray.I3Bool(True)
            frame["IsCascade"]=icetray.I3Bool(False)
            print("UpMu event!", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
            return True
        else:
            return False
tray.Add(printy, streams=[icetray.I3Frame.Physics])

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
        if p.time < timeWindow.start or p.time > timeWindow.stop:
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
tray.AddModule(collectStats, "CollectStats")

def add_primary(frame):
    if "I3MCTree" in frame:
        def sanitize(particle):
            if particle is None:
                return dataclasses.I3Particle()
            else:
                return particle

        mcTree = frame["I3MCTree"]
        primary = None
        neutrino = None
        for p in mcTree:
            if mcTree.depth(p) != 0: continue

            if p.is_neutrino:
                if neutrino is None or p.energy > neutrino.energy:
                    neutrino = p

            if primary is None or p.energy > primary.energy:
                primary = p
        if "MCPrimary" not in frame:
            frame["MCPrimary"] = sanitize(primary)
tray.Add(add_primary)

tray.AddModule('I3Writer', 'writer',
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics ,  icetray.I3Frame.Stream('M'),  icetray.I3Frame.Stream('S')],
        filename=options.OUTPUT)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

os.system('bzip2 -f %s' % options.OUTPUT)

del tray

stop_time = time.asctime()

print('Started:', start_time)
print('Ended:', stop_time)
