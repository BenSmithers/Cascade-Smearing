#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v1/icetray-start
#METAPROJECT /data/user/nwandkowsky/tarballs/icerec.V04-11-02

# ./L5_nugen.py -o nue_lowe.i3

from I3Tray import *
from icecube import icetray, dataio, dataclasses
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

infiles=[options.GCD]
infiles.append(options.INPUT)
outfile = options.OUTPUT

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', FilenameList=infiles)
print "Reading input file...", len(infiles)

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

def check_fits(frame):
    # check MonopodFit
    if frame["L4MonopodFit"].pos.z<600 and frame["L4MonopodFit"].pos.z>-600 and numpy.sqrt(frame["L4MonopodFit"].pos.x*frame["L4MonopodFit"].pos.x+frame["L4MonopodFit"].pos.y*frame["L4MonopodFit"].pos.y)<600:
        return True
    else:
        return False
tray.Add(check_fits)

def final_filter(frame):
    del frame["IsUpgoingMuon"]
    del frame["IsCascade"]
    del frame["IsHese"]

    ###########
    ##  HESE
    ###########
    if (frame["L4VetoLayer0"].value+frame["L4VetoLayer1"].value==0 and frame["HomogenizedQTot"].value>6000.) or frame["IsHESE_ck"].value==True:
        print "HESE event! ",frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
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
        if frame["CascadeFilter"].value==True and (frame["L4VetoTrackOfflineVetoCharge"].value<2 and \
                inside_volume_mono_offline and inside_volume_milli_offline) and \
                (numpy.log10(frame["L4MonopodFit"].energy)-numpy.cos(frame["L4MonopodFit"].dir.zenith)>2.5 or \
                numpy.cos(frame["L4MonopodFit"].dir.zenith)<0.5):
                    frame["IsCascade"]=icetray.I3Bool(True)
            frame["IsUpgoingMuon"]=icetray.I3Bool(False)
            print "CASCADE (cascade): ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,frame['HomogenizedQTot'].value#, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
            return True
        else:
            frame["IsCascade"]=icetray.I3Bool(False)
    else:
        if (frame["MuonFilter"].value==True and frame["L4VetoTrackMilliOfflineVetoCharge"].value<2 and \
                frame["L4VetoTrackOfflineVetoCharge"].value<2 and inside_volume_milli_offline and inside_volume_mono_offline):# or\
                #(frame["L4VetoTrackMilliSplitVetoCharge"].value<2 and inside_volume_milli_split):
            frame["IsCascade"]=icetray.I3Bool(True)
            frame["IsUpgoingMuon"]=icetray.I3Bool(False)
            print "CASCADE (track): ", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id,frame['HomogenizedQTot'].value#, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
            return True
        else:
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
            and frame["TrackFit"].dir.zenith>1.5 and frame["MuonFilter"].value==True and frame["OnlineL2Filter"].value==True:
                frame["IsUpgoingMuon"]=icetray.I3Bool(True)
        print "Upgoing Muon event!", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id
    else:
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)
    return frame["IsUpgoingMuon"].value==True or frame["IsHese"].value==True or frame["IsCascade"].value==True
tray.Add(final_filter)

def remove_coincident(frame):
    distance = numpy.sqrt( (frame["SPEFit4_split"].pos.x-frame["SplineMPE_split"].pos.x)**2+(frame["SPEFit4_split"].pos.y-frame["SplineMPE_split"].pos.y)**2 +(frame["SPEFit4_split"].pos.z-frame["SplineMPE_split"].pos.z)**2)
    frame["SPEFit4Distance"]=dataclasses.I3Double(distance)
    distance = numpy.sqrt( (frame["SPEFitSingle_split"].pos.x-frame["SplineMPE_split"].pos.x)**2+(frame["SPEFitSingle_split"].pos.y-frame["SplineMPE_split"].pos.y)**2 +(frame["SPEFitSingle_split"].pos.z-frame["SplineMPE_split"].pos.z)**2)
    frame["SPEFitSingleDistance"]=dataclasses.I3Double(distance)
    not_coinc = frame["SPEFit4Distance"].value<200 and frame["SPEFitSingleDistance"].value<200
    if frame["IsUpgoingMuon"].value==True and not_coinc==False:
        return False

    if frame["L4MonopodFit_AvgDistQ"].value>150. and frame["TrackFit_AvgDistQ"].value>110:
        return False
    else:
        return True
tray.Add(remove_coincident)

from icecube import VHESelfVeto
tray.AddModule('VertexInFiducialVolume', 'vertex_inside')
tray.AddModule('TauGeneratesMuon', 'muon_from_tau')

def check_type(frame):
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

    if "MCPrimary" in frame:
        del frame["IsCascade_true"]
        del frame["IsTrack_true"]
        if frame["MCPrimary"].type==12 or frame["MCPrimary"].type==-12:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)

        if (frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14) and frame["I3MCWeightDict"]["InteractionType"]==2.0:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)
        elif (frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14) and "VertexOutgoingLepton" in frame:
            if(frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14) and frame["I3MCWeightDict"]["InteractionType"]==1.0 and (frame["VertexOutgoingLepton"].type==13 or frame["VertexOutgoingLepton"].type==-13):
                frame["IsCascade_true"]= icetray.I3Bool(False)
                frame["IsTrack_true"]= icetray.I3Bool(True)
            else:
                frame["IsCascade_true"]= icetray.I3Bool(True)
                frame["IsTrack_true"]= icetray.I3Bool(False)
        elif (frame["MCPrimary"].type==14 or frame["MCPrimary"].type==-14):
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)

        if (frame["MCPrimary"].type==16 or frame["MCPrimary"].type==-16) and frame["I3MCWeightDict"]["InteractionType"]==2.0:
            frame["IsCascade_true"]= icetray.I3Bool(True)
            frame["IsTrack_true"]= icetray.I3Bool(False)
        elif  (frame["MCPrimary"].type==16 or frame["MCPrimary"].type==-16) and frame["I3MCWeightDict"]["InteractionType"]==1.0:
            frame["IsCascade_true"]= icetray.I3Bool(False)
            frame["IsTrack_true"]= icetray.I3Bool(True)

    if frame["IsHese"].value==False:
        del frame["IsCascade_reco"]
        del frame["IsTrack_reco"]
        if frame["L4StartingTrackHLCSplitVetoCharge"].value>2.5:
            frame["IsTrack_reco"]=icetray.I3Bool(True)
            frame["IsCascade_reco"]=icetray.I3Bool(False)
        else:
            frame["IsTrack_reco"]=icetray.I3Bool(False)
            frame["IsCascade_reco"]=icetray.I3Bool(True)
    else:
        if frame["L4StartingTrackHLCSplitVetoCharge"].value<1.:
            del frame["IsCascade_reco"]
            del frame["IsTrack_reco"]
            frame["IsTrack_reco"]=icetray.I3Bool(False)
            frame["IsCascade_reco"]=icetray.I3Bool(True)
        elif frame["L4StartingTrackHLCSplitVetoCharge"].value>100.:
            del frame["IsCascade_reco"]
            del frame["IsTrack_reco"]
            frame["IsTrack_reco"]=icetray.I3Bool(True)
            frame["IsCascade_reco"]=icetray.I3Bool(False)
tray.Add(check_type)

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

tray.Add(millipede.MonopodFit, 'L5MonopodAmpFit', Seed='L5CscdLlh',
        PhotonsPerBin=-1, **millipede_config)
tray.Add(millipede.MonopodFit, 'L5MonopodPreFit', Seed='L5MonopodAmpFit',
        PhotonsPerBin=5, **millipede_config)
tray.AddSegment(millipede.MonopodFit, 'L5MonopodFit4', Seed='L5MonopodPreFit',
        PhotonsPerBin=5, Iterations=4, **millipede_config)

tray.AddModule('I3Writer', 'writer',
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[  icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
        filename=outfile)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

os.system('bzip2 -f %s' % outfile)

del tray

stop_time = time.asctime()

print 'Started:', start_time
print 'Ended:', stop_time
