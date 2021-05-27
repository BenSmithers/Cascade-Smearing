#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/combo/py3-v4.1.1/

##!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v1/icetray-start
##METAPROJECT /data/user/nwandkowsky/tarballs/icerec.V04-11-02

# ./L4.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz -i /data/ana/Cscd/StartingEvents/NuGen/NuE/medium_energy/IC86_2011/l3_new/1/l3_00000001.i3.bz2 -o 2011_L4_nugen.i3

# ./L4.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2012.56063_V1.i3.gz -i 2012_nugen_L3.i3.bz2 -o 2012_nugen_L4.i3

# ./L4.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2013.56429_V1.i3.gz -i 2013_nugen_L3.i3.bz2 -o 2013_nugen_L4.i3

from I3Tray import *
from icecube import icetray, phys_services, dataio, dataclasses
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
parser.add_option("-r", "--runnumber", type="int", default=1, dest="RUNNUMBER", help="The run number for this simulation")
# parse cmd line args, bail out if anything is not understood
(options,args) = parser.parse_args()

outfile = options.OUTPUT

tray = I3Tray()

tray.AddModule('I3Reader', 'reader', FilenameList=[options.GCD, options.INPUT])
print("Reading input file...", options.INPUT)

def fix_header(frame):
    header = frame['I3EventHeader']
    frame['I3EventHeader_original']=header
    header.run_id = int(options.RUNNUMBER)
    del frame['I3EventHeader']
    frame['I3EventHeader'] = header
tray.Add(fix_header)

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

############################
## DOWNGOING TRACK VETO
############################
from icecube import linefit, improvedLinefit, lilliput
tray.AddSegment( improvedLinefit.simple, 'LineFit_offline', inputResponse = 'TWSRTOfflinePulses', fitName = 'LineFit_offline' )
tray.AddSegment( lilliput.I3SinglePandelFitter, 'SPEFitSingle_offline', pulses = 'TWSRTOfflinePulses', seeds = ['LineFit_offline'])
tray.AddSegment( lilliput.I3IterativePandelFitter, 'SPEFit2_offline', pulses = 'TWSRTOfflinePulses', n_iterations = 2, seeds = [ 'SPEFitSingle_offline'])

from icecube import DomTools
tray.AddModule('I3TimeWindowCleaning<I3RecoPulse>', 'TimeWindowHLC',
        InputResponse='OfflinePulsesHLC', OutputResponse='TWOfflinePulsesHLC',
        TimeWindow = 6000*I3Units.ns)

from icecube.CascadeL3_IC79.level3.reco import AxialCscdLlhFitter
tray.Add(AxialCscdLlhFitter, "L4CascadeLlh", Pulses="TWOfflinePulsesHLC", Seed="SPEFit2_offline")

from icecube.photonics_service import I3PhotoSplineService
base = os.path.expandvars('$I3_DATA/photon-tables/splines/ems_mie_z20_a10.%s.fits')
pxs = I3PhotoSplineService(base % "abs", base % "prob", 0)

from icecube.wavedeform import AddMissingTimeWindow
#tray.AddModule(AddMissingTimeWindow, 'OfflinePulsesTimeRange', Pulses='OfflinePulses')

from icecube import  millipede
exclusions = tray.AddSegment(millipede.HighEnergyExclusions, 'monopod_DOM_exclusions',
        Pulses='OfflinePulses',
        CalibrationErrata='CalibrationErrata',
        ExcludeSaturatedDOMs='OfflineInIceCalibrationErrata',
        ExcludeBrightDOMs='BrightDOMs',
        ExcludeDeepCore=False,
        BadDomsList='BadDomsList',
        SaturationWindows='SaturationTimes'
        )
config = dict(Pulses="OfflinePulses", PartialExclusion=False, CascadePhotonicsService=pxs, BadDOMs=exclusions, Parametrization='HalfSphere', MinTimeWidth=15)
tray.Add(millipede.MonopodFit, 'L4MonopodFit', Seed="L4CascadeLlh", BinSigma=2, **config)

tray.AddModule('I3LCPulseCleaning', 'lcclean2',
        Input='SplitPulses',
        OutputHLC='SplitPulsesHLC',
        OutputSLC='SplitPulsesSLC')

from simsegments import TrackVetoes
tray.Add(TrackVetoes, 'Offline',pulses="OfflinePulses",Vertex='L4MonopodFit')
tray.Add(TrackVetoes, 'Split',pulses="SplitPulses",Vertex='L4MonopodFit')


################################
## ALL THE TRACK FIT GOODIES
################################

def selectFinalFit(frame):
    fits=[ "SplineMPE_split", "MPEFit_split","SPEFit4_split","SPEFitSingle_split","ImprovedLineFit_split"]
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
            if(frame.Has(fitName+"FitParams")):
                frame.Put(resultName+"FitParams",frame.Get(fitName+"FitParams"))
            break
    if(not frame.Has(resultName+"Fallback")):
        frame.Put(resultName+"Fallback",icetray.I3Int(len(fits)))

#fetch masks or pulses transparently
def getRecoPulses(frame,name):
    pulses = frame[name]
    if pulses.__class__ == dataclasses.I3RecoPulseSeriesMapMask:
        pulses = pulses.apply(frame)
    return pulses

from icecube import improvedLinefit, lilliput
tray.AddSegment(improvedLinefit.simple,"ImprovedLineFit_split",
        inputResponse="SplitPulses",
        fitName="ImprovedLineFit_split"
        )

tray.AddSegment(lilliput.I3SinglePandelFitter,"SPEFitSingle_split",
        domllh="SPE1st",
        pulses="SplitPulses",
        seeds=["ImprovedLineFit_split"]
        )

tray.AddSegment(lilliput.I3IterativePandelFitter,"SPEFit4_split",
        domllh="SPE1st",
        n_iterations=4,
        pulses="SplitPulses",
        seeds=["SPEFitSingle_split","ImprovedLineFit_split"]
        )

def check_fits(frame):
    # check SPEFitSingle
    if frame["SPEFitSingle_split"].pos.z<500 and frame["SPEFitSingle_split"].pos.z>-500 and numpy.sqrt(frame["SPEFitSingle_split"].pos.x*frame["SPEFitSingle_split"].pos.x+frame["SPEFitSingle_split"].pos.y*frame["SPEFitSingle_split"].pos.y)<550:
        spe_is_inside = True
    else:
        spe_is_inside = False
        if  frame['HomogenizedQTot'].value<6000:
            return False
    # check SPEFit4
    if frame["SPEFit4_split"].pos.z<500 and frame["SPEFit4_split"].pos.z>-500 and numpy.sqrt(frame["SPEFit4_split"].pos.x*frame["SPEFit4_split"].pos.x+frame["SPEFit4_split"].pos.y*frame["SPEFit4_split"].pos.y)<550:
        spe_is_inside = True
    else:
        spe_is_inside = False
        if  frame['HomogenizedQTot'].value<6000:
            return False
    # check MonopodFit
    if frame["L4MonopodFit"].pos.z<500 and frame["L4MonopodFit"].pos.z>-500 and numpy.sqrt(frame["L4MonopodFit"].pos.x*frame["L4MonopodFit"].pos.x+frame["L4MonopodFit"].pos.y*frame["L4MonopodFit"].pos.y)<550:
        monopod_is_inside = True
    else:
        monopod_is_inside = False
        if  frame['HomogenizedQTot'].value<6000:
            return False
    if frame["L4MonopodFit"].energy == 0.0 or numpy.isnan(frame["L4MonopodFit"].energy):
        #print "Bad reco event..."
        return False
tray.Add(check_fits)

tray.AddSegment(lilliput.I3SinglePandelFitter,"MPEFit_split",
        domllh="MPE",
        pulses="SplitPulses",
        seeds=["SPEFit4_split","SPEFitSingle_split","ImprovedLineFit_split"]
        )

load('mue')
tray.AddModule("muex", "muex_angular4",
        Pulses="SplitPulses",
        rectrk="",
        result="MuEXAngular4",
        lcspan=0,
        repeat=4,
        usempe=True,
        detail=False,
        energy=False,
        icedir=os.path.expandvars("$I3_BUILD/mue/resources/ice/mie"))

# spline MPE

infmuonprobsplinepath = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits"
infmuonampsplinepath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'

from icecube import photonics_service
spline_mie=photonics_service.I3PhotoSplineService(infmuonampsplinepath,infmuonprobsplinepath, 4)
llh="MPE"

tray.AddService("I3GulliverMinuitFactory", "Minuit", Algorithm="SIMPLEX", MaxIterations=1000, Tolerance=0.01)

tray.AddService("I3SimpleParametrizationFactory", "SimpleTrack",
        StepX = 20*I3Units.m,
        StepY = 20*I3Units.m,
        StepZ = 20*I3Units.m,
        StepZenith = 0.1*I3Units.radian,
        StepAzimuth= 0.2*I3Units.radian,
        BoundsX = [-2000*I3Units.m, 2000*I3Units.m],
        BoundsY = [-2000*I3Units.m, 2000*I3Units.m],
        BoundsZ = [-2000*I3Units.m, 2000*I3Units.m])

tray.AddService("I3BasicSeedServiceFactory", "SplineSeedMPE",
        FirstGuesses=["MuEXAngular4"])

from icecube import spline_reco
tray.AddService("I3SplineRecoLikelihoodFactory","LLHSplineMPE",
        PhotonicsService=spline_mie,
        Pulses="SplitPulses",
        Likelihood=llh,
        NoiseRate=10*I3Units.hertz)

tray.AddModule("I3SimpleFitter", "SplineMPE_split",
        SeedService="SplineSeedMPE",
        Parametrization="SimpleTrack",
        LogLikelihood="LLHSplineMPE",
        Minimizer="Minuit")

tray.AddModule(selectFinalFit,"FinalFit")

tray.AddModule("muex", "muex_differential",
        Pulses = "SplitPulses",
        rectrk = "TrackFit",
        result = "SplineMPEMuEXDifferential",
        detail = True, # differential
        energy = True,
        lcspan = 0,
        icedir=os.path.expandvars("$I3_BUILD/mue/resources/ice/mie"))

def CleanDeepCore(frame):
    mask = dataclasses.I3RecoPulseSeriesMapMask(frame, "SplitPulses", lambda om, idx, pulse: om.string <= 78)
    frame["SplitPulses_NoDC"] = mask
tray.AddModule(CleanDeepCore, 'nodc_again')

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


tray.AddModule(ComputeChargeWeightedDist,"CCWD_Track",Pulses="SplitPulses_NoDC",Track="TrackFit")
tray.AddModule(ComputeChargeWeightedDist,"CCWD_Cascade",Pulses="SplitPulses_NoDC",Track="L4MonopodFit")

def remove_coinc(frame):
    if frame["L4MonopodFit_AvgDistQ"].value>150. and frame["TrackFit_AvgDistQ"].value>110:
        return False
    else:
        return True
tray.Add(remove_coinc)

# MILLIPEDE

muon_service = photonics_service.I3PhotoSplineService('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits', '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits', 0, 400.0)
cascade_service_mie = photonics_service.I3PhotoSplineService('/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/LowEnergyCorrectedCascades_z20_a10_HE.abs.fits', '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/LowEnergyCorrectedCascades_z20_a10_HE.prob.fits', 0.0, 400.0)

exclusionList = \
        tray.AddSegment(millipede.HighEnergyExclusions, 'millipede_DOM_exclusions',
                Pulses = "SplitPulses",
                ExcludeDeepCore='DeepCoreDOMs',
                ExcludeSaturatedDOMs='SaturatedDOMs',
                ExcludeBrightDOMs='BrightDOMs_milli',
                BadDomsList='BadDomsList',
                CalibrationErrata='CalibrationErrata',
                SaturationWindows='SaturationTimes'
                )

common_millipede_params = dict(
        ## params for I3MillipedeBase<ConditionalModule>
        MuonPhotonicsService=muon_service,
        CascadePhotonicsService=cascade_service_mie,
        ShowerRegularization=1e-9,
        PhotonsPerBin=15,
        DOMEfficiency=0.99,
        ExcludedDOMs=exclusionList,
        PartialExclusion=False,
        ReadoutWindow='SplitPulsesTimeRange',
        Pulses='SplitPulses',
        )

tray.AddModule(AddMissingTimeWindow, 'SplitPulsesTimeRange', Pulses='SplitPulses')
tray.AddModule('MuMillipede', 'mumillipede',
        ## params for MuMillipede
        SeedTrack="TrackFit",
        MuonSpacing=0.*I3Units.m,
        ShowerSpacing=2.5*I3Units.m,
        Output="MuonLosses",
        ## common parameters for all millipede modules
        **common_millipede_params
        )

global i
i=1
def count(frame):
    if "MuonLosses" in frame:
        global i
        loss_vec = frame["MuonLosses"]
        for i in range(0,len(loss_vec)):
            if loss_vec[i].energy>1. and loss_vec[i].pos.z<500 and loss_vec[i].pos.z>-500 and math.sqrt(loss_vec[i].pos.x*loss_vec[i].pos.x+loss_vec[i].pos.y*loss_vec[i].pos.y)<550:
                frame["MillipedeFirstLoss"] = loss_vec[i]
                break
        if 'MillipedeFirstLoss' not in frame:
            frame["MillipedeFirstLoss"] = loss_vec[0]
        i = i+1
    else:
        part = dataclasses.I3Particle()
        part.pos = dataclasses.I3Position(1950.,1950.,1950.)
        frame["MillipedeFirstLoss"]=part

tray.AddModule(count,"mycounter")

from icecube.CascadeL3_IC79.level4.veto import VetoMarginCalculator
tray.AddModule(VetoMarginCalculator, 'MillipedeBorderDistance', Vertex="MillipedeFirstLoss")

def millipede_loss(frame):
    if "MuonLosses" in frame:
        loss_vec = frame["MuonLosses"]
        deposited_energy = 0
        for i in range(0,len(loss_vec)):
            deposited_energy = deposited_energy + loss_vec[i].energy
    else:
        deposited_energy = 0.
    frame["MillipedeDepositedEnergy"] = dataclasses.I3Double(deposited_energy)
tray.AddModule(millipede_loss,"milli_loss")

tray.Add(TrackVetoes, 'MilliSplit', pulses="SplitPulses",Vertex='MillipedeFirstLoss')
tray.Add(TrackVetoes, 'MilliOffline', pulses="OfflinePulses",Vertex='MillipedeFirstLoss')

def final_filter(frame):
    #if 'MuonWeight_GaisserH4a' not in frame:
    #    frame['MuonWeight_GaisserH4a']=dataclasses.I3Double(0.)

    ###########
    ##  HESE
    ###########
    if (frame["L4VetoLayer0"].value+frame["L4VetoLayer1"].value<3 and frame["HomogenizedQTot"].value>6000.) or frame["IsHESE_ck"].value==True:
        frame["IsHese"]=icetray.I3Bool(True)
        frame["IsCascade"]=icetray.I3Bool(False)
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)
        print(frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, " HESE") #, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
        return True
    else:
        frame["IsHese"]=icetray.I3Bool(False)

    ###########
    ##  CASCADES
    ###########

    if frame["L4VetoTrackOfflineVetoCharge"].value<2 or frame["L4VetoTrackSplitVetoCharge"].value<2  or frame["L4VetoTrackMilliOfflineVetoCharge"].value<2 or frame["L4VetoTrackMilliSplitVetoCharge"].value<2:
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)
        frame["IsCascade"]=icetray.I3Bool(True)
        print(frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id, " CASCADE",frame['HomogenizedQTot'].value)#, frame['MuonWeight_GaisserH4a'].value*3600*24*340/160000
        return True
    else:
        frame["IsCascade"]=icetray.I3Bool(False)


    ###########
    ##  UPMU
    ###########
    if ( frame['L4UpgoingTrackOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackOfflineVetoCharge'].value > \
            frame['L4VetoTrackOfflineVetoCharge'].value and frame['L4UpgoingTrackOfflineVetoChannels'].value > 3 ) or \
            ( frame['L4UpgoingTrackSplitVetoCharge'].value > 10 and frame['L4UpgoingTrackSplitVetoCharge'].value > \
            frame['L4VetoTrackSplitVetoCharge'].value and frame['L4UpgoingTrackSplitVetoChannels'].value > 3 ) or \
            ( frame['L4UpgoingTrackMilliSplitVetoCharge'].value > 10 and frame['L4UpgoingTrackMilliSplitVetoCharge'].value > \
            frame['L4VetoTrackMilliSplitVetoCharge'].value and frame['L4UpgoingTrackMilliSplitVetoChannels'].value > 3 or \
            ( frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > 10 and frame['L4UpgoingTrackMilliOfflineVetoCharge'].value > \
            frame['L4VetoTrackMilliOfflineVetoCharge'].value and frame['L4UpgoingTrackMilliOfflineVetoChannels'].value > 3 )):
        frame["IsUpgoingMuon"]=icetray.I3Bool(True)
        print("Upgoing Muon event!", frame["I3EventHeader"].run_id, frame["I3EventHeader"].event_id)
    else:
        frame["IsUpgoingMuon"]=icetray.I3Bool(False)

    return frame["IsUpgoingMuon"].value==True or frame["IsHese"].value==True or frame["IsCascade"].value==True
tray.Add(final_filter)

def check_type(frame):
    # neutrino types: 12 nue, 14 numu, 16 nutau
    # interaction types: 1 CC, 2 NC, 3 GR
    if "MCPrimary" in frame:
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
tray.Add(check_type)

def discriminate(frame):
    # if frame["IsUpgoingMuon"].value==True and frame["TrackFit"].dir.zenith<1.6:
    #    return False
    if frame["IsUpgoingMuon"].value==True:
        frame["IsCascade_reco"]=icetray.I3Bool(False)
        frame["IsTrack_reco"]=icetray.I3Bool(True)
        return True

    if frame["L4MonopodFit"].energy<=1e4 and (frame["L4StartingTrackHLCOfflineVetoCharge"].value>10 or (frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>5 and frame["L4StartingTrackHLCOfflineVetoCharge"].value>5) ):
        frame["IsCascade_reco"]=icetray.I3Bool(False)
        frame["IsTrack_reco"]=icetray.I3Bool(True)
        return True

    if frame["L4MonopodFit"].energy<=1e4:
        frame["IsCascade_reco"]=icetray.I3Bool(True)
        frame["IsTrack_reco"]=icetray.I3Bool(False)
        return True
    # print "Is track reco"

    if frame["L4MonopodFit"].energy>1e4 and frame["HomogenizedQTot"].value<6000. and (frame["L4MonopodFit_AvgDistQ"].value-25)*0.9<frame["TrackFit_AvgDistQ"].value and frame["L4StartingTrackHLCMilliOfflineVetoCharge"].value>1000 and frame["L4StartingTrackHLCOfflineVetoCharge"].value>3:
        frame["IsCascade_reco"]=icetray.I3Bool(False)
        frame["IsTrack_reco"]=icetray.I3Bool(True)
        return True

    if frame["L4MonopodFit"].energy>1e4 and frame["HomogenizedQTot"].value<6000.:
        frame["IsCascade_reco"]=icetray.I3Bool(True)
        frame["IsTrack_reco"]=icetray.I3Bool(False)
        return True

    if frame["HomogenizedQTot"].value>6000. and (frame["L4MonopodFit_AvgDistQ"].value-25)*0.9<frame["TrackFit_AvgDistQ"].value:
        frame["IsCascade_reco"]=icetray.I3Bool(True)
        frame["IsTrack_reco"]=icetray.I3Bool(False)
    else:
        frame["IsCascade_reco"]=icetray.I3Bool(False)
        frame["IsTrack_reco"]=icetray.I3Bool(True)

tray.Add(discriminate)


tray.Add('Delete', "final_delete", Keys=['Cuts', 'InIceRawData', 'PoleEHESummaryPulseInfo', 'CascadeFilt_CscdLlh', 'CascadeFilt_LFVel', 'CascadeFilt_ToiVal', 'CleanTriggerHierarchy_IT', 'ClusterCleaningExcludedStations', 'IceTopRawData', 'JEBEventInfo', 'OfflineIceTopHLCPulseInfo', 'OfflineIceTopHLCTankPulses', 'OfflineIceTopHLCVEMPulses', 'OfflineIceTopSLCVEMPulses', 'TankPulseMergerExcludedStations'])


tray.AddModule('I3Writer', 'writer',
        DropOrphanStreams=[icetray.I3Frame.DAQ],
        Streams=[ icetray.I3Frame.Stream('S'), icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('M')],
        filename=outfile)

tray.AddModule('TrashCan', 'thecan')

tray.Execute()
tray.Finish()

os.system('bzip2 -f %s' % outfile)

del tray

stop_time = time.asctime()

print('Started:', start_time)
print('Ended:', stop_time)
