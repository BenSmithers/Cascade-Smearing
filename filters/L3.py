#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT gsiftp://gridftp-users.icecube.wisc.edu/data/user/nwandkowsky/tarballs/icerec.V05-01-07

# ./L3.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz -i /data/ana/Cscd/StartingEvents/NuGen_new/NuE/medium_energy/IC86_flasher_p1=0.3_p2=0.0/l2/1/l2_00000001.i3.zst -o L3_nugen.i3.zst

# ./L3.py -g /cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_2016.57531_V0.i3.gz -i L2_muongun.i3.zst -o L3_muongun.i3.zst

# ./L3.py -g /data/exp/IceCube/2011/filtered/level2pass2/0823/Run00118595/Level2pass2_IC86.2011_data_Run00118595_0823_4_20_GCD.i3.zst -i /data/exp/IceCube/2011/filtered/level2pass2/0823/Run00118595/Level2pass2_IC86.2011_data_Run00118595_Subrun00000000_00000000.i3.zst -o L3_data.i3.zst

from I3Tray import *
from icecube import icetray, phys_services, dataio, dataclasses
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

infile = "gsiftp://gridftp-users.icecube.wisc.edu"+options.INPUT
outfile = "gsiftp://gridftp-users.icecube.wisc.edu"+options.OUTPUT

tray = I3Tray()

tray.AddSegment(dataio.I3Reader, 'reader', FilenameList=[options.GCD, infile])

#####################################
## INITIAL CLEANING
#####################################

def streamcut(frame):
    return frame["I3EventHeader"].sub_event_stream!='IceTopSplit' and frame["I3EventHeader"].sub_event_stream!='NullSplit'
tray.Add(streamcut)

tray.Add('Delete','initial_clean',Keys=['BeaconLaunches',
    'CVMultiplicity', 'CVStatistics', 'CalibratedSLC',
    'CascadeContainmentTagging_L2', 'CascadeContainmentTagging_Singles_L2',
    'CascadeDipoleFit_L2', 'CascadeDipoleFit_L2Params', 'CascadeFillRatio_L2',
    'CascadeImprovedLineFit_L2', 'CascadeImprovedLineFit_L2Params', 'CascadeLast_IC_Singles_L2',
    'CascadeLast_IC_Singles_L2Params', 'CascadeLast_L2', 'CascadeLast_L2Params',
    'CascadeLineFitSplit1_L2', 'CascadeLineFitSplit1_L2Params', 'CascadeLineFitSplit2_L2', 'CascadeLineFitSplit2_L2Params',
    'CascadeLineFit_L2', 'CascadeLineFit_L2Params', 'CascadeLlhVertexFitSplit1_L2', 'CascadeLlhVertexFitSplit1_L2Params',
    'CascadeLlhVertexFitSplit2_L2', 'CascadeLlhVertexFitSplit2_L2Params',
    'CascadeLlhVertexFit_IC_Singles_L2', 'CascadeLlhVertexFit_IC_Singles_L2Params',
    'CascadeLlhVertexFit_L2', 'CascadeLlhVertexFit_L2Params',
    'CascadeSplitPulses_L21', 'CascadeSplitPulses_L22',
    'CascadeToISplit1_L2', 'CascadeToISplit1_L2Params', 'CascadeToISplit2_L2', 'CascadeToISplit2_L2Params',
    'CorsikaMoonMJD', 'CorsikaSunMJD',
    'CscdL2_Topo1stPulse_HLC0', 'CscdL2_Topo1stPulse_HLCSplitCount', 'CscdL2_Topo_HLCFirstPulses',
    'EHEATWDPortiaPulseSRT', 'EHEATWDPulseSeriesSRT', 'EHEBestPortiaPulseSRT',
    'EHEDSTShieldParameters_ImpLF', 'EHEDSTShieldParameters_SPE12', 'EHEFADCPortiaPulseSRT',
    'EHEFADCPulseSeriesSRT', 'EHEOpheliaBTWSRT', 'EHEOpheliaParticleBTWSRT', 'EHEOpheliaParticleSRT',
    'EHEOpheliaParticleSRT_ImpLF', 'EHEOpheliaSRT', 'EHEOpheliaSRT_ImpLF', 'EHEPortiaEventSummary',
    'EHEPortiaEventSummarySRT',
    'Estres_CausalQTot', 'Estres_Homogenized_QTot',
    'FilterMask_NullSplit0',
    'FiniteRecoCuts', 'FiniteRecoFit', 'FiniteRecoLlh', 'Homogenized_QTot', 'HuberFit',
    'InIceRawData', 'LargestOMKey', 'LineFit', 'LineFitEHE', 'LineFitEHEParams', 'LineFitParams',
    'MPEFit', 'MPEFitCharacteristics', 'MPEFitCramerRaoParams', 'MPEFitFitParams', 'MPEFitMuEX',
    'OfflinePulsesHLC', 'OfflinePulsesSLC',
    'OnlineL2_BayesianFit', 'OnlineL2_BayesianFitFitParams',
    'OnlineL2_BestFit', 'OnlineL2_BestFitFitParams', 'OnlineL2_BestFit_CramerRaoParams',
    'OnlineL2_BestFit_CramerRao_cr_azimuth', 'OnlineL2_BestFit_CramerRao_cr_zenith',
    'OnlineL2_BestFit_DirectHitsA', 'OnlineL2_BestFit_DirectHitsB', 'OnlineL2_BestFit_DirectHitsC',
    'OnlineL2_BestFit_DirectHitsD', 'OnlineL2_BestFit_DirectHitsE', 'OnlineL2_BestFit_MuEx', 'OnlineL2_BestFit_MuEx_r',
    'OnlineL2_BestFit_Name', 'OnlineL2_CleanedMuonPulses', 'OnlineL2_HitMultiplicityValues', 'OnlineL2_HitMultiplicityValuesIC',
    'OnlineL2_HitStatisticsValues', 'OnlineL2_HitStatisticsValuesIC', 'OnlineL2_MPEFit', 'OnlineL2_MPEFitFitParams',
    'OnlineL2_SPE2itFit', 'OnlineL2_SPE2itFitFitParams', 'OnlineL2_SplineMPE', 'OnlineL2_SplineMPEFitParams',
    'OnlineL2_SplineMPE_Characteristics', 'OnlineL2_SplineMPE_CharacteristicsIC', 'OnlineL2_SplineMPE_CharacteristicsNoRCut',
    'OnlineL2_SplineMPE_CharacteristicsNoRCutIC', 'OnlineL2_SplineMPE_CramerRaoParams',
    'OnlineL2_SplineMPE_CramerRao_cr_azimuth', 'OnlineL2_SplineMPE_CramerRao_cr_zenith',
    'OnlineL2_SplineMPE_DirectHitsA', 'OnlineL2_SplineMPE_DirectHitsB',
    'OnlineL2_SplineMPE_DirectHitsC', 'OnlineL2_SplineMPE_DirectHitsD', 'OnlineL2_SplineMPE_DirectHitsE',
    'OnlineL2_SplineMPE_DirectHitsICA', 'OnlineL2_SplineMPE_DirectHitsICB', 'OnlineL2_SplineMPE_DirectHitsICC',
    'OnlineL2_SplineMPE_DirectHitsICD', 'OnlineL2_SplineMPE_DirectHitsICE', 'OnlineL2_SplineMPE_MuE',
    'OnlineL2_SplineMPE_MuEx', 'OnlineL2_SplineMPE_MuEx_r', 'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_MuEres',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_Muon', 'OnlineL2_SplineMPE_TruncatedEnergy_AllBINS_Neutrino',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_MuEres', 'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_AllDOMS_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_BINS_MuEres',
    'OnlineL2_SplineMPE_TruncatedEnergy_BINS_Muon', 'OnlineL2_SplineMPE_TruncatedEnergy_BINS_Neutrino',
    'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_MuEres', 'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_DOMS_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Muon',
    'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_Neutrino', 'OnlineL2_SplineMPE_TruncatedEnergy_ORIG_dEdX',
    'OnlineL2_SplitGeo1_BayesianFit', 'OnlineL2_SplitGeo1_BayesianFitFitParams',
    'OnlineL2_SplitGeo1_Linefit', 'OnlineL2_SplitGeo1_SPE2itFit', 'OnlineL2_SplitGeo1_SPE2itFitFitParams',
'OnlineL2_SplitGeo2_BayesianFit', 'OnlineL2_SplitGeo2_BayesianFitFitParams', 'OnlineL2_SplitGeo2_Linefit',
    'OnlineL2_SplitGeo2_SPE2itFit', 'OnlineL2_SplitGeo2_SPE2itFitFitParams', 'OnlineL2_SplitTime1_BayesianFit',
    'OnlineL2_SplitTime1_BayesianFitFitParams', 'OnlineL2_SplitTime1_Linefit',
    'OnlineL2_SplitTime1_SPE2itFit', 'OnlineL2_SplitTime1_SPE2itFitFitParams',
    'OnlineL2_SplitTime2_BayesianFit', 'OnlineL2_SplitTime2_BayesianFitFitParams',
    'OnlineL2_SplitTime2_Linefit', 'OnlineL2_SplitTime2_SPE2itFit', 'OnlineL2_SplitTime2_SPE2itFitFitParams',
    'PassedAnyFilter', 'PassedConventional', 'PassedKeepSuperDSTOnly',
    'PoleCascadeFilter_CscdLlh', 'PoleEHESummaryPulseInfo', 'PoleMuonLinefit', 'PoleMuonLlhFit', 'PoleMuonLlhFitFitParams',
    'RTTWOfflinePulses_FR_WIMP', 'SPEFit12EHE', 'SPEFit12EHEFitParams',
    'SPEFit2', 'SPEFit2Characteristics', 'SPEFit2CramerRaoParams', 'SPEFit2FitParams', 'SPEFit2MuEX_FSS',
    'SPEFitSingle', 'SPEFitSingleEHE', 'SPEFitSingleEHEFitParams', 'SPEFitSingleFitParams', 'SRTInIcePulses',
    'SRTInIcePulses_IC_Singles_L2CleanedKeys', 'SRTInIcePulses_WODCCleanedKeys',
    'SplineMPE_Estres', 'SplineMPE_EstresFitParams',
    'TWOfflinePulsesHLC', 'TWOfflinePulses_FR_WIMP', 'TWOfflinePulses_FR_WIMPTimeRange',
  'UncleanedInIcePulsesTimeRange', 'WIMPrecoTopoSplitSplitCount',
  'splittedDOMMap','splittedDOMMapSRT', 'I3MCTree_preMuonProp'    ])

def filter_cut(frame):
    NPassedFilters = 0

    if frame["FilterMask"]["CascadeFilter_13"].condition_passed:
        NPassedFilters = NPassedFilters+1

    if frame["FilterMask"]["HighQFilter_17"].condition_passed:
        NPassedFilters = NPassedFilters+1

    if frame["FilterMask"]["MuonFilter_13"].condition_passed:
        NPassedFilters = NPassedFilters+1

    if frame["FilterMask"]["MESEFilter_15"].condition_passed:
        NPassedFilters = NPassedFilters+1

    if frame["FilterMask"]["HESEFilter_15"].condition_passed:
        NPassedFilters = NPassedFilters+1

    return NPassedFilters>0
tray.Add(filter_cut, streams=[icetray.I3Frame.Physics])

#####################################
## TOPOLOGICAL SPLITTER
#####################################

from icecube.STTools.seededRT.configuration_services import I3DOMLinkSeededRTConfigurationService
seededRTConfig = I3DOMLinkSeededRTConfigurationService(
        ic_ic_RTRadius              = 150.0*I3Units.m,
        ic_ic_RTTime                = 1000.0*I3Units.ns,
        treat_string_36_as_deepcore = False,
        useDustlayerCorrection      = False,
        allowSelfCoincidence        = True
        )

tray.AddModule('I3SeededRTCleaning_RecoPulseMask_Module', 'seededrt',
        InputHitSeriesMapName  = "SplitInIcePulses",
        OutputHitSeriesMapName = "SplitSRTInIcePulses",
        STConfigService        = seededRTConfig,
        SeedProcedure          = 'HLCCoreHits',
        NHitsThreshold         = 2,
        MaxNIterations         = 3,
        Streams                = [icetray.I3Frame.Physics]
        )

from icecube.icetray import I3PacketModule
class SplitToQ(I3PacketModule):
    '''
    Take a Null split  (A single P frame for each Q)
    and fold down to a single Q frame.
    '''
    def __init__(self, context):
        I3PacketModule.__init__(self, context, icetray.I3Frame.DAQ)
      self.AddOutBox('OutBox')
    def Configure(self):
        pass
    def FramePacket(self, frames):
        # Purge all non-native items from these frames, otherwise you loose GCD
      #frames[0].purge()
      #frames[1].purge()
      # Merge the single P frame into the Q frame
        if len(frames)>1:
            frames[0].merge(frames[1])
            #print frames[1]
          # Set all stop to be DAQ for all frame objects in Q frame
            if frames[1].Has("I3TriggerHierarchy"):
                frames[0].change_stream("I3TriggerHierarchy",icetray.I3Frame.DAQ)
            if frames[1].Has('FilterMask'):
                frames[0].change_stream('FilterMask',icetray.I3Frame.DAQ)
            if frames[1].Has('HESE_CausalQTot'):
                frames[0].change_stream('HESE_CausalQTot',icetray.I3Frame.DAQ)
            if frames[1].Has('HESE_VHESelfVeto'):
                frames[0].change_stream('HESE_VHESelfVeto',icetray.I3Frame.DAQ)
            if frames[1].Has('I3EventHeader'):
                frames[0].change_stream('I3EventHeader',icetray.I3Frame.DAQ)
            if frames[1].Has('SplitInIceDSTPulses'):
                frames[0].change_stream('SplitInIceDSTPulses',icetray.I3Frame.DAQ)
            if frames[1].Has('SplitInIceDSTPulsesTimeRange'):
                frames[0].change_stream('SplitInIceDSTPulsesTimeRange',icetray.I3Frame.DAQ)
            if frames[1].Has('SplitInIcePulses'):
                frames[0].change_stream('SplitInIcePulses',icetray.I3Frame.DAQ)
            if frames[1].Has('SplitInIcePulsesTimeRange'):
                frames[0].change_stream('SplitInIcePulsesTimeRange',icetray.I3Frame.DAQ)
            if frames[1].Has('SplitSRTInIcePulses'):
                frames[0].change_stream('SplitSRTInIcePulses',icetray.I3Frame.DAQ)
            if frames[1].Has('TriggerSplitterLaunchWindow'):
                frames[0].change_stream('TriggerSplitterLaunchWindow',icetray.I3Frame.DAQ)
          # Push *JUST* the q frame.
        self.PushFrame(frames[0])

tray.AddModule(SplitToQ, 'movetoq')

from icecube import TopologicalSplitter
tray.AddModule("I3TopologicalSplitter","TopologicalSplit",
        SubEventStreamName='topological_split',
        InputName="SplitSRTInIcePulses",
        OutputName="TSInIcePulses",
        Multiplicity=4,
        TimeWindow=4000*I3Units.ns,
        TimeCone=800*I3Units.ns,
        SaveSplitCount=True)

#tray.Add('Dump')
#####################################
## LAYER VETO, HESE
#####################################
from icecube import VHESelfVeto
tray.Add('HomogenizedQTot', 'qtot_topo', Output='HomogenizedQTot_toposplit', Pulses='TSInIcePulses')

def remove_small_splits(frame):
    if frame['HomogenizedQTot_toposplit'].value<90.:
        return False
tray.Add(remove_small_splits)

from icecube import DomTools
tray.Add('I3LCPulseCleaning', 'cleaning', OutputHLC='SplitInIcePulsesHLC', OutputSLC='', Input='SplitInIcePulses')

tray.Add('HomogenizedQTot', 'qtot', Output='HomogenizedQTot', Pulses='SplitInIcePulsesHLC')

def CleanDeepCore(frame):
    pulsemap = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'SplitInIcePulsesHLC')
    mask = dataclasses.I3RecoPulseSeriesMapMask(frame, 'SplitInIcePulsesHLC', lambda om, idx, pulse: om.string <= 78)
    frame['SplitInIcePulsesHLC_NoDC'] = mask
tray.Add(CleanDeepCore, 'nodc')

from icecube.CascadeL3_IC79.level4.veto import LayerVeto
tray.AddSegment(LayerVeto, "VetoLayer", Pulses='SplitInIcePulsesHLC_NoDC')

def layer_veto_cut(frame):
    if "HESE_CausalQTot" in frame:
        if frame["HESE_CausalQTot"].value>6000. and frame["HESE_VHESelfVeto"].value==False:
            print "HESE (ck) event!", frame["I3EventHeader"].event_id
            frame["IsHESE_ck"]=icetray.I3Bool(True)
        else:
            frame["IsHESE_ck"]=icetray.I3Bool(False)
    else:
        frame["IsHESE_ck"]=icetray.I3Bool(False)

    nstring = len(set(om.string for om in dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'SplitInIcePulsesHLC_NoDC').keys()))
    layer_veto_charge = frame['VetoLayer0'].value + frame['VetoLayer1'].value
    if nstring>3 and layer_veto_charge<3 and frame['HomogenizedQTot'].value > 6e3:
        frame["IsHESE_jvs"]=icetray.I3Bool(True)
        print "HESE (jvs) event!", frame["I3EventHeader"].event_id
        return True
    elif nstring>3 and layer_veto_charge==0 and frame['HomogenizedQTot'].value > 100.:
        frame["IsHESE_jvs"]=icetray.I3Bool(False)
        print "Potential MESE event!", layer_veto_charge, frame['HomogenizedQTot'].value, frame['HomogenizedQTot_toposplit'].value, frame["I3EventHeader"].event_id
        return True
    elif frame["IsHESE_ck"].value==True:
        return True
    else:
        return False
tray.Add(layer_veto_cut)

#####################################
## ALL THE CASCADE FIT GOODIES
#####################################

tray.AddModule( 'I3TimeWindowCleaning<I3RecoPulse>', 'TWC_HLC',
        InputResponse = 'SplitInIcePulsesHLC', # ! Use pulse series
        OutputResponse = 'TWSplitInIcePulsesHLC', # ! Name of cleaned pulse series
        TimeWindow = 6000 * I3Units.ns
        )

from icecube import clast, cscd_llh
tray.AddModule('I3CLastModule', 'CascadeLast',
        Name = "CascadeLast",
        InputReadout = 'TWSplitInIcePulsesHLC'
        )

tray.AddModule( 'I3CscdLlhModule', 'CascadeLlh',
        InputType = 'RecoPulse', # ! Use reco pulses
        RecoSeries = 'TWSplitInIcePulsesHLC', # ! Name of input pulse series
        FirstLE = True, # Default
        SeedWithOrigin = False, # Default
        SeedKey = "CascadeLast", # ! Seed fit - CLast reco
        MinHits = 8, # ! Require 8 hits
        AmpWeightPower = 0.0, # Default
        ResultName = "CascadeLlhVertexFit", # ! Name of fit result
        Minimizer = 'Powell', # ! Set the minimizer to use
        PDF = 'UPandel', # ! Set the pdf to use
        ParamT = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamX = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamY = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        ParamZ = '1.0, 0.0, 0.0, false',   # ! Setup parameters
        )
from icecube.photonics_service import I3PhotoSplineService
base = os.path.expandvars('$I3_DATA/photon-tables/splines/ems_mie_z20_a10.%s.fits')
pxs = I3PhotoSplineService(base % "abs", base % "prob", 0)

from icecube.millipede import HighEnergyExclusions, MonopodFit
exclusions = tray.AddSegment(HighEnergyExclusions, Pulses="SplitInIcePulses",ExcludeDeepCore=False,BadDomsList='BadDomsList')

millipede_config = dict(Pulses="SplitInIcePulses", PartialExclusion=True, CascadePhotonicsService=pxs, BadDOMs=exclusions, Parametrization='HalfSphere')

tray.AddSegment(MonopodFit, 'MonopodAmpFit', Seed="CascadeLlhVertexFit", PhotonsPerBin=-1, **millipede_config)
tray.AddSegment(MonopodFit, 'MonopodPreFit', Seed="MonopodAmpFit", PhotonsPerBin=5, **millipede_config)
tray.AddSegment(MonopodFit, 'MonopodFit4',Seed='MonopodPreFit', PhotonsPerBin=5,
        Iterations=4, DOMEfficiency=0.99, MintimeWidth=15, **millipede_config)

#####################################
## INCOMING TRACK VETO W/O CUTS
#####################################

from icecube.CascadeL3_IC79.level4.veto import ACausalHitSelector, TrackVeto, LayerVeto, VetoMarginCalculator
cleaned_pulses = tray.AddSegment(ACausalHitSelector, 'ACausalSplitInIcePulses_mono', pulses='SplitInIcePulses', vertex="MonopodFit4")
tray.AddSegment(TrackVeto, 'VetoTrack_mono', TimeWindow=(-15, 1000), Pulses=cleaned_pulses, Radius=100,
        Direction='down', TrackType=dataclasses.I3Particle.StoppingTrack, Vertex="MonopodFit4")

tray.AddSegment(TrackVeto, 'UpgoingTrack_mono', TimeWindow=(-30, 500), Pulses=cleaned_pulses, Radius=100,
        Direction='up', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MonopodFit4")

cleaned_pulses_hlc = tray.AddSegment(ACausalHitSelector, 'ACausalSplitInIcePulsesHLC_mono', pulses='SplitInIcePulsesHLC', vertex="MonopodFit4")
tray.AddSegment(TrackVeto, 'StartingTrackHLC_mono', TimeWindow=(-30, 500), Pulses=cleaned_pulses_hlc, Radius=100,
        Direction='all', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MonopodFit4")

tray.AddSegment(TrackVeto, 'UpgoingTrackHLC_mono', TimeWindow=(-30, 500), Pulses=cleaned_pulses_hlc, Radius=100,
        Direction='up', TrackType=dataclasses.I3Particle.StartingTrack, NSide=8, Vertex="MonopodFit4")

tray.AddModule(VetoMarginCalculator, 'VetoTrackMargin_mono', Vertex="MonopodFit4", Geometry="I3Geometry")

################################
## ALL THE TRACK FIT GOODIES
################################

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

# spline MPE
infmuonprobsplinepath = "/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_prob_z20a10_V2.fits"
infmuonampsplinepath = '/cvmfs/icecube.opensciencegrid.org/data/photon-tables/splines/InfBareMu_mie_abs_z20a10_V2.fits'
spline=I3PhotoSplineService(infmuonampsplinepath,infmuonprobsplinepath, 4)
llh="MPE"

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

# millipede
tray.AddModule("MuMillipede", "millipede",
        MuonPhotonicsService=None,
        CascadePhotonicsService=pxs,
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
                if first_loss_found==False:
                    frame["MillipedeFirstLoss"] = losses[i]
                    first_loss_found=True
    else:
        deposited_energy = 0.

    if "MillipedeFirstLoss" not in frame:
        frame["MillipedeFirstLoss"] = frame["MonopodFit4"]
    frame["MillipedeDepositedEnergy"] = dataclasses.I3Double(deposited_energy)
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

def get_track_length(frame):
    if "Millipede_SplineMPE" in frame:
        losses = frame["Millipede_SplineMPE"]
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
    else:
        frame["TrackLength"] = dataclasses.I3Double(0.)
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
tray.AddModule(ComputeChargeWeightedDist,"CCWD_Cascade",Pulses="TWTSInIcePulsesHLC_NoDC",Track="MonopodFit4")

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
tray.AddModule(collectStats, "CollectStats")

#tray.Add('Dump')
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
