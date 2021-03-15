#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Metaprojects/icerec.h5/build


##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.XLevel/build

import numpy as np
import tables
import argparse
from icecube.icetray import I3Frame
import math, glob, sys, os
from os.path import expandvars
from icecube import MuonInjector
from optparse import OptionParser
import numpy as np
from icecube import tableio, hdfwriter
from I3Tray import *
from icecube import icetray, dataio, dataclasses
from icecube.common_variables import direct_hits
from icecube.common_variables import hit_multiplicity
from icecube.icetray import I3Units
from icecube.common_variables import track_characteristics
from icecube import phys_services

from icecube import VHESelfVeto, DomTools

from icecube import weighting
from icecube.weighting.weighting import from_simprod
from icecube.weighting.fluxes import GaisserH3a
from icecube import truncated_energy
load( "libtruncated_energy" )
load("libmue")

default_xs_location="/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/jobs/files/xs_iso/"

parser = OptionParser(usage='%prog [OPTIONS] FILE1 [FILE2, FILE3, ...]',
    description='Concatenate tables, optionally converting between supported table formats.')

parser.add_option('--dsdxdy_nu_CC', 	default=default_xs_location+"dsdxdy_nu_CC_iso.fits")
parser.add_option('--dsdxdy_nubar_CC', 	default=default_xs_location+"dsdxdy_nubar_CC_iso.fits")
parser.add_option('--dsdxdy_nu_NC', 	default=default_xs_location+"dsdxdy_nu_NC_iso.fits")
parser.add_option('--dsdxdy_nubar_NC', 	default=default_xs_location+"dsdxdy_nubar_NC_iso.fits")
parser.add_option('--Flux', 			default='/data/ana/SterileNeutrino/IC86/HighEnergy/Analysis/scripts/jobs/Resources/GolemFit/resources/FluxOscCalculator/atmospheric_0_0.000000_0.000000_0.000000_0.000000_0.000000_0.000000.hdf5')

parser.add_option("-g","--gcdfile",dest="gcdfile",type="string",
		default = "/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz")

parser.add_option("-i","--infile", dest="infile", type="string",
		default = "/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Ares/IC86.AVG/XLevel/domeff_0.97/00001-01000/XLevel_00_11_00828.i3.zst")

parser.add_option("-o","--outfile",dest="outfile",type="string",
		default = "/data/user/saxani/test_h5.h5")

parser.add_option("-c","--configuration",dest="configuration",type="string",
		default = "/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/Nominal/Gen/00001-01000/Generation_data.lic")

parser.add_option("--update_golden",type="string",default = "False")
parser.add_option("--update_diamond",type="string",default = "False")
parser.add_option("--update_platinum",type="string",default = "True")
parser.add_option("--update_titanium",type="string",default = "False")
parser.add_option("--update_hese",  type="string",default = "False")

parser.add_option("--update_main",type="string",default = "False")
parser.add_option("--update_cuts",type="string",default = "False")
parser.add_option("--update_control",type="string",default = "False")
parser.add_option("--update_charges",type="string",default = "False")
parser.add_option("--update_EZ",type="string",default = "False")
parser.add_option("--update_lite",type="string",default = "True")
parser.add_option("--update_everything",type="string",default = "False")

options,args = parser.parse_args()
gcdfile 	= options.gcdfile
infile		= options.infile
outfile 	= options.outfile
infiles 	= [gcdfile,infile]

update_hese     = options.update_hese
update_diamond  = options.update_diamond
update_golden   = options.update_golden
update_platinum = options.update_platinum

update_main     = options.update_main
update_cuts     = options.update_cuts
update_control  = options.update_control
update_EZ       = options.update_EZ
update_lite     = options.update_lite
update_charges  = options.update_charges
update_everything = options.update_everything

print(infiles)
check_infile = dataio.I3File(infiles[-1])
frame = check_infile.pop_daq()
corsika_file = False
muon_gun_file = False
if 'Multisim' in infiles[-1]:
       print('Looks like a MultiSim file')
       Multisim = True
else:
       Multisim = False

if frame.Has('I3MCTree'):
    print('Infile is MC')
    isMC = True
    if 'CORSIKA' in infiles[-1]:
        if 'muon_gun' in  infiles[-1]:
            corsika_file = True
            muon_gun_file = True
            print('muon_gun INFO --------------------------------------------------------')
            from icecube import MuonGun

        else:
            print('CORSIKA')
            corsika_file = True
            corsika_set = int(infiles[-1].split('/')[7])
            flux = GaisserH3a()
            generator = from_simprod(int(corsika_set))
            cache = weighting.SimprodNormalizations()
            cache.refresh(corsika_set)
            #for dataset in (corsika_set):
            #	cache.refresh(str(dataset))
   
 
else:
	print('Infile is Data')
	isMC = False
	run_num = int(infile.split('Run')[-1].split('_')[0])
	print(run_num)
	#live_time = float(infile.split('LT')[-1].split('_')[0])/1000.
                
if isMC:
    if corsika_file:
        print('CORSIKA file')
        if muon_gun_file:
            print('Actuallly, it si a muon gun file')
            #generator = harvest_generators(infiles[1:])
    else:
        dsdxdy_nu_CC 	= options.dsdxdy_nu_CC
        dsdxdy_nubar_CC = options.dsdxdy_nubar_CC
        dsdxdy_nu_NC 	= options.dsdxdy_nu_NC
        dsdxdy_nubar_NC = options.dsdxdy_nubar_NC
        
        Flux = options.Flux
        print(Flux)

def harvest_generators(infiles):
    """
    Harvest serialized generator configurations from a set of I3 files.
    """
    #from icecube.icetray.i3logging import log_info as log
    generator = None
    for infile in infiles:
        f_num = infile.split('/')[-2]
        f_set = infile.split('/')[-4]
        f_name = infile.split('/')[-1].replace('XLevel','photon').replace('.zst','.bz2')
        fname = '/data/ana/Cscd/StartingEvents/MuonGun/'+f_set+'/photon/'+f_num+'/'+f_name

        f = dataio.I3File(fname)
        fr = f.pop_frame(icetray.I3Frame.Stream('S'))
        #fr = f.pop_frame()
        f.close()
        if fr is not None:
            for k in fr.keys():
                v = fr[k]
                if isinstance(v, MuonGun.GenerationProbability):
                    #log('%s: found "%s" (%s)' % (fname, k, type(v).__name__), unit="MuonGun")
                    if generator is None:
                        generator = v
                    else:
                        generator += v
    return generator


def get_TE(tray):
	tray.AddService( "I3PhotonicsServiceFactory", "PhotonicsServiceMu",
		PhotonicsTopLevelDirectory  ="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/",
		DriverFileDirectory         ="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/SPICEMie/driverfiles",
		PhotonicsLevel2DriverFile   = "mu_photorec.list", 
		PhotonicsTableSelection     = 2,
		ServiceName                 = "PhotonicsServiceMu" )

	tray.AddModule( "I3TruncatedEnergy",
		RecoPulsesName          = "SRTInIcePulses",
		RecoParticleName        = "TrackFit",
		ResultParticleName      = 'TruncatedEnergy',
		I3PhotonicsServiceName  = 'PhotonicsServiceMu',
		UseRDE                  = True)

def hese_stuff(tray):
        pulses_split = 'InIcePulses'
        tray.AddModule('I3LCPulseCleaning', 'cleaning', OutputHLC='HLCPulses', OutputSLC='', Input=pulses_split)
        
        tray.AddModule('VHESelfVeto', 'selfveto',
            Pulses              = 'HLCPulses',
            Geometry            = "I3Geometry",
            TopBoundaryWidth    = 90.,#args.veto_top_boundary_width,
            BottomBoundaryWidth = 10.,#args.veto_bottom_boundary_width,
            DustLayer           = -135.,#args.veto_dust_layer,
            DustLayerWidth      = 90.,#args.veto_dust_layer_width,
            TimeWindow          = 3000.,#args.veto_time_window * I3Units.ns,
            VertexThreshold     = 250.,#args.veto_vertex_threshold,
            VetoThreshold       = 3.)
        #store_parameter(tray, 'VHESelfVeto')
        tray.AddModule('HomogenizedQTot', 'qtot_causal', Pulses=pulses_split, Output='CausalQTot', VertexTime='VHESelfVetoVertexTime')


def get_variables(frame):
    if not "I3EventHeader" in frame:
        print("Frame Failed: No I3EventHeader")
        return False
    #event_id = str(frame["I3EventHeader"].run_id)+ '_' + str(frame["I3EventHeader"].event_id) + '_'+str(frame["I3EventHeader"].sub_event_id)
    #print(float(str(frame["I3EventHeader"].run_id)+'.'+str(frame["I3EventHeader"].event_id)))
    event_id = float(str(frame["I3EventHeader"].run_id)+'.'+str(frame["I3EventHeader"].event_id))
    #print(event_id)
    frame.Put("EventID"		,dataclasses.I3Double(frame["I3EventHeader"].event_id))
    frame.Put("RunID"		,dataclasses.I3Double(frame["I3EventHeader"].run_id))
    if isMC:
        if corsika_file:
            if muon_gun_file:
                #weights = frame['MuonWeight'].value
                #frame.Put("weights"         ,dataclasses.I3Double(weights))
                #print(weights)
                #print(frame['weights'])
                #frame.Put("weights"			,dataclasses.I3Double(weights))
                frame.Put("FinalStateX"		,dataclasses.I3Double(NaN))
                frame.Put("FinalStateY"		,dataclasses.I3Double(NaN))
                frame.Put("TotalColumnDepth",dataclasses.I3Double(NaN))
                frame.Put("ImpactParameter"	,dataclasses.I3Double(NaN))
                frame.Put("NuEnergy"		,dataclasses.I3Double(NaN))
                frame.Put("NuZenith"		,dataclasses.I3Double(NaN))
                frame.Put("NuAzimuth"		,dataclasses.I3Double(NaN))
                frame.Put("MuonEnergy"		,dataclasses.I3Double(NaN))
                frame.Put("MuonZenith"		,dataclasses.I3Double(NaN))
                frame.Put("MuonZ"			,dataclasses.I3Double(NaN))
                frame.Put("MuonAzimuth"		,dataclasses.I3Double(NaN))
                frame.Put("HadronEnergy"	,dataclasses.I3Double(NaN))
                frame.Put("HadronZenith"	,dataclasses.I3Double(NaN))
                frame.Put("TrueMuExEnergy"	,dataclasses.I3Double(NaN))
                frame.Put("InjectionRadius" ,dataclasses.I3Double(NaN)) 
                frame.Put("TrueMuExZenith"	,dataclasses.I3Double(NaN))
                frame.Put("VHESelfVeto"  ,dataclasses.I3Double(NaN))
                if update_hese == 'True':
                    frame.Put("CausalQTot"  ,dataclasses.I3Double(NaN))

            else:
                energy = frame['MCPrimaries'][0].energy
                ptype = frame['MCPrimaries'][0].type
                weights = flux(energy, ptype)/generator(energy, ptype)
                #print(weights)
                if weights < 4.5e-06:
                    return False
                frame.Put("FinalStateX"		,dataclasses.I3Double(NaN))
                frame.Put("FinalStateY"		,dataclasses.I3Double(NaN))
                frame.Put("TotalColumnDepth",dataclasses.I3Double(NaN))
                frame.Put("ImpactParameter"	,dataclasses.I3Double(NaN))
                frame.Put("NuEnergy"		,dataclasses.I3Double(NaN))
                frame.Put("NuZenith"		,dataclasses.I3Double(NaN))
                frame.Put("NuAzimuth"		,dataclasses.I3Double(NaN))
                muon_energy 	=  dataclasses.get_most_energetic_muon(frame['I3MCTree']).energy
                muon_zenith 	=  dataclasses.get_most_energetic_muon(frame['I3MCTree']).dir.zenith
                muon_azimuth 	=  dataclasses.get_most_energetic_muon(frame['I3MCTree']).dir.azimuth
                muon_Z 			=  dataclasses.get_most_energetic_muon(frame['I3MCTree']).pos.z
                #print(muon_energy)
                frame.Put("MuonEnergy"		,dataclasses.I3Double(muon_energy))
                frame.Put("MuonZenith"		,dataclasses.I3Double(muon_zenith))
                frame.Put("MuonZ"			,dataclasses.I3Double(muon_Z))
                frame.Put("MuonAzimuth"		,dataclasses.I3Double(muon_azimuth))
                frame.Put("HadronEnergy"	,dataclasses.I3Double(NaN))
                frame.Put("HadronZenith"	,dataclasses.I3Double(NaN))
                frame.Put("weights"			,dataclasses.I3Double(weights))
                frame.Put("TrueMuExEnergy"	,dataclasses.I3Double(NaN))
                frame.Put("InjectionRadius" ,dataclasses.I3Double(NaN)) 
                frame.Put("TrueMuExZenith"	,dataclasses.I3Double(NaN))

       
            frame.Put("weights"			,dataclasses.I3Double(weights))


            # HESE related values
            ''' 
            # 1 means that the event was vetoed or non-existant. 0 Means it is a good hese event.
            if "VHESelfVeto" not in frame:
                frame.Put("VHESelfVeto"        ,dataclasses.I3Double(1))
            else:
                if frame["VHESelfVeto"].value:
                    del frame["VHESelfVeto"]
                    frame.Put("VHESelfVeto"        ,dataclasses.I3Double(1))
                else:
                    del frame["VHESelfVeto"]
                    frame.Put("VHESelfVeto"        ,dataclasses.I3Double(0))
           
            #HESE
            causalqtot = frame["CausalQTot"].value
            del frame["CausalQTot"]
            frame.Put("CausalQTot"        ,dataclasses.I3Double(causalqtot))
            '''
    
            #if frame.Stop==icetray.I3Frame.Physics:
            
            #print(frame['I3MCTree'])

            I3MCTree = frame['I3MCTree']
            muon = I3MCTree[1]
            hadron = I3MCTree[2]

            x = -muon.pos.x
            y = -muon.pos.y
            z = -muon.pos.z
            theta = muon.dir.zenith
            phi   = muon.dir.azimuth
            r     = 100000.
            x2 = r * np.sin(theta) * np.cos(phi)
            y2 = r * np.sin(theta) * np.sin(phi)
            z2 = r * np.cos(theta)
            p2 = np.asarray([x2,y2,z2])
            p1 = np.asarray([0,0,0])
            p  = np.asarray([x,y,z])
            dist = (p-p1)-np.dot((p-p1),(p2-p1))*(p2-p1)/np.linalg.norm(p2-p1)**2
            injection_radius = np.linalg.norm(dist)

            frame.Put("InjectionRadius" ,dataclasses.I3Double(injection_radius)) 
            frame.Put("MuonEnergy"		,dataclasses.I3Double(I3MCTree[1].energy))
            frame.Put("MuonZenith"		,dataclasses.I3Double(I3MCTree[1].dir.zenith))
            frame.Put("MuonZ"			,dataclasses.I3Double(I3MCTree[1].dir.z))
            frame.Put("MuonAzimuth"		,dataclasses.I3Double(I3MCTree[1].dir.azimuth))
            frame.Put("HadronEnergy"	,dataclasses.I3Double(I3MCTree[2].energy))
            frame.Put("HadronZenith"	,dataclasses.I3Double(I3MCTree[2].dir.zenith))
            #frame.Put("TrueMuExEnergy"	,dataclasses.I3Double(frame['MuExTrue'].energy))
            #frame.Put("TrueMuExZenith"	,dataclasses.I3Double(frame['MuExTrue'].dir.zenith))
            frame.Put("PionWeight"		,dataclasses.I3Double(pion_weight))
            frame.Put("KaonWeight"		,dataclasses.I3Double(kaon_weight))

    else:
        frame.Put("FinalStateX"		,dataclasses.I3Double(NaN))
        frame.Put("FinalStateY"		,dataclasses.I3Double(NaN))
        frame.Put("TotalColumnDepth",dataclasses.I3Double(NaN))
        frame.Put("ImpactParameter"	,dataclasses.I3Double(NaN))
        frame.Put("NuEnergy"		,dataclasses.I3Double(NaN))
        frame.Put("NuZenith"		,dataclasses.I3Double(NaN))
        frame.Put("NuAzimuth"		,dataclasses.I3Double(NaN))
        frame.Put("MuonEnergy"		,dataclasses.I3Double(NaN))
        frame.Put("MuonZenith"		,dataclasses.I3Double(NaN))
        frame.Put("MuonZ"			,dataclasses.I3Double(NaN))
        frame.Put("MuonAzimuth"		,dataclasses.I3Double(NaN))
        frame.Put("HadronEnergy"	,dataclasses.I3Double(NaN))
        frame.Put("HadronZenith"	,dataclasses.I3Double(NaN))
        weights = 1.0
        frame.Put("weights"			,dataclasses.I3Double(weights))
        frame.Put("TrueMuExEnergy"	,dataclasses.I3Double(NaN))
        frame.Put("InjectionRadius" ,dataclasses.I3Double(NaN)) 
        frame.Put("TrueMuExZenith"	,dataclasses.I3Double(NaN))
        frame.Put("RunNum"		,dataclasses.I3Double(run_num))

    frame.Put("MuExEnergy"		,dataclasses.I3Double(frame['MuEx'].energy))
    frame.Put("MuExZenith"		,dataclasses.I3Double(frame['MuEx'].dir.zenith))
    frame.Put("MuExAzimuth"		,dataclasses.I3Double(frame['MuEx'].dir.azimuth))
    frame.Put("MuExCOGR"		,dataclasses.I3Double(frame['MuEx'].pos.r))
    frame.Put("MuExCOGZ"		,dataclasses.I3Double(frame['MuEx'].pos.z))
    frame.Put("SplitGeo1DirN"		,dataclasses.I3Double(frame['SPEFit4SplitGeo1_dh'].n_dir_doms))
    frame.Put("SplitGeo2DirN"		,dataclasses.I3Double(frame['SPEFit4SplitGeo2_dh'].n_dir_doms))
    frame.Put("SplitTime1DirN"		,dataclasses.I3Double(frame['SPEFit4SplitTime1_dh'].n_dir_doms))
    frame.Put("SplitTime2DirN"		,dataclasses.I3Double(frame['SPEFit4SplitTime2_dh'].n_dir_doms))
    frame.Put("TrackFitDirN"		        ,dataclasses.I3Double(frame['TrackFit_dh'].n_dir_doms))
    
    #frame.Put("SplitGeo1NChan"		,dataclasses.I3Double(frame['SPEFit4SplitGeo1_hm'].n_hit_doms))
    #frame.Put("SplitGeo2NChan"		,dataclasses.I3Double(frame['SPEFit4SplitGeo2_hm'].n_hit_doms))
    #frame.Put("SplitTime1NChan"		,dataclasses.I3Double(frame['SPEFit4SplitTime1_hm'].n_hit_doms))
    #frame.Put("SplitTime2NChan"		,dataclasses.I3Double(frame['SPEFit4SplitTime2_hm'].n_hit_doms))

    #frame.Put("MuExZenith"		,dataclasses.I3Double(frame['MuExCompat'].dir.zenith))
    #frame.Put("MuExAzimuth"		,dataclasses.I3Double(frame['MuExCompat'].dir.azimuth))
    #frame.Put("MuExCOGR"		,dataclasses.I3Double(frame['MuExCompat'].pos.r))
    #frame.Put("MuExCOGZ"		,dataclasses.I3Double(frame['MuExCompat'].pos.z))
   
    #print(frame['MuExCompat'].energy,frame['MuEx'].energy)
    if Multisim == True:
        MultisimAmplitudes = frame['MultisimAmplitudes'][:7]
        MultisimPhases = frame['MultisimPhases'][:7]
        del frame['MultisimAmplitudes']
        del frame['MultisimPhases']
        del frame['MultisimModes']
        #print(MultisimAmplitudes)
        frame.Put('MultisimAmplitudes' ,dataclasses.I3VectorDouble(MultisimAmplitudes))
        frame.Put('MultisimPhases'     ,dataclasses.I3VectorDouble(MultisimPhases))
        #frame.Put('MultisimModes'      ,dataclasses.I3VectorDouble(MultisimModes))

    pulse_series 	= dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'TTPulses_NoDC')
    geo 			= frame.Get('I3Geometry')
    om_geo 			= geo.omgeo
    Qtot 			= 0
    SrQ 			= 0
    SzQ				= 0
    ASzQ			= 0
    ASrQ			= 0
    Sr				= 0
    Sz				= 0
    SrQN			= 0
    SzQN			= 0
    NChan = len(pulse_series)
    charges 		= []
    depth 			= []
    depth_charges 	= []
    for om, pulses in pulse_series:
        om_pos = om_geo[om].position
        om_z = om_pos[2]
        om_r = np.sqrt(om_pos[0]**2+om_pos[1]**2+om_pos[2]**2)
        Qom = 0
        for pulse in pulses:
            Qom += pulse.charge
        charges.append(Qom)
        depth.append(om_z)
        depth_charges.append(Qom)
        SrQ += om_r * Qom
        SzQ += om_z * Qom
        SrQN += om_r * Qom * len(pulses)
        SzQN += om_z * Qom * len(pulses)
        Qtot += Qom
        Sr += om_r
        Sz += om_z
    if int(Qtot) == 0:
        ASzQ = NaN
        ASrQ = NaN
    else:
        ASzQ = SzQ/Qtot
        ASrQ = SrQ/Qtot

    if not muon_gun_file:
        charge_counts, bin_edges   = np.histogram(charges, bins = 60, range = (0,3), weights = np.ones(len(charges))*weights)
        charge_err   , bin_edges   = np.histogram(charges, bins = 60, range = (0,3), weights = np.ones(len(charges))*(weights**2))
        frame.Put('charges'		,dataclasses.I3VectorDouble(charge_counts))
        frame.Put('charges_sw2' ,dataclasses.I3VectorDouble(charge_err))

        depth_counts, bin_edges = np.histogram(depth, bins = 60, range = (-500,500),weights = depth_charges)
    frame.Put('SrQ_NoDC'	,dataclasses.I3Double(SrQ))
    frame.Put('SzQ_NoDC'	,dataclasses.I3Double(SzQ))
    frame.Put('SrQN_NoDC'	,dataclasses.I3Double(SrQN))
    frame.Put('SzQN_NoDC'	,dataclasses.I3Double(SzQN))
    frame.Put('ASrQ_NoDC'	,dataclasses.I3Double(ASrQ))
    frame.Put('ASzQ_NoDC'	,dataclasses.I3Double(ASzQ))
    frame.Put('NChan_NoDC'	,dataclasses.I3Double(NChan))
    frame.Put('Qtot_NoDC'	,dataclasses.I3Double(Qtot))
    frame.Put('Sz_NoDC'		,dataclasses.I3Double(Sz))
    frame.Put('Sr_NoDC'		,dataclasses.I3Double(Sr))
    #frame.Put('depth'		,dataclasses.I3VectorDouble(depth_counts))
    frame.Put('event_id'	,dataclasses.I3Double(event_id))

    pulse_series 	= dataclasses.I3RecoPulseSeriesMap.from_frame(frame,'TTPulses')
    Qtot 			= 0
    SrQ 			= 0
    SzQ				= 0
    ASzQ			= 0
    ASrQ			= 0
    Sr				= 0
    Sz				= 0
    SrQN			= 0
    SzQN			= 0
    NChan = len(pulse_series)
    charges 		= []
    depth 			= []
    depth_charges 	= []
    for om, pulses in pulse_series:
        om_pos = om_geo[om].position
        om_z = om_pos[2]
        om_r = np.sqrt(om_pos[0]**2+om_pos[1]**2+om_pos[2]**2)
        Qom = 0
        for pulse in pulses:
            Qom += pulse.charge
        SrQ += om_r * Qom
        SzQ += om_z * Qom
        SrQN += om_r * Qom * len(pulses)
        SzQN += om_z * Qom * len(pulses)
        Qtot += Qom
        Sr += om_r
        Sz += om_z
    if int(Qtot) == 0:
        ASzQ = NaN
        ASrQ = NaN
    else:
        ASzQ = SzQ/Qtot
        ASrQ = SrQ/Qtot
    frame.Put('SrQ'		,dataclasses.I3Double(SrQ))
    frame.Put('SzQ'		,dataclasses.I3Double(SzQ))
    frame.Put('SrQN'	,dataclasses.I3Double(SrQN))
    frame.Put('SzQN'	,dataclasses.I3Double(SzQN))
    frame.Put('ASrQ'	,dataclasses.I3Double(ASrQ))
    frame.Put('ASzQ'	,dataclasses.I3Double(ASzQ))
    frame.Put('NChan'	,dataclasses.I3Double(NChan))
    frame.Put('Qtot'	,dataclasses.I3Double(Qtot))
    frame.Put('Sz'		,dataclasses.I3Double(Sz))
    frame.Put('Sr'		,dataclasses.I3Double(Sr))

    deepCoreStrings = range(79,87)
    TTpulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses')
    track=frame['TrackFit']

    AvgDistQ_NoDC = 0
    AvgDistQ_DC = 0
    Qtot_NoDC = 0
    Qtot_DC = 0

    SzQ_NoDC = 0
    SzQ_DC = 0

    SrQ_NoDC = 0
    SrQ_DC = 0

    for omkey,pulse in TTpulses:
        if omkey[0] in deepCoreStrings:
            om_pos = om_geo[omkey].position
            Dist = phys_services.I3Calculator.closest_approach_distance(track,om_pos)
            Z = om_pos[2]
            R = np.sqrt(om_pos[0]**2+om_pos[1]**2+om_pos[2]**2)
            Qdom = 0
            for q in pulse:
                Qdom+=q.charge
            Qtot_DC += Qdom
            AvgDistQ_DC += Dist * Qdom
            SzQ_DC += Z*Qdom
            SrQ_DC += R*Qdom
        else:
            om_pos=om_geo[omkey].position
            Dist=phys_services.I3Calculator.closest_approach_distance(track,om_pos)
            Z = om_pos[2]
            R = np.sqrt(om_pos[0]**2+om_pos[1]**2+om_pos[2]**2)
            Qdom=0
            for q in pulse:
                Qdom+=q.charge
            Qtot_NoDC += Qdom
            AvgDistQ_NoDC += Dist*Qdom
            SzQ_NoDC += Z*Qdom
            SrQ_NoDC += R*Qdom

    if(Qtot_DC == 0):
        COGz_DC = NaN
        COGr_DC = NaN
        AvgDistQ_DC=NaN
    else:
        COGz_DC = SzQ_DC/Qtot_DC
        COGr_DC = SrQ_DC/Qtot_DC
        AvgDistQ_DC /= Qtot_DC

    if(Qtot_NoDC == 0):
        COGz_NoDC = NaN
        COGr_NoDC = NaN
        AvgDistQ_NoDC=NaN
    else:
        COGz_NoDC = SzQ_NoDC/Qtot_NoDC
        COGr_NoDC = SrQ_NoDC/Qtot_NoDC
        AvgDistQ_NoDC /= Qtot_NoDC

    frame.Put("AvgDistQ_DC"		,dataclasses.I3Double(AvgDistQ_DC))
    frame.Put("AvgDistQ_NoDC"	,dataclasses.I3Double(AvgDistQ_NoDC))
    frame.Put("COGz_DC"			,dataclasses.I3Double(COGz_DC))
    frame.Put("COGz_NoDC"		,dataclasses.I3Double(COGz_NoDC))
    frame.Put("COGr_DC"			,dataclasses.I3Double(COGr_DC))
    frame.Put("COGr_NoDC"		,dataclasses.I3Double(COGr_NoDC))

    '''
    TTpulses = dataclasses.I3RecoPulseSeriesMap.from_frame(frame, 'TTPulses')
    track 	= frame['TrackFit']
    phi 	= track.dir.zenith
    theta   = track.dir.azimuth
    x  	  	= track.pos.x
    y		= track.pos.y
    z		= track.pos.z

    e_x = np.sin(phi)*np.cos(theta)
    e_y = np.sin(phi)*np.sin(theta)
    e_z = np.cos(phi)

    x1 = x + 10000 * e_x
    y1 = y + 10000 * e_y
    z1 = z + 10000 * e_z
    p1 = np.asarray([x1,y1,z1])

    x2 = x - 10000 * e_x
    y2 = y - 10000 * e_y
    z2 = z - 10000 * e_z
    p2 = np.asarray([x2,y2,z2])

    AvgDistQ = 0
    Qtot = 0
    for om, pulses in TTpulses:
        om_pos = om_geo[om].position
        p0 = np.asarray([om_pos[0],om_pos[1],om_pos[2]])

        dist = np.linalg.norm(np.cross(p0-p1,p0-p2))/np.linalg.norm(p2-p1)
        Qdom = 0
        for q in pulse:
            Qdom+=q.charge
        Qtot += Qdom
        AvgDistQ += dist*Qdom
    AvgDistQ = AvgDistQ/Qtot
    frame.Put("my_AvgDistQ"		,dataclasses.I3Double(AvgDistQ))
    '''
    return True




def get_hese(frame):
    veto 		= frame.Get("VHESelfVeto").value
    causalqtot 	= frame.Get("CausalQTot").value
    #print(veto,causalqtot)
    if int(veto) == 1:
        return False
    if causalqtot < 6000.:
        return False
    return True

def remove_hese(frame):
    veto 		= frame.Get("VHESelfVeto").value
    causalqtot 	= frame.Get("CausalQTot").value

    if int(veto) == 0:# Throw out events that would make it into HESE
        if causalqtot >= 6000.:
            return False
    return True



def get_platinum(frame):
    if np.logical_or(get_diamond(frame),get_golden(frame)):
        #print('keeping event')
        return True
    else:
        return False


def get_lite(frame):
	keysDict ={
			#'CausalQTot'   			:frame['CausalQTot'].value,
			#'VHESelfVeto'           :frame['VHESelfVeto'].value,
			'PionWeight'   			:frame['PionWeight'].value,
			'KaonWeight'   			:frame['KaonWeight'].value,
			'weights'   			:frame['weights'].value,
			#'oneweight'   			:frame['oneweight'].value,
			'MuExEnergy'   			:frame['MuExEnergy'].value,
			'MuExZenith'     		:frame['MuExZenith'].value,
			'MuExAzimuth'     		:frame['MuExAzimuth'].value,
			'NuEnergy'    			:frame['NuEnergy'].value,
			'NuZenith'  			:frame['NuZenith'].value,
			'NuAzimuth'             :frame['NuAzimuth'].value,
			'FinalStateX'	       	:frame['FinalStateX'].value,
			'FinalStateY'          	:frame['FinalStateY'].value,
			'TotalColumnDepth'     	:frame['TotalColumnDepth'].value,
			'ImpactParameter'     	:frame['ImpactParameter'].value,
			'PrimaryType'     	    :frame['PrimaryType'].value,
			'FinalType0'     		:frame['FinalType0'].value,
			'FinalType1'     		:frame['FinalType1'].value,
			}
	frame.Put('keysDict',dataclasses.I3MapStringDouble(keysDict))
	return True


def get_control(frame):
    frame.Put('C_MuExAzimuth'     ,dataclasses.I3Double(frame['MuExAzimuth'].value))
    frame.Put('C_DirL'     ,dataclasses.I3Double(frame.Get("TrackFit_dh").dir_track_length))
    frame.Put('C_DirS'     ,dataclasses.I3Double(frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness))
    frame.Put('C_DirNDoms'     ,dataclasses.I3Double(frame.Get("TrackFit_dh").n_dir_doms))
    frame.Put('C_RLogL'     ,dataclasses.I3Double(frame['RLogL'].value))
    frame.Put('C_BayesLLHR'     ,dataclasses.I3Double(frame['BayesLLHR'].value))
    frame.Put('C_ParaboloidSigma'     ,dataclasses.I3Double(frame['CorrectedParaboloidSigma'].value))
    frame.Put('C_AvgDistQ'     ,dataclasses.I3Double(frame['TrackFit_AvgDistQ'].value))
    frame.Put('C_Separation'     ,dataclasses.I3Double(frame['trackfit_separation'].value))
    frame.Put('C_NChan'     ,dataclasses.I3Double(frame['NChan'].value))
    frame.Put('C_Qtot'     ,dataclasses.I3Double(frame['Qtot'].value))
    frame.Put('C_Overburden'     ,dataclasses.I3Double(frame['Overburden'].value))
    frame.Put('C_DirL_NoDC'     ,dataclasses.I3Double(frame.Get("TrackFit_NoDC_dh").dir_track_length))
    frame.Put('C_DirS_NoDC'     ,dataclasses.I3Double(frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness))
    frame.Put('C_OuterCharge'     ,dataclasses.I3Double(frame['TTPulses_total_outer_charge'].value))
    frame.Put('C_NChan_NoDC'     ,dataclasses.I3Double(frame['NChan_NoDC'].value))
    frame.Put('C_Qtot_NoDC'     ,dataclasses.I3Double(frame['Qtot_NoDC'].value))
    frame.Put('C_MuExCOGZ'     ,dataclasses.I3Double(frame['MuExCOGZ'].value))
    frame.Put('C_MuExCOGR'     ,dataclasses.I3Double(frame['MuExCOGR'].value))
    frame.Put('C_SplitGeo1DirN'     ,dataclasses.I3Double(frame['SplitGeo2DirN'].value))
    frame.Put('C_SplitGeo2DirN'     ,dataclasses.I3Double(frame['SplitGeo2DirN'].value))
    frame.Put('C_SplitTime1DirN'     ,dataclasses.I3Double(frame['SplitTime1DirN'].value))
    frame.Put('C_SplitTime2DirN'     ,dataclasses.I3Double(frame['SplitTime2DirN'].value))
    frame.Put('C_TrackFitDirN'     ,dataclasses.I3Double(frame['TrackFitDirN'].value))
    return True

def get_cuts(frame):
	keysDict ={
			#'CausalQTot'   			:frame['CausalQTot'].value,
			#'VHESelfVeto'           :frame['VHESelfVeto'].value,
			'EventID'   			:frame['EventID'].value,
			'RunID'   			    :frame['RunID'].value,
			'MuonZenith'    		:frame['MuonZenith'].value,
			'MuonEnergy' 			:frame['MuonEnergy'].value,
			'MuonAzimuth'			:frame['MuonAzimuth'].value,
			'weights'		        :frame['weights'].value,
			'MuExCOGZ'    			:frame['MuExCOGZ'].value,
			'MuExCOGR'    			:frame['MuExCOGR'].value,
			'MuExAzimuth'  			:frame['MuExAzimuth'].value,
			'MuExEnergy'   			:frame['MuExEnergy'].value,
			'MuExZenith'   			:frame['MuExZenith'].value,
			'NuEnergy'    			:frame['NuEnergy'].value,
			'NuZenith'  			:frame['NuZenith'].value,
			'NuAzimuth'             :frame['NuAzimuth'].value,
			'PrimaryType'     	    :frame['PrimaryType'].value,
			'oneweight'   			:frame['oneweight'].value,
			'DirL'					:frame.Get("TrackFit_dh").dir_track_length,
			'DirS'					:frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness,
			'DirL_NoDC'				:frame.Get("TrackFit_NoDC_dh").dir_track_length,
			'DirS_NoDC'    			:frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness,
			'DirNStrings'			:frame.Get("TrackFit_dh").n_dir_strings,
			'DirNPulses'			:frame.Get("TrackFit_dh").n_dir_pulses,
			'DirNDoms'				:frame.Get("TrackFit_dh").n_dir_doms,
			'DirNStrings_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_strings,
			'DirNPulses_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_pulses,
			'DirNDoms_NoDC'				:frame.Get("TrackFit_NoDC_dh").n_dir_doms,
			'OuterCharge'         	:frame['TTPulses_total_outer_charge'].value,
			'RLogL'               	:frame['RLogL'].value,
			'BayesLLHR'           	:frame['BayesLLHR'].value,
			'ParaboloidSigma'    	:frame['CorrectedParaboloidSigma'].value,
			'AvgDistQ'             	:frame['TrackFit_AvgDistQ'].value,
			'Separation'           	:frame['trackfit_separation'].value,
			'NChan_NoDC'			:frame['NChan_NoDC'].value,
			'Qtot_NoDC' 			:frame['Qtot_NoDC'].value,
			'NChan' 				:frame['NChan'].value,
			'Qtot' 					:frame['Qtot'].value,
			'Overburden'           	:frame['Overburden'].value,
			'InjectionRadius'     	:frame['InjectionRadius'].value,
            'AvgDistQ_DC'           :frame['AvgDistQ_DC'].value,
            'AvgDistQ_NoDC'         :frame['AvgDistQ_NoDC'].value,
            'COGz_DC'               :frame['COGz_DC'].value,
            'COGz_NoDC'             :frame['COGz_NoDC'].value,
            'COGr_DC'               :frame['COGr_DC'].value,
            'COGr_NoDC'             :frame['COGr_NoDC'].value,
	}
	frame.Put('keysDict',dataclasses.I3MapStringDouble(keysDict))
	return True

def get_main(frame):
	keysDict ={
			#'CausalQTot'   			:frame['CausalQTot'].value,
			#'VHESelfVeto'           :frame['VHESelfVeto'].value,
			'EventID'   			:frame['EventID'].value,
			'RunID'   			    :frame['RunID'].value,
			'MuonZenith'    		:frame['MuonZenith'].value,
			'MuonEnergy' 			:frame['MuonEnergy'].value,
			'MuonAzimuth'			:frame['MuonAzimuth'].value,
			'weights'		        :frame['weights'].value,
			'MuExCOGZ'    			:frame['MuExCOGZ'].value,
			'MuExCOGR'    			:frame['MuExCOGR'].value,
			'MuExAzimuth'  			:frame['MuExAzimuth'].value,
			'MuExEnergy'   			:frame['MuExEnergy'].value,
			'MuExZenith'   			:frame['MuExZenith'].value,
			'NuEnergy'    			:frame['NuEnergy'].value,
			'NuZenith'  			:frame['NuZenith'].value,
			'DirL'					:frame.Get("TrackFit_dh").dir_track_length,
			'DirS'					:frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness,
			'DirL_NoDC'				:frame.Get("TrackFit_NoDC_dh").dir_track_length,
			'DirS_NoDC'    			:frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness,
			'DirNStrings'			:frame.Get("TrackFit_dh").n_dir_strings,
			'DirNPulses'			:frame.Get("TrackFit_dh").n_dir_pulses,
			'DirNDoms'				:frame.Get("TrackFit_dh").n_dir_doms,
			'DirNStrings_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_strings,
			'DirNPulses_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_pulses,
			'DirNDoms_NoDC'				:frame.Get("TrackFit_NoDC_dh").n_dir_doms,
			'OuterCharge'         	:frame['TTPulses_total_outer_charge'].value,
			'RLogL'               	:frame['RLogL'].value,
			'BayesLLHR'           	:frame['BayesLLHR'].value,
			'ParaboloidSigma'    	:frame['CorrectedParaboloidSigma'].value,
			'AvgDistQ'             	:frame['TrackFit_AvgDistQ'].value,
			'Separation'           	:frame['trackfit_separation'].value,
			'NChan_NoDC'			:frame['NChan_NoDC'].value,
			'Qtot_NoDC' 			:frame['Qtot_NoDC'].value,
			'NChan' 				:frame['NChan'].value,
			'Qtot' 					:frame['Qtot'].value,
			'Overburden'           	:frame['Overburden'].value,
			'InjectionRadius'     	:frame['InjectionRadius'].value,
            'AvgDistQ_DC'           :frame['AvgDistQ_DC'].value,
            'AvgDistQ_NoDC'         :frame['AvgDistQ_NoDC'].value,
            'COGz_DC'               :frame['COGz_DC'].value,
            'COGz_NoDC'             :frame['COGz_NoDC'].value,
            'COGr_DC'               :frame['COGr_DC'].value,
            'COGr_NoDC'             :frame['COGr_NoDC'].value,
	}
	frame.Put('keysDict',dataclasses.I3MapStringDouble(keysDict))
	return True

def get_full(frame):
	keysDict ={
			'MuonZenith'    		:frame['MuonZenith'].value,
			'MuonEnergy' 			:frame['MuonEnergy'].value,
			'MuonAzimuth'			:frame['MuonAzimuth'].value,
			'TrueMuExEnergy'	    :frame['TrueMuExEnergy'].value,
			'TrueMuExZenith'    	:frame['TrueMuExZenith'].value,
			'weights'		        :frame['weights'].value,
			'MuExCOGZ'    			:frame['MuExCOGZ'].value,
			'MuExCOGR'    			:frame['MuExCOGR'].value,
			'MuExAzimuth'  			:frame['MuExAzimuth'].value,
			'MuExEnergy'   			:frame['MuExEnergy'].value,
			'MuExZenith'   			:frame['MuExZenith'].value,
			'DirL'					:frame.Get("TrackFit_dh").dir_track_length,
			'DirS'					:frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness,
			'DirL_NoDC'				:frame.Get("TrackFit_NoDC_dh").dir_track_length,
			'DirS_NoDC'    			:frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness,
			'DirNStrings'			:frame.Get("TrackFit_dh").n_dir_strings,
			'DirNPulses'			:frame.Get("TrackFit_dh").n_dir_pulses,
			'DirNDoms'				:frame.Get("TrackFit_dh").n_dir_doms,
			'DirNStrings_NoDC'		:frame.Get("TrackFit_NoDC_dh").n_dir_strings,
			'DirNPulses_NoDC'		:frame.Get("TrackFit_NoDC_dh").n_dir_pulses,
			'DirNDoms_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_doms,
			'OuterCharge'         	:frame['TTPulses_total_outer_charge'].value,
			'RLogL'               	:frame['RLogL'].value,
			'BayesLLHR'           	:frame['BayesLLHR'].value,
			'ParaboloidSigma'    	:frame['CorrectedParaboloidSigma'].value,
			'AvgDistQ'             	:frame['TrackFit_AvgDistQ'].value,
			'Separation'           	:frame['trackfit_separation'].value,
			'Overburden'           	:frame['Overburden'].value,
			'NChan_NoDC'			:frame['NChan_NoDC'].value,
			'Qtot_NoDC' 			:frame['Qtot_NoDC'].value,
			'NChan' 				:frame['NChan'].value,
			'Qtot' 					:frame['Qtot'].value,
			'Sz_NoDC' 				:frame['Sz_NoDC'].value,
			'Sr_NoDC' 				:frame['Sr_NoDC'].value,
			'SrQ' 					:frame['SrQ'].value,
			'SzQ' 					:frame['SzQ'].value,
			'SrQN' 					:frame['SrQN'].value,
			'SzQN' 					:frame['SzQN'].value,
			'ASrQ' 					:frame['ASrQ'].value,
			'ASzQ'					:frame['ASzQ'].value,
			'Sz' 					:frame['Sz'].value,
			'Sr' 					:frame['Sr'].value,
			'InjectionRadius'     	:frame['InjectionRadius'].value,
	}
	frame.Put('keysDict',dataclasses.I3MapStringDouble(keysDict))
	return True

def get_everything(frame):
	keysDict ={
			#'RunNum'				:frame['RunNum'].value,
			'AvgDistQ_DC'			:frame['AvgDistQ_DC'].value,
			'AvgDistQ_NoDC'			:frame['AvgDistQ_NoDC'].value,
			'COGz_DC'				:frame['COGz_DC'].value,
			'COGz_NoDC'				:frame['COGz_NoDC'].value,
			'COGr_DC'				:frame['COGr_DC'].value,
			'COGr_NoDC'				:frame['COGr_NoDC'].value,
			'MuonZenith'    		:frame['MuonZenith'].value,
			'MuonEnergy' 			:frame['MuonEnergy'].value,
			'MuonAzimuth'			:frame['MuonAzimuth'].value,
			'NuEnergy'	    		:frame['NuEnergy'].value,
			'NuZenith'     			:frame['NuZenith'].value,
			'NuAzimuth'				:frame['NuAzimuth'].value,
			'TrueMuExEnergy'	    :frame['TrueMuExEnergy'].value,
			'TrueMuExZenith'	    :frame['TrueMuExZenith'].value,
			'weights'		        :frame['weights'].value,
			'MuExCOGZ'    			:frame['MuExCOGZ'].value,
			'MuExCOGR'    			:frame['MuExCOGR'].value,
			'MuExAzimuth'  			:frame['MuExAzimuth'].value,
			'MuExEnergy'   			:frame['MuExEnergy'].value,
			'MuExZenith'   			:frame['MuExZenith'].value,
			'DirNDStrings'			:frame['dir_N_D_strings'].value,
			'DirNDDoms'				:frame['dir_N_D_doms'].value,
			'DirNDPulses'			:frame['dir_N_D_pulses'].value,
			'DirLD'					:frame['dir_L_D'].value,
			'DirSD'					:frame['dir_S_D'].value,
			'DirNStrings'			:frame.Get("TrackFit_dh").n_dir_strings,
			'DirNPulses'			:frame.Get("TrackFit_dh").n_dir_pulses,
			'DirNDoms'				:frame.Get("TrackFit_dh").n_dir_doms,
			'DirL'					:frame.Get("TrackFit_dh").dir_track_length,
			'DirS'					:frame.Get("TrackFit_dh").dir_track_hit_distribution_smoothness,
			'DirNStrings_NoDC'		:frame.Get("TrackFit_NoDC_dh").n_dir_strings,
			'DirNPulses_NoDC'		:frame.Get("TrackFit_NoDC_dh").n_dir_pulses,
			'DirNDOMs_NoDC'			:frame.Get("TrackFit_NoDC_dh").n_dir_doms,
			'DirL_NoDC'				:frame.Get("TrackFit_NoDC_dh").dir_track_length,
			'DirS_NoDC'    			:frame.Get("TrackFit_NoDC_dh").dir_track_hit_distribution_smoothness,
			'OuterCharge'         	:frame['TTPulses_total_outer_charge'].value,
			'RLogL'               	:frame['RLogL'].value,
			'BayesLLHR'           	:frame['BayesLLHR'].value,
			'ParaboloidSigma'    	:frame['CorrectedParaboloidSigma'].value,
			'AvgDistQ'             	:frame['TrackFit_AvgDistQ'].value,
			'Separation'           	:frame['trackfit_separation'].value,
			'Overburden'           	:frame['Overburden'].value,
			'SrQ_NoDC' 				:frame['SrQ_NoDC'].value,
			'SzQ_NoDC' 				:frame['SzQ_NoDC'].value,
			'SrQN_NoDC' 			:frame['SrQN_NoDC'].value,
			'SzQN_NoDC' 			:frame['SzQN_NoDC'].value,
			'ASrQ_NoDC' 			:frame['ASrQ_NoDC'].value,
			'ASzQ_NoDC' 			:frame['ASzQ_NoDC'].value,
			'NChan_NoDC'			:frame['NChan_NoDC'].value,
			'Qtot_NoDC' 			:frame['Qtot_NoDC'].value,
			'Sz_NoDC' 				:frame['Sz_NoDC'].value,
			'Sr_NoDC' 				:frame['Sr_NoDC'].value,
			'SrQ' 					:frame['SrQ'].value,
			'SzQ' 					:frame['SzQ'].value,
			'SrQN' 					:frame['SrQN'].value,
			'SzQN' 					:frame['SzQN'].value,
			'ASrQ' 					:frame['ASrQ'].value,
			'ASzQ'					:frame['ASzQ'].value,
			'NChan' 				:frame['NChan'].value,
			'Qtot' 					:frame['Qtot'].value,
			'Sz' 					:frame['Sz'].value,
			'Sr' 					:frame['Sr'].value,
			'InjectionRadius'     	:frame['InjectionRadius'].value,
            "SplitGeo1DirN"         :frame['SplitGeo2DirN'].value,
            "SplitGeo2DirN"         :frame['SplitGeo2DirN'].value,
            'SplitTime1DirN'        :frame['SplitTime1DirN'].value,
            'SplitTime2DirN'        :frame['SplitTime2DirN'].value,
            'TrackFitDirN'          :frame['TrackFitDirN'].value,
            'AvgDistQ_DC'           :frame['AvgDistQ_DC'].value,
            'AvgDistQ_NoDC'         :frame['AvgDistQ_NoDC'].value,
            'COGz_DC'               :frame['COGz_DC'].value,
            'COGz_NoDC'             :frame['COGz_NoDC'].value,
            'COGr_DC'               :frame['COGr_DC'].value,
            'COGr_NoDC'             :frame['COGr_NoDC'].value,
	}
	frame.Put('keysDict',dataclasses.I3MapStringDouble(keysDict))
	return True

# We need:
# 1. No Cuts for developing diamond selection
# 2. golden_EZ for easy plotting
# 3. golden_lite for sterilizer
# 4. golden full for comparing cuts


if isMC == True:
    if corsika_file:
        if muon_gun_file:
            print('yah')
            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_platinum)
            tray.Add(get_main)
            
            outputKeys = ['keysDict'] 
            h5file = outfile.replace('.h5','_main_platinum.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()
            
            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_main)
            outputKeys = ['keysDict']
            h5file = outfile.replace('.h5','_main.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()

        else:
            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_platinum)
            tray.Add(get_main)
            
            outputKeys = ['keysDict'] 
            h5file = outfile.replace('.h5','_main_platinum.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()

            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_diamond)
            tray.Add(get_main)
            outputKeys = ['keysDict'] 
            h5file = outfile.replace('.h5','_main_diamond.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()

            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_golden)
            tray.Add(get_main)
            outputKeys = ['keysDict'] 
            h5file = outfile.replace('.h5','_main_golden.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()

            tray = I3Tray()
            tray.Add("I3Reader",FilenameList = infiles)
            tray.Add(get_variables)
            tray.Add(get_main)
            outputKeys = ['keysDict']
            h5file = outfile.replace('.h5','_main.h5')
            print('Saving: '+h5file)
            tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                    ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                    ("SubEventStreams",["TTrigger"]),
                    ("keys",outputKeys)
                    )
            tray.Execute()
            tray.Finish()

    else:
        
        if update_platinum =='True':
            print('Updating Platinum')
            if update_control == 'True':
                tray = I3Tray()
                tray.Add("I3Reader",FilenameList = infiles)
                tray.Add(get_variables)
                tray.Add(get_platinum)
                tray.Add(get_control)
                outputKeys = ['oneweight','MuExEnergy','MuExZenith','MuExAzimuth','NuEnergy','NuZenith','NuAzimuth',
                                'FinalStateX','FinalStateY','TotalColumnDepth','ImpactParameter','PrimaryType','FinalType0','FinalType1',
                                'C_MuExAzimuth','C_DirL','C_DirS','C_DirNDoms','C_RLogL','C_BayesLLHR','C_ParaboloidSigma','C_AvgDistQ',
                                'C_Separation','C_NChan','C_Qtot','C_Overburden',
                                'C_DirL_NoDC','C_DirS_NoDC','C_OuterCharge','C_NChan_NoDC',
                                'C_Qtot_NoDC','C_MuExCOGZ','C_MuExCOGR','C_SplitGeo1DirN','C_SplitGeo2DirN',
                                'C_SplitTime1DirN','C_SplitTime2DirN','C_TrackFitDirN']
                h5file = outfile.replace('.h5','_control_platinum.h5')
                print('Saving: '+h5file)
                tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                        ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                        ("SubEventStreams",["TTrigger"]),
                        ("keys",outputKeys)
                        )
                tray.Execute()
                tray.Finish()

            if update_lite == 'True':
                tray = I3Tray()
                tray.Add("I3Reader",FilenameList = infiles)
                tray.Add(get_variables)
                tray.Add(get_platinum)
                tray.Add(get_lite)
                outputKeys = ['oneweight','MuExEnergy','MuExZenith','MuExAzimuth','NuEnergy','NuZenith','NuAzimuth','FinalStateX','FinalStateY','TotalColumnDepth','ImpactParameter','PrimaryType','FinalType0','FinalType1']
                h5file = outfile.replace('.h5','_lite_platinum.h5')
                print('Saving: '+h5file)
                tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                        ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                        ("SubEventStreams",["TTrigger"]),
                        ("keys",outputKeys)
                        )
                tray.Execute()
                tray.Finish()

			
elif isMC != True:
    if update_platinum == 'True':
        tray = I3Tray()
        tray.Add("I3Reader",FilenameList = infiles)
        tray.Add(get_variables)
        tray.Add(get_platinum)
        tray.Add(get_everything)
        outputKeys = ['keysDict','charges','MuExEnergy','MuExZenith','MuExAzimuth','MuEx']
        h5file = outfile.replace('.h5','_everything_platinum.h5')
        print('Saving: '+h5file)
        tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                ("SubEventStreams",["TTrigger"]),
                ("keys",outputKeys)
                )
        i3streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')]
        i3file = outfile.replace('.h5','_everything_platinum.i3.bz2')
        tray.AddModule("I3Writer","i3writer")(
                ("Filename",i3file),
                ("Streams",i3streams),
                ('DropOrphanStreams',[icetray.I3Frame.DAQ]),
                )

        tray.Execute()
        tray.Finish()

    if update_diamond == 'True':
        tray = I3Tray()
        tray.Add("I3Reader",FilenameList = infiles)
        tray.Add(get_variables)
        tray.Add(get_diamond)
        tray.Add(get_everything)
        outputKeys = ['keysDict','charges','MuExEnergy','MuExZenith','MuExAzimuth','MuEx']
        h5file = outfile.replace('.h5','_everything_diamond.h5')
        print('Saving: '+h5file)
        tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                ("SubEventStreams",["TTrigger"]),
                ("keys",outputKeys)
                )
        i3streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')]
        i3file = outfile.replace('.h5','_everything_diamond.i3.bz2')
        tray.AddModule("I3Writer","i3writer")(
                ("Filename",i3file),
                ("Streams",i3streams),
                ('DropOrphanStreams',[icetray.I3Frame.DAQ]),
                )

        tray.Execute()
        tray.Finish()
    if update_golden =='True':
        tray = I3Tray()
        tray.Add("I3Reader",FilenameList = infiles)
        tray.Add(get_variables)
        tray.Add(get_golden)
        tray.Add(get_everything)
        outputKeys = ['keysDict','charges','MuExEnergy','MuExZenith','MuExAzimuth','MuEx']
        h5file = outfile.replace('.h5','_everything_golden.h5')
        print('Saving: '+h5file)
        tray.AddModule(tableio.I3TableWriter, "hdfwriter")(
                ("tableservice",hdfwriter.I3HDFTableService(h5file)),
                ("SubEventStreams",["TTrigger"]),
                ("keys",outputKeys)
                )
        i3streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics, icetray.I3Frame.Stream('S'), icetray.I3Frame.Stream('M')]
        i3file = outfile.replace('.h5','_everything_golden.i3.bz2')
        tray.AddModule("I3Writer","i3writer")(
                ("Filename",i3file),
                ("Streams",i3streams),
                ('DropOrphanStreams',[icetray.I3Frame.DAQ]),
                )

        tray.Execute()
        tray.Finish()

