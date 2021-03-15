#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.trunk/build/

##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.trunk.WaveDeformUpdate/build

### The simulation L1 processing script ###
from I3Tray import *
from icecube import icetray, dataclasses, dataio, filter_tools, trigger_sim
from icecube import phys_services
from icecube.filterscripts import filter_globals
from icecube.filterscripts.all_filters import OnlineFilter
from icecube.phys_services.which_split import which_split
import os, sys, time
import subprocess
from math import log10, cos, radians
from optparse import OptionParser
from os.path import expandvars

def make_parser():
	"""Make the argument parser"""
	from optparse import OptionParser
	parser = OptionParser()
	
	parser.add_option("-i", "--infile", action="store",
		type="string", default="/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/Noise/Ares/IC86.AVG/Det/domeff_1.27/00001-01000/Ares_Det_00_22_00651.i3.zst", dest="infile",
		help="Input i3 file(s)  (use comma separated list for multiple files)")
    
	parser.add_option("-g", "--gcdfile", action="store",
		type="string", default="/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz", dest="gcdfile",
		help="GCD file for input i3 file")
    
	parser.add_option("-o", "--outfile", action="store",
		type="string", default="/data/user/saxani/test_L1.i3.zst", dest="outfile",
		help="Output i3 file")
    
	parser.add_option("--move",
		type="string", default="False", dest="move",
		help="Make a temp file, then move it?")
    
	parser.add_option("--osg", action="store",
		type="string", default="False", dest="osg",
		help="Do you want to run on the OSG??")
    
	parser.add_option("-n", "--num", action="store",
		type="int", default=0, dest="num",
		help="Number of frames to process")
    
	parser.add_option("--qify", action="store_true",
		default=False, dest="qify",
		help="Apply QConverter, use if file is P frame only")
    
	parser.add_option("--MinBiasPrescale", action="store",
		type="int", default=None, dest="MinBiasPrescale",
		help="Set the Min Bias prescale to something other than default")
    
	parser.add_option("--photonicsdir", action="store",
		type="string", default="/cvmfs/icecube.opensciencegrid.org/data/photon-tables",
		dest="photonicsdir", help="Directory with photonics tables")
    
	parser.add_option("--disable-gfu", action="store_false",
		default=True, dest="GFU",
		help="Do not run GFU filter")
	return parser


def main(options, stats={}):
    """The main L1 script"""

    ########################################
    gcdfile = options['gcdfile']
    infile  = options['infile']
    outfile = options['outfile']
    move    = str(options['move'])
    print(move)
    osg     = options['osg']
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
        #outfile_temp = str(os.getcwd()+'/'+outfile.split('/')[-1])
    else:
        outfile_temp = '/data/ana/SterileNeutrino/IC86/HighEnergy/MC/scripts/temp/'+outfile.split('/')[-1]
        outfile = outfile

    infiles = [gcdfile, infile]
    print(infiles)
    tray = I3Tray()
    tray.Add(dataio.I3Reader, "reader", FilenameList=infiles)

    # run online filters
    online_kwargs = {}
    if options['photonicsdir']:
        online_kwargs.update({
            'SplineRecoAmplitudeTable': os.path.join(options['photonicsdir'],'splines','InfBareMu_mie_abs_z20a10.fits'),
            'SplineRecoTimingTable': os.path.join(options['photonicsdir'],'splines','InfBareMu_mie_prob_z20a10.fits'),
            'hese_followup_base_GCD_filename': gcdfile,
            #'base_GCD_filename': gcdfile,
        })
    if options['GFU'] is not None:
        online_kwargs['gfu_enabled'] = options['GFU']
    tray.AddSegment(OnlineFilter, "OnlineFilter",
                    decode=False, simulation=True,
                    vemcal_enabled=False,
                    #alert_followup_omit_GCD_diff = True,
                    **online_kwargs
                    )

    # make random service
    seed = os.getpid()
    filter_mask_randoms = phys_services.I3GSLRandomService(seed)

    # override MinBias Prescale
    filterconfigs = filter_globals.filter_pairs + filter_globals.sdst_pairs
    if options['MinBiasPrescale']:
        for i,filtertuple in enumerate(filterconfigs):
            if filtertuple[0] == filter_globals.FilterMinBias:
                del filterconfigs[i]
                filterconfigs.append((filtertuple[0],options['MinBiasPrescale']))
                break
    print(filterconfigs)

    # Generate filter Masks for all P frames
    tray.AddModule(filter_tools.FilterMaskMaker, "MakeFilterMasks",
                   OutputMaskName = filter_globals.filter_mask,
                   FilterConfigs = filterconfigs,
                   RandomService = filter_mask_randoms)

    # Merge the FilterMasks
    tray.AddModule("OrPframeFilterMasks", "make_q_filtermask",
                   InputName = filter_globals.filter_mask,
                   OutputName = filter_globals.qfilter_mask)


    #Q+P frame specific keep module needs to go first, as KeepFromSubstram
    #will rename things, let's rename post keep.
    def is_Q(frame):
        return frame.Stop==frame.DAQ

    simulation_keeps = [
            'BackgroundI3MCTree',
            'BackgroundI3MCTreePEcounts',
            'BackgroundI3MCPESeriesMap',
            'BackgroundI3MCTree_preMuonProp',
            'BackgroundMMCTrackList',
            'BeaconLaunches',
            'CorsikaInteractionHeight',
            'CorsikaWeightMap',
            'EventProperties',
            'GenerationSpec',
            'I3LinearizedMCTree',
            'I3MCTree',
            'I3MCTreePEcounts',
            'I3MCTree_preMuonProp',
            'I3MCPESeriesMap',
            'I3MCPulseSeriesMap',
            'I3MCPulseSeriesMapParticleIDMap',
            'I3MCWeightDict',
            'IceXModelNumber',#Modified May 16th.
            'MultisimModes',#Modified May 16th.
            'MultisimPhases',#Modified May 16th.
            'MultisimAmplitudes',#Modified May 16th.
            'LeptonInjectorProperties',
            'MCHitSeriesMap',
            'MCPrimary',
            'MCPrimaryInfo',
            'MMCTrackList',
            'PolyplopiaInfo',
            'PolyplopiaPrimary',
            'RNGState',
            'SignalI3MCPEs',
            'SimTrimmer', # for SimTrimmer flag
            'TimeShift', # the time shift amount
            'WIMP_params', # Wimp-sim
           ]

    keep_before_merge = filter_globals.q_frame_keeps + [
                            'InIceDSTPulses', # keep DST pulse masks
                            'IceTopDSTPulses',
                            'LeptonInjectorProperties',
                            'CalibratedWaveformRange', # keep calibration info
                            'UncleanedInIcePulsesTimeRange',
                            'SplitUncleanedInIcePulses',
                            'SplitUncleanedInIcePulsesTimeRange',
                            'SplitUncleanedInIceDSTPulsesTimeRange',
                            'CalibrationErrata',
                            'SaturationWindows',
                            'InIceRawData', # keep raw data for now
                            'IceTopRawData',
                           ] + simulation_keeps

    tray.AddModule("Keep", "keep_before_merge",
                   keys = keep_before_merge,
                   If=is_Q
                   )

    ## second set of prekeeps, conditional on filter content, based on newly created Qfiltermask
    #Determine if we should apply harsh keep for events that failed to pass any filter
    ##  Note: excluding the sdst_streams entries

    tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckAll",
                   FilterNameList = filter_globals.filter_streams,
                   FilterResultName = filter_globals.qfilter_mask,
                   DecisionName = "PassedAnyFilter",
                   DiscardEvents = False,
                   Streams = [icetray.I3Frame.DAQ]
                   )
    def do_save_just_superdst(frame):
        if frame.Has("PassedAnyFilter"):
            if not frame["PassedAnyFilter"].value:
                return True    #  <- Event failed to pass any filter.
            else:
                return False # <- Event passed some filter

        else:
            print("Failed to find key frame Bool!!")
            return False

    keep_only_superdsts = filter_globals.keep_nofilterpass+[
                             'PassedAnyFilter',
                                'LeptonInjectorProperties',
                             'InIceDSTPulses',
                             'IceTopDSTPulses',
                             'SplitUncleanedInIcePulses',
                             'SplitUncleanedInIcePulsesTimeRange',
                             'SplitUncleanedInIceDSTPulsesTimeRange',
                             'RNGState',
                             ] + simulation_keeps
    tray.AddModule("Keep", "KeepOnlySuperDSTs",
                   keys = keep_only_superdsts,
                   If = do_save_just_superdst
                   )

    ## Now clean up the events that not even the SuperDST filters passed on.
    tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheckSDST",
                   FilterNameList = filter_globals.sdst_streams,
                   FilterResultName = filter_globals.qfilter_mask,
                   DecisionName = "PassedKeepSuperDSTOnly",
                   DiscardEvents = False,
                   Streams = [icetray.I3Frame.DAQ]
                   )

    def dont_save_superdst(frame):
        if frame.Has("PassedKeepSuperDSTOnly") and frame.Has("PassedAnyFilter"):
            if frame["PassedAnyFilter"].value:
                return False  #  <- these passed a regular filter, keeper
            elif not frame["PassedKeepSuperDSTOnly"].value:
                return True    #  <- Event failed to pass SDST filter.
            else:
                return False # <- Event passed some  SDST filter
        else:
            print("Failed to find key frame Bool!!")
            return False

    tray.AddModule("Keep", "KeepOnlyDSTs",
                   keys = filter_globals.keep_dst_only
                          + ["PassedAnyFilter","PassedKeepSuperDSTOnly",
                             filter_globals.eventheader]+simulation_keeps,
                   If = dont_save_superdst
                   )


    tray.Add('Dump')
    ## Frames should now contain only what is needed.  now flatten, write/send to server
    ## Squish P frames back to single Q frame, one for each split:
    tray.AddModule("KeepFromSubstream","null_stream",
                   StreamName = filter_globals.NullSplitter,
                   KeepKeys = filter_globals.null_split_keeps,
                   )

    in_ice_keeps = filter_globals.inice_split_keeps + filter_globals.onlinel2filter_keeps
    in_ice_keeps = in_ice_keeps + ['I3EventHeader',
                                    'LeptonInjectorProperties',
                                   'SplitUncleanedInIcePulses',
                                   'TriggerSplitterLaunchWindow',
                                   'I3TriggerHierarchy',
                                    'SplitUncleanedInIceDSTPulsesTimeRange',#Added by Kotoyo
                                   'GCFilter_GCFilterMJD']
    tray.AddModule("Keep", "inice_keeps",
                   keys = in_ice_keeps,
                   If = which_split(split_name=filter_globals.InIceSplitter),
                   )


    tray.AddModule("KeepFromSubstream","icetop_split_stream",
                   StreamName = filter_globals.IceTopSplitter,
                   KeepKeys = filter_globals.icetop_split_keeps,
                   )

    # Apply small keep list (SuperDST/SmallTrig/DST/FilterMask for non-filter passers
    # Remove I3DAQData object for events not passing one of the 'filters_keeping_allraw'
    tray.AddModule("I3IcePickModule<FilterMaskFilter>","filterMaskCheck",
                   FilterNameList = filter_globals.filters_keeping_allraw,
                   FilterResultName = filter_globals.qfilter_mask,
                   DecisionName = "PassedConventional",
                   DiscardEvents = False,
                   Streams = [icetray.I3Frame.DAQ]
                   )

    ## Clean out the Raw Data when not passing conventional filter
    def I3RawDataCleaner(frame):
        if not (('PassedConventional' in frame and
                 frame['PassedConventional'].value == True) or
                ('SimTrimmer' in frame and
                 frame['SimTrimmer'].value == True)
               ):
            frame.Delete('InIceRawData')
            frame.Delete('IceTopRawData')

    tray.AddModule(I3RawDataCleaner,"CleanErrataForConventional",
                   Streams=[icetray.I3Frame.DAQ])

    ###################################################################
    ########### WRITE STUFF                  ##########################
    ###################################################################

    # clean up special case for bz2 files to get around empty file bug
    '''	
    bzip2_files = []
    if options['outfile'] and options['outfile'].endswith('.bz2'):
        options['outfile'] = options['outfile'][:-4]
        bzip2_files.append(options['outfile'])
    '''
    # Write the physics and DAQ frames
    i3streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics,icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation]
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

    '''
    tray.AddModule("I3Writer", "EventWriter",
                   filename=options['outfile'],
                   Streams=[icetray.I3Frame.DAQ,icetray.I3Frame.Physics,
                            icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation]
                   )
    '''
    tray.AddModule("TrashCan","byebye")

    if options['num'] >= 0:
        tray.Execute(options['num'])
    else:
        tray.Execute()

    tray.Finish()

    # print more CPU usage info. than speicifed by default
    tray.PrintUsage(fraction=1.0)
    for entry in tray.Usage():
        stats[entry.key()] = entry.data().usertime
    '''
    if bzip2_files:
        # now do bzip2
        if os.path.exists('/usr/bin/bzip2'):
           subprocess.check_call(['/usr/bin/bzip2', '-f']+bzip2_files)
        elif os.path.exists('/bin/bzip2'):
           subprocess.check_call(['/bin/bzip2', '-f']+bzip2_files)
        else:
           raise Exception('Cannot find bzip2')
    '''
    # clean up forcefully in case we're running this in a loop
    if move == 'True':
        if osg == "True":
            copy_to_NPX(outfile)
        else:
            os.system(str("mv "+str(outfile_temp) + " " +str(outfile)))
    else:
        if osg =="True":
            print('You need move option to be True if using the OSG')
            pass

    del tray

### iceprod stuff ###
try:
	from iceprod.modules import ipmodule
except ImportError, e:
	print('Module iceprod.modules not found. Will not define IceProd Class')
else:
	Level1SimulationFilter = ipmodule.FromOptionParser(make_parser(),main)
### end iceprod stuff ###


if __name__ == '__main__':
	# run as script from the command line
	# get parsed args
	parser = make_parser()
	(options,args) = parser.parse_args()

	# convert to dictionary
	opts = vars(options)
	#opts['infile'] = opts['infile'].split(',')

	# call main function
	main(opts)

