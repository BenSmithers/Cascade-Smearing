#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py2-v2/icetray-start
#METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.trunk/build/
##METAPROJECT /data/ana/SterileNeutrino/IC86/HighEnergy/MC/Metaprojects/icerec.trunk.WaveDeformUpdate/build

### The main L2 processing script ###

import os
import sys
import subprocess

from I3Tray import *
from icecube import icetray, dataclasses, dataio 
from icecube.icetray import I3PacketModule, I3Units

from icecube.filterscripts.offlineL2.level2_all_filters import OfflineFilter
from icecube.filterscripts.offlineL2 import SpecialWriter

def make_parser():
	"""Make the argument parser"""
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-s","--simulation", action="store_true",
		default=True, dest="mc", help="Mark as simulation (MC)")
    
	parser.add_option("-i", "--infile", action="store",
		type="string", default="/data/ana/SterileNeutrino/IC86/HighEnergy/MC/Systematics/Noise/IC86.AVG/L1/domeff_1.27/05001-06000/Noise_L1_00_22_05326.i3.zst", dest="infile",
		help="Input i3 file(s)  (use comma separated list for multiple files)")
    
	parser.add_option("-g", "--gcdfile", action="store",
		type="string", default="/data/ana/SterileNeutrino/IC86/HighEnergy/SPE_Templates/SPE_harvesting/SPE_GCD/GCD_Files_for_JC/GCD_923_OWD_50ns_v3_Pass2/GeoCalibDetectorStatus_IC86.AVG_Pass2.i3.gz", dest="gcdfile",
		help="GCD file for input i3 file")
    
	parser.add_option("-o", "--outfile", action="store",
		type="string", default="/data/user/saxani/test_L2.i3.zst", dest="outfile",
		help="Output i3 file")
    
	parser.add_option("--move",
		type="string", default="False", dest="move",
		help="Make a temp file, then move it?")
    
	parser.add_option("--osg", action="store",
		type="string", default="False", dest="osg",
		help="Do you want to run on the OSG??")   
    
	parser.add_option("-n", "--num", action="store",
		type="int", default=-1, dest="num",
		help="Number of frames to process")
    
	parser.add_option("--dstfile", action="store",
		type="string", default=None, dest="dstfile",
		help="DST root file (should be .root)")
    
	parser.add_option("--gapsfile", action="store",
		type="string", default=None, dest="gapsfile",
		help="gaps text file (should be .txt)")
    
	parser.add_option("--icetopoutput", action="store",
		type="string", default=None, dest="icetopoutput",
		help="Output IceTop file")
    
	parser.add_option("--eheoutput", action="store",
		type="string", default=None, dest="eheoutput",
		help="Output EHE i3 file")
    
	parser.add_option("--slopoutput", action="store",
		type="string", default=None, dest="slopoutput",
		help="Output SLOP file")
    
	parser.add_option("--rootoutput", action="store",
		type="string", default=None, dest="rootoutput",
		help="Output root file")
    
	parser.add_option("--photonicsdir", action="store",
		type="string", default="/cvmfs/icecube.opensciencegrid.org/data/photon-tables/",
		dest="photonicsdir", help="Directory with photonics tables")
    
	return parser

def main(options, stats={}):
    """The main L2 processing script"""
    tray = I3Tray()

    '''
    # make list of input files from GCD and infile
    if isinstance(options['infile'],list):
        infiles = [options['gcdfile']]
        infiles.extend(options['infile'])
    else:
        infiles = [options['gcdfile'], options['infile']]
    print('infiles: ',infiles)
    '''
    ########################################                                                                                                                                                
    gcdfile = options['gcdfile']                                                                                                                                                            
    infile = options['infile']                                                                                                                                                              
    outfile = options['outfile']                                                                                                                                                            
    move = str(options['move'])                                                                                                                                                             
    osg = options['osg'] 

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


    # test access to input and output files
    # read input files
    tray.Add( dataio.I3Reader, "Reader",
        Filenamelist = infiles
    )

    tray.AddSegment(OfflineFilter, "OfflineFilter",
        dstfile=options['dstfile'],
        mc=options['mc'],
        doNotQify=options['mc'],
        photonicsdir=options['photonicsdir']
        )

    ###################################################################
    ########### WRITE STUFF                  ##########################
    ###################################################################

    '''	
    # clean up special case for bz2 files to get around empty file bug
    bzip2	_files = []
    if options['outfile'] and options['outfile'].endswith('.bz2'):
        options['outfile'] = options['outfile'][:-4]
        bzip2_files.append(options['outfile'])
    if options['dstfile'] and options['dstfile'].endswith('.bz2'):
        options['dstfile'] = options['dstfile'][:-4]
        bzip2_files.append(options['dstfile'])
    if options['icetopoutput'] and options['icetopoutput'].endswith('.bz2'):
        options['icetopoutput'] = options['icetopoutput'][:-4]
        bzip2_files.append(options['icetopoutput'])
    if options['eheoutput'] and options['eheoutput'].endswith('.bz2'):
        options['eheoutput'] = options['eheoutput'][:-4]
        bzip2_files.append(options['eheoutput'])
    if options['slopoutput'] and options['slopoutput'].endswith('.bz2'):
        options['slopoutput'] = options['slopoutput'][:-4]
        bzip2_files.append(options['slopoutput'])
    '''
    i3streams = [icetray.I3Frame.DAQ,icetray.I3Frame.Physics,icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation]
    if move == "True":
        tray.AddModule("I3Writer","i3writer")(
            ("Filename",outfile_temp),
            ("Streams",i3streams),
            ("DropOrphanStreams",[icetray.I3Frame.DAQ])
            )
    else:
        tray.AddModule("I3Writer","i3writer")(
            ("Filename",outfile),
            ("Streams",i3streams),
            ("DropOrphanStreams",[icetray.I3Frame.DAQ])
            )

    #tray.AddModule( "I3Writer", "EventWriter" ,
    #    Filename = options['outfile'],
    #    Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics,
    #               icetray.I3Frame.TrayInfo, icetray.I3Frame.Simulation],
    #    DropOrphanStreams = [icetray.I3Frame.DAQ],
    #)

    # special outputs
    ''' 
    if options['eheoutput']:
        tray.AddSegment(SpecialWriter.EHEWriter, "write_ehe",
            Filename=options['eheoutput']
        )
    if options['slopoutput']:
        tray.AddSegment(SpecialWriter.SLOPWriter, "write_slop",
            Filename=options['slopoutput']
        )
    if options['rootoutput']:
        tray.AddSegment(SpecialWriter.RootWriter, "write_root",
            Filename = options['rootoutput']
        )
    if options['gapsfile']:
        tray.AddSegment(SpecialWriter.GapsWriter, "write_gaps",
            Filename = options['gapsfile'],
            MinGapTime = 1
        )
    if options['icetopoutput']:
        # this needs to be the last output
        tray.AddModule('Delete','deleteicetopextrakeys',
            keys=['InIceRawData','CleanInIceRawData']
        )
        tray.AddSegment(SpecialWriter.IceTopWriter, "write_icetop",
            Filename=options['icetopoutput']
        )
    '''
    #tray.AddModule("TrashCan", "Bye")

    if int(options['num']) >= 0:
        tray.Execute(options['num'])
    else:
        tray.Execute()

    tray.Finish()

    # print more CPU usage info. than speicifed by default
    '''
    tray.PrintUsage(fraction=1.0) 
    for entry in tray.Usage():
        stats[entry.key()] = entry.data().usertime
    '''
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
	Level2Filter = ipmodule.FromOptionParser(make_parser(),main)
### end iceprod stuff ###


if __name__ == '__main__':
	# run as script from the command line
	# get parsed args
	parser = make_parser()
	(options,args) = parser.parse_args()
	opts = {}
	# convert to dictionary
	for name in parser.defaults:
		value = getattr(options,name)
		if name == 'infile' and ',' in value:
			value = value.split(',') # split into multiple inputs
		opts[name] = value

	# call main function
	main(opts)
