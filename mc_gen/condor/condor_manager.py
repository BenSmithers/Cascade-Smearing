#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4/icetray-start
#METAPROJECT /data/user/bsmithers/metaprojects/snobo/py3-v4.0.1/

import os
import warnings

# Benjamin Smithers
# benjamin.smithers@mavs.uta.edu
# bsmithers@mail.icecube.wisc.edu

# this script can be used to create a series of dags you can submit
# it produces individual ones to be run with individual scripts

# TODO:
#       add argparse like interface for better UX
#       allow the output files to be saved somewhere that's not my personal directory 

# to someone who might want to use this:
#   + you may want to change which gcd file is specified in the 'write_dag_chain' function to one that's not in my /data/user folder
#   + you will need to change the 'workpath' variable so somewhere you have write permission
#   + depending on what kind of generation you want to run, you might have to change the hardcoded arguments in the 'write_dag_chain' function
# Not sure where these arguments come from? You may want to check the scripts that are called like
#   /path/to/script/1-process_gen.py --help

def make_condor_submit(execPath, name, gpus, memory):
    """
    writes the job.submit file to go with the dag. 
    execPath is which executable we will call, this writes the submit file with name "$name.submit" 
    will call for $gpus GPUs and $memory megabytes of memory
    """
       
    # just open a file, write in the file, and close the file...
    condor_file = open("{}.submit".format(name), 'w' )
    condor_file.write("executable = {}\n".format(execPath))
    condor_file.write("arguments  = $(arg)\n" )
    condor_file.write("\n")
    
    condor_file.write("should_transfer_files = true\n")
    condor_file.write("when_to_transfer_output = ON_EXIT\n")
    condor_file.write("getenv = False\n")
    condor_file.write("output = /scratch/bsmithers/{}.$(cluster).out\n".format(name))
    condor_file.write("error  = /scratch/bsmithers/{}.$(cluster).err\n".format(name))
    condor_file.write("log    = /scratch/bsmithers/{}.$(cluster).log\n".format(name))
    condor_file.write("\n")
    condor_file.write("request_memory = {}\n".format(memory))
    condor_file.write("request_gpus   = {}\n".format(gpus))
    condor_file.write("queue 1\n")
    condor_file.close()


# yes there are better ways of doing this.
# no, I don't care
# 
# to get the i'th digit, you shift the decimal point over a few orders of magnitude,
#   round the number (lobbing off the stuff to the right of the dot),
#   and find the modulus of 10, therefore lobbing off everything before the ones place
def get_5d_string( number ):
    """
    takes an integer < 100000, returns a string of the number that's the 5 characters long.
    If the number is < 10000, it puts zeros at the front to make the string 5 chars long
    """
    thing = ""
    for i in range(5):
        thing+= str( int(number/( 10**(4-i))) % 10 )
    return(thing)

def write_dag_chain(execPath, dataPath, base_name, nJobs, nEvents):
    """
    writes a dag file that will generate a sample all the way from 
    generation -> particle propagation -> photon propagation -> detector sim -> (more eventually)

    Arguments:
    1. folder to executables
    2. folder to save data
    3. base name for output data
    4. number of jobs
    5. number of events per job (default 2500)

    Several arguments are hard-coded in. These are all explicitly written for SnowSuite scripts, though it might be easy to modify to suit other needs  
    """
    
    # check to make sure these are the right data types and the folders actually exist
    if not isinstance(nJobs, int):
        raise Exception("nJobs should be type {}, received type {}".format(int, type(nJobs)))
    if not isinstance(nEvents, int):
        raise Exception("nEvents should be type {}, received type {}".format(int, type(nEvents)))
    if not os.path.exists(dataPath):
        raise Exception("Folder {} does not exist".format(dataPath))
    if not os.path.exists(execPath):
        raise Exception("Executable folder {} does not exist.".format(execPath))
    if not isinstance(base_name, str):
        raise Exception("base_name should be type {}, received {}".format(str, type(base_name)))
    
    # make sure my gcd file is still there
    gcd_file    = "/data/user/bsmithers/calibration/GeoCalibDetectorStatus_AVG_55697-57531_PASS2_SPE_wNoiseUpdate.i3.gz"  
    if not os.path.isfile(gcd_file):
        raise Exception("Could not find gcd file at {}".format(gcd_file))

    seedroot    = 5000 # doesn't really matter what this is, but it gives control over the simulation
    
    #arguments are provided as a list ["generic arguments", "input file"]
    # the input thingy may be provided as none, in which case there are no inputs
    generator_name          = "Agen"
    generator               = "{}/process-Gen.py".format(execPath)
    
    # TODO  may need to be frequently changed! 
    #
    #
    # 
    generation_arguments    = ["-i 2 --Emin 100 --Emax 1000000 -t Electron -c NC -n {}".format(nEvents), None, True]
    #
    #
    # TODO
   
    corsika_name = "BCorsika"
    corsika = "{}/gen-CORS.py".format(execPath)
    corsika_args = ["-n {} -g{} ".format(nEvents*15, gcd_file), None, True]

    poly_name = "CPolyplop"
    poly = "{}/process-Polyplopia.py".format(execPath)
    poly_arguments         =["--mctype LeptonInjector --MCTreeName I3MCTree_preMuonProp --OutputMCTreeName I3MCTree --TimeWindow 40 ", generator_name,True]

    config_file = "/home/bsmithers/software_dev/cascade/mc_gen/Snowstorm_FullSystematics.yml"

    phot_name              = "Bphoton"
    phot           = "{}/process-SnowStorm.py".format(execPath)
    phot_arguments         = ["-g {} -c {} --domoversizefactor {} --events-per-model 100".format(gcd_file, config_file, 1.0), generator_name, True]
    
    det_name              = "Cdetector"
    det       = "{}/process-Det.py".format(execPath)
    det_arguments         = [" ", phot_name, True]
    
    l1_name              = "Dlevel1"
    l1       = "{}/process-L1.py".format(execPath)
    l1_arguments         = [" ", det_name, False]


    # default I3 file suffix! 
    suffix = ".i3"
    print("Using file suffix: {}".format(suffix))
    
    
    # dictionary for the names of each level of the MC generation and their arguments 
    level_arguments =  {base_name+generator_name : generation_arguments, 
                        #base_name+corsika_name : corsika_args,
                        #base_name+poly_name : poly_arguments,
                        base_name+phot_name : phot_arguments,
                        base_name+det_name : det_arguments,
                        base_name+l1_name : l1_arguments}
    
    # make sure each of the necessary executables exist and are executable. If they are not executable just issue a warning! 
    for executable in [generator,poly, phot, det, l1]:
        if not os.path.isfile(executable):
            raise Exception("Executable doesn't exist at {}".format(executable))
        # the specified executable is not an executable, the job will fail unless this is changed. Warn the user. 
        if not os.access(executable, os.X_OK):
            warnings.warn("WARNING! CANNOT EXECUTE {}".format(executable))
    
    
    dag_file = open("chain_dag.dag",'w')
    
    
    
    # use that dictionary, get the keys (which will be the data file names too)
    file_names = [jobname for jobname in level_arguments.keys() ]
    file_names.sort()
    for outname in file_names:
        for job in range(nJobs):
            # JOB thing00002 thing.submit
            dag_file.write("JOB  {}{}  {}.submit\n".format(outname, get_5d_string(job), outname ))
    
    # write the arguments from each dictionary entry into the dag
    for outname in file_names:
        for job in range(nJobs):
            # -o /path/to/output/thingout_00002.suffix 
            outfile = "-o {}{}_{}{}".format(dataPath,outname,get_5d_string(job),suffix) 
            
            if level_arguments[outname][1] is not None: 
                # gen is the only level that doesn't take an input
                # -i /path/to/output/thingin_000002.suffix
                infile = "-i {}{}_{}{}".format(dataPath, base_name+level_arguments[outname][1], get_5d_string(job), suffix)
            else:
                infile = ""
            
            if level_arguments[outname][2]:
                args = "-s {} ".format(seedroot+job)
            else:
                args = ""
            args +=  " {} {} {}".format(level_arguments[outname][0] , outfile, infile ) 
            
            #VARS thing00002 -arg1 A -arg2 B -arg3 C [...]
            dag_file.write('VARS  {}{}  arg="{}"\n'.format(outname, get_5d_string(job), args)) 
    
    #write the family dag tree
    for job in range(nJobs):
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[0], get_5d_string(job), file_names[1],get_5d_string(job)))
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[1], get_5d_string(job), file_names[2],get_5d_string(job)))
        dag_file.write('PARENT {}{} CHILD {}{}\n'.format(file_names[2], get_5d_string(job), file_names[3],get_5d_string(job)))
        
    dag_file.close()
    
    # make the .submit files
    make_condor_submit(generator,  file_names[0], gpus=0, memory=1000)
   # make_condor_submit(corsika,  file_names[1], gpus=0, memory=1000)
    #make_condor_submit(poly,  file_names[2], gpus=0, memory=1000)
    make_condor_submit(phot, file_names[1], gpus=1, memory=4000)
    make_condor_submit(det, file_names[2], gpus=0, memory=2000)
    make_condor_submit(l1, file_names[3], gpus=0, memory=2000)


nJobs = 2

root_dir = "/home/bsmithers/software_dev/cascade/mc_gen/"
workpath = "/data/user/bsmithers/data/nc_mc_cascades/"

# exec dir | data output | name | jobs | events per job
write_dag_chain(root_dir, workpath, "cascade", nJobs, 2000)

