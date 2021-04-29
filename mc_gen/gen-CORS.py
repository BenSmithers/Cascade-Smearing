#!/bin/sh /cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/icetray-start
#METAPROJECT /data/user/bsmithers/icetray/install/

"""
This script is used to generate CORSIKA air showers and the resulting muons. 
"""

import time
import logging
logging.basicConfig()
rootLogger = logging.getLogger('')
rootLogger.setLevel(logging.INFO)
from optparse import OptionParser

import sys
from shutil import rmtree

if __name__=="__main__":
    #Create parser to accept arguments from user
    parser = OptionParser()
    parser.add_option("-s","--seed",dest="seed",default=1,type="int",
                      help="Require a unique number under 4e9 for each job.")
    parser.add_option("-o","--outfile",dest="outfile",type="string",
                      help="Outfile name",
                      default = "test")
    parser.add_option("--Emin", dest="Emin",type="float",
                      help="Min energy to simulate [GeV]",
                      default = 600)
    parser.add_option("--Emax", dest="Emax",
                      type="float",
                      help="Maximum energy to smulate [GeV]",
                      default = 1e11)
    parser.add_option("-g","--gcd", dest="gcd",
                      type="string",
                      help="GCD File to use",
                      default="/cvmfs/icecube.opensciencegrid.org/data/GCD/GeoCalibDetectorStatus_IC86.55697_corrected_V2.i3.gz")
    parser.add_option("-r","--runno", dest="runno",
                     type="int",
                     help="Enter the run number",
                     default=1)
    parser.add_option("-n","--nEvents", dest="nEvents",
                      type="int",
                      help="Number of events",
                      default = 10000)
    options, args = parser.parse_args()

    print("Parsed?")
    outfile = options.outfile
    e_min = options.Emin
    e_max = options.Emax
    runno = options.runno
    seed  = options.seed
    n_events = options.nEvents
    gcd = options.gcd

    print(outfile)

    from icecube.simprod.modules import Corsika5ComponentGenerator as c5
    #from icecube.simprod.modules import CorsikaGenerator as c5
    import os

    #del options
    #del args
    #sys.argv = [sys.argv[0]]
    summary_file =os.path.join( os.path.dirname(outfile), "summary_{}_.json".format(runno))

    cors = c5()
    cors.SetParameter("eprimarymin",e_min)
    cors.SetParameter("eprimarymax", e_max )
    cors.SetParameter("outputfile", outfile)
    cors.SetParameter("gcdfile", gcd)
    cors.SetParameter("nshowers", n_events)
    cors.SetParameter("seed", seed)
    cors.SetParameter("runnum",runno)
    cors.SetParameter("summaryfile", summary_file)
    cors.SetParameter("cvmfs", "/cvmfs/icecube.opensciencegrid.org/py3-v4.1.1/")
    cors.SetParser(parser)
    print("Executing")
    stats = {}
    print("This is my outputfile {}".format(cors.outputfile))
    cors.Execute({})

    # clean up the leftovers 

    if os.path.exists("bsmithers_icesoft"):
        rmtree("bsmithers_icesoft")
    if os.path.exists("dcors1"):
        rmtree("dcors1")
    if os.path.exists("outputCors5Comp"):
        rmtree("outputCors5Comp")
    
    stupid_log_files = glob("CORSIKA000000*")
    for f in stupid_log_files:
        os.remove(f)
