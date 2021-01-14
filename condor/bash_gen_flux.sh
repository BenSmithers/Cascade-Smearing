#!/usr/bin/bash

mkdir ${_CONDOR_SCRATCH_DIR}/data
cp /data/user/bsmithers/GolemFit/sources/nuSQuIDS/data/xsections/csms_square.h5 ${_CONDOR_SCRATCH_DIR}/data/.

source ~/.bashrc

GOLEM
# export HDF5_USE_FILE_LOCKING=false 
python3 /home/bsmithers/software_dev/cascade/generate_fluxes.py $@
