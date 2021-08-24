#!/usr/bin/bash

mkdir ${_CONDOR_SCRATCH_DIR}/data
cp /data/user/bsmithers/GolemFit/sources/nuSQuIDS/data/xsections/csms_square.h5 ${_CONDOR_SCRATCH_DIR}/data/.
cp /data/user/bsmithers/NuFSGenMC_nominal.dat ${_CONDOR_SCRATCH_DIR}/data/.
echo "Copied file to ${_CONDOR_SCRATCH_DIR}/data/"
ls ${_CONDOR_SCRATCH_DIR}/data/

source ~/.bashrc

GOLEM
# export HDF5_USE_FILE_LOCKING=false 
python3 /home/bsmithers/software_dev/cascade/sensitivity/generate_all_integrated_fluxes.py $@
rm -r ${_CONDOR_SCRATCH_DIR}/data/
