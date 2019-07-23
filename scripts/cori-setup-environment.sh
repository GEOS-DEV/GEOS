#!/bin/bash

# export environment variables
export CRAYPE_LINK_TYPE=dynamic
export HDF5_USE_FILE_LOCKING=FALSE
export XTPE_LINK_TYPE=dynamic

# load new cmake
module load cmake/3.11.4

# load new gcc to support c++14
module load gcc/7.3.0
# reload intel + cray with new gcc
module load intel/19.0.0.117
module load craype/2.5.16