#!/bin/bash

export CRAYPE_LINK_TYPE=dynamic
export HDF5_USE_FILE_LOCKING=FALSE
export XTPE_LINK_TYPE=dynamic

# load new cmake
module load cmake/3.14.4

module load git/2.21.0
module load git-lfs/2.8.0

# # load new gcc to support c++14
# module load gcc/8.3.0
# # reload intel + cray with new gcc
# module load intel/19.0.3.199
# module swap craype/2.6.2
# module swap cray-libsci/19.02.1

# module load -f gcc/7.3.0
# # reload intel + cray with new gcc
# module load -f intel/19.0.0.117
# module load -f craype/2.5.18
