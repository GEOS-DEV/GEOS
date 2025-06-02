include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/tioga-cce-19-rocm-6.4.0.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/amdgpu-base.cmake)

# MPI
set(MPI_HOME /opt/cray/pe/mpich/8.1.33.1/ofi/crayclang/18.0 CACHE PATH "")
