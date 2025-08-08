include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/dane-gcc-12.cmake)
set(GEOS_TPL_DIR "${GEOS_TPL_ROOT_DIR}/2024-12-06/install-${CONFIG_NAME}-release" CACHE PATH "Override path" FORCE)

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic CACHE PATH "")

# ATS
set(ATS_ARGUMENTS "--machine slurm112"  CACHE STRING "")

include(${CMAKE_CURRENT_LIST_DIR}/llnl-cpu-base.cmake)
