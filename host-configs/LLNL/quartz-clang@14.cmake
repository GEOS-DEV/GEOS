set( CONFIG "clang@14")

#build-quartz-clang@14-release
string(TOLOWER ${CMAKE_BUILD_TYPE} CMAKE_BUILD_TYPE_LOWER)
set( GEOS_TPL_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../tplInstall-quartz-${CONFIG}-${CMAKE_BUILD_TYPE_LOWER} CACHE PATH " Path to root thirdpartylibrary install ")
set( GEOS_TPL_DIR ${GEOS_TPL_ROOT_DIR} CACHE PATH " GEOS_TPL_DIR is used within the tpls.cmake file ")

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/LLNL/quartz-clang@14.cmake)

# C++
# The "-march=native -mtune=native" which LvArray adds breaks the PVT package.
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "" FORCE)

# Fortran
set(CMAKE_Fortran_COMPILER /usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gfortran CACHE PATH "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")

# MPI
set(MPI_HOME /usr/tce/packages/mvapich2/mvapich2-2.3.6-clang-14.0.6-magic CACHE PATH "")

# Set some MPM compiler options 
set(MPMPROFILER_OPTION 1) 	# 0 off, 1 something, 2 something else.  doesn't work yet
set(MPM_OPTION 0)      		# 0 default, 1 something, 2 something else.  doesn't work yet

configure_file(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/mpmConfigFiles/config.h.in ${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/physicsSolvers/solidMechanics/mpmConfigFiles/config.h @ONLY)
include(${CMAKE_CURRENT_LIST_DIR}/quartz-base.cmake)