#######################################
#
# Pangea2 / intlog - Intel 19 build
#
# Load modules in this order:
#   1) cmake/3.23.2
#   2) gcc/8.3.0
#   3) openmpi/3.1.0
#   4) python/3.9.13
#
########################################


set( CONFIG_NAME "pangea2-intlog-gcc@8.3.0" CACHE PATH "" ) 

#######################################
# Compilers
#######################################

set( CMAKE_HOME "/data_local/sw/gcc/RHEL7/8.3.0" CACHE PATH "" )
set( CMAKE_C_COMPILER       "${CMAKE_HOME}/bin/gcc" CACHE PATH "" )
set( CMAKE_CXX_COMPILER     "${CMAKE_HOME}/bin/g++" CACHE PATH "" )
set( CMAKE_Fortran_COMPILER "${CMAKE_HOME}/bin/gfortran" CACHE PATH "" )

set( ENABLE_FORTRAN OFF CACHE BOOL "" FORCE )

set( ENABLE_MPI ON CACHE BOOL "" FORCE )
set( MPI_HOME "/data_local/sw/OpenMPI/RHEL7/4.0.3/gcc/8.3.0" CACHE PATH "" )
set( MPI_C_COMPILER       "${MPI_HOME}/bin/mpicc"   CACHE PATH "" )
set( MPI_CXX_COMPILER     "${MPI_HOME}/bin/mpicxx"  CACHE PATH "" )
set( MPI_Fortran_COMPILER "${MPI_HOME}/bin/mpifort" CACHE PATH "" )
set( MPIEXEC              "${MPI_HOME}/bin/mpirun" CACHE PATH "" )
set( MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "" )

set( ENABLE_MKL ON CACHE BOOL "" )
set( MKL_ROOT /data_local/sw/intel/RHEL7/compilers_and_libraries_2019.3.199/linux/mkl )
set( MKL_INCLUDE_DIRS ${MKL_ROOT}/include CACHE STRING "" )
set( MKL_LIBRARIES ${MKL_ROOT}/lib/intel64/libmkl_intel_lp64.so
                   ${MKL_ROOT}/lib/intel64/libmkl_sequential.so
                   ${MKL_ROOT}/lib/intel64/libmkl_core.so
                   CACHE STRING "" )

set( GEOSX_BUILD_OBJ_LIBS    OFF CACHE BOOL "" )
set( GEOSX_BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE )

#######################################
# GENERAL TPLs
#######################################
set( GEOSX_TPL_DIR "/data/pau901/VDG_GSI/MelvinREY/WorkEnv/GEOSX_repos/thirdPartyLibs/install-pangea2-intlog-gcc@8.3.0-release" CACHE PATH "" FORCE )

set( ENABLE_PYGEOSX           ON CACHE BOOL "" )
set( ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )
set( ENABLE_VTK               ON CACHE BOOL "" FORCE )
set( ENABLE_PVTPackage        ON CACHE BOOL "" FORCE )
set( ENABLE_MATHPRESSO        ON CACHE BOOL "" FORCE )
set( ENABLE_BENCHMARKS        OFF CACHE BOOL "" FORCE )

#######################################
# Solvers dependencies
#######################################
set( ENABLE_HYPRE ON CACHE BOOL "" FORCE )
set( ENABLE_PETSC OFF CACHE BOOL "" FORCE )

#######################################
# Docs & schema generation
#######################################
set( ENABLE_XML_UPDATES ON CACHE BOOL "" FORCE ) 
set( ENABLE_DOCS        ON CACHE BOOL "" FORCE )
set( ENABLE_DOXYGEN     OFF CACHE BOOL "" FORCE )
set( ENABLE_SPHINX      ON CACHE BOOL "" FORCE )
set( ENABLE_UNCRUSTIFY  ON CACHE BOOL "" FORCE )

#######################################
# RAJA/CHAI SETUP
#######################################
option( RAJA_ENABLE_TBB "" OFF )
option( ENABLE_CALIPER  "Enables CALIPER" ON )
set( ENABLE_CUDA        "OFF"       CACHE PATH "" FORCE )
set( CHAI_BUILD_TYPE    "cpu-no-rm" CACHE PATH "" FORCE )
set( CHAI_ARGS          ""          CACHE PATH "" FORCE )
set( ENABLE_OPENMP      "OFF"       CACHE PATH "" FORCE )
set( RAJA_ENABLE_OPENMP "OFF"       CACHE PATH "" FORCE )





include( ${CMAKE_CURRENT_LIST_DIR}/../tpls.cmake )
