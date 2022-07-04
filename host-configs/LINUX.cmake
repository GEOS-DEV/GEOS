# file : LINUX.cmake
# detect host and name the configuration file
site_name(HOST_NAME)
set(CONFIG_NAME "LINUX" CACHE PATH "")
message( "CONFIG_NAME = ${CONFIG_NAME}" )
# set paths to C, C++, and Fortran compilers. Note that while GEOSX_env does not contain any Fortran code,
# some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
set(CMAKE_C_COMPILER "/usr/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/bin/gfortran" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -std=c++2a -march=native -mtune=native" CACHE STRING "")
set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -DNDEBUG -march=native -mtune=native" CACHE STRING "")
# enable MPI and set paths to compilers and executable.
# Note that the MPI compilers are wrappers around standard serial compilers.
# Therefore, the MPI compilers must wrap the appropriate serial compilers specified
# in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/usr/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/usr/bin/mpif90" CACHE PATH "")
set(MPIEXEC "/usr/bin/mpirun" CACHE PATH "")
# disable CUDA and OpenMP
set(CUDA_ENABLED "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)
# enable PAMELA and PVTPackage
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)
set(ENABLE_HDF5 ON CACHE BOOL "" FORCE)
set(ENABLE_CONDUIT ON CACHE BOOL "" FORCE)
set(ENABLE_SILO ON CACHE BOOL "" FORCE)
set(ENABLE_HYPRE ON CACHE BOOL "" FORCE)
set(ENABLE_UNCRUSTIFY ON CACHE BOOL "" FORCE)
# path to the installed TPL's directory
#set(GEOSX_TPL_DIR "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/" CACHE PATH "")
set(HDF5_DIR         "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/hdf5" CACHE PATH "")

set(CONDUIT_DIR      "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/conduit" CACHE PATH "")

set(SILO_DIR         "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/silo" CACHE PATH "")

set(PUGIXML_DIR      "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/pugixml" CACHE PATH "")

set(RAJA_DIR         "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/raja" CACHE PATH "")

set(CHAI_DIR         "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/chai" CACHE PATH "")

set(UMPIRE_DIR       "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/chai" CACHE PATH "")

set(ADIAK_DIR        "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/adiak" CACHE PATH "")

set(ASTYLE_DIR       "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/astyle" CACHE PATH "")

set(AXOM_DIR         "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/axon" CACHE PATH "")

set(CALIPER_DIR      "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/caliper" CACHE PATH "")

set(MATHPRESSO_DIR   "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/mathpresso" CACHE PATH "")

set(METIS_DIR        "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/metis" CACHE PATH "")

set(PARMETIS_DIR     "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/parmetis" CACHE PATH "")

set(SUPERLU_DIST_DIR "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/superlu_dist" CACHE PATH "")

set(SUITESPARSE_DIR  "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/suitesparse" CACHE PATH "")

set(HYPRE_DIR        "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/hypre" CACHE PATH "")

set(VTK_DIR          "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/vtk" CACHE PATH "")

set(UNCRUSTIFY_DIR   "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/uncrustify" CACHE PATH "")


set(FMT_DIR   "/home/m3d/aurelie/codes/thirdPartyLibs/install-LINUX-release/fmt" CACHE PATH "")


# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )
# disable TRILINOS
set(ENABLE_TRILINOS "OFF" CACHE PATH "" FORCE )
