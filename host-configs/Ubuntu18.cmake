# file: your-platform.cmake

# detect host and name the configuration file
site_name(HOST_NAME)
set(CONFIG_NAME "your-platform" CACHE PATH "")
message( "CONFIG_NAME = ${CONFIG_NAME}" )

# set paths to C, C++, and Fortran compilers. Note that while GEOSX does not contain any Fortran code,
# some of the third-party libraries do contain Fortran code. Thus a Fortran compiler must be specified.
set(CMAKE_C_COMPILER "/usr/local/bin/gcc" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/local/bin/g++" CACHE PATH "")
set(CMAKE_Fortran_COMPILER "/usr/local/bin/gfortran" CACHE PATH "")
set(ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)

# enable MPI and set paths to compilers and executable.
# Note that the MPI compilers are wrappers around standard serial compilers.
# Therefore, the MPI compilers must wrap the appropriate serial compilers specified
# in CMAKE_C_COMPILER, CMAKE_CXX_COMPILER, and CMAKE_Fortran_COMPILER.
set(ENABLE_MPI ON CACHE PATH "")
set(MPI_C_COMPILER "/home/test/mpich-install/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/home/test/mpich-install/bin/mpicxx" CACHE PATH "")
set(MPI_Fortran_COMPILER "/home/test/mpich-install/bin/mpifort" CACHE PATH "")
set(MPIEXEC "/home/test/mpich-install/bin/mpirun" CACHE PATH "")

# disable CUDA and OpenMP
set(CUDA_ENABLED "OFF" CACHE PATH "" FORCE)
set(ENABLE_OPENMP "OFF" CACHE PATH "" FORCE)

# enable PAMELA and PVTPackage
set(ENABLE_PAMELA ON CACHE BOOL "" FORCE)
set(ENABLE_PVTPackage ON CACHE BOOL "" FORCE)

# enable tests
set(ENABLE_GTEST_DEATH_TESTS ON CACHE BOOL "" FORCE )


set(HDF5_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/hdf5" CACHE PATH "")
set(CONDUIT_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/conduit" CACHE PATH "")
set(SILO_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/silo" CACHE PATH "")
set(PUGIXML_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/pugixml" CACHE PATH "")
set(RAJA_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/raja" CACHE PATH "")
set(UMPIRE_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/chai/share/umpire" CACHE PATH "")
set(CHAI_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/chai" CACHE PATH "")
set(VTK_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/vtk" CACHE PATH "")
set(MATHPRESSO_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/mathpresso" CACHE PATH "")
set(METIS_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/metis" CACHE PATH "")
set(PARMETIS_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/parmetis" CACHE PATH "")
set(SUPERLU_DIST_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/superlu_dist" CACHE PATH "")
set(SUITESPARSE_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/suitesparse" CACHE PATH "")
set(HYPRE_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/hypre" CACHE PATH "")
set(TRILINOS_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/trilinos" CACHE PATH "")
set(FMT_DIR "/home/test/GEOSX1/thirdPartyLibs/install-Ubuntu18-release/fmt" CACHE PATH "")






