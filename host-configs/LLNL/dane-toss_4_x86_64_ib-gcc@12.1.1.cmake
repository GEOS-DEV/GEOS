#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@12.1.1/tip5n7rqr2qzxyu675utpsa2hwmlcf7u
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/g++" CACHE PATH "")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -DNDEBUG" CACHE STRING "")

set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "")

#--------------------------------------------------------------------------------
# CMake Standard
#--------------------------------------------------------------------------------

set(BLT_CXX_STD "c++17" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicxx" CACHE PATH "")

set(MPIEXEC "srun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# ROCm/HIP
#--------------------------------------------------------------------------------

set(ENABLE_HIP OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/chai-2025.12.0-gj4cge5skn5pmcq7h4rd7vg47zygtx63" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/raja-2025.12.0-7tz3a6vkc7e642t7zmpxexyi5lqsfafg" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/umpire-2025.12.0-qxwp33yhikgs2fkl2vdeobbmal6bqw6n" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/camp-2025.12.0-v3got3dbxmrizv3nk7m6gftc37aen4ym" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/caliper-2.14.0-lwnj4cxlrnclq564deuu7tljugndvtnt" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/adiak-0.4.0-6lz2sfblhg3haec5wghtcahttyei5h3v" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/hdf5-1.14.6-7n343tz7mztayupsigyaj7chzbsltayo" CACHE PATH "")

set(ENABLE_ADIOS2 ON CACHE BOOL "")

set(ADIOS2_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/adios2-2.12.1-ohdg2rgdnmmxbiipnwx72yhwergwn7ek" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/conduit-0.9.5-j25t5t5kvxjhbt3vzvrig6ejqyjuy3pl" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/silo-4.12.0-wo3vlynxc3jkwkunrlikk3645pg76gog" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/pugixml-1.13-vnpomoy3m4sh7zsnz5uqfn2hp6zh2mzw" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/vtk-9.4.2-mu24yirl5ddhfv3qzhhbwilxxsjn47dp" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/fmt-11.0.2-fbj372jbpfsajq5bzeyhdj4ouz4pgbzm" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/metis-5.1.0-pxv4jh6tetihokv5u7keagd7cytygua7" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/parmetis-4.0.3-rbzv62lqv7bqm54x5ofwvlgw43wjcalo" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/scotch-7.0.8-eqsq7le3z5srqx2ezsez4yvwu2vq5n6g" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/superlu-dist-9.2.1-n2ffm6w6nh6s7lc3puw6uiffldvtuto4" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/suite-sparse-5.10.1-d6clmff6msrxrgoecsx6one2cbju6a4t" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/trilinos-16.1.0-bhndzpfd5ytpinku5xk54weqiprrc7mq" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-l3nuajfedgl7ova5u2fheltf6euxupvq" CACHE PATH "")

set(ENABLE_HYPREDRV OFF CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-6fr7ry746mawo3hg35sap7j43nfijbzi" CACHE PATH "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "")

set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/" CACHE PATH "")

set(Python3_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/python3" CACHE PATH "")

set(ENABLE_PYGEOSX ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(SPHINX_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/doxygen-1.8.20-lbqhlrg3msfivysu2rfxq674siosbmob/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-yqfbyrbgcnjjpkfaijkdtsejzjtpp4ig/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-12_tpls/gcc-12.1.1/mathpresso-geos-iz4kph2ijg7inklauzcits2sivybrttp" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

