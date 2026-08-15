#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@13.3.1/usup2wgpduknmwrwuols33qbmrmhz5dq
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-13.3.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-13.3.1-magic/bin/g++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-13.3.1-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-13.3.1-magic/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/chai-2025.12.0-qbuudn4fk6m7wftxulklbk5cw7ozacjf" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/raja-2025.12.0-g37z6m2jmvyqb7qiwwa37hom3uyp2ypp" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/umpire-2025.12.0-hqtbmsfujzpjxj3kl3jk3bmt3a7k5ety" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/camp-2025.12.0-3q3quidchcepoysqupvlwla577e3smyt" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/caliper-2.14.0-xixabnkvywrgd6a7v766yr74mobs73tx" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/adiak-0.4.0-u3s7u7hqxlxok3vyus63lfao77t3wvdq" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/hdf5-1.14.6-ptzjlsdhytvyf5dxp6js6iugmpfd6pqc" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/conduit-0.9.5-6ockgir74o2gncru2tfhabauss3leyqb" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/silo-4.12.0-zwey5kdj32q7nfz6o6gbgay5mln6qozf" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/pugixml-1.13-vykyowhrjriqfjysczgzijjvo6fbimew" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/vtk-9.4.2-xllx5twzha2typ6vabecuwrdb7dt5vdn" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/fmt-11.0.2-g27pu3mcspwrp46iuo32tkyrrxrhzdrk" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/metis-5.1.0-7lnhkbv2fdcbtcrv5ubv7og3l6cor54b" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/parmetis-4.0.3-3mayiemdzkvic4myatatjp37gx7wmk3k" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/scotch-7.0.8-smdrzzop27f3svwxzwq63aoaol7wwyzu" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/superlu-dist-9.2.1-lkn6ca4m4oejawk3pauulv7qmeiruuhq" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/suite-sparse-5.10.1-avjp6ihuavwtlx4xgtotehbz2noze2fh" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/trilinos-16.1.0-ylg344vhxpyzdj5lb4y4l7lwls3iqrba" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-obgrfae6tjeywyz5xufz77nm44qexv55" CACHE PATH "")

set(ENABLE_HYPREDRV OFF CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-s75vosid4wxvttywirk75zdaqxdb5zw7" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/doxygen-1.8.20-x3azuxonudnfy6nnvzjyovjkjfncioq2/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-wp72eon5hcohsatqntmggcspmd3ibryy/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-12/dane-gcc-13_tpls/gcc-13.3.1/mathpresso-geos-fzt67tuymmzyrl5lhcprtffdhnanu4xq" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

