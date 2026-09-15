#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@13.3.1/kifr3hiw7m2flijsqjc6frr3sazhluhr
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

set(BLT_CXX_STD "c++20" CACHE STRING "")
#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-13.3.1-magic/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-13.3.1-magic/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "srun" CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
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
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/chai-2026.07.0-erllnggecomf3hqxdwxp6miwkz476u3p" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/raja-2026.07.0-pmazaf7uf4dnzti7vxd244qh3fju6fmk" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/umpire-2026.07.1-r4onlnlhnw5onrffszfbavamfoyx6df7" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/camp-2026.07.1-5i3dqmdhnizbf76tgoojjazi77szs2h6" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/caliper-2.14.0-ea5aqg3g2x6fgs3josqj2i3s6cjugcot" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/adiak-0.4.0-ln3qxbtrzibrmsoj2dcnedtjktoeavug" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/hdf5-1.14.6-xoitq72at2yzpsyyfckgr53ggbyohkj5" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/conduit-0.9.5-465lxwaml5s6yhdboptdazqjjkc7c54r" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/silo-4.12.0-ypq6aylz7s4vhsaadexftwo7y4dlwkc4" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/pugixml-1.13-3gjh2dlkrvqjpmtlvmppujdsniiofkeo" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/vtk-9.7.0-5zxp24mu63fko5hxxqolnhhly2q2tjzp" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/fmt-12.1.0-xo5ymndyuse2nczmykaaypwcuwsrsyk5" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/metis-5.1.0-7zomiygwzikvuccil6c7ihjpltbivt5y" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/parmetis-4.0.3-k6cyzwra45uns7fpv5br7vvh7m6t42ra" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/scotch-7.0.8-o3o2e4l775tdp6qupdb4jzzrwyf3v35y" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/superlu-dist-9.2.1-ekjkj6lhitkd5uxjvhsmhb2t3v6vlbif" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/suite-sparse-5.10.1-26fl4m4timzwuen7wffieaiuoryo7tqi" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-ok3wifjeigmvu56jxotepszhoybppeyw" CACHE PATH "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-cd3urpkqkk2nqirveyapkzgqc7p5s7au" CACHE PATH "")
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
set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/doxygen-1.8.20-p3q5ard6nqqhdkxt5ddwfzuex5logv7z/bin/doxygen" CACHE PATH "")
#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")
set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-sxgqpjvzb6fue3ghevunym3kqg2olnri/bin/uncrustify" CACHE PATH "")
#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-13_tpls/gcc-13.3.1/mathpresso-geos-x4u2ycetkh5sj7kl4cza2pxtvda4gr6r" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")
