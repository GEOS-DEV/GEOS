#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@14.0.6/o76uk6nkkf5espsolkoeq4aanxuyafiq
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-14.0.6-magic/bin/clang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-14.0.6-magic/bin/clang++" CACHE PATH "")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native -mtune=native" CACHE STRING "")
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
set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "srun" CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")
#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA ON CACHE BOOL "")
set(CMAKE_CUDA_STANDARD "20" CACHE STRING "")
set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-12.6.0" CACHE PATH "")
set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")
set(CMAKE_CUDA_ARCHITECTURES "90" CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")
set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")
#--------------------------------------------------------------------------------
# ROCm/HIP
#--------------------------------------------------------------------------------

set(ENABLE_HIP OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/chai-2026.07.0-ncjkh2hbgzyeudhjmpqs7xvfbjyt7dk5" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/raja-2026.07.0-z2p727ohmjfmwlpwu52mtg3dccfthczi" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/umpire-2026.07.1-vtmsjk54euoib5fg72id36l42bi3ffnb" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/camp-2026.07.1-oxknvdxcwzyjcfwbbvaiiz2halwbzeyt" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/caliper-2.14.0-5rgrh6rzvsxs7rp3fbiey753drmsezcv" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/adiak-0.4.0-dbgiv6uwizl66u76tgktgqrsg7c4ucpa" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hdf5-1.14.6-vfjzdpmm7ombhjdvwvmnn5w53twxowpg" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/conduit-0.9.5-5ayszdj34sbd45onwl36aqmem7h7wwqi" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/silo-4.12.0-hvlrfpptyilahrs75yfxz3bc2ssmqoza" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/pugixml-1.13-k3cwmmpsh6booq32x2baxgpqxdfz7zlq" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/vtk-9.7.0-efw6p6yu6rdixhflf7wzksy6hahzwgl3" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/fmt-12.1.0-y23hj65ucudnadoh2ebjwwsmec4st2wu" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/metis-5.1.0-leqxeahwgia3lo4ibl354h2nvq5x2lxe" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/parmetis-4.0.3-3pgxn2lv3bb3wdsaaxay3uqnvrlo64id" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/scotch-7.0.8-ldqysttu2hgyxwvohfnppdm4hqxppcq2" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/superlu-dist-9.2.1-tysxw3ozi5yb2xztlli22e3fqt4jvzkw" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/suite-sparse-5.10.1-wrk3dc7zkozobil6vj5qz4jahy5srx4u" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hypre-git.f1374fb6182c9e730abaa82f865a87b36f9ad50a_master-mofmhaxtb57fugykarlwennewrhmtmup" CACHE PATH "")
set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hypredrive-git.98989bef31e865d738c1b678ee058bb0e5dd635e_master-d2kqrseau6gjhrcaxetspogfhxdzph4f" CACHE PATH "")
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

set(ENABLE_DOCS OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")
set(ENABLE_SPHINX OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/mathpresso-geos-vk4hgtaedmc754nikl57m7damhndcyve" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")
