#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: llvm-amdgpu@6.4.3/ttoqi52mggogog43da5bbe6cdmups5s3
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-8.3.1/cmake-3.24.2-ywx52e32uh6gkxzuyubpwkulzgdvxyh6/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/opt/rocm-6.4.3/llvm/bin/amdclang" CACHE PATH "")
set(CMAKE_CXX_COMPILER "/opt/rocm-6.4.3/llvm/bin/amdclang++" CACHE PATH "")
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
set(MPI_C_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicxx" CACHE PATH "")
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

set(ENABLE_HIP ON CACHE BOOL "")
set(CMAKE_HIP_STANDARD "20" CACHE STRING "")
set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.3/bin/hipcc" CACHE PATH "")
set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")
set(GPU_TARGETS "gfx942" CACHE STRING "")
set(AMDGPU_TARGETS "gfx942" CACHE STRING "")
set(ROCM_PATH "/opt/rocm-6.4.3" CACHE PATH "")
#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/chai-2026.07.0-bhdqmyan3vjecugnknsqpyyl4jvc7axb" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/raja-2026.07.0-naplhhdojltk7zaazpdiycdzrc3lgzfc" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/umpire-2026.07.1-enbbftzynvwitnpbcpwbg4rvm2jmibd5" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/camp-2026.07.1-afsrg42mxsrd4tixk5btnv5lelv3vobo" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/caliper-2.14.0-6lzeq264s6suyo6smnytk66jrl4ehhll" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/adiak-0.4.0-yuqe5ylmsnuy2iifh232fv6r6e323jsk" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hdf5-1.14.6-sadc6b7korndde4egvmwkh635uttq4to" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/conduit-0.9.5-2xhuqhsv4nfzfuwoevtd3hezl2weyslf" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/silo-4.12.0-nukovk6iwlltcuphtha5v3e2xbyzi4o6" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/pugixml-1.13-zzs25jilf4nmsol5yrrbdfqlze4hk2oi" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/vtk-9.7.0-76tdcqygh2hhk2ycttgbiexcca7dzsvn" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/fmt-12.1.0-tziq2pevpj5nfqv66dlpt6rxx7haovbk" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/metis-5.1.0-e2xbrr2ghhfa6smpiehrxnzz5rv4apb3" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/parmetis-4.0.3-sluwoamkj5m7sfyrpc3ax2c5bs5yi4oo" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/scotch-7.0.8-qi73qny4dr3bmomkbxk6y7bvpdp4ttyb" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/superlu-dist-9.2.1-ko3qhelgleaejtxri4uw457ta2iv4f2h" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/suite-sparse-5.10.1-cethsrccixxuyxbri3qee4cwjvs4zhem" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypre-git.f1374fb6182c9e730abaa82f865a87b36f9ad50a_master-z5j7kwxg6zy67hrvgizpotcmcnjez4z3" CACHE PATH "")
set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypredrive-git.98989bef31e865d738c1b678ee058bb0e5dd635e_master-3uzl2fr565l5fpbubqofxofrk2pfifus" CACHE PATH "")
set(ENABLE_PETSC OFF CACHE BOOL "")
set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "")
set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")
#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/tce/packages/python/python-3.9.12" CACHE PATH "")
set(Python3_EXECUTABLE "/usr/tce/packages/python/python-3.9.12/bin/python3" CACHE PATH "")
set(ENABLE_PYGEOSX OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(ENABLE_DOCS OFF CACHE BOOL "")
set(ENABLE_DOXYGEN OFF CACHE BOOL "")
set(ENABLE_SPHINX OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")
set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-xgy6bvgcz5pvi6bmeuwnzgmwacbgxpxn/bin/uncrustify" CACHE PATH "")
#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/mathpresso-geos-h2s2s3flqivtos7pz3nvjpmbx3gt3lmd" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
