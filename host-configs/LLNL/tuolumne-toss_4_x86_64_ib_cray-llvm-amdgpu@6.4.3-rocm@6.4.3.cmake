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

set(BLT_CXX_STD "c++17" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicxx" CACHE PATH "")

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

set(ENABLE_HIP ON CACHE BOOL "")

set(CMAKE_HIP_STANDARD "17" CACHE STRING "")

set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.3/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")

set(ROCM_PATH "/opt/rocm-6.4.3" CACHE PATH "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/chai-2025.12.0-lsam7xrrn2xoyutixbrodqfptgmxmedr" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/raja-2025.12.0-seu5ffy4sb266huku6vnuhkmcaa7dxvo" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/umpire-2025.12.0-fblzq75z5dodltuknm4y5n5fovsoukvf" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/camp-2025.12.0-j2avlxjynm2ox4xrchozqfx56tf2jjl4" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/caliper-2.14.0-bztfgt3etletndi2b7yl3rxgitgbhyoz" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/adiak-0.4.0-vco7ybffav5qujuwixxgfxmmaffnkho3" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hdf5-1.14.6-ombvlnozyrresuh7fqv5invm4yqvusfz" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/conduit-0.9.5-usraryaoqwgfnck5m2s3z7kjv7qxkkbi" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/silo-4.12.0-vnjbftwc546varqz44lqvwlim6eywwca" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/pugixml-1.13-zzs25jilf4nmsol5yrrbdfqlze4hk2oi" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/vtk-9.4.2-zq655fcbq5bdwlj5sshzqyt2pgdg225l" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/fmt-11.2.0-f2suv2dia6zrhi4ncyt44yq3yvpptoao" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/metis-5.1.0-e2xbrr2ghhfa6smpiehrxnzz5rv4apb3" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/parmetis-4.0.3-qaocm7xjf636zgs7egy5wfsh7mvjncal" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/scotch-7.0.8-4ebpprbctqec3la3h4g3vjre77lwsrot" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/superlu-dist-9.2.1-n2grbg575iyxklxtaiszyajjqu7qba5b" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/suite-sparse-5.10.1-vqhne5l2mkrmzrl5ulupdvy64zlczcth" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-tabxldomvumrv6qzsubntlg65e3viit3" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-xptsclqele2ekwezfremc3jrysyjgzxh" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-xgy6bvgcz5pvi6bmeuwnzgmwacbgxpxn/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/mathpresso-geos-h2s2s3flqivtos7pz3nvjpmbx3gt3lmd" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

