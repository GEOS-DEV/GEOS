#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: llvm-amdgpu@6.4.2/bttqptgb4i44beu2uawq5hj2dlqmogvr
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-8.3.1/cmake-3.29.2-pj7wkdymawrxifvtbcwpu7mrpzuqy6wk/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/opt/rocm-6.4.2/llvm/bin/amdclang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/opt/rocm-6.4.2/llvm/bin/amdclang++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-9.0.1-rocmcc-6.4.2/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-9.0.1-rocmcc-6.4.2/bin/mpicxx" CACHE PATH "")

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

set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.2/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")

set(ROCM_PATH "/opt/rocm-6.4.2" CACHE PATH "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-ba3iqmvoqdxiqub76d57mxbzxk25zshl" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-qusctgvbt4inmtlui2df7ak5nudfq22p" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-r3kevsc4fo63hqywbuldmygfcvqqht42" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-w7tycuxwwygdldmukzqkozgcobn26ukl" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-3gr4uhu77yfxrmfcdbdaaml6lc7issdm" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/adiak-0.4.0-h34szgtj6jr4ytoa7dq4eeyt6tjopjkf" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/hdf5-1.12.1-icrqfpdr4ym3njhmk27f5ptqnqfi5uks" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-6nwlumanhowitptuwkfnquimin7u57eg" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/silo-4.11.1-bsd-4jhwtwuyfeez7oxzl37psk3oshry6dsm" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/pugixml-1.13-ug5xbnmta4mhtus2py2r5u7iiqv4kqcd" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/vtk-9.4.2-pok4gvwtgi4xhnz57mddo7a2azsbzmv2" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/fmt-10.0.0-dizkjxylsab7lxjm74yyij2wdbbva2p5" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/metis-5.1.0-f7oiyqytqqkvgxuvnzjnw7bfgfpcn5nj" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/parmetis-4.0.3-jvnwbcvgpkgwd3e6haovagcetjenpnnn" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/scotch-7.0.8-f4xm7ddgopligqtfk7vcbiyp5wci6kvg" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-nzrmoolujjo4ty7lbprykmikovz5x4po" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/suite-sparse-5.10.1-vyn45j5myamsgtijhp2f4f6g5b3fsrww" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-3eifpyatzxdbpf6sutvte32vzcb7l76f" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "")

set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/tce/packages/python/python-3.12.2" CACHE PATH "")

set(Python3_EXECUTABLE "/usr/tce/packages/python/python-3.12.2/bin/python3" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-eiztnbtnwd36zq63yzgifpdyama4lxyz/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/mathpresso-geos-ywc6gduc4oh4agkxixjx3al3bcaxgsjt" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

