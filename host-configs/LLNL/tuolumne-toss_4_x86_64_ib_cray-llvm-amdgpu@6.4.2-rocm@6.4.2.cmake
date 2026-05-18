#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: llvm-amdgpu@6.4.2/5a4xtwpgpfug54gw5nc6p4vqsf2364ee
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

# CMake native HIP support rejects the hipcc wrapper as CMAKE_HIP_COMPILER.
# Use ROCm Clang directly and provide the ROCm root for HIP package discovery.
set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.2/llvm/bin/amdclang++" CACHE PATH "")
set(CMAKE_HIP_COMPILER_ROCM_ROOT "/opt/rocm-6.4.2" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")

set(ROCM_PATH "/opt/rocm-6.4.2" CACHE PATH "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-oxvdirkozvm4x66dax5vrsb33exiqfum" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-rdalbu5th2w2cfmm5jsyinousas446ec" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-lrjilgouiothrpcbkvhw5lktce5g2u7t" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-cs5ijhqciobi5ckjaycopg4i4ykuti4j" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-l7xuxz6wxaibkphqevuecd6az4juoq2m" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/adiak-0.4.0-zvp7ufgwijgcae7rhgsjfr74enxlqaot" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/hdf5-1.14.6-ep2obi6q5fhmkuacdopelfze7pzxplik" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-74uuxtwgtjx3s5tlrizirhivpiuqdy5l" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/silo-4.12.0-g76pwu4p32ugpy3ogtigdc6pjdkomuc5" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/pugixml-1.13-yp4oimttnqplvruzsbhdsin3uzxnb4tv" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/vtk-9.4.2-dp2xijiidlvhcuuxesuyxqeubnwlfo4h" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/fmt-10.0.0-f6w3axllcfw3rxb7yo2ldxyuj22ajzub" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/metis-5.1.0-f2zrc2lk64xewglc3r7tnxtfydft567f" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/parmetis-4.0.3-hth7hi7njixzvdus44tgxyf7yvvofq6p" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/scotch-7.0.8-e3opnl7aps7meqctzvpyt7nhaneaawhw" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/superlu-dist-9.2.1-wyp6a7kynbdjmpjqd2ls5e7qtqpl5ihq" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/suite-sparse-5.10.1-u4pyorgxqcum2tteducx4ixleaagwzta" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/hypre-git.804c217767f8e5371c6db34328e5945e40fbcb5c_master-6qi6snkrzkptyynjw3ja2iyd6ftekvdc" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-wovfv3aqsidh37ab57fghaclh3knba72/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/tuolumne-llvm-amdgpu-6.4.2-rocm-6.4.2_tpls/llvm-amdgpu-6.4.2/mathpresso-geos-rphdrakff73drj445fhjlgfscqzhrbxm" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

