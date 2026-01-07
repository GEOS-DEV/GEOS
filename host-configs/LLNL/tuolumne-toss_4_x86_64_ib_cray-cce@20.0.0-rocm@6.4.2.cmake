#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: cce@20.0.0/jgeq6ldcfe7mr2jh7ubrodkfeslpjo7v
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-8.3.1/cmake-3.29.2-pj7wkdymawrxifvtbcwpu7mrpzuqy6wk/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/cce-tce/cce-20.0.0/bin/craycc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/cce-tce/cce-20.0.0/bin/crayCC" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-9.0.1-rocmcc-6.4.2-cce-20.0.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/cray-mpich-tce/cray-mpich-9.0.1-rocmcc-6.4.2-cce-20.0.0/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-2xaqayak6p4lymu5vh2vhyp73ykjb4fa" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-6aqiptzuvooxe7hhy565vvth4a5avy77" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-q4rcpcn7l2aysxfkiiwmtlk3nh6brh3r" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-6fj3io27pd6zqd5lw2rr5stajmth5fqr" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-exwybomvhukg2t34cmbzqivgeuoyvngh" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/adiak-0.4.0-anij2zxqxqowlhq73j7yb6bwaz5ysuac" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/hdf5-1.12.1-ghw6mf5cqce7zwgoc6plvrgj2hhkwtjy" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-p5d7nvdubcgoagkmeo544qqnsbyem4es" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/silo-4.11.1-bsd-7h53ibxcoactx4irohzocdv22jpxj2dr" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/pugixml-1.13-nxlztejlbuggihqi6padoxn3qp4ffurq" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/vtk-9.4.2-cz4wehsanworogxcsf7nrwu56oleapvd" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/fmt-10.0.0-mm2fkoftthjyxxlgb4jmalaz7umb75ro" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/metis-5.1.0-6b5cm45inumnfp32i6pnedt4orytni4t" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/parmetis-4.0.3-urqmxwomewy43rj5v3uoiocht563f75v" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/scotch-7.0.8-ayxlcpfy5d5u4oxdlhz5thom2s7lnwqe" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-7tmrsc5dxaxc3xu2k4ztabgn6pcjgv6t" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/suite-sparse-5.10.1-4cqu4zutkycg43klwzjdqsrhvlawhvoq" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-qxcxi4buscxkwiygmerfuzv3tk2hahq5" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-hqx2qehtetfyrovl7nanjnylrb5yeq5u/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/tuolumne-cce-20-rocm-6.4.2_tpls/cce-20.0.0/mathpresso-geos-rwhiq5pnbjlloviwbihn2p5el64ihoqs" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

