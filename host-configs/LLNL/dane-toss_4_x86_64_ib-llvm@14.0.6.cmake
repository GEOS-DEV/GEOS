#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@14.0.6/4vnsgoicqdfa2y3y4fimnxwn4spjab4h
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

set(BLT_CXX_STD "c++17" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-14.0.6-magic/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/chai-2025.12.0-urbcgmtlhlhrr6vzvuhqw3tujm7lipyx" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/raja-2025.12.0-kdz2sziciypij2bhyuvhxbxnyberi235" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/umpire-2025.12.0-oe6irwruzadu5coxjogtqmwfahukhyqq" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/camp-2025.12.0-7x4s3ppkngxucz77r6xjaar5pomko27n" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/caliper-2.14.0-4mecvrawwuimridt2vcmcvq3xrttryik" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/adiak-0.4.0-gmrdaggcchlixqxe6ngceealawcebbb7" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/hdf5-1.14.6-nmv7kx5odiuczkygx3vrwmsur3fpnspk" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/conduit-0.9.5-yxn2ze2pyp7leldjzzrtorjwsmzsfoym" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/silo-4.12.0-jpyl4bi46pjlb3vtmrcm2ievdaxico6f" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/pugixml-1.13-qxmbxcyf6esso4yjb5377fkaq5suab6v" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/vtk-9.4.2-wtc5dpvsaqj2bsn2hv26h4tq6kn3la6u" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/fmt-11.2.0-dezyl57k6wjfhfuu6hmtgjqzxyiq7vsw" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/metis-5.1.0-jqfs2lt52bwl4zwjughckbhrvs54hqpm" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/parmetis-4.0.3-ccsm27ijmnhatpsods57uv3nk4x5a3lo" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/scotch-7.0.8-vkrey3rhteittppbr7q3ddrwftxtqiwt" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/superlu-dist-9.2.1-75nhlzs236zvqllhleo7ssnmkpwt64xd" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/suite-sparse-5.10.1-spiteyzsnusnpctx737runf2qdcvxndh" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-kgpct6wzeai5orhj2g6cyqsclp3gack2" CACHE PATH "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-iq6jnjrhetz4v7dpy66x5jti3wvfoghy" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/doxygen-1.8.20-kqb5aukvkkn5rkwpjdszechtcthzuv7o/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-clappbh25vxeky7k4lr3qurtcrjkwjkk/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-14_tpls/llvm-14.0.6/mathpresso-geos-sn363eudk74prv5hmgok5kchbo3l243l" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

