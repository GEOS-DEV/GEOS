#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@13.3.1/vp4lpvkssphj4dgtfp75wttn622bvigr
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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-ezuzan2hxi3thik3yquycfp6px3i6p5p" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-7lze2byq355dy3o3xkaltscep6ykfe3e" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-53jbw7cbuvlp6ndacj6i4ke2pssjw3j2" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-oq3gmigoepeyftgmgpauh2bubj2mlahp" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-oolcb7mkc7hgoqmah3zlajcbi6wanc2x" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/adiak-0.4.0-sromwlpvz44dh2m5zjkro4oapef3qbc6" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/hdf5-1.12.1-opvgp2rb5xjnz5as4nkixvudjvq57mty" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-flxhl7p67o53jmexzfg3yfgbqb5ggffg" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/silo-4.11.1-bsd-ojvxcphuoyzdhzivwccmwncxz2duap6m" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/pugixml-1.13-xrny6uxyzgaqjifkwudtcxwqsqhc7e43" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/vtk-9.4.2-j2ldzmygouwz3z3l2a6eapiuwbsqlxs2" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/fmt-10.0.0-3j3bjq4wygq2sox46ywfy66mtxqhgaa5" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/metis-5.1.0-tiuancpjex2f66gnwpnzcfmxzgqvbn3v" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/parmetis-4.0.3-fqn5coe374oz4gclepbycoyv3f2d6bay" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/scotch-7.0.8-sjzijgtweyd3ryn5xuyrudqgebkzvvpr" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-2rlaabejsnrd6wpizwd6kulpflawymcc" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/suite-sparse-5.10.1-apg2yj6qsdbnlhx2643xw6phtuoulg5f" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/trilinos-16.1.0-coa7tsuza4smus6e2t7lbtwhw7clujr2" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-sujc3abytrzh3wey2jpx4kk627tkzlmo" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/doxygen-1.8.20-ka4fk3pa4mkdi7yll3utv3j4yyeumxai/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-qstnvriew3watzlvi3nxwq5ztr5atqts/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-gcc-13_tpls/gcc-13.3.1/mathpresso-geos-umxnfdpmnvgrzbbo5s7gedzjrdla6jcu" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

