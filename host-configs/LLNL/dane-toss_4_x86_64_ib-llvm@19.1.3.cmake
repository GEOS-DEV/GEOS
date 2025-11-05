#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@19.1.3/kmyrajzdclojf3cplhwuj4qub6tx2gcn
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-19.1.3-magic/bin/clang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-19.1.3-magic/bin/clang++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-7vcryfp3qie5hyu6wywgtkc2ylsicztq" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-qowzuw4upjvqnkex5qv7bpxtepanq3gk" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-2xhxuag4uw5khgqagck2kupru42zsqi6" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-jarvyqms2om4vs3gtbkcitgc47nd4kv3" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-h7yk2n54vt4wu4yochiptv5lcwqorltj" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/adiak-0.4.0-t37gghn6ffprh7tfvtyyfgpjfzqg2wuk" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/hdf5-1.12.1-77ah5nvdk4vsp565diomznfciplzgsh5" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-nhtyjmf2lect7mbqaykfrepxud7ye63q" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/silo-4.11.1-bsd-u54hi3o5ghejm5yu7di2gkuxzsqecf54" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/pugixml-1.13-djt6ncbwewwueb53zx22fpcei7tpnu64" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/vtk-9.4.2-4i6sfrxmezqtfql3t5ptie6n74j22hje" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/fmt-10.0.0-wjbhfzxwrphvpdj3a4x5s7axhvxib7w5" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/metis-5.1.0-6reqfcfwkwyh6omphj4jurjr3ppsrw7v" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/parmetis-4.0.3-4cac2qzz62pe5nlaprqgdwuteyjl4vms" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/scotch-7.0.8-mwvdig3nq2noghgnplnx7tjaaz3rzej2" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-l7s6ytczis2wioi2nl3gbfp2iazg3mbv" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/suite-sparse-5.10.1-tgukpvq65zy7bmb2ymm3ulrbkkgph7ji" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/trilinos-16.1.0-27k7hi32qyyo6x62t642cksqjqgny54d" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-pg3a3fmztrjxe3hintcizpc5d2usw6q5" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/doxygen-1.8.20-gasqqxe42cf7tfauydyfrema57c2wugc/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-wyy27wd7hurnl6jgiadguxggm2prusac/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/dane-clang-19_tpls/llvm-19.1.3/mathpresso-geos-qt2kiicwiqndzanst73ngmwt4r2umxms" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

