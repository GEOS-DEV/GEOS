#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: blueos_3_ppc64le_ib_p9
# Compiler Spec: clang@=13.0.1
# CMake executable path: /usr/tce/packages/cmake/cmake-3.29.2/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/clang/clang-13.0.1-gcc-8.3.1/bin/clang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/clang/clang-13.0.1-gcc-8.3.1/bin/clang++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-13.0.1-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-clang-13.0.1-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

set(ENABLE_WRAP_ALL_TESTS_WITH_MPIEXEC ON CACHE BOOL "")

set(MPIEXEC "lrun" CACHE PATH "")

set(MPIEXEC_NUMPROC_FLAG "-n" CACHE PATH "")

#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA ON CACHE BOOL "")

set(CMAKE_CUDA_STANDARD "17" CACHE PATH "")

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-12.0.0" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_ARCHITECTURES "70" CACHE STRING "")

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -arch sm_70" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG -Xcompiler -DNDEBUG -Xcompiler -O3 -Xcompiler -mcpu=powerpc64le -Xcompiler -mtune=powerpc64le" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_RELWITHDEBINFO "-g -lineinfo ${CMAKE_CUDA_FLAGS_RELEASE}" CACHE STRING "")

set(CMAKE_CUDA_FLAGS_DEBUG "-g -G -O0 -Xcompiler -O0" CACHE STRING "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/chai-git.df7741f1dbbdc5fff5f7d626151fdf1904e62b19_develop-nmmjzxisecb5fpqfnvsm7q25riyc4v6n" CACHE PATH "")

set(RAJA_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/raja-git.4d7fcba55ebc7cb972b7cc9f6778b48e43792ea1_develop-54zp6bxl3yepw3nh3awk455q22vakj2u" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/umpire-git.abd729f40064175e999a83d11d6b073dac4c01d2_develop-ztjsbkdojxcrgqnqrkvm7wq5oz4lm7cw" CACHE PATH "")

set(CAMP_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/camp-git.0f07de4240c42e0b38a8d872a20440cb4b33d9f5_main-fs7yeokq6dq3pql5hndzczn36ykkfhon" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-4gu3btnlno2sbv56vhd3ux2uvgoupkgi" CACHE PATH "")

set(adiak_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/adiak-0.4.0-hqwcvthhjo5zs4b3u52glnlcxk5bp34n/lib/cmake/adiak" CACHE PATH "")

set(HDF5_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/hdf5-1.12.1-dycb4qq65jod3vakmo4jfcknwbba4qp5" CACHE PATH "")

set(CONDUIT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-7nrioisimbmm3mqcl7uenprhkz36og6u" CACHE PATH "")

set(SILO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/silo-4.11.1-bsd-bdktrlhcoynsflkrcgumryqbtu45gezb" CACHE PATH "")

set(PUGIXML_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/pugixml-1.13-dn3wo74qjiclno55tao7i6q3qsmmtbmm" CACHE PATH "")

set(VTK_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/vtk-9.3.1-hjvim2htouza4aw6laupkeluggn7vzzl" CACHE PATH "")

set(FMT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/fmt-10.0.0-u3hormd7jprgxpemsdwy67f5wdd3rcks" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_ESSL ON CACHE BOOL "")

set(ESSL_INCLUDE_DIRS "/usr/tcetmp/packages/essl/essl-6.3.0.2/include" CACHE PATH "")

set(ESSL_LIBRARIES /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/libessl.so
                   /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/liblapackforessl.so
                   /usr/tcetmp/packages/essl/essl-6.3.0.2/lib64/liblapackforessl_.so CACHE STRING "")

set(FORTRAN_MANGLE_NO_UNDERSCORE ON CACHE BOOL "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/metis-5.1.0-yk3zg6ifqxkue25iqr3waysqbkzjbiim" CACHE PATH "")

set(PARMETIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/parmetis-4.0.3-hfu52rbhm3xhhbi33unvxrsydvvdqejd" CACHE PATH "")

set(SCOTCH_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/scotch-7.0.3-v3r45jclgkkywzyft37lqqoguhperyxm" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-2zru2u6fegs5bsz44wckbxcaxuwwdivw" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/suite-sparse-5.10.1-rxlyv2tiek7hl3ln7ps7ykguwog3tax4" CACHE PATH "")

set(TRILINOS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/trilinos-16.0.0-wo5v5rma22soustvdev64v5thopnxpct" CACHE PATH "")

set(HYPRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/hypre-git.21e5953ddc6daaa24699236108866afa597a415c_2.32.0-git.33-sh2e22ra27h65kqhgn4fikhhbg6tjqjt" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python/" CACHE PATH "")

set(Python3_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/lassen-gcc-python/python/bin/python3" CACHE PATH "")

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
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-03-21/lassen-clang-13-cuda-12_tpls/clang-13.0.1/mathpresso-geos-z4ebytw2o2tmltrko7jwbko5x6jvadod" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--ats jsrun_omp --ats jsrun_bind=packed" CACHE STRING "")

