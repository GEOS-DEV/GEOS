#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@12.1.1/yrpaj5xpatcxrcwfoq54ecujlpt2jmnf
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/g++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicxx" CACHE PATH "")

set(MPIEXEC "srun" CACHE PATH "")

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

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-12.6.0" CACHE PATH "")

set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")

set(CMAKE_CUDA_ARCHITECTURES "90" CACHE STRING "")

set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations -arch sm_90" CACHE STRING "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-wzjalvcalpkeyxe5siy7gwqrbtnnyxxz" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-jbhkrbgaxfn5mzaggta4rk7wv5p2chgk" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-zt4bjb6zfhp2tbd25epz5e42jwcgycaa" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-wczbiqynrets26vxqrsjudbufzvk43px" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-rlix6c6ragjlpu4zz4ianbzm3mocd2xj" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/adiak-0.4.0-jaobjsovvsr5cxhvur4vhbjih46cnyxn" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hdf5-1.14.6-asvxotlx2q56he4mdfdtjcnb4jq62zic" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-kytlb7oa3gohj7endcxxge7e7mt7hvdq" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/silo-4.12.0-aqrzg22nh5byo7vajr7y6iqtegjbfb5e" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/pugixml-1.13-tpluwqvpaqeb2shk27wuse7syhub4dau" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/vtk-9.4.2-hqjv4sjl75rmbodofocrwfcen7o3v2yx" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/fmt-10.0.0-6u2gjlwcibsgjei5ulun3tyfpvyxpfxw" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/metis-5.1.0-cmykznrnhtc5biethufzl6wzrihlpp64" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/parmetis-4.0.3-bhxwfuyf3dlsmdw3slllqqpxftodvhfc" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/scotch-7.0.8-7g27ocdgn4ub55aic7ucj7yuywi2e5vp" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/superlu-dist-9.2.1-giy6hfrd4pyhfrax4ugfdv4ogz2566qu" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/suite-sparse-5.10.1-zt627tebpvlpsxqqxq6szn37r3canc5b" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/trilinos-16.1.0-xaabteedwoghcxz37axn6hbri23e66ey" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hypre-git.804c217767f8e5371c6db34328e5945e40fbcb5c_master-ox2zixiagpynjqojuuzr23z465of6lgi" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

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
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/mathpresso-geos-sh3rclmkawcdazdhnwwncursinf4a534" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

