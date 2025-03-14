#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: blueos_3_ppc64le_ib_p9
# Compiler Spec: gcc@=8.3.1
# CMake executable path: /usr/tce/packages/cmake/cmake-3.29.2/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-8.3.1/bin/g++" CACHE PATH "")

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

set(MPI_C_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-8.3.1/bin/mpicxx" CACHE PATH "")

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

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-11.8.0" CACHE PATH "")

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

set(CHAI_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/chai-git.df7741f1dbbdc5fff5f7d626151fdf1904e62b19_develop-jvnmbz3nfrjym77m2m2ddxfgyigf2hqv" CACHE PATH "")

set(RAJA_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/raja-git.4d7fcba55ebc7cb972b7cc9f6778b48e43792ea1_develop-ktim3j27la76tmlr6ryj2vjfllltmz4o" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/umpire-git.abd729f40064175e999a83d11d6b073dac4c01d2_develop-xmi5nbcsff7dz76ljlmr5catjxa4hzjn" CACHE PATH "")

set(CAMP_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/camp-git.0f07de4240c42e0b38a8d872a20440cb4b33d9f5_main-d4ymmapzrvlzw6jwjiyenqd5l2sy7pg7" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-rw6xo4fv4xn4oapuiypwwn7cz3m73naz" CACHE PATH "")

set(adiak_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/adiak-0.4.0-n3yugcwz5o3ggvi2nrwybgmfattinedh/lib/cmake/adiak" CACHE PATH "")

set(HDF5_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/hdf5-1.12.1-obzanrc4jpjjzno25h4iohq6f3hr2j27" CACHE PATH "")

set(CONDUIT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-ubpigj3aj4jsqbi7j7bbibmg2ryo7xoj" CACHE PATH "")

set(SILO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/silo-4.11.1-bsd-roheyk37frcivqqo4l3a4ldq7d5hx6st" CACHE PATH "")

set(PUGIXML_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/pugixml-1.13-s3vycldwlzcj5s6pua2hyghvdoqm6hu3" CACHE PATH "")

set(VTK_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/vtk-9.3.1-xoeo4iubpe3qwprc7o2j5lglbp4y6m7t" CACHE PATH "")

set(FMT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/fmt-10.0.0-6jiaf3rryjoxxkptvy2jcrnoajdpe6mv" CACHE PATH "")

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

set(METIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/metis-5.1.0-nn7ngjgltv2ijx4kswtnoeyfgwbq24k2" CACHE PATH "")

set(PARMETIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/parmetis-4.0.3-chbckvoyrbustqebbsxhrrw7e3ngjx43" CACHE PATH "")

set(SCOTCH_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/scotch-7.0.3-4n72uptvwu74yw3nplp2dbyhxt5if6jf" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-co27fkcshr4pdskjot7sa3bfhccymaws" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/suite-sparse-5.10.1-rgp2xi4ihgkbjm2rzgspizeler7w6otm" CACHE PATH "")

set(TRILINOS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/trilinos-16.0.0-k35ukkbitl64gzvdfe2dlvt5sg6r7t67" CACHE PATH "")

set(HYPRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/hypre-git.c893886d15eb57e87dd36efec23693ece3ddc88e_2.32.0-git.4-mesu7comouhq5aiqrdoj6p5hqyyu2t2f" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-02-28/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/mathpresso-geos-tf3bai7si63pwb6nuqqkrliyiipl7z2m" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ATS_ARGUMENTS "--ats jsrun_omp --ats jsrun_bind=packed" CACHE STRING "")

