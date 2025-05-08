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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-7wnlaxzcihkuizymdnnci56uhbaqvdgh" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-nnv4plp73ftllcxgst6szuok76v45g5q" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-jdqm56pwdndks4nqjlvdcbm5f67d5ymf" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-leabuhdy46cwgzrdgy6euqvfzi3sdv6j" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-4o2hij6pbcavyav2djiifne5yoh6judq" CACHE PATH "")

set(adiak_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/adiak-0.4.0-jrqmoii2lzu2imzhwsazd3rujttkgqws/lib/cmake/adiak" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/hdf5-1.12.1-oedmzeyclljshhp6wvflvydyfelztffk" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-v5rwkgv56zdmh23lwyuynu3fdurelckt" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/silo-4.11.1-bsd-wglxufqf4vmt5zsxt37py277yq3g7nqx" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/pugixml-1.13-wx4umca2stjej7clxk3rfbzb2mpdtt5g" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/vtk-9.4.2-jvz4a6eg5ei3ux2mbfxiflo243pg7tt2" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/fmt-10.0.0-mtehyjpjzuugvbuslhqthupvhvwnw4bu" CACHE PATH "")

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

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/metis-5.1.0-4skjvcwya7vkf6kd7q2no2difomuuc7n" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/parmetis-4.0.3-n6yetnfnev3yu2zhcxt3nzsipieufbcl" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/scotch-7.0.3-r6l2zc25vdmkegz7jkf6kxo57trtm32m" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-g5y5xuks34y2uiwjqobhjj2zbdnocfwm" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/suite-sparse-5.10.1-fmerfb2ylan4xzq6oshxe7627xcebbuc" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/trilinos-16.1.0-l7f2twq6lqvaql2ir4iwqub6budkh4if" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/hypre-git.be52325a3ed8923fb93af348b1262ecfe44ab5d2_2.33.0-git.12-jo7kcwt5ayn5qdygl5xmk532qloiprbu" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/lassen-gcc-8-cuda-11_tpls/gcc-8.3.1/mathpresso-geos-wykn2hvzqjjlwjn2j6nh37e5ofsontwq" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--ats jsrun_omp --ats jsrun_bind=packed" CACHE STRING "")

