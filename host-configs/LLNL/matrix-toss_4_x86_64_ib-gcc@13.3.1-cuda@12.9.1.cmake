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

set(ENABLE_CUDA ON CACHE BOOL "")

set(CMAKE_CUDA_STANDARD "17" CACHE PATH "")

set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-12.9.1" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-jhlxehgtcgsppbs4xj4rwtb3baravkjz" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-56jftnohdf7x3uerkkyh4dpq2bwljqjw" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-wfspghqnsxpiqhrf53dzibwfvwfhkn2h" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-nuit5llrewiawfafdtfy65fbt5vps53f" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-oolcb7mkc7hgoqmah3zlajcbi6wanc2x" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/adiak-0.4.0-sromwlpvz44dh2m5zjkro4oapef3qbc6" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/hdf5-1.12.1-opvgp2rb5xjnz5as4nkixvudjvq57mty" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-flxhl7p67o53jmexzfg3yfgbqb5ggffg" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/silo-4.11.1-bsd-ojvxcphuoyzdhzivwccmwncxz2duap6m" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/pugixml-1.13-xrny6uxyzgaqjifkwudtcxwqsqhc7e43" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/vtk-9.4.2-j2ldzmygouwz3z3l2a6eapiuwbsqlxs2" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/fmt-10.0.0-3j3bjq4wygq2sox46ywfy66mtxqhgaa5" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/metis-5.1.0-tiuancpjex2f66gnwpnzcfmxzgqvbn3v" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/parmetis-4.0.3-fqn5coe374oz4gclepbycoyv3f2d6bay" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/scotch-7.0.8-sjzijgtweyd3ryn5xuyrudqgebkzvvpr" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-2rlaabejsnrd6wpizwd6kulpflawymcc" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/suite-sparse-5.10.1-apg2yj6qsdbnlhx2643xw6phtuoulg5f" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/trilinos-16.1.0-coa7tsuza4smus6e2t7lbtwhw7clujr2" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-d2io6ghkh74z772roqcm7w7vsrsk4uj7" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-10-30/matrix-gcc-13-cuda-12.9_tpls/gcc-13.3.1/mathpresso-geos-umxnfdpmnvgrzbbo5s7gedzjrdla6jcu" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

