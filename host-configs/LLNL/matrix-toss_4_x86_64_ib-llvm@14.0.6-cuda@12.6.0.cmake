#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@14.0.6/ul5bul3z5rs5ykers5jxmtuyaovscvqr
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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-f53bhm6scvibecphxaahannyb67ne4gn" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-36pk2ixifqzhhqrloqbvev475eroq2ns" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-j2evmfqon42guvjuz5zxlqht523xgl7k" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-eqycl5hesixa2mzurda2uw4rqpu3xjho" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-p5aeqjrvvexyrdal4epo7it2fbdgzzxk" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/adiak-0.4.0-7fjc7gjmpgx4jytzryqlsg2anuv27xcg" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hdf5-1.14.6-n4d653mlykuikmfab3ss6ruvejw24dml" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-l4n3oasabottshhvnmpynvnr5qv4m5ym" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/silo-4.12.0-m72pwkottbfattkmfy5et4y2zexnt3hv" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/pugixml-1.13-qnpx3zh4q3eouffejz63xhyrdjvjkitl" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/vtk-9.4.2-fqjvu6cakp4mkhr7s5ypn5f2ftbka5hg" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/fmt-10.0.0-qbvletrbia2c5guh6kzewyt6i5ua4q22" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/metis-5.1.0-y7tsdri3aun762o4qjxji66y6ud4yjsk" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/parmetis-4.0.3-uvoynbg75vl2fda5gupkf572btgexrw7" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/scotch-7.0.8-5utttqumhe2j2kb2swrnxt6iebuectv4" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/superlu-dist-9.2.1-2da4io7hvinici6gompzakkl2o2mqzcx" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/suite-sparse-5.10.1-dccejparafjkgc4jjp2bl3hwhvjz4xtg" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/trilinos-16.1.0-oouxqanvjkbleahqbwwkslj55tpjygym" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hypre-git.804c217767f8e5371c6db34328e5945e40fbcb5c_master-dotoyonkjyzkvzh74aojacvlzuwblxbt" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-01-28/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/mathpresso-geos-4npqxs2paek425irgmzf2s5bniqoqkg2" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

