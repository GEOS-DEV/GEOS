#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@19.1.3/ae3qod26lywpyu252dwhujh5dusql3mp
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

set(BLT_CXX_STD "c++20" CACHE STRING "")
#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "srun" CACHE STRING "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")
#--------------------------------------------------------------------------------
# OpenMP
#--------------------------------------------------------------------------------

set(ENABLE_OPENMP ON CACHE BOOL "")
#--------------------------------------------------------------------------------
# Cuda
#--------------------------------------------------------------------------------

set(ENABLE_CUDA ON CACHE BOOL "")
set(CMAKE_CUDA_STANDARD "20" CACHE STRING "")
set(CUDA_TOOLKIT_ROOT_DIR "/usr/tce/packages/cuda/cuda-12.9.1" CACHE PATH "")
set(CMAKE_CUDA_COMPILER "${CUDA_TOOLKIT_ROOT_DIR}/bin/nvcc" CACHE PATH "")
set(CMAKE_CUDA_ARCHITECTURES "90" CACHE STRING "")
set(CMAKE_CUDA_FLAGS "-restrict --expt-extended-lambda -Werror cross-execution-space-call,reorder,deprecated-declarations" CACHE STRING "")
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
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/chai-2026.07.0-dsac5qbpyupfo2dr5ttjgmbiefxn7b3k" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/raja-2026.07.0-enesfqspzuazpp7guqyoc3ceb2x5eedf" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/umpire-2026.07.1-22rzxotzft3h2mxlbekteiveorbywt6r" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/camp-2026.07.1-u3vorarioz2p7pquj3rhepapaitrpplt" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/caliper-2.14.0-ca4cdk2vf73oij5gt3gyetizhuhmeh4u" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/adiak-0.4.0-mtsop6dnpc3vb7aejlrhngdxmwk3mufc" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/hdf5-1.14.6-y2impgpjstmtreab63cq7e3a2lga6ytp" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/conduit-0.9.5-3if4nxrqxrzv5r5hndoaqb7pb5nhjg5v" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/silo-4.12.0-3kyu3a3anmija4kgqeuvjcw6bb2gwacc" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/pugixml-1.13-xks5inlgbhx6d5cxjuhv73zcgihwx5dh" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/vtk-9.7.0-w3suyneqpeekwwbpzzajtuihcvpyjhsi" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/fmt-12.1.0-twujbdwfgpvuj3jmjhnnitsxuleoy5fl" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/metis-5.1.0-4tacinv42rn7o7gqpdoheisytvza4llm" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/parmetis-4.0.3-4rrg4at5mrqqifeio4s25vwqwx5wthab" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/scotch-7.0.8-ssll6ps7tobxvrm4zzubyk5yvihgxwzn" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/superlu-dist-9.2.1-vmldzue2ikjwvzcsj6gmusk4bdd5elic" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/suite-sparse-5.10.1-j2h5wmk2fwnhhvo7ceqrgo7o2lrbhh6k" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/hypre-git.f1374fb6182c9e730abaa82f865a87b36f9ad50a_master-yaoeizw5orew2nkiinp6culhnd64ya6v" CACHE PATH "")
set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/hypredrive-git.98989bef31e865d738c1b678ee058bb0e5dd635e_master-lu3mzdo5vzueh5onmzdkbqmif7zckclk" CACHE PATH "")
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
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/matrix-llvm-19-cuda-12.9_tpls/llvm-19.1.3/mathpresso-geos-m5g4ftkfck3mm3e66xdckii3dy4skjtm" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")
