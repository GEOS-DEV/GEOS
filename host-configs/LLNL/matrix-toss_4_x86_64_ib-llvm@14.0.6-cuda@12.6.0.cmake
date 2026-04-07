#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@14.0.6/6hueebr6kxrx6r2je7lfkrj7gk6sdghn
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

set(CMAKE_CUDA_STANDARD "17" CACHE STRING "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/chai-2025.12.0-434reoatceiqmp5nhsvqhmqqh5v66ekl" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/raja-2025.12.0-ftze7ut2hidtld6s6557vwo7kqd7m2d2" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/umpire-2025.12.0-zvalmcwps5dsita5uhmu5g3zzwdbux7x" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/camp-2025.12.0-hzukgfyf3dys4u6ehe6tjpq7r2cwpjgb" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/caliper-2.14.0-oahv5uuriqn2vyqbzespzkrfsghxjfew" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/adiak-0.4.0-iitjnfmanurcrmlyau4qm7kpep6q6vxt" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hdf5-1.14.6-ottiokh52c2uvkujg6s4yc63afy24e7d" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/conduit-0.9.5-y5dwa2wdrmelzcpyx7yzayckzo3gpgrm" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/silo-4.12.0-24vfiuq62fnn5ui6uqdhvouz65zytin2" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/pugixml-1.13-rd24lxieweuz24w7ycbaryfyk6lbznmb" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/vtk-9.4.2-fnpua5soyrcqbkrqnnh4zqt6yn5hem5e" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/fmt-11.0.2-voxs6pm2ttqbwt3gszcbvcb2xmugiosj" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/metis-5.1.0-7ultypmkrm5nnur3o6ntl557ut4zpi3a" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/parmetis-4.0.3-l7zdwdi3icsmd2ei7lmnytikc4kitsrf" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/scotch-7.0.8-kfvnwp3qofbqvm4duccjctx23zsjfyj3" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/superlu-dist-9.2.1-qyykp4rzu22h7hzag5xubujl2oxl6iqu" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/suite-sparse-5.10.1-6s5wuk2jv66hevvxm4tjnrxmad6cuqw3" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/trilinos-16.1.0-ugaxqluqwxpbnupdzgkjqjcjx42k5d4n" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/hypre-git.8b0093306228fef1b92384d9face7fbe5a63b460_master-rp25dcvyczbf7kusq2gn52o2fqufn5eo" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

set(ENABLE_HYPREDRV OFF CACHE BOOL "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-04-03/matrix-llvm-14-cuda-12.6_tpls/llvm-14.0.6/mathpresso-geos-kxwkykrmexwzjykwccwn3p4hl54nsmym" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

