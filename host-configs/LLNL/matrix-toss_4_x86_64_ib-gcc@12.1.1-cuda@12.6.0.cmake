#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@=12.1.1
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

set(MPIEXEC_EXECUTABLE "/usr/bin/srun" CACHE PATH "")
set(MPIEXEC_NUMPROC_FLAG "-n" CACHE STRING "")

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
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-ae5pyj2jzvvrdicofvbx2kpw5mhtamdp" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-wvx5iofamobz7sse4y2vye4yg6lk2fjn" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-jhsr6nwtyfvb2i5ud62e6yccktxhgc7h" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-afnsvlgatskbbvxyp7dx7neopqftwp24" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-xolotrctq2ml5kcg4dminlshg6vzagor" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/adiak-0.4.0-bz6cshrbl4zpowg6uvhdo7ki7ss2tgjr" CACHE PATH "")

set(ZLIB_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/zlib-1.3.1-54oviesqjo6upo65ry5fitqz5omfw57n" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hdf5-1.12.1-oghlibs4dsel6eqrwu2xp2gvg4pjz2md" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-37y7d7y2e735gaelsdl2ub7itkak6a77" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/silo-4.11.1-bsd-c4eskt44znon3uq54ozahcwpq6y4uusl" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/pugixml-1.13-hsfh2l3cmzc7piwx4ex3azvz5tk75tfs" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/vtk-9.4.2-ygyhmgwabjr7arolnf46wopotcbmlsl5" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/fmt-10.0.0-e7nurcfb4jtbpu75wmofg7hdurzjaf5n" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2022.1.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_intel_lp64.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_core.so
                  /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-13.3.1/llvm-19.1.3-gy2lu5xbi4csr2k47emlajzfs5mlsd4g/bin/../lib/x86_64-unknown-linux-gnu/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/metis-5.1.0-b2tyegtpk5zhjutrvdsgukfjuzfl35bg" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/parmetis-4.0.3-wc7b6rp7ap6mtmu2zuucl7dwd3d72giq" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/scotch-7.0.8-wipiko44zee5banzpq4q4eororx5uycl" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-czyec6xb66h4ls7pbli6iv7efc7na6kp" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/suite-sparse-5.10.1-4giwaxmvci2yc7uo27hxwwpdhdb2pcuq" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/trilinos-16.1.0-4pgm7gjhs3ksivasgfzkquyipnxr5tgl" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-bbkm7xoooezhtxgfmg35khi3mkzhjwwz" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

set(ENABLE_PETSC OFF CACHE BOOL "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/mathpresso-geos-pmas5poag5egtkkogmkxzbowtlsnhznm" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

