#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: clang@=14.0.6
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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-rpztyhc3d6ar4qd7pvhgiurzhbkp42j3" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-aezavc2a3xo65mlwe7f43ruklbm65rr7" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-gsdr3wacbaatethtrelqsmulthkihzly" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-dwd6adgrnz7rwqbrwsnl2xb3q56dmhgt" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-jgcljxcmnwpmu2w3e6fwepntqqjt274h" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/adiak-0.4.0-shzoi42twk7k6gbnkjwsbdddzu4rhrkw" CACHE PATH "")

set(ZLIB_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/zlib-1.3.1-2p2eqwjlvvz42sc2q5ml2yw2eymz5kog" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/hdf5-1.12.1-4sbm456runiwk6j7qazogadcz435xvrg" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-5fv2jm2puzfkclt4k7jc6nnvbfyzezim" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/silo-4.11.1-bsd-jwyj3qy5zssnecqucxc3jxnotuyxyh6w" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/pugixml-1.13-octpifpder5dpo6tjewt2n3vl2n2vlok" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/vtk-9.4.2-l4y7h3mnnwnagx5uk5u5l2qokoaqyd3l" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/fmt-10.0.0-4gzwrffoa2pl6v7uwyynj5yegriur7ik" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2022.1.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_intel_lp64.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_core.so
                  /lib/../lib64/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/metis-5.1.0-ar4qtsbmnbmtlul5qx3pofbmf3jakucm" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/parmetis-4.0.3-xvx5rgufkonp2talkmpnzp452pgebq55" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/scotch-7.0.8-gdgkqvornryacilixrkvn7sm2po74jsj" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-fqieoarboqxtgg5zhahckdj4l54unqq3" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/suite-sparse-5.10.1-ccwpkt522s5m3xwuvcj2dn4ctqykfshz" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/trilinos-16.1.0-cnv2mfchuntcex5yhk2fhrlyxeu7xexo" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-eqounuyoubqxwe4gnx4id352f5h5v7iz" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/matrix-clang-14-cuda-12.6_tpls/clang-14.0.6/mathpresso-geos-5kkbozv2oghijthzn3etg66evaegymsx" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

