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

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-f4wgpqutnikc45ymbfntdnbwcbpxp6ww" CACHE PATH "")

set(RAJA_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-dn7yccm4x5q2kz6zl7pinkmtmksf3yoc" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-k6gdhiuktlcu5tfk4vuhzyqhj57ew7it" CACHE PATH "")

set(CAMP_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-e2mec4ypeahwfmjk2ntrbonnjlo6cahm" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-d6jzlhx4ofqairetauzhiw7rcocb3lj5" CACHE PATH "")

set(adiak_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/adiak-0.4.0-vxgcukydyqrm6ciuxrw4leflputw6r5e/lib/cmake/adiak" CACHE PATH "")

set(HDF5_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/hdf5-1.12.1-yqvru4fil7emvnajdmu3tq736rlzesza" CACHE PATH "")

set(CONDUIT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-sowkah5jcivfrcna3lbmxxhntokwin6n" CACHE PATH "")

set(SILO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/silo-4.11.1-bsd-dmyx6w72zap3cqtrtrdoaggki6dbdywv" CACHE PATH "")

set(PUGIXML_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/pugixml-1.13-qdosppqlevwbzkbkr25tycjlts6dxuoq" CACHE PATH "")

set(VTK_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/vtk-9.3.1-sgjxp6i2a5hlvzux3ufnmvrr3zoa22su" CACHE PATH "")

set(FMT_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/fmt-10.0.0-f2lzpypwoyhiracsh44luu46ezmch72t" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(ENABLE_MKL ON CACHE BOOL "")

set(MKL_INCLUDE_DIRS "/usr/tce/packages/mkl/mkl-2022.1.0/include" CACHE PATH "")

set(MKL_LIBRARIES /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gf_lp64.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_gnu_thread.so
                  /usr/tce/packages/mkl/mkl-2022.1.0/mkl/2022.1.0/lib/intel64/libmkl_core.so
                  /lib/../lib64/libomp.so
                  /lib64/libpthread.so
                  /lib64/libm.so
                  /lib64/libdl.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/metis-5.1.0-2bggpuqbh4sjefm2hisblqpofiek27dg" CACHE PATH "")

set(PARMETIS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/parmetis-4.0.3-dy34gwdh5lk4bmzrucykto667hgj7njl" CACHE PATH "")

set(SCOTCH_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/scotch-7.0.3-vquie4qxvhrhcbsaagdkjlwrvxyasj7r" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-uf3dqwltl5532qvkrcrse2nzk32drtv5" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/suite-sparse-5.10.1-kgaikvqvpf5quqvhu4bfm3tbixwcoyfo" CACHE PATH "")

set(TRILINOS_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/trilinos-16.0.0-kt3gwm2xmyhg3t2z4pfwx3w75wkbu3ga" CACHE PATH "")

set(HYPRE_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/hypre-git.21e5953ddc6daaa24699236108866afa597a415c_2.32.0-git.33-4ssxay255dsjyqvkgcfcpgsbp5mktmic" CACHE PATH "")

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

set(SPHINX_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/python/quartz-gcc-python/python/bin/sphinx-build" CACHE PATH "")

set(DOXYGEN_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/doxygen-1.8.20-xj3wmf2ztwfreg3ao4zdjyjlami7lrb4/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-3tyetsqkpeonbslmn2iithvvvyfhozc6/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/gapps/GEOSX/thirdPartyLibs/2025-04-03_spack/ruby-clang-14_tpls/clang-14.0.6/mathpresso-geos-pxokzqriqhovt3ggd37a3ocrhzp3p7ho" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

