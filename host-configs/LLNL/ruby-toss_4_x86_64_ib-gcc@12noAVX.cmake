#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@=12noAVX
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-10.3.1/cmake-3.26.3-nz532rvfpaf5lf74zxmplgiobuhol7lu/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/gcc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/gcc/gcc-12.1.1-magic/bin/g++" CACHE PATH "")

set(CMAKE_CXX_FLAGS "-march=x86-64-v2 -mno-avx512f" CACHE PATH "")

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

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-44liirgqt77y5jx6chuabzw6v5rpj22r" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-jkwf235s5cewqgnoe6me7gv4hty5hfuw" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-f5nd3amxl5vi32gpxuebjy3hza3i46xr" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-ozknqlixdd5t2fggyzzahyp2uwgtpqxe" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-fpt5fqjl57zcfp7juq2723exiqbiafue" CACHE PATH "")

set(adiak_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/adiak-0.4.0-ic6ousqruhancnroudaj2ramn26nnzxt/lib/cmake/adiak" CACHE PATH "")

set(ZLIB_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/zlib-1.3.1-yed2ynyxzfkhibckigaw4ev7r4urixwx" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/hdf5-1.12.1-gde7rzb5axcjosct7jnzil5mqznuphym" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-xwnrnqd7jpyvg7e6ak3n2lrjb2teflv3" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/silo-4.11.1-bsd-qgi3yuajumhnem5jyfxudgaxymsxsry2" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/pugixml-1.13-czxwgb3p5ogoeiqftrp3bqahzr6lgr5n" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/vtk-9.4.2-74r4y5vhyb24cm4ljesn3zmj4t7a2rhg" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/fmt-10.0.0-2mwtki3p7acfmyns4d27f6pvtwgbxdaw" CACHE PATH "")

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

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/metis-5.1.0-hemq6k3gyqbppgsvok3b2kcrtx24oglo" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/parmetis-4.0.3-jbpzgcnoj3trxe3jjjbxxqjxgm4666h2" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/scotch-7.0.3-frgphccdsmcv57rmsmizi7kzn7ggyazg" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-kutepd3wey2t6z72ln7tcyyeclxpfj4m" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/suite-sparse-5.10.1-i6e7qqrkbistqhfir3lmlhw5oc5ghwoh" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/trilinos-16.1.0-ummkoh33hrslvibofwc6m4xktisxltn2" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/hypre-git.be52325a3ed8923fb93af348b1262ecfe44ab5d2_2.33.0-git.12-wytsytvye27jyj64tncsyavp56yfdafb" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/doxygen-1.8.20-gimiza4fo3e3l6hk7is2i6ldsx5tzjdz/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-eecq2jjtthgyzmdi3sk5enloz6xw5bzw/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-05-05/ruby-gcc-12noAVX_tpls/gcc-12noAVX/mathpresso-geos-vqn7cjtawyzaitoivbgkxqbclcdhog2j" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

