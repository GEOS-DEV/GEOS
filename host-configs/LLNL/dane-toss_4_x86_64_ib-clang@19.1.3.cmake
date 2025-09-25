#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: clang@=19.1.3
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

set(BLT_CXX_STD "c++17" CACHE STRING "")

#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3-magic/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3-magic/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-ggrnk5rsnvzuqf6om3cdjajchpxihrza" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-cs5kvvqnh2wl2jrf7occuti2yhiwvbuq" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-slqtula7bnjuvddhzgsxaklgvrusqlov" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-tjmotdymah2kgcepcscw42qteury5vg7" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-lq3stlxw4jc6gkrjo7ltgyezhuzcyks7" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/adiak-0.4.0-5idtmidayajcjfk2xioba4mfgbuy22sy" CACHE PATH "")

set(ZLIB_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/zlib-1.3.1-fyqxsc4kj5j5ktc6k5jzljordczfrpq7" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/hdf5-1.12.1-bjb5ar6qmtqmmomz5b47br2hrqzdy2lm" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-gsc2vgc6tdehpvrfuvvdodyzo664veoa" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/silo-4.11.1-bsd-opb2yktncpevjf3o5vypxpvt5d2ttp4m" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/pugixml-1.13-e242tpuklhfyhmeqxsdnn3ymg27o3dpl" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/vtk-9.4.2-ue42puketxfzrudwvtle6a5x6ppvow7j" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/fmt-10.0.0-nqkdz2rhyvtyjt3yygsewulw326ekno6" CACHE PATH "")

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

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/metis-5.1.0-k5jawftu22mvnwcdb7o4vgiy5tseipgm" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/parmetis-4.0.3-l4nr3kqb3q3zwou7miaz43vq2qkjmnta" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/scotch-7.0.8-wzdshazr5pw5cw2445atxbqnlbq2ptpe" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/superlu-dist-git.0f6efc377df2440c235452d13d28d2c717f832a1_6.3.0-git.8-g4utcupo4mxwlczbkzxj7tgf4htl7onm" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/suite-sparse-5.10.1-u3pgd5dblz7nio5s4jbfxsgfqsr5anmq" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/trilinos-16.1.0-ztes25cijyufrb25fnbo4um473wpnjvj" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/hypre-git.907a2d07b64fe47bdde4540c54665c83ced83a2c_2.33.0-git.20-xeixfkejyg2zvxgjkizecladcennahsv" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/doxygen-1.8.20-pybghjqgvozfobh4n2gqpstqaqyqtr2u/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-2xgcoibttoazu3bdbk2gj3atqc5mnale/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2025-09-23/dane-clang-19_tpls/clang-19.1.3/mathpresso-geos-z4ylqrbgcoqyeyoafiwsifmvc5za4xhs" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm56" CACHE STRING "")

