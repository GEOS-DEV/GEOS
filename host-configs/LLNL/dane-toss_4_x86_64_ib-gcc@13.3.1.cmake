#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@13.3.1/uy2jpfm5ainwsww2bboqa727eyowk7c2
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

set(ENABLE_CUDA OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# ROCm/HIP
#--------------------------------------------------------------------------------

set(ENABLE_HIP OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/chai-git.4b9060b18b9bec1167026cfb3132bd540c4bd56b_develop-c6kv3hjgayugw3m54whnvqrwv763et3n" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/raja-git.1d70abf171474d331f1409908bdf1b1c3fe19222_develop-iucoqkktaz627ezgtlyd2cv7bjd5ivdy" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/umpire-git.1ed0669c57f041baa1f1070693991c3a7a43e7ee_develop-d7dgxnh7hujciq5utyw4ylg2rb2ropxk" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/camp-git.ee0a3069a7ae72da8bcea63c06260fad34901d43_main-tcbkhejrohvuizsbmqewihnki7ko2fdt" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/caliper-git.287b7f3ad2d12f520aad04268d44f353cd05403c_2.12.0-4mfajxjfzt45oamvc5bhigyinly5v3oo" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/adiak-0.4.0-bvqjdhele3xlnjovgesbiuiykjuavcwm" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/hdf5-1.14.6-l557c2ijv4nstgnlm3lsqnouiy4476w2" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/conduit-git.ad86e316ad56a75c099d30ca5ce75cff275b5924_develop-y4ngi2h5jrzoxnin3i5u2yk7p3zszfzp" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/silo-4.12.0-pl4e72jltht7z2ffufy774uq5wuptyhz" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/pugixml-1.13-uekbyefkladmhu3txy7fsrphjqqwpc3d" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/vtk-9.4.2-csi6mfw6xz6cdc6fkvedj6wejaev3m2j" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/fmt-10.0.0-oyro5yd5ntu3yp3tjoc2q6mvg2xaoy4s" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/metis-5.1.0-gub5ytd4ni6besdgq4a7v7634da4lrum" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/parmetis-4.0.3-l2457z6c5si54kdzwidospl5xqyshgnj" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/scotch-7.0.8-p3qtu42xktbvvp3sz4dksh2u6fggq6lw" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/superlu-dist-9.2.1-l7j4zhv6q6vpp2bzzfzzaamh26xoxto4" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/suite-sparse-5.10.1-cg225lh2z6cwornathi4obwcud4x2xgd" CACHE PATH "")

set(TRILINOS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/trilinos-16.1.0-x7lb6xh2vn3nazrtg53e26cggt3k3ie6" CACHE PATH "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/hypre-git.8b0093306228fef1b92384d9face7fbe5a63b460_master-bhenjbz47e5h4j4e5mm5mk4n7ncsevv2" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/doxygen-1.8.20-llexzql2nksyaxm2nn3cwqwwskgbcvmj/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-qn72owjuzr7toqjrykjju7dp3xi4qtwr/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# addr2line
#--------------------------------------------------------------------------------

set(ENABLE_ADDR2LINE ON CACHE BOOL "")

set(ADDR2LINE_EXEC  "/usr/bin/addr2line" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-02-24/dane-gcc-13_tpls/gcc-13.3.1/mathpresso-geos-xc5mfpfv55njpzic7bgfnn2ld55k6k4m" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

