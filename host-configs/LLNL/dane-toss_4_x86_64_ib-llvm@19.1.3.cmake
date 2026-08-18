#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: llvm@19.1.3/hjgi63xrm7o26tcgiguvstjhzvyeda3s
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

set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-clang-19.1.3/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/chai-2025.12.0-xlhj7sclyhnsiulvp6qc5c4ezlfbaiew" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/raja-2025.12.0-bkzlpnqnixwi2ogm2k3lf2randqmyicg" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/umpire-2025.12.0-7zrjp32yg6qqgl2smkzf6dcmciskh37v" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/camp-2025.12.0-e76op5cdiskl624xxxvsyf4hmwwuun6r" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/caliper-2.14.0-xgp2z6wwfczkqvtjwtbczwuy4d5b4wew" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/adiak-0.4.0-l4rn2eei2l5js3s6vt32a7mrh4dbbojy" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/hdf5-1.14.6-uls2v2o6f6hp46g6bvgsi6sn2xqndym2" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/conduit-0.9.5-u6iyb6fyjwfma4tlrdyqm5iufqeqypux" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/silo-4.12.0-4aojvzxfzbmb7yhggqizphrhk6q3n7yh" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/pugixml-1.13-llxtretjci2iarvgcsdcjvayyxyo4s7s" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/vtk-9.4.2-dtgbaevs2f6pnae4dmn5s4lfarmegbuu" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/fmt-11.2.0-y4etnmq2kvkm5ko6ex5vfvtnlvbtwagb" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/metis-5.1.0-sn2fyabl2rntxmxlqgmazrfdtj6fujob" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/parmetis-4.0.3-mpeju2kfj72ijsmamenmmudki7zn4qvu" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/scotch-7.0.8-ckpvu2ubc5nazyigjqiawgzb3aiswobf" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/superlu-dist-9.2.1-ekotqlk77zanouekc22otu4ppgicvrwc" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/suite-sparse-5.10.1-k5aes4hkcsouibigmhupy6vqsw7nzdcc" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-nvg2iom2kunatpesimhnxhovsx4uwkjc" CACHE PATH "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-2oghowp44bkklt4ur3a354olqg3amyw6" CACHE PATH "")

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

set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/gcc-13.3.1/doxygen-1.8.20-h6tfpgng4j7bbntgjqqwaaq42etylkqt/bin/doxygen" CACHE PATH "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-lpj2ebpivfqzdje35npcui5riil54cun/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/dane-llvm-19_tpls/llvm-19.1.3/mathpresso-geos-cxrldz2fh27jcxzlb74mqpd3tvacb5uz" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

