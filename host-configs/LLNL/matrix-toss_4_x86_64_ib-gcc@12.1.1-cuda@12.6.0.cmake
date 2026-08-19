#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@12.1.1/juhzhguo5cx4jdi4eeq2i5tmtrhqwpvv
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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/chai-2025.12.0-tcvodmrdymi3a5sdbd26mzb4ltpxedsg" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/raja-2025.12.0-tysowueoi52e6dldv7rc72uyj5rw6hev" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/umpire-2025.12.0-gkjrtdopalcxq5qrs5innyiwl75ulnlb" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/camp-2025.12.0-pb46uhuxkjs52ckqqnzxj5u7klke3tk2" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/caliper-2.14.0-mds3ilcapuq6d6xhaoez2fhqcm74gl2k" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/adiak-0.4.0-w24kil5cdhfoghlnvjsduqe3v6cgrmxw" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hdf5-1.14.6-kp5nz5niclj56qdwea7wt766dfvubldk" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/conduit-0.9.5-wkwt2pbolfwfpyfzc2xanjuysw6bahij" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/silo-4.12.0-6p6aqqnpqjtsyzofqlma5sjkyvzlcaht" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/pugixml-1.13-l6ugbcq26b2424xde5w66vrabv4zcgu5" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/vtk-9.4.2-276jg4hcm7affh7zas5vfvo6ouu6hpsm" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/fmt-11.2.0-uidxs2v3pdxvvrbrihwui3shhmcumdsj" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/metis-5.1.0-lnhdc7ibpdw7c3grhslxc56hphtxcsmc" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/parmetis-4.0.3-y55zlqz7hkjej5xsl66teevhkt3quqmz" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/scotch-7.0.8-sq2iimq7dhjmbptm276nxup2cplkfsc4" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/superlu-dist-9.2.1-ym374pfv4ubu4mkv43pxw5uxlpboaetv" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/suite-sparse-5.10.1-hsiexjmfd2ouykgq2ibkyuem2nk4uwzc" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-2d7w7adfjv7gecq744oe3fseh6l5nmpf" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "CUDA" CACHE STRING "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-fxyi34uxlhdu6nav3zzmx4pzkllny4tc" CACHE PATH "")

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

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/matrix-gcc-12-cuda-12.6_tpls/gcc-12.1.1/mathpresso-geos-5hxcbnh5k2suv2evhmqzibqafzdsi3x5" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")

