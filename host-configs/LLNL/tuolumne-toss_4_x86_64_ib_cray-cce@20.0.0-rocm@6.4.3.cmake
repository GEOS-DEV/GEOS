#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: cce@20.0.0/3lsh34k6vtdv3cgqje3mxi5tudetifrc
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-8.3.1/cmake-3.24.2-ywx52e32uh6gkxzuyubpwkulzgdvxyh6/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/usr/tce/packages/cce-tce/cce-20.0.0/bin/craycc" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/usr/tce/packages/cce-tce/cce-20.0.0/bin/crayCC" CACHE PATH "")

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

set(MPI_C_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/cray/20.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/cray/20.0/bin/mpicxx" CACHE PATH "")

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

set(ENABLE_HIP ON CACHE BOOL "")

set(CMAKE_HIP_STANDARD "17" CACHE STRING "")

set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.3/bin/hipcc" CACHE PATH "")

set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")

set(ROCM_PATH "/opt/rocm-6.4.3" CACHE PATH "")

#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/chai-2025.12.0-xnyu3aqqd7benfr7fkwlates66mruuth" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/raja-2025.12.0-amgdvrnkwl2oh5lpddgvyupteldgl6ip" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/umpire-2025.12.0-l7ni3umhaf7f4rwoy4qhjfyrnye3oit6" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/camp-2025.12.0-p2srwdveimzi7csrjw7o3u4ryxpaz5wm" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/caliper-2.14.0-xfkj7r2zaqfuuugu2wdkppydyr3fwrph" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/adiak-0.4.0-6rim24h5wdcrhoaxgpw5mldcufzkfqxb" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hdf5-1.14.6-b55srpmbowvfsql65n2652uh5xj3n3b5" CACHE PATH "")

set(ENABLE_ADIOS2 ON CACHE BOOL "")

set(ADIOS2_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/adios2-2.12.1-2ji4ttcd6rwaubb4v63jbd6isdhtuls5" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/conduit-0.9.5-74psmq2xghvj4rfuqymomncufbxjxb7b" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/silo-4.12.0-dvjbdqytwnsndxfskcpkdodgj6zr7k5n" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/pugixml-1.13-o2jle7fb6idvyr32ibcsvm552yq5xqbe" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/vtk-9.4.2-ztzvsmp6mlmtyk5o3rpqbgn7phg7ervu" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/fmt-11.0.2-x7b5m6zdgekxl3mmjyampxc2pmt4lgur" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/metis-5.1.0-oo45sjxzdfx4e4z3qgrmnudzvwfboi7n" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/parmetis-4.0.3-onkwfgqwdwkund4tmpvzo3gif3cewf7y" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/scotch-7.0.8-wphygatd54wcwziiuvi63jzkzx3pvrp7" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/superlu-dist-9.2.1-t3fmcjdkmgtbi3il5fji6t4mihzrxzni" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/suite-sparse-5.10.1-x3vbz2scbieouiqnnqxoirfsvrx7fj5t" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypre-git.8b0093306228fef1b92384d9face7fbe5a63b460_master-bwdj6fuiqwuyd25pegpk3s6zu6bfzvju" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypredrive-git.d062c48cc8d26a9bee92f33f3c588d071c230d30_master-kruuawah55xa3lyknsq6e7e6qyow7fok" CACHE PATH "")

set(ENABLE_PETSC OFF CACHE BOOL "")

set(ENABLE_CALIPER_HYPRE ON CACHE BOOL "")

set(GEOS_LA_INTERFACE "Hypre" CACHE STRING "")

#--------------------------------------------------------------------------------
# Python
#--------------------------------------------------------------------------------

set(Python3_ROOT_DIR "/usr/tce/packages/python/python-3.9.12" CACHE PATH "")

set(Python3_EXECUTABLE "/usr/tce/packages/python/python-3.9.12/bin/python3" CACHE PATH "")

set(ENABLE_PYGEOSX OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Documentation
#--------------------------------------------------------------------------------

set(ENABLE_DOCS OFF CACHE BOOL "")

set(ENABLE_DOXYGEN OFF CACHE BOOL "")

set(ENABLE_SPHINX OFF CACHE BOOL "")

#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-6hvpggvwvhbhd322pnqivylw2ywdkhaq/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/mathpresso-geos-w3ssu5ocek533pxesneyk3urbxhnhgmw" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

