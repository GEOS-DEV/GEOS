#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: cce@20.0.0/jsxhkuinfphvs3f63jhiol4t24cdwomd
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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/chai-2025.12.0-4vm3j6ldec5ff7ynrjvlreivoqko4dud" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/raja-2025.12.0-76bbsdqjimlfobyfgzmihb74rss57rdz" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/umpire-2025.12.0-l5csty6zd7kni55knksm5rxkpub2ak5v" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/camp-2025.12.0-4q4op6apqykoaou3joqiep4pfj4kk65p" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/caliper-2.14.0-ddbvcxro7xql2luo4wtbsaidlzlhyb7m" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/adiak-0.4.0-zv47sfcwzudtrk2ivemdyajlckq4pon5" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hdf5-1.14.6-3ibotnphfnhuygyzokkhmwiqpmeald7s" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/conduit-0.9.5-bvumr65f4dzktb5iobjhnpzv6ixwmyac" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/silo-4.12.0-6237ggfbykukb3h2cky6krkm4ch2jymd" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/pugixml-1.13-e3pmghwewcexjn2kz6s3szgy4wpkqdos" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/vtk-9.4.2-rbtq744shnttfnsowdpmpvmasz3d4f6o" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/fmt-11.2.0-f3ryh7fjo7ivavi64qeysr2dacahptbi" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/metis-5.1.0-xo5ilqtsu22vqo3ap7ofbfnpbljsd3xc" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/parmetis-4.0.3-mpxoxzshxm43om6nglgeoyvuvc2ducof" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/scotch-7.0.8-fqkfrn3fe3m565ybbjjueu4caabho6i4" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/superlu-dist-9.2.1-pixwzad2cf64k5pd6f4gjedaa2hwsito" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/suite-sparse-5.10.1-fwkhklsjinduguawym7mi2k7go7knilw" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-bvg6pjtba4fwdwgcwqajirjum6k33adf" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-3yyxtxatnqfxnohp7olwyezxuawukrsh" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-tn3ncjjuhqn3ol3jac5qisybktur24c2/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-17/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/mathpresso-geos-4cc3g6a2qymg3l7ps2nadn5rrbjh6gsn" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

