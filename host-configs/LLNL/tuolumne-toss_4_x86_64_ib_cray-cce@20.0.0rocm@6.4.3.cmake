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

set(BLT_CXX_STD "c++20" CACHE STRING "")
#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/cray/20.0/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/cray/20.0/bin/mpicxx" CACHE PATH "")
set(MPIEXEC "srun" CACHE STRING "")
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
# ROCm/HIP
#--------------------------------------------------------------------------------

set(ENABLE_HIP ON CACHE BOOL "")
set(CMAKE_HIP_STANDARD "20" CACHE STRING "")
set(CMAKE_HIP_COMPILER "/opt/rocm-6.4.3/bin/hipcc" CACHE PATH "")
set(CMAKE_HIP_ARCHITECTURES "gfx942" CACHE STRING "")
set(GPU_TARGETS "gfx942" CACHE STRING "")
set(AMDGPU_TARGETS "gfx942" CACHE STRING "")
set(ROCM_PATH "/opt/rocm-6.4.3" CACHE PATH "")
#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/chai-2026.07.0-bq4ryit6nm7bvbxfvhfqjbzqfdbw7gph" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/raja-2026.07.0-f3jtrixuf73bf3i3425en36vmppvewq7" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/umpire-2026.07.1-3nmb55tc2co4pwx2rjci4ljygr64espm" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/camp-2026.07.1-bk433mseis7t76udur56vhihauvdagv7" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/caliper-2.14.0-kbcanpm7bzcc4777gpqku2jl36ljou6v" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/adiak-0.4.0-wlig6ibpw5hvzp4ntfbhgk2w3wffvwhd" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hdf5-1.14.6-klbzlt347ci4gzv5yygc7bbe4fjr2cxi" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/conduit-0.9.5-jzuyvsppsgw5rhnqtn2vimjygenxcogk" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/silo-4.12.0-27btgh5efwcmt7i6nmhdle3c7rkq6plx" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/pugixml-1.13-e3pmghwewcexjn2kz6s3szgy4wpkqdos" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/vtk-9.7.0-xecjel6gtqmgwe5rndm65oz7bdr3zddh" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/fmt-12.1.0-5rnwcvpwh3zz4rrlbumqopuvgthhehbt" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/metis-5.1.0-xo5ilqtsu22vqo3ap7ofbfnpbljsd3xc" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/parmetis-4.0.3-l3jd56x35efw4qywynw4r2zak4wj6jww" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/scotch-7.0.8-2npx2xrhwj4auha6ucxmrmjpeux56o4u" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/superlu-dist-9.2.1-qvndf4hzwau2unfteuvvdrzike42xkap" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/suite-sparse-5.10.1-v64mai2jmeey7jm6l4psipm5ljvuc3lp" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypre-git.f1374fb6182c9e730abaa82f865a87b36f9ad50a_master-nmvgpeziaqadsgh6vmdj2jo73sxqdqmi" CACHE PATH "")
set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/hypredrive-git.98989bef31e865d738c1b678ee058bb0e5dd635e_master-gvzyslft3kysy2prv2ubfwsajuy6npkm" CACHE PATH "")
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
set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-tn3ncjjuhqn3ol3jac5qisybktur24c2/bin/uncrustify" CACHE PATH "")
#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/tuolumne-cce-20-rocm-6.4.3_tpls/cce-20.0.0/mathpresso-geos-4cc3g6a2qymg3l7ps2nadn5rrbjh6gsn" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
