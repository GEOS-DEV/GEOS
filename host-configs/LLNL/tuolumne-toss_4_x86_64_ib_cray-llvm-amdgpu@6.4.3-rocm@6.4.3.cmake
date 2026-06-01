#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib_cray
# Compiler Spec: llvm-amdgpu@6.4.3/saipjvhhes22pxj5rnkjw46erthbkkln
# CMake executable path: /usr/tce/backend/installations/linux-rhel8-x86_64/gcc-8.3.1/cmake-3.24.2-ywx52e32uh6gkxzuyubpwkulzgdvxyh6/bin/cmake
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Compilers
#--------------------------------------------------------------------------------

set(CMAKE_C_COMPILER "/opt/rocm-6.4.3/llvm/bin/amdclang" CACHE PATH "")

set(CMAKE_CXX_COMPILER "/opt/rocm-6.4.3/llvm/bin/amdclang++" CACHE PATH "")

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

set(MPI_C_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicc" CACHE PATH "")

set(MPI_CXX_COMPILER "/opt/cray/pe/mpich/9.0.1/ofi/amd/6.0/bin/mpicxx" CACHE PATH "")

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

set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/chai-2025.12.0-bui2qjacbspelk7pwk3ap474padtsriy" CACHE PATH "")

set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/raja-2025.12.0-phw2uaei64dmmcxcle5ct3jsso262eny" CACHE PATH "")

set(ENABLE_UMPIRE ON CACHE BOOL "")

set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/umpire-2025.12.0-7nyetx7q5rljh7z2bcchdutu7ugkxpas" CACHE PATH "")

set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/camp-2025.12.0-n4t7sxktutnfffuwc2e6z2hya4dbxiti" CACHE PATH "")

#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")

set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/caliper-2.14.0-nchdbsxkkfwgrezu5azi4iaiwnjx6zy3" CACHE PATH "")

set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/adiak-0.4.0-x7wwsvhjq2ehtw4r7jvzjecdpcxiwfl7" CACHE PATH "")

set(ZLIB_DIR "/usr" CACHE PATH "")

set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hdf5-1.14.6-qqfx3gxy555optxhwmffsjtwuo7wpwlo" CACHE PATH "")

set(ENABLE_ADIOS2 ON CACHE BOOL "")

set(ADIOS2_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/adios2-2.12.1-xzmxlump2er7gcy65dpdxbaof4dkczm2" CACHE PATH "")

set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/conduit-0.9.5-qqbspejhqzbjqxr7antnfuedgxehjysl" CACHE PATH "")

set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/silo-4.12.0-ale723hmqzp6vwhbqp3r7xnslpj73pww" CACHE PATH "")

set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/pugixml-1.13-yo2g6mdyjdcz3saiq6cewsmilagincwn" CACHE PATH "")

set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/vtk-9.4.2-7rwq4ausmeiweksccifubmdxgb6kmgjw" CACHE PATH "")

set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/fmt-11.0.2-vjwb72ugeomwqz4kcgt64id6mnkgqzhp" CACHE PATH "")

#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/liblapack.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/metis-5.1.0-fah7yojtygjosrxaywxuhy6arel2a5o3" CACHE PATH "")

set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/parmetis-4.0.3-ztyzko3dlrcb5mfhkvxdwvgoyyq7wbl3" CACHE PATH "")

set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/scotch-7.0.8-5t6uibtbtbfdbkisch7rgxvy44pgw6pm" CACHE PATH "")

set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/superlu-dist-9.2.1-m66d72gxlkboj7tdkfmuxww2auhzraii" CACHE PATH "")

set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/suite-sparse-5.10.1-iwxia6d76g2r6dwnto2iajidfws2j3ts" CACHE PATH "")

set(ENABLE_TRILINOS OFF CACHE BOOL "")

set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypre-git.8b0093306228fef1b92384d9face7fbe5a63b460_master-jiobpmxaz3fp5dp6wt5tyqewmzv3a3q6" CACHE PATH "")

set(ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "")

set(ENABLE_HYPREDRV ON CACHE BOOL "")

set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/hypredrive-git.d062c48cc8d26a9bee92f33f3c588d071c230d30_master-wvj777kfsqo3vihtwzulpt5ovdrh32ne" CACHE PATH "")

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

set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-b7rjjhdbf6yqcn5br6rzvsngyvb2i4ha/bin/uncrustify" CACHE PATH "")

#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")

set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-05-31/tuolumne-llvm-amdgpu-6.4.3_tpls/llvm-amdgpu-6.4.3/mathpresso-geos-nolryyeptxklghofeiszkkyyhd6don35" CACHE PATH "")

set(ENABLE_XML_UPDATES ON CACHE BOOL "")

set(ENABLE_GRPC OFF CACHE BOOL "")

set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")

