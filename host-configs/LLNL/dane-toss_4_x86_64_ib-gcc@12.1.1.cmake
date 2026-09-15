#################################################################################
# Generated host-config - Edit at own risk!
#################################################################################
#--------------------------------------------------------------------------------
# SYS_TYPE: toss_4_x86_64_ib
# Compiler Spec: gcc@12.1.1/fsawo3vnk2dqu4cbg7rziuvd4qxm6vfs
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

set(BLT_CXX_STD "c++20" CACHE STRING "")
#--------------------------------------------------------------------------------
# MPI
#--------------------------------------------------------------------------------

set(ENABLE_MPI ON CACHE BOOL "")
set(MPI_C_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicc" CACHE PATH "")
set(MPI_CXX_COMPILER "/usr/tce/packages/mvapich2/mvapich2-2.3.7-gcc-12.1.1-magic/bin/mpicxx" CACHE PATH "")
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

set(ENABLE_HIP OFF CACHE BOOL "")
#--------------------------------------------------------------------------------
# Performance Portability TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CHAI ON CACHE BOOL "")
set(CHAI_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/chai-2026.07.0-arw7jdudvycraljfktabld4sb2of6o3e" CACHE PATH "")
set(RAJA_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/raja-2026.07.0-wrium6b6635mflhweziytv5prwjkcwsy" CACHE PATH "")
set(ENABLE_UMPIRE ON CACHE BOOL "")
set(UMPIRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/umpire-2026.07.1-4q4nq74i6ailhner73qwqw6d4r2jhai2" CACHE PATH "")
set(CAMP_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/camp-2026.07.1-hx7cyxvstw43hq53oyydcqkq2yewljat" CACHE PATH "")
#--------------------------------------------------------------------------------
# IO TPLs
#--------------------------------------------------------------------------------

set(ENABLE_CALIPER ON CACHE BOOL "")
set(CALIPER_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/caliper-2.14.0-wafnkyyvt3d6bpcgybrvxz5ihymutmm2" CACHE PATH "")
set(ADIAK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/adiak-0.4.0-burrsvbtmzww3jfh2konwvbt3xzctddz" CACHE PATH "")
set(ZLIB_DIR "/usr" CACHE PATH "")
set(HDF5_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/hdf5-1.14.6-4z7uyydq22vfjhcu3te3p4a7lpzkvs62" CACHE PATH "")
set(CONDUIT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/conduit-0.9.5-oebbjwcq2klhalpbbsvesgpba6egkxsb" CACHE PATH "")
set(SILO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/silo-4.12.0-mwsrjsnnd4tqnukagvjj7ujc7djwk6tc" CACHE PATH "")
set(PUGIXML_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/pugixml-1.13-5qhtoeyb5yfogjqvjnl5dpuxtl6dllkk" CACHE PATH "")
set(VTK_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/vtk-9.7.0-23skvwl7s3xxsehpxqj7nsytdi53ujrp" CACHE PATH "")
set(FMT_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/fmt-12.1.0-gtzebp3hnmzxf2kvqucafutphmhnnjod" CACHE PATH "")
#--------------------------------------------------------------------------------
# System Math Libraries
#--------------------------------------------------------------------------------

set(BLAS_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

set(LAPACK_LIBRARIES /usr/lib64/libopenblas.so CACHE STRING "")

#--------------------------------------------------------------------------------
# Math TPLs
#--------------------------------------------------------------------------------

set(METIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/metis-5.1.0-znhibjpryqyq2pp5kjsx7e6xwvnezemr" CACHE PATH "")
set(PARMETIS_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/parmetis-4.0.3-imhjqrmgovranp6v6ebenueucxoa6sw4" CACHE PATH "")
set(SCOTCH_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/scotch-7.0.8-fp4wqafwsmtktisjoo5ahgkd3mpybt7d" CACHE PATH "")
set(SUPERLU_DIST_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/superlu-dist-9.2.1-7jgnu7ut74ng5hd6wt2kg3xb43xhvxh4" CACHE PATH "")
set(SUITESPARSE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/suite-sparse-5.10.1-fm6efx2e5yjn23talubodykqxvxc25xu" CACHE PATH "")
set(ENABLE_TRILINOS OFF CACHE BOOL "")
set(HYPRE_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/hypre-git.9fbaf60dc9435e71ff5af984f1e12e2bf8be6ad8_master-gyqf5o522eztxjwd5bo6jfrj4kqxwonl" CACHE PATH "")
set(ENABLE_HYPREDRV ON CACHE BOOL "")
set(HYPREDRV_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/hypredrive-git.4eb4f1b126332844feaaf941e32ae5dc125e5bdc_master-izaqyrygfxrcyv3lzm7fblwfrwctkqjx" CACHE PATH "")
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
set(DOXYGEN_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/doxygen-1.8.20-6wkrdmcf4ghqjrtt6hr3cd7cdxjz7whn/bin/doxygen" CACHE PATH "")
#--------------------------------------------------------------------------------
# Development tools
#--------------------------------------------------------------------------------

set(ENABLE_UNCRUSTIFY ON CACHE BOOL "")
set(UNCRUSTIFY_EXECUTABLE "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/uncrustify-git.401a4098bce9dcc47e024987403f2d59d9ba7bd2_0.70.1-git.319-s4vgg4xqfynoelpw7b7k7od5yv65y2jq/bin/uncrustify" CACHE PATH "")
#--------------------------------------------------------------------------------
# Other
#--------------------------------------------------------------------------------

set(ENABLE_MATHPRESSO ON CACHE BOOL "")
set(MATHPRESSO_DIR "/usr/WS1/GEOS/GEOSX/TPLs_2026-08-25/dane-gcc-12_tpls/gcc-12.1.1/mathpresso-geos-exvyulsbndfdllbs46fw6ga6vric6ki7" CACHE PATH "")
set(ENABLE_XML_UPDATES ON CACHE BOOL "")
set(ENABLE_GRPC OFF CACHE BOOL "")
set(GEOS_BUILD_SHARED_LIBS ON CACHE BOOL "")
set(ATS_ARGUMENTS "--machine slurm112" CACHE STRING "")
