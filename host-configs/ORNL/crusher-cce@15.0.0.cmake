set(CCE_VERSION 15.0.0)

include(${CMAKE_CURRENT_LIST_DIR}/../../src/coreComponents/LvArray/host-configs/ORNL/crusher-cce@${CCE_VERSION}.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/crusher-base.cmake)

set( CONDUIT_DIR "${GEOSX_TPL_DIR}/conduit" CACHE PATH "" )
set( HDF5_DIR "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )

# set(ENABLE_SILO FALSE CACHE BOOL "" )
set(ENABLE_SILO TRUE CACHE BOOL "" )
set(SILO_DIR "${GEOSX_TPL_DIR}/silo" CACHE PATH "" )

set(ENABLE_VTK FALSE CACHE BOOL "" )
# set(ENABLE_VTK TRUE CACHE BOOL "" )
# set(VTK_DIR "${GEOSX_TPL_DIR}/vtk" CACHE PATH "" )
# # these are needed due to how vtk searches for dependencies.. hahahahahaha
# set(GLEW_ROOT "${GEOSX_TPL_DIR}/glew" CACHE PATH "" )
# set(OSMESA_ROOT "/sw/crusher/spack-envs/base/opt/linux-sles15-x86_64/gcc-7.5.0/mesa-21.2.1-cbfrg4yfizu6ouirpoql7j7l62emid5c/" CACHE PATH "" )
# set(HDF5_ROOT "${GEOSX_TPL_DIR}/hdf5" CACHE PATH "" )
# set(NetCDF_ROOT "${GEOSX_TPL_DIR}/netcdf-c" CACHE PATH "" )
# set(THEORA_ROOT "${GEOSX_TPL_DIR}/libtheora" CACHE PATH "" )
# set(OGG_ROOT "${GEOSX_TPL_DIR}/libogg" CACHE PATH "" )
# set(JsonCpp_ROOT "${GEOSX_TPL_DIR}/jsoncpp" CACHE PATH "" )
# set(GL2PS_ROOT "${GEOSX_TPL_DIR}/gl2ps" CACHE PATH "" )
# set(PNG_ROOT "${GEOSX_TPL_DIR}/libpng" CACHE PATH "" )
# set(pugixml_DIR "${GEOSX_TPL_DIR}/pugixml/lib64/cmake/pugixml" CACHE PATH "" )
# set(SQLite3_ROOT "${GEOSX_TPL_DIR}/sqlite" CACHE PATH "" )
# set(LibPROJ_ROOT "${GEOSX_TPL_DIR}/proj" CACHE PATH "" )
# set(X11_ROOT "${GEOSX_TPL_DIR}/libx11" CACHE PATH "" )
# set(EXPAT_ROOT "${GEOSX_TPL_DIR}/expat" CACHE PATH "" )
# set(double-conversion_ROOT "${GEOSX_TPL_DIR}/double-conversion" CACHE PATH "" )
# set(LZ4_ROOT "${GEOSX_TPL_DIR}/lz4" CACHE PATH "" )
# set(utf8cpp_ROOT "${GEOSX_TPL_DIR}/utf8cpp" CACHE PATH "" )
# set(Eigen3_ROOT "${GEOSX_TPL_DIR}/eigen" CACHE PATH "" )
# set(JPEG_ROOT "${GEOSX_TPL_DIR}/libjpeg-turbo" CACHE PATH "" )
# set(TIFF_ROOT "${GEOSX_TPL_DIR}/libtiff" CACHE PATH "" )
# set(Freetype_ROOT "${GEOSX_TPL_DIR}/freetype" CACHE PATH "" )
# set(SEACASIoss_ROOT "${GEOSX_TPL_DIR}/seacas" CACHE PATH "" )
# set(CGNS_ROOT "${GEOSX_TPL_DIR}/cgns" CACHE PATH "" )

set( BLAS_DIR "/opt/rocm-${ROCM_VERSION}/" CACHE PATH "" )

set( PUGIXML_DIR "${GEOSX_TPL_DIR}/pugixml" CACHE PATH "" )
set( FMT_DIR "${GEOSX_TPL_DIR}/fmt" CACHE PATH "" )
set( SUITESPARSE_DIR "${GEOSX_TPL_DIR}/suite-sparse" CACHE PATH "" )

# HYPRE options
# set( ENABLE_HYPRE_DEVICE "CPU" CACHE STRING "" )
set( ENABLE_HYPRE_DEVICE "HIP" CACHE STRING "" )

# set( ENABLE_HYPRE_MIXINT FALSE CACHE STRING "" )
set( ENABLE_HYPRE_MIXINT TRUE CACHE STRING "" )

# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-18-ga592bbd12_cce-15.0.0_rel" CACHE PATH "" ) # CPU SERIAL (WORKS)
# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-200-ge907ce401_cce-15.0.0_mixint_rel" CACHE PATH "" ) # cpu serial mixint (working depending on solver config)

# set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/hypre_v2.27.0-201-g6ad0909ec_cce-15.0.0_rocm-5.4.0_mixint_umpire-2022.3.0_caliper_rel" CACHE PATH "" ) # rocm, no unified memory, int64/int32 (working dependeing on solver config)
set( HYPRE_DIR "/gpfs/alpine/geo127/world-shared/hypre/install-amd-gpu/" CACHE PATH "" )
set( ENABLE_CALIPER ON CACHE BOOL "" FORCE )
set( CALIPER_DIR "${GEOSX_TPL_DIR}/caliper" CACHE PATH "" )
