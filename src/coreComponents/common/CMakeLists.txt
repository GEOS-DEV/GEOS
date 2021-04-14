
#
# Specify all headers
#
set(common_headers
    ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
    BufferAllocator.hpp
    DataLayouts.hpp
    DataTypes.hpp
    FieldSpecificationOps.hpp
    GeosxMacros.hpp
    initializeEnvironment.hpp
    Logger.hpp
    MpiWrapper.hpp
    Path.hpp
    GEOS_RAJA_Interface.hpp
    Stopwatch.hpp
    Tensor.hpp
    TimingMacros.hpp
   )

#
# Specify all sources
#
set(common_sources
    BufferAllocator.cpp
    DataTypes.cpp
    initializeEnvironment.cpp
    Logger.cpp
    MpiWrapper.cpp
    Path.cpp
   )

set( dependencyList lvarray pugixml RAJA chai conduit::conduit )

if ( ENABLE_MKL )
  set( dependencyList ${dependencyList} mkl )
endif()

if( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif( )

if( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif( )

if( ENABLE_MPI )
  set( dependencyList ${dependencyList} mpi )
endif()

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

blt_add_library( NAME                  common
                 SOURCES               ${common_sources}
                 HEADERS               ${common_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${GEOSX_BUILD_OBJ_LIBS}
#                 SHARED TRUE
               )

target_include_directories( common PUBLIC ${CMAKE_BINARY_DIR}/include )
target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks(PREFIX common )

add_subdirectory( unitTests )
