
#
# Specify all headers
#
set(common_headers
    ${CMAKE_BINARY_DIR}/include/common/GeosxConfig.hpp
    BufferAllocator.hpp
    DataTypes.hpp
    Path.hpp
    GeosxMacros.hpp
    Stopwatch.hpp
    TimingMacros.hpp
    Logger.hpp
    DataLayouts.hpp
    Tensor.hpp
   )

#
# Specify all sources
#
set(common_sources
    BufferAllocator.cpp
    DataTypes.cpp
    Logger.cpp
    Path.cpp
   )

set( dependencyList lvarray pugixml RAJA chai conduit::conduit )

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
  set( dependencyList ${dependencyList} caliper )
endif()

blt_add_library( NAME                  common
                 SOURCES               ${common_sources}
                 HEADERS               ${common_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
#                 SHARED TRUE
               )

target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents )

geosx_add_code_checks(PREFIX common )

add_subdirectory( unitTests )
