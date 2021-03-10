#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    BufferOpsDevice.hpp
    ExecutableGroup.hpp
    Group.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    MappedVector.hpp
    ObjectCatalog.hpp
    ReferenceWrapper.hpp
    RestartFlags.hpp
    ConduitRestart.hpp
    wrapperHelpers.hpp
    Wrapper.hpp
    WrapperBase.hpp
    xmlWrapper.hpp
    )

#
# Specify all sources
#
set(dataRepository_sources
    BufferOpsDevice.cpp
    ExecutableGroup.cpp
    Group.cpp
    ConduitRestart.cpp
    Wrapper.cpp
    WrapperBase.cpp
    xmlWrapper.cpp
    )

set( extraComponentsLinkList "")

if( ENABLE_MPI )
  set( extraComponentsLinkList ${extraComponentsLinkList} mpi )
endif()

if( BUILD_OBJ_LIBS )
  set (dependencyList common lvarray codingUtilities pugixml fileIO)
else()
  set( dependencyList 
       common 
       mpiCommunications
       managers
       codingUtilities)
endif()

if ( ENABLE_OPENMP )
  set( dependencyList ${dependencyList} openmp )
endif()

if ( ENABLE_CUDA )
  set( dependencyList ${dependencyList} cuda )
endif()

set( dependencyList ${dependencyList} common conduit )

if( ENABLE_PYGEOSX )
  set( dependencyList ${dependencyList} pylvarray )
endif()

if( ENABLE_CALIPER )
  set( dependencyList ${dependencyList} caliper adiak )
endif()

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            ${dependencyList}
                 OBJECT                ${buildAsObj}
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
