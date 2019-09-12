
#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    ExecutableGroup.hpp
    Group.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    MappedVector.hpp
    ObjectCatalog.hpp 
    ReferenceWrapper.hpp
    RestartFlags.hpp
    ConduitRestart.hpp
    WrapperHelpers.hpp
    Wrapper.hpp
    WrapperBase.hpp
    xmlWrapper.hpp
    )

#
# Specify all sources
#
set(dataRepository_sources
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



set( dependencyList common )

if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit )
endif()

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            ${dependencyList}
                 #OBJECT                TRUE
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)

geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
