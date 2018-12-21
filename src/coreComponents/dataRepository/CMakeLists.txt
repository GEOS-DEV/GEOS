
#
# Specify all headers
#
set(dataRepository_headers
    BufferOps_inline.hpp
    BufferOps.hpp
    InputFlags.hpp
    KeyIndexT.hpp
    ManagedGroup.hpp
    ExecutableGroup.hpp
    MappedVector.hpp
    ReferenceWrapper.hpp
    SidreWrapper.hpp
    ViewWrapper.hpp
    ViewWrapperBase.hpp
    RestartFlags.hpp )

#
# Specify all sources
#
set(dataRepository_sources
    ManagedGroup.cpp
    SidreWrapper.cpp
    ViewWrapper.cpp
    ViewWrapperBase.cpp )
    


set( extraComponentsListLink "")


if( ENABLE_MPI )
  set( extraComponentsListLink ${extraComponentsListLink} mpi )
endif()




set( dependencyList common cxx-utilities codingUtilities pugixml fileIO )
if( CONDUIT_FOUND )
  set( dependencyList ${dependencyList} conduit )
endif()

if( ATK_FOUND )
  set( dependencyList ${dependencyList} sidre )
endif()

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            ${dependencyList}
               )

target_include_directories( dataRepository PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks(PREFIX dataRepository )

add_subdirectory( unitTests )
