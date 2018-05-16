
#
# Specify all headers
#
set(dataRepository_headers
    CommBufferOps.hpp
    CommBufferOps_inline.hpp
    KeyIndexT.hpp
    ManagedGroup.hpp
    MappedVector.hpp    
    SidreWrapper.hpp
    ViewWrapper.hpp
    ViewWrapperBase.hpp
    Buffer.hpp
    RestartFlags.hpp )

#
# Specify all sources
#
set(dataRepository_sources
    CommBufferOps.cpp
    ManagedGroup.cpp
    SidreWrapper.cpp
    ViewWrapper.cpp
    ViewWrapperBase.cpp )
    


set( extraComponentsListLink "")


if( ENABLE_MPI )
  set( extraComponentsListLink ${extraComponentsListLink} mpi )
endif()

if( USE_ATK )
  set( extraComponentsListLink ${extraComponentsListLink} )
endif()

#message( "extraComponentsListLink = ${extraComponentsListLink}" )
#message("adding geosx_core library")
#
#set( dependencyList common cxx-utilities  )
#set( dependencyList2 pugixml )
#foreach( lib ${dependencyList} )
#    set( dependencyList2 ${dependencyList2} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )
#endforeach()
#message( "dependencyList2 = ${dependencyList2}" )

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            common cxx-utilities pugixml
                 SHARED FALSE
               )

target_include_directories( common PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)

geosx_add_code_checks(PREFIX dataRepository )