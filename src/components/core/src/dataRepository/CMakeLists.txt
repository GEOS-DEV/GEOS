
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
message("adding geosx_core library")

blt_add_library( NAME                  dataRepository
                 SOURCES               ${dataRepository_sources}
                 HEADERS               ${dataRepository_headers}
                 DEPENDS_ON            common cxx-utilities pugixml
                                       ${extraComponentsListLink} 
                 SHARED FALSE
               )
                 
geosx_add_code_checks(PREFIX dataRepository )