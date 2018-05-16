
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )




set( coreLibs 
     codingUtilities
     common
     constitutive
     dataRepository
     fileIO
     finiteElement
     math
     mesh
     managers
     meshUtilities
     MPI_Communications
     physicsSolvers
     python
     systemSolverInterface
     )


set( coreLibLinkLine "")
foreach( lib ${coreLibs} )
    add_subdirectory(${lib})
    if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
        set( coreLibLinkLine ${coreLibLinkLine} $<TARGET_OBJECTS:${lib}> )
    endif()
endforeach()



# Python config
set( extraComponentsListLink "")
#if( ENABLE_PYTHON )
#  set( extraComponentsListLink ${extraComponentsListLink} python_interp )
#endif()

if( ENABLE_OPENMP)
  set( extraComponentsListLink ${extraComponentsListLink} openmp )
endif()

if( ENABLE_MPI )
  set( extraComponentsListLink ${extraComponentsListLink} mpi )
endif()

if( USE_ATK )
  set( extraComponentsListLink ${extraComponentsListLink} )
endif()

#message( "extraComponentsListLink = ${extraComponentsListLink}" )
message("adding geosx_core library")

#blt_add_library( NAME                  geosx_core
#                 SOURCES               ${geosx_core_sources}
#                 HEADERS               ${geosx_headers}
#                 DEPENDS_ON            ${coreLibs}
#                                       ${extraComponentsListLink} 
#                                       ${thirdPartyLibs}
#                 SHARED TRUE
#               )


message( "coreLibLinkLine = ${coreLibLinkLine}" )
add_library ( geosx_core SHARED ${coreLibLinkLine} )
target_link_libraries(geosx_core ${coreLibLinkLine})

geosx_add_code_checks(PREFIX core )

