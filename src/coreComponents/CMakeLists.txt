set( subdirs
     math
     cxx-utilities
     common
     linearAlgebra
     codingUtilities
     fileIO
     dataRepository
     managers
     constitutive
     finiteElement
     finiteVolume
     mesh
     meshUtilities
     mpiCommunications
     physicsSolvers
     wells
     python
     )

#set( subdirs
#     codingUtilities
#     common
#     constitutive
#     dataRepository
#     fileIO
#     finiteElement
#     finiteVolume
#     math
#     mesh
#     managers
#     meshUtilities
#     MPI_Communications
#     physicsSolvers
#     wells
#     python
#     linearAlgebra
#     cxx-utilities
#     )

set( coreLibLinkLine "" PARENT_SCOPE)
set( geosx_core_include_dirs "" )
set( geosx_core_defines "" )


set( buildAsObj FALSE)
if ( BUILD_OBJ_LIBS)
  set(buildAsObj TRUE)
endif()

foreach( lib ${subdirs} )
    add_subdirectory(${lib})
    if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python") )
      list( APPEND coreLibs ${lib})
    endif()
endforeach()

# Python config
set( extraComponentsLinkList "" )
#if( ENABLE_PYTHON )
#  set( extraComponentsLinkList ${extraComponentsLinkList} python_interp )
#endif()

if ( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

if ( ENABLE_CUDA )
  set( extraComponentsLinkList ${extraComponentsLinkList} cuda )
endif()

if( ENABLE_MPI )
  set( extraComponentsLinkList ${extraComponentsLinkList} mpi )
endif()

if( GEOSX_USE_ATK )
  set( extraComponentsLinkList ${extraComponentsLinkList} axom )
endif()


if( BUILD_OBJ_LIBS )
  blt_add_library ( NAME geosx_core
                    SOURCES dummy.cpp
                    DEPENDS_ON  ${coreLibs}
                    SHARED TRUE
                  )
else()
  blt_combine_static_libraries( NAME          geosx_core
                                SOURCE_LIBS   ${coreLibs}
                                LIB_TYPE SHARED
                                LINK_PREPEND  ${GEOSX_LINK_PREPEND_FLAG}
                                LINK_POSTPEND ${GEOSX_LINK_POSTPEND_FLAG}
                              )
endif()

#set( geosx_core_libs geosx_core CACHE PATH "" FORCE)

target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks( PREFIX core )

set( enforceUncrustifyLibs
     common
     dataRepository
     cxx-utilities )

if( ENABLE_UNCRUSTIFY )
  foreach( lib ${enforceUncrustifyLibs} )
    add_test( NAME testUncrustifyCheck_${lib}
              COMMAND ${CMAKE_MAKE_PROGRAM} ${lib}_uncrustify_check
              WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            )
  endforeach()
endif()

add_subdirectory( unitTests )
