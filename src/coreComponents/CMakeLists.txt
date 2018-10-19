

set( subdirs
     codingUtilities
     common
     constitutive
     dataRepository
     fileIO
     finiteElement
     finiteVolume
     math
     mesh
     managers
     meshUtilities
     MPI_Communications
     physicsSolvers
     python
     systemSolverInterface
     linearAlgebraInterface
     cxx-utilities
     )

set( coreLibLinkLine "" PARENT_SCOPE)
set( geosx_core_include_dirs "" )
set( geosx_core_defines "" )

foreach( lib ${subdirs} )
    add_subdirectory(${lib})
    if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
      list( APPEND coreLibs ${lib} )
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

if( GEOSX_USE_ATK )
  set( extraComponentsListLink ${extraComponentsListLink} )
endif()

blt_combine_static_libraries( NAME          geosx_core
                              SOURCE_LIBS   ${coreLibs}
                              LIB_TYPE SHARED
                              LINK_PREPEND  ${GEOSX_LINK_PREPEND_FLAG} 
                              LINK_POSTPEND ${GEOSX_LINK_POSTPEND_FLAG}
                            )

#set( geosx_core_libs geosx_core CACHE PATH "" FORCE) 

target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks(PREFIX core )


add_subdirectory( unitTests )
