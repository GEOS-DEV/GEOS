# Specific ordering to build object libraries
set( subdirs
     math
     cxx-utilities
     common
     linearAlgebraInterface
     codingUtilities
     fileIO
     dataRepository
     managers
     constitutive
     finiteElement
     finiteVolume
     mesh
     meshUtilities
     MPI_Communications
     physicsSolvers
     wells
     python
     )

set( coreLibLinkLine "" PARENT_SCOPE)
set( geosx_core_include_dirs "" )
set( geosx_core_defines "" )

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
  set( extraComponentsLinkList ${extraComponentsLinkList} )
endif()


blt_add_library ( NAME geosx_core
                  SOURCES dummy.cpp
                  DEPENDS_ON  ${coreLibs}
                  SHARED TRUE
                )

#set( geosx_core_libs geosx_core CACHE PATH "" FORCE) 

target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks(PREFIX core )

add_test( NAME testUncrustifyCheck
          COMMAND ${CMAKE_MAKE_PROGRAM} common_uncrustify_check dataRepository_uncrustify_check
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          )

add_subdirectory( unitTests )
