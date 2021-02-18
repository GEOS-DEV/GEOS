set( subdirs
     math
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
     python
     virtualElement
     )


set( coreLibLinkLine "" PARENT_SCOPE)
set( geosx_core_include_dirs "" )
set( geosx_core_defines "" )


set( buildAsObj FALSE)
if ( BUILD_OBJ_LIBS)
  set(buildAsObj TRUE)
endif()

list( APPEND coreLibs lvarray)
add_subdirectory(LvArray)

if( ENABLE_PYGEOSX )
  list( APPEND coreLibs pylvarray )
endif()

foreach( lib ${subdirs} )
  if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
    list( APPEND coreLibs ${lib})
  endif()
endforeach()

if( NOT GEOSX_BUILD_SHARED_LIBS )
  set( geosx_core_list "" )
  foreach( lib ${coreLibs} )
    list( APPEND geosx_core_list ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG})
  endforeach()
  set (geosx_core_libs "${geosx_core_list}" CACHE INTERNAL "")
endif()

foreach( lib ${subdirs} )
    add_subdirectory(${lib})
endforeach()

set( extraComponentsLinkList "" )

if ( ENABLE_OPENMP )
  set( extraComponentsLinkList ${extraComponentsLinkList} openmp )
endif()

if ( ENABLE_CUDA )
  set( extraComponentsLinkList ${extraComponentsLinkList} cuda )
endif()

if( ENABLE_MPI )
  set( extraComponentsLinkList ${extraComponentsLinkList} mpi )
endif()

if( GEOSX_BUILD_SHARED_LIBS )
  set( GEOSX_BUILD_SHARED_BOOL TRUE CACHE BOOL "" )
  set( CSL_TYPE SHARED )
else()
  set( GEOSX_BUILD_SHARED_BOOL FALSE CACHE BOOL "" )
  set( CSL_TYPE STATIC )
endif()

if( BUILD_OBJ_LIBS )
  blt_add_library ( NAME geosx_core
                    SOURCES dummy.cpp
                    DEPENDS_ON  ${coreLibs}
                    SHARED ${GEOSX_BUILD_SHARED_BOOL}
                  )
else()
  blt_combine_static_libraries( NAME          geosx_core
                                SOURCE_LIBS   ${coreLibs}
                                LIB_TYPE ${CSL_TYPE}
                                LINK_PREPEND  ${GEOSX_LINK_PREPEND_FLAG}
                                LINK_POSTPEND ${GEOSX_LINK_POSTPEND_FLAG}
                              )
endif()

# To install the library with the runtime path used during the building
set_target_properties( geosx_core PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE )

# To install the shared library
install(TARGETS geosx_core LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)

target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/coreComponents)


geosx_add_code_checks( PREFIX coreComponents 
                       EXCLUDES cmake )

if( ENABLE_UNCRUSTIFY )
  add_test( NAME testUncrustifyCheck
            COMMAND sh -c "${CMAKE_MAKE_PROGRAM} coreComponents_uncrustify_check main_uncrustify_check 2>&1 >/dev/null"
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
endif()

add_subdirectory( unitTests )
