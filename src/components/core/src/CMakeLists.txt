

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


set( coreLibLinkLine "" PARENT_SCOPE)
foreach( lib ${coreLibs} )
    add_subdirectory(${lib})
    if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
      if( ${BUILD_SHARED_LIBS} OR ${ENABLE_SHARED_LIBS} )
        set( coreLibLinkLine ${coreLibLinkLine} ${lib} )      
      else()
        set( coreLibLinkLine ${coreLibLinkLine} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )      
      endif()
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


set( geosx_core ${coreLibLinkLine} ${extraComponentsListLink} CACHE PATH "" FORCE)
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

#message( "coreLibLinkLine = ${coreLibLinkLine}" )
#add_library ( geosx_core SHARED )
#target_link_libraries(geosx_core ${coreLibLinkLine})

#execute_process(COMMAND ${CMAKE_AR} rcs ${TARGET_LOCATION} ${all_object_locations}
#                WORKING_DIRECTORY ${scratch_dir} RESULT_VARIABLE exe_result)
#  if(NOT "${exe_result}" STREQUAL "0")
#    message(FATAL_ERROR "COMBINE_UNIX_ARCHIVES: archive construction process failed exe_result='${exe_result}'")
#endif()

geosx_add_code_checks(PREFIX core )

