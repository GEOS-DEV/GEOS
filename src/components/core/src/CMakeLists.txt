

set( subdirs
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
set( geosx_core_include_dirs "" )
set( geosx_core_defines "" )

set( coreLibs "" )
foreach( lib ${subdirs} )
    add_subdirectory(${lib})
    if( NOT ( ${lib} STREQUAL "math" OR ${lib} STREQUAL "python" ) )
      list( APPEND coreLibs ${lib} )
      if( ${BUILD_SHARED_LIBS} OR ${ENABLE_SHARED_LIBS} )
        set( coreLibLinkLine ${coreLibLinkLine} ${lib} )      
      else()
        set( coreLibLinkLine ${coreLibLinkLine} ${GEOSX_LINK_PREPEND_FLAG} ${lib} ${GEOSX_LINK_POSTPEND_FLAG} )      


        get_target_property( liblist ${lib} INCLUDE_DIRECTORIES )
        message( "INCLUDE_DIRECTORIES for ${lib}: ${liblist}")
        foreach( libdir ${liblist} )
          if ( NOT ${libdir} IN_LIST geosx_core_include_dirs )
            list( APPEND geosx_core_include_dirs ${libdir} )
          endif()
        endforeach()
        
        get_target_property( definelist ${lib} INTERFACE_COMPILE_DEFINITIONS )
        if( NOT ${definelist} MATCHES "NOTFOUND" )
        foreach( definesymbol ${definelist} )
          if ( NOT ${definesymbol} IN_LIST geosx_core_defines )
            list( APPEND geosx_core_defines ${definesymbol} )
          endif()
        endforeach()
        endif()
        
        
      endif()
    endif()
endforeach()
#message( "geosx_core INCLUDE_DIRECTORIES = ${geosx_core_include_dirs}" )
#message( "geosx_core INTERFACE_COMPILE_DEFINITIONS = ${geosx_core_defines}" )



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


if( TRUE )
  message("adding geosx_core library as a single dynamic lib")
  add_library ( geosx_core SHARED )
  target_link_libraries(geosx_core PRIVATE ${coreLibLinkLine})
  set( geosx_core_libs geosx_core CACHE PATH "" FORCE)

  if( FALSE )
    target_include_directories( geosx_core SYSTEM PUBLIC ${geosx_core_include_dirs} )
    target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src)
    target_include_directories( geosx_core PUBLIC ${CMAKE_SOURCE_DIR}/components/cxx-utilities/src/src)
    target_include_directories( geosx_core SYSTEM PUBLIC ${CMAKE_SOURCE_DIR}/thirdparty/optionparser/src)
    target_compile_definitions( geosx_core PUBLIC ${geosx_core_defines})
  else()
    blt_register_library( NAME geosx_core
                          DEPENDS_ON ${coreLibs}
                          LIBRARIES geosx_core
                        )
  endif()                
else()
  message("adding geosx_core library as a collection of static libs")
  set( geosx_core ${coreLibLinkLine} ${extraComponentsListLink} CACHE PATH "" FORCE)
  set( geosx_core_libs ${geosx_core} CACHE PATH "" FORCE)
  
endif()

geosx_add_code_checks(PREFIX core )

