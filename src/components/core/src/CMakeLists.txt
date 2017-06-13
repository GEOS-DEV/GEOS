
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )
    
add_subdirectory(codingUtilities)
add_subdirectory(common)
add_subdirectory(constitutive)
add_subdirectory(dataRepository)
add_subdirectory(fileIO)
add_subdirectory(finiteElement)
#add_subdirectory(legacy)
add_subdirectory(math)
add_subdirectory(schema)
add_subdirectory(managers)
add_subdirectory(MeshUtilities)
add_subdirectory(MPI_Communications)
add_subdirectory(PhysicsSolvers)
add_subdirectory(python)


# Python config
set( extraComponentsListLink "")
if( ENABLE_PYTHON )
  set( extraComponentsListLink ${extraComponentsListLink} python_interp )
endif()


if( ENABLE_MATHPRESSO )
  set( extraComponentsListLink ${extraComponentsListLink} mathpresso )
endif()


message("adding geosx_core library")
blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers}"
                 HEADERS_OUTPUT_SUBDIR
                     geosx_core
                 DEPENDS_ON
                     cxx-utilities 
                 SYSTEM_DEPENDS_ON
                     slic sidre conduit hdf5 pugixml ${extraComponentsListLink} ${thirdPartyLibs}
                 )
                 
#target_include_directories( geosx_core PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
#target_include_directories( geosx_core PUBLIC ${CMAKE_CURRENT_LIST_DIR}/legacy )
#target_include_directories( geosx_core PUBLIC /Users/rrsettgast/Codes/geos/gpaclib/include )

blt_add_executable(NAME geosx_corex
                   SOURCES main/main.cpp
                   DEPENDS_ON ${GEOSX_LINK_PREPEND_FLAG} geosx_core ${GEOSX_LINK_POSTPEND_FLAG}
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja chai cxx-utilities ${extraComponentsListLink}
                   )
#target_include_directories( geosx_corex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} )
                   
#set_target_properties(core  PROPERTIES COMPILE_FLAGS )
#install(TARGETS core DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
                   
