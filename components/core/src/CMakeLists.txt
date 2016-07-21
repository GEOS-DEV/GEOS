
#
# Specify all headers
#
set(geosx_headers "" )

#
# Specify all sources
#

set(geosx_core_sources "" )
    
add_subdirectory(dataRepository)
add_subdirectory(codingUtilities)
add_subdirectory(PhysicsSolvers)


blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_CURRENT_SOURCE_DIR})

blt_add_library( NAME
                     geosx_core
                 SOURCES
                     "${geosx_core_sources}"
                 HEADERS
                     "${geosx_headers}"
                 HEADERS_OUTPUT_SUBDIR
                     geosx_core
                 SYSTEM_DEPENDS_ON
                     slic sidre conduit hdf5 raja
                 )

blt_add_executable(NAME geosx_corex
                   SOURCES main/main.cpp
                   DEPENDS_ON ${GEOSX_LINK_PREPEND_FLAG} geosx_core ${GEOSX_LINK_POSTPEND_FLAG}
                   SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja 
                   )
                   
#set_target_properties(core  PROPERTIES COMPILE_FLAGS )
#install(TARGETS core DESTINATION lib EXPORT ${PROJECT_NAME}-targets)
                   