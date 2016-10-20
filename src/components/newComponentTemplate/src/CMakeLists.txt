
#
# Specify all headers
#
set(componentHeaders "NewComponent.hpp" )

#
# Specify all sources
#

set(componentSources "NewComponent.cpp" )

blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

message("adding newComponentTemplate library")
blt_add_library( NAME             newComponentTemplate
                 SOURCES          ${componentSources}
                 HEADERS          ${componentHeaders}
                 DEPENDS_ON        geosx_core cxx-utilities
                 SYSTEM_DEPENDS_ON slic sidre conduit hdf5 RAJA chai
                 )

target_include_directories( newComponentTemplate PUBLIC ${CMAKE_CURRENT_LIST_DIR} )
target_include_directories( newComponentTemplate PUBLIC ${CMAKE_SOURCE_DIR}/components/core/src )
