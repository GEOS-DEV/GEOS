
#
# Specify all headers
#
set(componentHeaders "NewComponent.hpp" )

#
# Specify all sources
#

set(componentSources "NewComponent.cpp" )
    
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

message("CMAKE_CURRENT_SOURCE_DIR 2 = ${CMAKE_CURRENT_SOURCE_DIR}")
message("CMAKE_SOURCE_DIR 2 = ${CMAKE_SOURCE_DIR}")

blt_add_library( NAME
                     componentName
                 SOURCES
                     "${componentSources}"
                 HEADERS
                     "${componentHeaders}"
                 DEPENDS_ON
                     geosx_core raja
                 )
                 