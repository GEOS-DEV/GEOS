
#
# Specify all headers
#
set(componentHeaders "NewComponent.hpp" )

#
# Specify all sources
#

set(componentSources "NewComponent.cpp" )
    
blt_append_custom_compiler_flag(FLAGS_VAR CMAKE_CXX_FLAGS DEFAULT -I${CMAKE_SOURCE_DIR}/components/core/src)

blt_add_library( NAME             componentName
                 SOURCES          ${componentSources}
                 HEADERS          ${componentHeaders}
                 DEPENDS_ON        geosx_core
                 SYSTEM_DEPENDS_ON slic sidre conduit hdf5 raja
                 )
                 