
#
# Specify all headers
#
set(componentHeaders "NewComponent.hpp" )

#
# Specify all sources
#

set(componentSources "NewComponent.cpp" )
    

blt_add_library( NAME
                     componentName
                 SOURCES
                     "${componentSources}"
                 HEADERS
                     "${componentHeaders}"
                 DEPENDS_ON
                     geosx_core
                 )
                 